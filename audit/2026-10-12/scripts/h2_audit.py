"""AUDIT reproduction of the H2 benchmark (not project code).

Model: H2, STO-3G, neutral singlet, nonrelativistic Born-Oppenheimer electronic Hamiltonian
projected onto the 2 spatial (4 spin) orbitals of the minimal basis. "Exact" below means exact
diagonalization of THIS finite 4-qubit Hamiltonian (= FCI in STO-3G), not the exact H2 energy.

Independent routes to the same number:
  (1) PySCF FCI (pyscf.fci)                       (2) analytic 2x2 CI in {|g^2>, |u^2>}
  (3) Pauli-string JW Hamiltonian (own code)      (4) dense JW from explicit ladder matrices (own code)
Usage:  python h2_audit.py [out.json]
"""
import json
import platform
import sys
import time

import numpy as np
import scipy
import pyscf
from pyscf import ao2mo, fci, gto, scf

import jw

TOL_E = 1e-9  # Eh, agreement required between independent routes


def integrals(mol, mf):
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri_ao = mol.intor("int2e")  # full (pq|rs) in AO basis, no permutational packing
    eri_own = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, eri_ao, optimize=True)
    # in-core PySCF transform from the SCF's packed AO integrals (avoids HDF5 temp files)
    eri_pyscf = ao2mo.restore(1, ao2mo.full(mf._eri, C), C.shape[1])
    return h1, eri_own, eri_pyscf


def run_point(R, detail=False):
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {R!r}", basis="sto-3g", unit="Angstrom",
                charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.chkfile = None
    mf.conv_tol = 1e-12
    e_rhf = float(mf.kernel())
    e_nuc = float(mol.energy_nuc())
    h1, eri, eri_pyscf = integrals(mol, mf)

    # RHF energy rebuilt from MO integrals (orbital 0 doubly occupied)
    e_rhf_rebuilt = e_nuc + 2 * h1[0, 0] + eri[0, 0, 0, 0]

    # (1) PySCF FCI
    e_fci = float(fci.FCI(mf).kernel()[0])

    # (2) analytic 2x2 CI (singles do not couple: g vs u symmetry)
    H11 = e_rhf_rebuilt
    H22 = e_nuc + 2 * h1[1, 1] + eri[1, 1, 1, 1]
    H12 = eri[0, 1, 0, 1]
    e_ci2 = float(np.linalg.eigvalsh(np.array([[H11, H12], [H12, H22]]))[0])

    # (3)+(4) Jordan-Wigner qubit Hamiltonian
    Hp = jw.molecular_qubit_hamiltonian(e_nuc, h1, eri)
    nq = 4
    M_pauli = jw.dense_from_paulis(Hp, nq)
    M_direct = jw.dense_hamiltonian_direct(e_nuc, h1, eri)
    N, Sz = jw.number_and_sz(nq)
    w, v = np.linalg.eigh(M_pauli)
    n_exp = np.real(np.einsum("ik,ij,jk->k", v.conj(), N, v))
    sz_exp = np.real(np.einsum("ik,ij,jk->k", v.conj(), Sz, v))
    sel = (np.abs(n_exp - 2) < 1e-8) & (np.abs(sz_exp) < 1e-8)
    e_jw = float(w[sel].min())
    sector_min = {int(round(n)): float(w[np.abs(n_exp - n) < 1e-8].min()) for n in range(nq + 1)}

    out = {
        "R_angstrom": R, "E_nuc": e_nuc, "E_RHF": e_rhf, "RHF_converged": bool(mf.converged),
        "E_FCI_pyscf": e_fci, "E_CI2x2_analytic": e_ci2, "E_JW_exact_N2_Sz0": e_jw,
        "E_corr": e_fci - e_rhf, "global_min_sector_N": min(sector_min, key=sector_min.get),
        "checks": {
            "rhf_rebuilt_minus_pyscf": float(e_rhf_rebuilt - e_rhf),
            "eri_own_vs_ao2mo_maxabs": float(np.abs(eri - eri_pyscf).max()),
            "h1_symmetric_maxabs": float(np.abs(h1 - h1.T).max()),
            "eri_8fold_maxabs": float(max(np.abs(eri - eri.transpose(1, 0, 2, 3)).max(),
                                          np.abs(eri - eri.transpose(0, 1, 3, 2)).max(),
                                          np.abs(eri - eri.transpose(2, 3, 0, 1)).max())),
            "pauli_coeff_max_imag": float(max(abs(c.imag) for c in Hp.values())),
            "dense_pauli_vs_direct_maxabs": float(np.abs(M_pauli - M_direct).max()),
            "hermiticity_maxabs": float(np.abs(M_pauli - M_pauli.conj().T).max()),
            "commutator_H_N_maxabs": float(np.abs(M_pauli @ N - N @ M_pauli).max()),
            "commutator_H_Sz_maxabs": float(np.abs(M_pauli @ Sz - Sz @ M_pauli).max()),
            "fci_vs_jw": e_fci - e_jw, "fci_vs_ci2": e_fci - e_ci2,
            "variational_rhf_minus_fci": e_rhf - e_fci,
        },
    }
    if detail:
        out["mo_energies"] = mf.mo_energy.tolist()
        out["h1_mo"] = h1.tolist()
        out["eri_mo_unique"] = {"(00|00)": eri[0, 0, 0, 0], "(11|11)": eri[1, 1, 1, 1],
                                "(00|11)": eri[0, 0, 1, 1], "(01|01)": eri[0, 1, 0, 1]}
        out["n_qubits"] = nq
        out["n_pauli_terms_incl_identity"] = len(Hp)
        out["pauli_terms"] = {jw.pauli_label(k, nq): float(c.real)
                              for k, c in sorted(Hp.items(), key=lambda kv: -abs(kv[1]))}
        out["pauli_1norm_excl_identity"] = float(sum(abs(c) for k, c in Hp.items() if k != (0, 0)))
        out["lowest_eigenvalue_by_particle_number"] = sector_min
    return out


def main():
    t0 = time.time()
    res = {
        "what": "Independent audit reproduction of H2/STO-3G (RHF, FCI, JW qubit Hamiltonian)",
        "model": {"molecule": "H2", "basis": "STO-3G", "units": "Angstrom", "charge": 0,
                  "spin_multiplicity": 1, "method": "RHF; FCI = exact diagonalization within STO-3G",
                  "geometry": "H at (0,0,0) and (0,0,R)",
                  "spin_orbital_order": "interleaved: qubit 2p = alpha(p), 2p+1 = beta(p)",
                  "mapping": "Jordan-Wigner, |1> = occupied"},
        "software": {"python": sys.version.split()[0], "numpy": np.__version__,
                     "scipy": scipy.__version__, "pyscf": pyscf.__version__,
                     "platform": platform.platform()},
        "tolerance_Eh": TOL_E,
    }
    res["single_points"] = [run_point(0.7414, detail=True), run_point(0.735, detail=True)]
    grid = [round(0.30 + 0.05 * i, 2) for i in range(55)]  # 0.30 ... 3.00 A
    scan = [run_point(R) for R in grid]
    res["scan"] = scan
    # minima on a fine grid around equilibrium
    fine = np.arange(0.680, 0.760001, 0.0005)
    fine_pts = [run_point(float(round(R, 4))) for R in fine]
    res["minimum_RHF_A"] = min(fine_pts, key=lambda d: d["E_RHF"])["R_angstrom"]
    res["minimum_FCI_A"] = min(fine_pts, key=lambda d: d["E_FCI_pyscf"])["R_angstrom"]
    # dissociation limit reference: one H atom in STO-3G (1 electron: HF is exact in the basis)
    hatom = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    mfh = scf.UHF(hatom)
    mfh.chkfile = None
    res["E_two_H_atoms_STO3G"] = 2 * float(mfh.kernel())
    allpts = res["single_points"] + scan + fine_pts
    worst = {k: max(abs(p["checks"][k]) for p in allpts)
             for k in res["single_points"][0]["checks"] if k != "variational_rhf_minus_fci"}
    res["worst_case_checks"] = worst
    res["all_routes_agree"] = bool(max(worst["fci_vs_jw"], worst["fci_vs_ci2"]) < TOL_E)
    res["variational_bound_ok"] = bool(min(p["checks"]["variational_rhf_minus_fci"] for p in allpts) > -TOL_E)
    res["global_min_always_in_N2"] = bool(all(p["global_min_sector_N"] == 2 for p in allpts))
    res["n_points_checked"] = len(allpts)
    res["wall_time_s"] = round(time.time() - t0, 2)
    out = sys.argv[1] if len(sys.argv) > 1 else "h2_audit.json"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    sp = res["single_points"]
    print(f"python {res['software']['python']} numpy {np.__version__} scipy {scipy.__version__} pyscf {pyscf.__version__}")
    for p in sp:
        print(f"R={p['R_angstrom']:.4f} A  E_nuc={p['E_nuc']:.10f}  E_RHF={p['E_RHF']:.10f}  "
              f"E_FCI={p['E_FCI_pyscf']:.10f}  E_JW={p['E_JW_exact_N2_Sz0']:.10f}  "
              f"E_CI2={p['E_CI2x2_analytic']:.10f}  qubits={p['n_qubits']}  terms={p['n_pauli_terms_incl_identity']}")
    print("worst-case checks:", json.dumps(worst, indent=None))
    print("all routes agree <1e-9:", res["all_routes_agree"], "| variational bound ok:", res["variational_bound_ok"],
          "| global min always N=2:", res["global_min_always_in_N2"], "| points:", res["n_points_checked"])
    print(f"minimum RHF R={res['minimum_RHF_A']} A, minimum FCI R={res['minimum_FCI_A']} A, "
          f"2*E(H)={res['E_two_H_atoms_STO3G']:.10f}, FCI@3.0A={scan[-1]['E_FCI_pyscf']:.10f}, RHF@3.0A={scan[-1]['E_RHF']:.10f}")
    print("time", res["wall_time_s"], "s ->", out)


if __name__ == "__main__":
    main()
