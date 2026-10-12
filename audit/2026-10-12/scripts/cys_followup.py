"""AUDIT follow-ups on the lowest RHF/STO-3G cysteine conformer:
 (a) size of the RHF->UHF (triplet) instability reported by PySCF's stability analysis;
 (b) Jordan-Wigner Pauli-term counts vs number of spatial orbitals, polynomial fit validated on held-out
     sizes, then extrapolated to the full space (46 orbitals, 92 qubits) and frozen-core space (35, 70).
Usage (from the repository root, after cys_audit.py):
    python audit/2026-10-12/scripts/cys_followup.py audit/2026-10-12/results/cys_followup.json \
           [audit/2026-10-12/results/cysteine_rhf_sto3g_lowest.xyz]
The geometry defaults to cysteine_rhf_sto3g_lowest.xyz next to the output JSON.
"""
import json
import os
import sys
import time

import numpy as np
from pyscf import ao2mo, gto, lib, scf
from pyscf.scf import stability

import jw


def load_xyz(path):
    rows = [l.split() for l in open(path).read().splitlines()[2:] if l.strip()]
    return [(r[0], [float(v) for v in r[1:4]]) for r in rows]


def main(out, xyz):
    res = {"geometry_file": os.path.basename(xyz)}
    mol = gto.M(atom=load_xyz(xyz), basis="sto-3g", unit="Angstrom", verbose=0)
    mf = scf.RHF(mol)
    mf.chkfile = None
    mf.conv_tol = 1e-11
    e_rhf = float(mf.kernel())

    # (a) external stability: lowest eigenvalue of the RHF->UHF orbital Hessian, then follow it
    mf.verbose = 4
    buf = []
    mf.stdout = type("S", (), {"write": lambda self, s: buf.append(s), "flush": lambda self: None})()
    mo_u_guess, stable_ext = stability.rhf_external(mf, return_status=True)
    mf.verbose = 0
    eig_lines = [l.strip() for l in "".join(buf).splitlines() if "lowest" in l.lower() or "instab" in l.lower()]
    mfu = scf.UHF(mol)
    mfu.chkfile = None
    mfu.conv_tol = 1e-11
    occ = mfu.get_occ(mo_energy=(mf.mo_energy, mf.mo_energy), mo_coeff=(mf.mo_coeff, mf.mo_coeff))
    dm0 = mfu.make_rdm1(mo_u_guess, occ)
    e_uhf = float(mfu.kernel(dm0))
    # second pass of UHF internal stability to make sure we reached a minimum
    mo_i, _, st_i, _ = mfu.stability(return_status=True)
    if not st_i:
        e_uhf = float(mfu.kernel(mfu.make_rdm1(mo_i, mfu.mo_occ)))
    s2, mult = mfu.spin_square()
    res["external_stability"] = {"E_RHF": e_rhf, "E_UHF_broken_symmetry": e_uhf, "E_UHF_minus_RHF_Eh": e_uhf - e_rhf,
                                 "E_UHF_minus_RHF_mEh": (e_uhf - e_rhf) * 1e3, "UHF_<S^2>": float(s2),
                                 "UHF_converged": bool(mfu.converged), "pyscf_log": eig_lines}
    print(json.dumps(res["external_stability"], indent=1), flush=True)

    # (b) Pauli-term counts on sub-blocks of the canonical MO integrals (generic C1 integrals)
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.full(mf._eri, C), C.shape[1])
    nocc = mol.nelectron // 2
    series = {}
    for n in range(2, 13):
        lo = max(0, nocc - n // 2)
        idx = slice(lo, lo + n)
        t = time.time()
        Hp = jw.molecular_qubit_hamiltonian(0.0, h1[idx, idx], eri[idx, idx, idx, idx])
        series[n] = len(Hp)  # includes identity
        print(f"norb={n:2d} qubits={2 * n:2d} pauli_terms={len(Hp):6d} ({time.time() - t:.1f}s)", flush=True)
    res["pauli_terms_by_norb"] = series
    fit_n = list(range(2, 10))
    coef = np.polyfit(fit_n, [series[n] for n in fit_n], 4)
    pred = {n: float(np.polyval(coef, n)) for n in series}
    res["poly4_fit_on_norb_2_to_9"] = {"coefficients_high_to_low": coef.tolist(),
                                       "heldout_errors": {n: pred[n] - series[n] for n in (10, 11, 12)},
                                       "max_abs_error_all": max(abs(pred[n] - series[n]) for n in series)}
    res["extrapolated_pauli_terms"] = {"norb46_qubits92": float(np.polyval(coef, 46)),
                                       "norb35_qubits70_frozen_core": float(np.polyval(coef, 35))}
    print(json.dumps(res["poly4_fit_on_norb_2_to_9"], indent=None))
    print(json.dumps(res["extrapolated_pauli_terms"]))
    with open(out, "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    out_path = sys.argv[1] if len(sys.argv) > 1 else "cys_followup.json"
    xyz_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(os.path.abspath(out_path)), "cysteine_rhf_sto3g_lowest.xyz")
    main(out_path, xyz_path)
