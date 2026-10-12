"""AUDIT of the historical cysteine RHF/STO-3G claims (not project code).

Historical claim under test: 46 spatial basis functions, 92 spin orbitals, 64 electrons,
E(RHF) ~ -709.111974 Eh, "approximate cysteine geometry" (geometry file NOT available).

What this script establishes, independently of the missing file:
  1. Whether the counts follow from the formula C3H7NO2S in STO-3G.
  2. The RHF/STO-3G energy range of chemically reasonable neutral L-cysteine geometries
     (idealized starting points AND RHF/STO-3G-optimized conformers, SCF stability checked).
  3. What kinds of input error produce an energy shift of the size observed.
  4. Resource numbers for the quantum layer (FCI dimension, active spaces, JW Pauli counts, 1-norm).
Usage (from the repository root; needs PySCF and geomeTRIC):
    python audit/2026-10-12/scripts/cys_audit.py audit/2026-10-12/results/cys_audit.json
Writes, next to the output JSON: <name>_progress.log and cysteine_rhf_sto3g_lowest.xyz.
"""
import json
import math
import os
import platform
import shutil
import sys
import tempfile
import time

import numpy as np
import scipy
import geometric
import pyscf
from pyscf import ao2mo, gto, mcscf, mp, scf
from pyscf.geomopt import geometric_solver

import jw
from cys_geom import l_cysteine, min_nonbonded_distance

E_HIST = -709.111974
BASIS = "sto-3g"
CONV = {"convergence_energy": 1e-7, "convergence_grms": 1e-4, "convergence_gmax": 1.5e-4,
        "convergence_drms": 4e-4, "convergence_dmax": 6e-4}
EH2KCAL, EH2EV = 627.509474, 27.211386
T0 = time.time()
PROGRESS_LOG = None  # set in main(): progress lines are also appended to this file


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    print(line, flush=True)
    if PROGRESS_LOG:
        with open(PROGRESS_LOG, "a") as fh:
            fh.write(line + "\n")


def rhf(atoms, unit="Angstrom", basis=BASIS, charge=0, spin=0):
    mol = gto.M(atom=atoms, basis=basis, unit=unit, charge=charge, spin=spin, verbose=0)
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    mf.chkfile = None
    mf.conv_tol = 1e-10
    mf.max_cycle = 200
    e = float(mf.kernel())
    if not mf.converged:  # second-order fallback, recorded
        mf = mf.newton()
        e = float(mf.kernel())
    return mol, mf, e


def geometry_report(mol):
    xyz = mol.atom_coords(unit="Angstrom")
    sym = [mol.atom_symbol(i) for i in range(mol.natm)]
    D = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    # covalent radii (Angstrom) for a connectivity sanity check
    rcov = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}
    bonds = []
    for i in range(mol.natm):
        for j in range(i + 1, mol.natm):
            if D[i, j] < 1.25 * (rcov[sym[i]] + rcov[sym[j]]):
                bonds.append((f"{sym[i]}{i}", f"{sym[j]}{j}", round(float(D[i, j]), 4)))
    np.fill_diagonal(D, np.inf)
    return {"n_bonds_detected": len(bonds), "bonds": bonds, "min_distance_A": round(float(D.min()), 4)}


def orbital_report(mol, mf):
    nocc = mol.nelectron // 2
    e = mf.mo_energy
    S = mf.get_ovlp()
    C = mf.mo_coeff
    s_idx = [i for i, lab in enumerate(mol.ao_labels()) if lab.split()[1] == "S"]
    homo = C[:, nocc - 1]
    w_s = float(np.dot(homo[s_idx], (S @ homo)[s_idx]))
    return {"n_occupied": nocc, "HOMO_index_0based": nocc - 1,
            "E_HOMO_Eh": float(e[nocc - 1]), "E_LUMO_Eh": float(e[nocc]),
            "gap_Eh": float(e[nocc] - e[nocc - 1]), "gap_eV": float((e[nocc] - e[nocc - 1]) * EH2EV),
            "koopmans_IP_eV": float(-e[nocc - 1] * EH2EV),
            "HOMO_mulliken_weight_on_S": round(w_s, 3),
            "n_virtual": int(mol.nao - nocc)}


def stability(mf):
    try:
        out = mf.stability(internal=True, external=True, return_status=True)
        return {"internal_stable": bool(out[2]), "external_stable(RHF->UHF)": bool(out[3])}
    except Exception as exc:  # recorded, not hidden
        return {"error": repr(exc)}


def counts(mol):
    per_elem = {}
    for i in range(mol.natm):
        el = mol.atom_symbol(i)
        nao_atom = sum((mol.bas_angular(b) * 2 + 1) * mol.bas_nctr(b)
                       for b in range(mol.nbas) if mol.bas_atom(b) == i)
        per_elem.setdefault(el, {"atoms": 0, "functions_per_atom": nao_atom})
        per_elem[el]["atoms"] += 1
    nocc = mol.nelectron // 2
    core = {"H": 0, "C": 1, "N": 1, "O": 1, "S": 5}
    ncore = sum(core[mol.atom_symbol(i)] for i in range(mol.natm))
    return {"formula_atoms": {k: v["atoms"] for k, v in per_elem.items()},
            "electrons": mol.nelectron, "spatial_basis_functions": mol.nao,
            "spin_orbitals": 2 * mol.nao, "basis_functions_per_atom": {k: v["functions_per_atom"] for k, v in per_elem.items()},
            "occupied_spatial": nocc, "virtual_spatial": mol.nao - nocc,
            "core_spatial_orbitals(frozen-core convention)": ncore,
            "jw_qubits_full_space": 2 * mol.nao,
            "fci_determinants_full": str(math.comb(mol.nao, nocc) ** 2),
            "fci_determinants_full_log10": round(2 * math.log10(math.comb(mol.nao, nocc)), 2),
            "fci_determinants_frozen_core": str(math.comb(mol.nao - ncore, nocc - ncore) ** 2),
            "fci_determinants_frozen_core_log10": round(2 * math.log10(math.comb(mol.nao - ncore, nocc - ncore)), 2)}


def koridon_one_norm(h1, eri):
    """Pauli 1-norm (excluding identity) of the spin-summed molecular Hamiltonian
    (Koridon et al., Phys. Rev. Research 3, 033127 (2021)); validated below against explicit JW."""
    n = h1.shape[0]
    t = h1 + np.einsum("pqrr->pq", eri) - 0.5 * np.einsum("prrq->pq", eri)
    lam_t = np.abs(t).sum()
    lam_v = 0.25 * np.abs(eri).sum()
    anti = eri - eri.transpose(0, 3, 2, 1)  # (pq|rs) - (ps|rq)
    p, q, r, s = np.ogrid[:n, :n, :n, :n]
    mask = (p > r) & (s > q)
    lam_v += 0.5 * np.abs(anti[np.broadcast_to(mask, anti.shape)]).sum()
    return float(lam_t + lam_v)


def main(out):
    global PROGRESS_LOG
    out = os.path.abspath(out)
    outdir = os.path.dirname(out)
    PROGRESS_LOG = os.path.splitext(out)[0] + "_progress.log"
    open(PROGRESS_LOG, "w").close()
    res = {"what": "Audit of historical cysteine RHF/STO-3G claims", "historical_E_RHF": E_HIST,
           "software": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                        "pyscf": pyscf.__version__, "platform": platform.platform(),
                        "geomeTRIC": getattr(geometric, "__version__", "unknown")},
           "model": {"molecule": "neutral (non-zwitterionic) L-cysteine HSCH2CH(NH2)COOH", "charge": 0,
                     "spin_multiplicity": 1, "basis": "STO-3G (PySCF built-in)", "method": "RHF",
                     "units": "Angstrom", "scf_conv_tol": 1e-10, "geomopt": "geomeTRIC, tight criteria", "geomopt_criteria": CONV}}

    # ---- 1. counts
    mol0, mf0, _ = rhf(l_cysteine(-60, 180, 160))
    res["counts"] = counts(mol0)
    log(f"counts: {res['counts']['electrons']} e, {res['counts']['spatial_basis_functions']} AO, {res['counts']['spin_orbitals']} spin orbitals")

    # ---- 2. conformers: idealized SP -> optimization -> tight SP
    seeds = [(-60, 180, 160), (60, 180, 160), (180, 180, 160), (-60, 60, -20), (60, -60, -20), (180, 60, -20)]
    confs = []
    best = None
    cwd0 = os.getcwd()
    scratch = tempfile.mkdtemp(prefix="cys_geomopt_")
    os.chdir(scratch)  # keep any optimizer scratch files out of the repository
    for chi1, chi2, psi in seeds:
        atoms = l_cysteine(chi1, chi2, psi)
        mol, mf, e_start = rhf(atoms)
        t = time.time()
        ok, mol_eq = geometric_solver.kernel(mf, maxsteps=300, **CONV)
        t_opt = time.time() - t
        mol_eq, mf_eq, e_opt = rhf(
            [(mol_eq.atom_symbol(i), mol_eq.atom_coord(i, unit="Angstrom").tolist()) for i in range(mol_eq.natm)])
        g = mf_eq.nuc_grad_method().kernel()
        c = {"seed(chi1,chi2,psi)": [chi1, chi2, psi], "E_idealized_start": e_start,
             "start_min_distance_A": round(min_nonbonded_distance(atoms), 3),
             "opt_converged": bool(ok), "opt_time_s": round(t_opt, 1), "E_optimized": e_opt,
             "scf_converged": bool(mf_eq.converged), "max_abs_gradient_Eh_per_bohr": float(np.abs(g).max()),
             "geometry": geometry_report(mol_eq), "orbitals": orbital_report(mol_eq, mf_eq),
             "stability": stability(mf_eq),
             "xyz_angstrom": [[mol_eq.atom_symbol(i)] + [round(float(v), 8) for v in mol_eq.atom_coord(i, unit="Angstrom")]
                              for i in range(mol_eq.natm)]}
        confs.append(c)
        log(f"seed {chi1,chi2,psi}: start {e_start:.8f}  opt {e_opt:.10f}  conv={ok}  bonds={c['geometry']['n_bonds_detected']}  {t_opt:.0f}s")
        if best is None or e_opt < best[2]:
            best = (mol_eq, mf_eq, e_opt, c)
    os.chdir(cwd0)
    shutil.rmtree(scratch, ignore_errors=True)
    res["conformers"] = confs
    mol_b, mf_b, e_b, c_b = best
    res["lowest_conformer_E"] = e_b
    e_opts = [c["E_optimized"] for c in confs]
    e_starts = [c["E_idealized_start"] for c in confs]
    res["E_range_optimized"] = [min(e_opts), max(e_opts)]
    res["E_range_idealized_starts"] = [min(e_starts), max(e_starts)]
    gap = E_HIST - e_b
    res["historical_minus_lowest"] = {"Eh": gap, "kcal_mol": gap * EH2KCAL, "eV": gap * EH2EV}
    with open(os.path.join(outdir, "cysteine_rhf_sto3g_lowest.xyz"), "w") as f:
        f.write(f"{mol_b.natm}\nL-cysteine neutral, RHF/STO-3G optimized (audit), E={e_b:.10f} Eh, PySCF {pyscf.__version__}\n")
        for row in c_b["xyz_angstrom"]:
            f.write(f"{row[0]:2s} {row[1]:14.8f} {row[2]:14.8f} {row[3]:14.8f}\n")

    # ---- 2b. harmonic frequencies of the lowest conformer (is it a true minimum?)
    try:
        t = time.time()
        from pyscf.hessian import thermo
        hess = mf_b.Hessian().kernel()
        freq = thermo.harmonic_analysis(mol_b, hess)["freq_wavenumber"]
        fr = np.asarray(freq)
        n_imag = int(np.sum(np.iscomplex(fr) & (np.abs(fr.imag) > 1e-6)))
        res["lowest_conformer_frequencies"] = {"n_imaginary": n_imag,
                                               "lowest_cm-1": [str(np.round(v, 1)) for v in fr[:4]],
                                               "time_s": round(time.time() - t, 1)}
        log(f"frequencies: n_imag={n_imag} lowest {fr[:3]}")
    except Exception as exc:
        res["lowest_conformer_frequencies"] = {"error": repr(exc)}
        log(f"frequencies failed: {exc!r}")

    # ---- 3. which input errors give a shift of ~1.76 Eh?
    hyp = {}
    xyz = np.array([r[1:] for r in c_b["xyz_angstrom"]])
    sym = [r[0] for r in c_b["xyz_angstrom"]]
    as_atoms = lambda X: [(s, x.tolist()) for s, x in zip(sym, X)]
    _, mfh, e = rhf(as_atoms(xyz), unit="Bohr")
    hyp["angstrom_numbers_read_as_bohr(x0.529)"] = {"E": e, "dE": e - e_b, "converged": bool(mfh.converged)}
    _, mfh, e = rhf(as_atoms(xyz / 0.52917721), unit="Angstrom")
    hyp["bohr_numbers_read_as_angstrom(x1.890)"] = {"E": e, "dE": e - e_b, "converged": bool(mfh.converged)}
    rng = np.random.default_rng(20261012)
    noise = {}
    for sigma in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
        des = []
        for _ in range(5):
            _, mfh, e = rhf(as_atoms(xyz + rng.normal(0, sigma, xyz.shape)))
            des.append(e - e_b)
        noise[f"sigma_{sigma:.2f}A"] = {"dE_mean": float(np.mean(des)), "dE_min": float(np.min(des)), "dE_max": float(np.max(des))}
    hyp["random_cartesian_noise_5_seeds"] = noise
    # zwitterion single point: move carboxyl H onto N (NH3+ / COO-)
    iN, iCA, iHO, iO2 = 0, 3, 13, 12
    X = xyz.copy()
    u = X[iN] - X[iCA]
    u /= np.linalg.norm(u)
    X[iHO] = X[iN] + 1.03 * u
    _, mfh, e = rhf(as_atoms(X))
    hyp["zwitterion_single_point(H moved O->N)"] = {"E": e, "dE": e - e_b, "converged": bool(mfh.converged)}
    # one H atom too close: thiol H at 0.6 A from S
    X = xyz.copy()
    iS, iHG = 8, 9
    v = X[iHG] - X[iS]
    X[iHG] = X[iS] + 0.6 * v / np.linalg.norm(v)
    _, mfh, e = rhf(as_atoms(X))
    hyp["one_H_misplaced(S-H=0.60A)"] = {"E": e, "dE": e - e_b, "converged": bool(mfh.converged)}
    # other bases with 46 functions
    m6, mf6, e6 = rhf(as_atoms(xyz), basis="sto-6g")
    hyp["sto-6g_same_geometry"] = {"nao": m6.nao, "E": e6, "dE": e6 - e_b}
    # triplet at same geometry (would still be 64 electrons)
    _, mft, et = rhf(as_atoms(xyz), spin=2)
    hyp["ROHF_triplet_same_geometry"] = {"E": et, "dE": et - e_b, "converged": bool(mft.converged)}
    res["error_hypotheses_vs_lowest_conformer"] = hyp
    log("hypotheses done")

    # ---- 4. correlation scale, active spaces, qubit-Hamiltonian resources
    t = time.time()
    e_mp2_fc = float(mp.MP2(mf_b).set(frozen=res["counts"]["core_spatial_orbitals(frozen-core convention)"]).kernel()[0])
    e_mp2_ae = float(mp.MP2(mf_b).kernel()[0])
    res["MP2_correlation_full_space"] = {"frozen_core": e_mp2_fc, "all_electron": e_mp2_ae, "time_s": round(time.time() - t, 1)}
    log(f"MP2 corr fc={e_mp2_fc:.6f} ae={e_mp2_ae:.6f}")
    cas = []
    for n in (2, 4, 6, 8, 10, 12):
        t = time.time()
        mc = mcscf.CASCI(mf_b, n, n)
        mc.verbose = 0
        e_cas = float(mc.kernel()[0])
        t_cas = time.time() - t
        h1eff, ecore = mc.get_h1eff()
        h2 = ao2mo.restore(1, mc.get_h2eff(), n)
        row = {"CAS(e,o)": [n, n], "qubits": 2 * n, "determinants": math.comb(n, n // 2) ** 2,
               "E_CASCI": e_cas, "E_corr_captured": e_cas - e_b, "casci_time_s": round(t_cas, 2),
               "one_norm_koridon": koridon_one_norm(h1eff, h2)}
        if n <= 10:
            t = time.time()
            Hp = jw.molecular_qubit_hamiltonian(ecore, h1eff, h2)
            row["jw_pauli_terms_incl_identity"] = len(Hp)
            row["jw_one_norm_explicit"] = float(sum(abs(c) for k, c in Hp.items() if k != (0, 0)))
            row["jw_build_time_s"] = round(time.time() - t, 2)
            if n <= 4:  # cross-check JW ground state (N, Sz sector) against CASCI
                M = jw.dense_from_paulis(Hp, 2 * n)
                N, Sz = jw.number_and_sz(2 * n)
                w, v = np.linalg.eigh(M)
                ne = np.real(np.einsum("ik,ij,jk->k", v.conj(), N, v))
                sz = np.real(np.einsum("ik,ij,jk->k", v.conj(), Sz, v))
                sel = (np.abs(ne - n) < 1e-8) & (np.abs(sz) < 1e-8)
                row["jw_ground_minus_casci"] = float(w[sel].min() - e_cas)
        cas.append(row)
        log(f"CAS({n},{n}): E={e_cas:.8f} corr={e_cas - e_b:.6f} terms={row.get('jw_pauli_terms_incl_identity')} "
            f"lam_explicit={row.get('jw_one_norm_explicit')} lam_formula={row['one_norm_koridon']:.6f}")
    res["active_spaces"] = cas
    # full-space Pauli 1-norm via the (validated) formula
    C = mf_b.mo_coeff
    h1 = C.T @ mf_b.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.full(mf_b._eri, C), C.shape[1])
    res["full_space_one_norm"] = koridon_one_norm(h1, eri)
    res["full_space_nonzero_mo_eri(|v|>1e-10)_fraction"] = float(np.mean(np.abs(eri) > 1e-10))
    res["wall_time_s"] = round(time.time() - T0, 1)
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    log(f"done -> {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cys_audit.json")
