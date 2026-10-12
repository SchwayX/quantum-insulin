"""Diagnose a cysteine geometry and its RHF energy (audit tool; read-only).

    python audit/2026-10-12/scripts/cys_check.py original_geometry.xyz [--unit Angstrom|Bohr]
           [--basis sto-3g] [--charge 0] [--spin 0] [--reference audit/2026-10-12/results/cys_audit.json]

Accepts a standard XYZ file (natoms / comment / "El x y z" lines) or a plain list of "El x y z" lines,
e.g. the atom block pasted out of a PySCF script. Use --unit Bohr if the original script used unit='Bohr'.
Reports composition, electron/orbital counts, interatomic-distance and bond-length outliers, a units
heuristic, and the RHF energy compared with the historical value and the audit reference.
The reference (lowest RHF/STO-3G conformer found by the audit) is read from results/cys_audit.json;
if that file is missing, the documented audit value is used (lowest optimized conformer
-710.8691575333 Eh; optimized conformers span -710.86700 to -710.86916 Eh).
Historical claim under test: -709.111974 Eh.
"""
import argparse
import json
import os
from collections import Counter

import numpy as np

HIST = -709.111974
DEFAULT_REFERENCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "results", "cys_audit.json")
FALLBACK_LOWEST = -710.8691575333  # audit 2026-10-12; used only if the results file is unavailable
RCOV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}
TYPICAL = {  # approximate typical ranges (Angstrom) for neutral amino acids
    ("C", "H"): (1.06, 1.12), ("H", "N"): (0.99, 1.04), ("H", "O"): (0.95, 1.00), ("H", "S"): (1.31, 1.38),
    ("C", "C"): (1.47, 1.56), ("C", "N"): (1.43, 1.50), ("C", "O"): (1.18, 1.45), ("C", "S"): (1.78, 1.87)}


def read_geometry(path):
    rows = [l.split() for l in open(path) if l.strip()]
    if rows and len(rows[0]) == 1 and rows[0][0].isdigit():
        rows = rows[2:]  # XYZ header
    atoms = [(r[0].capitalize(), [float(v) for v in r[1:4]]) for r in rows if len(r) >= 4]
    return atoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geometry")
    ap.add_argument("--unit", default="Angstrom")
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=0, help="2S (number of unpaired electrons)")
    ap.add_argument("--reference", default=DEFAULT_REFERENCE, help="cys_audit.json holding the audit reference energy")
    a = ap.parse_args()
    try:
        with open(a.reference) as fh:
            audit_lowest = float(json.load(fh)["lowest_conformer_E"])
        reference_source = os.path.normpath(a.reference)
    except (OSError, KeyError, ValueError):
        audit_lowest = FALLBACK_LOWEST
        reference_source = "built-in fallback (audit 2026-10-12)"
    atoms = read_geometry(a.geometry)
    sym = [s for s, _ in atoms]
    X = np.array([x for _, x in atoms]) * (0.52917721 if a.unit.lower().startswith("b") else 1.0)
    rep = {"composition": dict(Counter(sym)), "expected_C3H7NO2S": dict(Counter("C C C H H H H H H H N O O S".split()))}
    rep["composition_ok"] = rep["composition"] == rep["expected_C3H7NO2S"]
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    np.fill_diagonal(D, np.inf)
    rep["min_distance_A"] = round(float(D.min()), 4)
    close, bonds, outliers = [], [], []
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            d = D[i, j]
            if d < 0.8:
                close.append((f"{sym[i]}{i}", f"{sym[j]}{j}", round(float(d), 3)))
            if d < 1.25 * (RCOV[sym[i]] + RCOV[sym[j]]):
                bonds.append((f"{sym[i]}{i}", f"{sym[j]}{j}", round(float(d), 3)))
                lo, hi = TYPICAL.get(tuple(sorted((sym[i], sym[j]))), (0, 9))
                if not lo <= d <= hi:
                    outliers.append((f"{sym[i]}{i}-{sym[j]}{j}", round(float(d), 3), (lo, hi)))
    rep["pairs_closer_than_0.8A"] = close
    rep["bonds_detected"] = len(bonds)
    rep["bonds_expected_for_cysteine"] = 13
    rep["bond_length_outliers"] = outliers
    heavy = [j for j, s in enumerate(sym) if s != "H"]
    h_nn = [float(D[i, heavy].min()) for i, s in enumerate(sym) if s == "H"]  # each H to nearest heavy atom
    rep["median_H_to_nearest_heavy_atom_A"] = round(float(np.median(h_nn)), 3) if h_nn else None
    if h_nn and np.median(h_nn) < 0.8:
        rep["units_hint"] = "X-H distances ~0.55 A: coordinates were probably Angstrom but read as Bohr"
    elif h_nn and np.median(h_nn) > 1.6:
        rep["units_hint"] = "X-H distances ~2 A: coordinates were probably Bohr but read as Angstrom"
    else:
        rep["units_hint"] = "no obvious unit error"
    from pyscf import gto, scf
    mol = gto.M(atom=[(s, x.tolist()) for s, x in zip(sym, X)], unit="Angstrom", basis=a.basis,
                charge=a.charge, spin=a.spin, verbose=0)
    mf = (scf.RHF if a.spin == 0 else scf.ROHF)(mol)
    mf.chkfile = None  # read-only tool: do not write checkpoint files
    mf.conv_tol = 1e-10
    e = float(mf.kernel())
    nocc = mol.nelectron // 2
    rep.update({"electrons": mol.nelectron, "spatial_basis_functions": mol.nao, "spin_orbitals": 2 * mol.nao,
                "E_RHF": e, "scf_converged": bool(mf.converged),
                "E_minus_historical": e - HIST, "audit_lowest_conformer_E": audit_lowest,
                "audit_reference_source": reference_source,
                "E_minus_audit_lowest_conformer": e - audit_lowest})
    if a.spin == 0:
        rep["HOMO_LUMO_gap_eV"] = float((mf.mo_energy[nocc] - mf.mo_energy[nocc - 1]) * 27.211386)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
