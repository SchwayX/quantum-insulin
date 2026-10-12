"""AUDIT: size scaling of electron, orbital and qubit counts from molecular formulas (PySCF basis definitions).
Counts only - no electronic-structure calculation is implied. Spherical (5d) functions, PySCF built-in basis sets."""
import json
import math
import sys

from pyscf import gto

Z = {"H": 1, "Li": 3, "C": 6, "N": 7, "O": 8, "S": 16}
CORE = {"H": 0, "Li": 1, "C": 1, "N": 1, "O": 1, "S": 5}
SYSTEMS = {
    "H2": {"H": 2},
    "LiH": {"Li": 1, "H": 1},
    "H2O": {"H": 2, "O": 1},
    "cysteine C3H7NO2S": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 1},
    "cystine C6H12N2O4S2 (disulfide)": {"C": 6, "H": 12, "N": 2, "O": 4, "S": 2},
    "human insulin C257H383N65O77S6": {"C": 257, "H": 383, "N": 65, "O": 77, "S": 6},
}
BASES = ["sto-3g", "6-31g*", "cc-pvdz"]


def nao_per_element(basis):
    out = {}
    for el in Z:
        spin = Z[el] % 2
        out[el] = gto.M(atom=f"{el} 0 0 0", basis=basis, spin=spin, verbose=0).nao
    return out


def main(out):
    table = {b: nao_per_element(b) for b in BASES}
    rows = []
    for name, f in SYSTEMS.items():
        ne = sum(Z[e] * n for e, n in f.items())
        row = {"system": name, "atoms": sum(f.values()), "electrons": ne,
               "frozen_core_orbitals": sum(CORE[e] * n for e, n in f.items())}
        for b in BASES:
            nao = sum(table[b][e] * n for e, n in f.items())
            row[f"{b}_spatial_orbitals"] = nao
            row[f"{b}_jw_qubits"] = 2 * nao
            row[f"{b}_eri_storage_8fold_GB"] = round(nao ** 4 / 8 * 8 / 1e9, 6)
        nocc = ne // 2
        n = row["sto-3g_spatial_orbitals"]
        row["sto-3g_fci_log10_determinants"] = round(2 * math.log10(math.comb(n, nocc)), 1) if ne % 2 == 0 else None
        rows.append(row)
    res = {"basis_functions_per_element": table, "systems": rows}
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(table))
    for r in rows:
        print(f"{r['system']:34s} atoms={r['atoms']:4d} e={r['electrons']:5d} | "
              + " | ".join(f"{b}: {r[b + '_spatial_orbitals']:5d} orb {r[b + '_jw_qubits']:6d} q {r[b + '_eri_storage_8fold_GB']:.3g} GB" for b in BASES)
              + f" | FCI(STO-3G) 10^{r['sto-3g_fci_log10_determinants']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "scaling.json")
