# Audit scripts and results (12 Oct 2026)

The findings are in [`AUDIT.md`](AUDIT.md). This folder holds what produced them: source code in [`scripts/`](scripts/) and generated output in [`results/`](results/). Nothing here came from the original project code, which was not available to the audit.

## Environment used, and how it differs from the reported project stack

| | This audit | Reported project stack |
|---|---|---|
| Python | 3.13.16 | 3.14.4 |
| NumPy / SciPy | 2.5.3 / 1.18.1 | 2.5.3 / 1.18.1 |
| PySCF | 2.14.0, built from git tag `v2.14.0` (c63a953) | 2.14.0 |
| XC library | xcfun (libxc unavailable; no DFT was used) | Unknown |
| geomeTRIC | git `e31d644` (18 May 2026) | Not used |
| h5py | **Stub**: any HDF5 access raises an error; all runs in memory with checkpoint files off | Real h5py |

PyPI was unreachable from the audit environment, so Python 3.14 packages could not be installed. For identical inputs, these differences should not change Hartree–Fock or FCI energies beyond about 10⁻⁹ Eh.

## Rerunning

Run from the repository root. In WSL2, keep the project under `/home/<user>/`, not `/mnt/c/`. To compare against the committed results without overwriting them, write to another folder:

```bash
mkdir -p /tmp/rerun
python audit/2026-10-12/scripts/env_report.py > /tmp/rerun/env_report.json      # records your environment (read-only)
python audit/2026-10-12/scripts/h2_audit.py     /tmp/rerun/h2_audit.json         # about 10 s
python audit/2026-10-12/scripts/scaling.py      /tmp/rerun/scaling.json          # a few seconds
python audit/2026-10-12/scripts/cys_audit.py    /tmp/rerun/cys_audit.json        # about 10 min; needs geomeTRIC
python audit/2026-10-12/scripts/cys_followup.py /tmp/rerun/cys_followup.json     # about 30 s; uses the xyz written by cys_audit.py
python audit/2026-10-12/scripts/key_numbers.py  /tmp/rerun                       # collects every quoted number
python audit/2026-10-12/scripts/cys_check.py    ORIGINAL_GEOMETRY.xyz            # diagnoses the historical -709.111974 Eh run
```

`cys_check.py` reads an XYZ file or plain `Element x y z` lines. Use `--unit Bohr` if the original script used Bohr.

## Files

| File | Purpose |
|---|---|
| `scripts/jw.py` | Minimal Jordan–Wigner mapping (interleaved spin orbitals: qubit 2p = α(p), 2p+1 = β(p)), with two independent dense constructions for cross-checks |
| `scripts/h2_audit.py` → `results/h2_audit.json` | H₂/STO-3G: RHF, PySCF FCI, analytic 2×2 CI and Jordan–Wigner exact diagonalization; Hermiticity, [H, N], [H, S_z] and integral symmetries; scan over 0.30–3.00 Å |
| `scripts/cys_geom.py` | Neutral L-cysteine built from textbook internal coordinates, with a chirality check |
| `scripts/cys_audit.py` → `results/cys_audit.json`, `results/cys_audit_progress.log` | Counts; six conformers optimized at RHF/STO-3G; SCF stability; orbitals; frequencies; deliberate-error tests; MP2; CASCI series; Pauli counts and 1-norms |
| `results/cysteine_rhf_sto3g_lowest.xyz` | Lowest RHF/STO-3G conformer found (no imaginary frequencies) |
| `scripts/cys_followup.py` → `results/cys_followup.json` | Size of the RHF→UHF instability; exact Pauli-term count formula and full-space extrapolation |
| `scripts/scaling.py` → `results/scaling.json` | Electrons, orbitals and qubits from H₂ to insulin (STO-3G, 6-31G*, cc-pVDZ) |
| `scripts/key_numbers.py` → `results/key_numbers.json` | Every number quoted in `AUDIT.md`, collected from the files above |
| `scripts/cys_check.py` | Read-only diagnostic for any cysteine geometry file: composition, bonds, unit errors, energy |
| `scripts/env_report.py` | Read-only record of the Python environment, packages, BLAS, CPU/RAM and WSL details |

## Changes since the first run

- **Bug fix in `cys_check.py`.** The first version compared energies against −710.8689081 Eh, the first conformer optimized, instead of the lowest minimum, −710.8691575 Eh. Its reported differences were therefore 0.25 mEh too small. It now reads the reference from `results/cys_audit.json`. A copy of the scripts handed over as a zip on 12 Oct 2026 contains the old version.
- **Layout.** Scripts moved into `scripts/` and outputs into `results/`. `cys_audit.py` now keeps optimizer scratch files in a temporary folder and records the installed geomeTRIC version instead of a hard-coded commit, and `key_numbers.py` replaces an ad-hoc extraction. All results were then regenerated from this layout; the reproducibility check is described at the end of [`AUDIT.md`](AUDIT.md#audit-environment-and-its-limits).
