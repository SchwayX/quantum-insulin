# Audit of historical results (12 October 2026)

This is a frozen record of the project's first audit, done before any implementation. The current status is kept in [`docs/roadmap.md`](../../docs/roadmap.md).

**How the audit was done.** The repository was empty, and the original files on the owner's WSL2 machine were not accessible. The audit therefore tests the historical *claims*, not the original code. Every number below was computed independently by the scripts in [`scripts/`](scripts/). The raw outputs are in [`results/`](results/), and [`results/key_numbers.json`](results/key_numbers.json) collects every quoted value. [`README.md`](README.md) explains how to rerun everything.

**Conventions.**
- Energies are in hartree (Eh): 1 Eh = 627.5095 kcal/mol = 27.2114 eV.
- Coordinates are in ångström.
- "Exact" means exact diagonalization within the stated finite basis. It never means the exact molecular energy.
- Qubit counts assume the Jordan–Wigner mapping: one qubit per spin orbital.

## A. Repository overview

- On 12 Oct 2026, `SchwayX/quantum-insulin` had no commits, branches or tags. `git ls-remote` returned no refs, and GitHub showed "This repository is empty".
- There was no source code, notebook, configuration, dependency file or test to inspect.
- The earlier work is on the owner's WSL2 machine. The audit could not reach it, and no copy existed in earlier conversations.

## B. Scientific and technical status at audit time

| Layer | Previously claimed | In repository | Status |
|---|---|---|---|
| A. Classical chemistry | Cysteine RHF/STO-3G | Nothing | Counts verified; energy not reproducible (section D, item 1) |
| B. Hamiltonians | "An initial H₂ Hamiltonian experiment" (no numbers given) | Nothing | Not assessable; reference rebuilt independently (section C) |
| C. Quantum algorithms | — | Nothing | Planned only |
| D. Reporting | — | Nothing | No README, environment specification, tests, result files or documentation |

The previous workflow (geometry → PySCF → RHF → printed numbers) kept no record of units, charge, spin multiplicity, convergence settings or software versions. That missing record, more than any single number, is why the earlier results cannot be relied on yet.

## C. Verified facts versus unverified claims

### Software environment

| Claim | Verdict | How it was checked |
|---|---|---|
| Linux under WSL2 | Not verified | The owner's machine was not accessible |
| Python 3.14.4 | The version exists; its installation was not verified. The newest 3.14 patch release on the audit date was 3.14.8 | CPython git tags |
| NumPy 2.5.3, SciPy 1.18.1 | Exist; the audit used exactly these | Installed in the audit environment |
| PySCF 2.14.0 | Exists, and is the newest release (18 Jul 2026); the audit used it | Git tag `v2.14.0`; PyPI release page |

### Historical calculations

| Claim | Verdict |
|---|---|
| H₂ Hamiltonian experiment | Not verifiable: no code or numbers |
| Approximate cysteine geometry | Not verifiable: no file, and no record of where it came from |
| 64 electrons | Correct for neutral C₃H₇NO₂S: 3×6 + 7×1 + 7 + 2×8 + 16 = 64 |
| 46 spatial basis functions | Correct for STO-3G: 3 C × 5 + 1 N × 5 + 2 O × 5 + 1 S × 9 + 7 H × 1 = 46. Not unique to STO-3G: STO-6G also gives 46 |
| 92 spin orbitals | Correct (2 × 46) |
| E(RHF) ≈ −709.111974 Eh | **Not reproducible** (see below) |

### Cysteine RHF/STO-3G energies

Model: neutral (non-zwitterionic) L-cysteine, HS–CH₂–CH(NH₂)–COOH, charge 0, singlet, RHF/STO-3G in PySCF 2.14.0, SCF convergence threshold 10⁻¹⁰ Eh. Six starting structures were built from textbook bond lengths and angles, varying the dihedrals χ1 (N–Cα–Cβ–S), χ2 (Cα–Cβ–S–H) and ψ′ (N–Cα–C=O). Each was then optimized with geomeTRIC using tight convergence criteria.

| χ1, χ2, ψ′ (°) | Unoptimized start (Eh) | Optimized (Eh) |
|---|---|---|
| −60, 180, 160 | −710.858637 | −710.8689081 |
| 60, 180, 160 | −710.857950 | −710.8669973 |
| 180, 180, 160 | −710.857598 | −710.8682875 |
| −60, 60, −20 | −710.858500 | −710.8684161 |
| 60, −60, −20 | −710.857978 | −710.8679093 |
| 180, 60, −20 | −710.859766 | **−710.8691575** |

- All six optimizations converged, with all 13 expected bonds intact and residual gradients of at most 1.4×10⁻⁵ Eh/bohr.
- The lowest structure has no imaginary vibrational frequencies, so it is a true minimum. Its coordinates are in [`results/cysteine_rhf_sto3g_lowest.xyz`](results/cysteine_rhf_sto3g_lowest.xyz).
- The six minima lie within 2.2 mEh of one another. Even the unoptimized starting structures are within 12 mEh of the lowest minimum.
- **The historical −709.111974 Eh is 1.757 Eh (1,103 kcal/mol; 47.8 eV) above the lowest minimum.** That is about 800 times the spread between conformers. The historical calculation cannot have been of a reasonable cysteine structure at this level of theory.

### H₂/STO-3G reference

Model: H₂, STO-3G, charge 0, singlet, with H atoms at (0, 0, 0) and (0, 0, R). The qubit Hamiltonian uses the Jordan–Wigner mapping with interleaved spin orbitals (qubit 2p = α(p), qubit 2p+1 = β(p); |1⟩ = occupied). Its identity term includes the nuclear repulsion.

| R (Å) | E_nuc (Eh) | E_RHF (Eh) | E_FCI (Eh) |
|---|---|---|---|
| 0.7414 | 0.7137539937 | −1.1166843871 | −1.1372701747 |
| 0.735 | 0.7199689944 | −1.1169989968 | −1.1373060358 |

Checks at 218 geometries (0.30–3.00 Å, plus a fine grid near the minimum):
- PySCF FCI, an analytic 2×2 CI, and the ground state of the qubit Hamiltonian in the N = 2, S_z = 0 sector agree to within 1.8×10⁻¹⁵ Eh.
- The qubit matrix built from Pauli strings matches one built directly from fermionic ladder-operator matrices to within 8.1×10⁻¹⁵ (largest element difference).
- The Hamiltonian is Hermitian and commutes exactly with the particle number N and with S_z. All Pauli coefficients are real.
- RHF ≥ FCI at every point, as the variational principle requires.
- The lowest eigenvalue over the whole 16-dimensional Fock space always lies in the N = 2 sector. That holds for H₂ but is not guaranteed for other molecules.
- On a 0.0005 Å grid, the minima are at R = 0.712 Å (RHF) and 0.735 Å (FCI).
- At 3.0 Å, RHF gives −0.6560 Eh while FCI gives −0.9336 Eh, close to the energy of two STO-3G hydrogen atoms (−0.9332 Eh). This is the known failure of RHF at dissociation.

At R = 0.7414 Å the mapped Hamiltonian has 4 qubits and 15 Pauli terms (identity included), with a Pauli 1-norm of 1.885 Eh (identity excluded):

| Pauli term(s) | Coefficient (Eh) |
|---|---|
| I | −0.098863969 |
| Z0, Z1 | +0.171197749 each |
| Z2, Z3 | −0.222785930 each |
| Z0 Z1 | +0.168622192 |
| Z2 Z3 | +0.174348442 |
| Z0 Z2, Z1 Z3 | +0.120544822 each |
| Z0 Z3, Z1 Z2 | +0.165867024 each |
| X0 Y1 Y2 X3, Y0 X1 X2 Y3 | +0.045322202 each |
| X0 X1 Y2 Y3, Y0 Y1 X2 X3 | −0.045322202 each |

## D. Critical errors and risks

**1. The historical cysteine energy does not correspond to cysteine at any reasonable geometry.** To see which kinds of input error produce a shift of this size, the audit introduced deliberate errors into the lowest structure and recomputed the energy:

| Deliberate error | ΔE vs lowest minimum (Eh) |
|---|---|
| Random Cartesian noise, σ = 0.05 Å (5 trials) | +0.05 to +0.10 |
| Random Cartesian noise, σ = 0.10 Å (5 trials) | +0.22 to +0.74 |
| Random Cartesian noise, σ = 0.20 Å (5 trials) | +0.57 to +1.99 |
| Random Cartesian noise, σ = 0.30 Å (5 trials) | +1.38 to +10.77 |
| Thiol H moved to 0.60 Å from S | +1.46 |
| Bohr coordinates read as Å (molecule expanded ×1.89) | +3.02 |
| Å coordinates read as Bohr (molecule shrunk ×0.53) | +17.15 |
| Zwitterion-like proton placement (single point) | +0.36 |
| ROHF triplet at the same geometry | +0.12 |
| STO-6G instead of STO-3G (also 46 functions) | −6.07 |

The +1.757 Eh gap fits a geometry with misplaced atoms or large coordinate errors. A whole-file units mix-up would give about +3 Eh or +17 Eh instead, and neither a different spin state nor a different 46-function basis fits. Only the original geometry file can confirm the cause; [`scripts/cys_check.py`](scripts/cys_check.py) diagnoses it. Anything derived from the historical run (orbital energies, the HOMO–LUMO gap, active-space choices) should be discarded.

**2. No record of inputs.** Results without their inputs, units, settings and software versions cannot be reproduced. From now on every result file should carry them; this audit's JSON files show the minimum.

**3. Full-size cysteine is not a realistic quantum target.**
- The full problem needs 92 qubits. Its configuration space has 5.75×10²² determinants (5.4×10¹⁸ with the 11 core orbitals frozen).
- The qubit Hamiltonian has 6,624,139 Pauli terms. For generic real integrals the count follows 1.5n⁴ − n³ + 2.5n² + 1 for n spatial orbitals. The formula was fitted on n = 2–9 and reproduced n = 10–12 exactly.
- The Pauli 1-norm is λ = 2,779 Eh. The standard sampling bound (λ/ε)² gives about 3×10¹² measurements for a single energy to ε = 1.6 mEh (about 1 kcal/mol), before measurement grouping or variance reduction.
- Active spaces are unavoidable, but small ones recover little correlation energy:

| Active space | Qubits | Determinants | Correlation recovered (Eh) | Pauli terms |
|---|---|---|---|---|
| CAS(2,2) | 4 | 4 | −0.0000005 | 27 |
| CAS(4,4) | 8 | 36 | −0.0006 | 361 |
| CAS(6,6) | 12 | 400 | −0.0021 | 1,819 |
| CAS(8,8) | 16 | 4,900 | −0.0227 | 5,793 |
| CAS(10,10) | 20 | 63,504 | −0.0307 | 14,251 |
| CAS(12,12) | 24 | 853,776 | −0.0342 | 29,737 (from the count formula) |

  For scale, frozen-core MP2 estimates the total correlation energy in this basis at −0.330 Eh. These are CASCI values with canonical RHF orbitals centred on the HOMO–LUMO gap. Better orbital choices (for example MP2 natural orbitals or CASSCF) should be tested before an active space is chosen.

**4. Interpretation traps.**
- **HOMO–LUMO gap.** The RHF/STO-3G gap is 15.5–15.7 eV across the six conformers. It is a difference of orbital energies, not an optical or fundamental gap, and minimal-basis virtual orbitals have little physical meaning. The HOMO is the sulfur lone pair (86.5–94.5% of its Mulliken population on S). Koopmans' theorem turns its energy into an ionization-energy estimate of 7.4–7.5 eV.
- **SCF stability.** Every conformer's RHF solution is internally stable, but each has an RHF → UHF instability. For the lowest conformer the broken-symmetry UHF solution is 1.14 mEh lower, with ⟨S²⟩ = 0.157. The spin-symmetric RHF determinant remains the appropriate reference for a closed-shell singlet; the instability should be recorded and rechecked in larger basis sets.
- **Chemistry.** In water, cysteine exists as a zwitterion, and all six cysteines in insulin are disulfide-bonded (A6–A11, A7–B7, A20–B19). A single gas-phase cysteine molecule says nothing about insulin chemistry or drug discovery.

**5. Traps for the quantum layers.**
- Exact diagonalization and VQE must be restricted to the correct particle-number and spin sector. This held for H₂ at every geometry, but nothing guarantees it for other molecules.
- Spin-orbital ordering changes the Pauli strings but not the energies, so it must be documented with every Hamiltonian.
- VQE on H₂ tests the plumbing. It is not evidence of capability or advantage.

**6. Dependencies.** Python 3.14 is a low risk for the classical stack, because PySCF 2.14.0 ships `py3-none` wheels that don't depend on the Python minor version. For the quantum layer, PyPI metadata on the audit date showed the following (not install-tested):

| Package | Latest version | Python 3.14 support |
|---|---|---|
| Qiskit | 2.5.2 | `cp310-abi3` wheel, usable on 3.14 |
| qiskit-aer | 0.17.2 | `cp314` wheel (added 4 Feb 2026) |
| PennyLane | 0.45.1 | Classifiers include 3.14 |
| ffsim | 0.0.84 | `cp311-abi3` wheel |
| cirq-core | 1.7.0 | Classifiers include 3.14 |
| OpenFermion | 1.8.1 | Classifiers stop at 3.13 |
| Qiskit Nature | 0.8.0 | Classifiers stop at 3.12 |
| openfermionpyscf | 0.5 | Last release December 2020 |

Under WSL2, the Linux VM gets 50% of the Windows RAM by default, and Microsoft recommends keeping project files in the Linux file system (for example `/home/<user>/`) rather than under `/mnt/c` for performance.

## E. Prioritized implementation plan

| # | Improvement | Why it matters | Effort |
|---|---|---|---|
| 1 | Commit the original files unchanged, plus the output of `scripts/env_report.py` | No historical result can be verified without its inputs | About 15 min (owner) |
| 2 | Run `scripts/cys_check.py` on the original geometry; mark −709.111974 Eh as superseded; adopt a documented geometry | The only reported number is wrong | Under 1 h |
| 3 | Milestone 1 skeleton: project layout, pinned environment, a results format that records every input setting and version, automated tests | Stops items 1–2 from recurring | About half a day |
| 4 | H₂ benchmark built on the validation prototyped here, with a bond-length scan and plot | Foundation for VQE | About half a day |
| 5 | A one-hour Python 3.14 trial, then one quantum SDK (or plain NumPy up to about 20 qubits); avoid openfermionpyscf | Avoids overlapping SDKs and dependency churn | 1–2 h |

## F. Compute requirements and likely limitations

Timings measured in the audit environment (2 CPU cores, 8 GB RAM):

| Calculation | Time |
|---|---|
| Whole H₂ audit (218 geometries) | 5–7 s |
| One cysteine RHF/STO-3G energy | Under 1 s |
| One cysteine geometry optimization | 45–80 s |
| Cysteine harmonic frequencies | About 30 s |
| Cysteine CASCI(12,12) | About 20 s |

Simulating qubits directly takes memory that doubles with every qubit. An ideal state vector needs 16 bytes × 2ⁿ: 28 qubits need 4 GiB and 30 need 16 GiB. A density matrix for noisy simulation needs 16 bytes × 4ⁿ: 14 qubits need 4 GiB. Under WSL2's default memory cap, about 28 qubits (ideal) and 12–14 qubits (noisy) are practical ceilings. Milestones 1–4 need no quantum hardware or credentials.

| System | Electrons | STO-3G orbitals (qubits) | cc-pVDZ orbitals (qubits) |
|---|---|---|---|
| H₂ | 2 | 2 (4) | 10 (20) |
| LiH | 4 | 6 (12) | 19 (38) |
| H₂O | 10 | 7 (14) | 24 (48) |
| Cysteine | 64 | 46 (92) | 137 (274) |
| Cystine (disulfide) | 126 | 90 (180) | 264 (528) |
| Human insulin, C₂₅₇H₃₈₃N₆₅O₇₇S₆ | 3,092 | 2,432 (4,864) | 7,609 (15,218) |

Storing just the cc-pVDZ two-electron integrals for insulin would take about 3.35 PB (n⁴/8 unique values × 8 bytes).

## G. First concrete task: capture the historical baseline

1. **Owner (about 15 minutes):** copy the original files, unchanged, onto a `historical-baseline` branch, and add two read-only outputs. Run this inside WSL2:

   ```bash
   git clone https://github.com/SchwayX/quantum-insulin.git ~/quantum-insulin
   cd ~/quantum-insulin
   git switch -c historical-baseline
   mkdir historical
   cp -r /path/to/original/files/* historical/
   python audit/2026-10-12/scripts/env_report.py > historical/env_report.json
   python audit/2026-10-12/scripts/cys_check.py historical/GEOMETRY_FILE > historical/cys_check.json
   git add historical
   git status        # the repository is public: check for keys, tokens or personal files first
   git commit -m "Add original files, unchanged"
   git push -u origin historical-baseline
   ```

   `cys_check.py` reads an XYZ file or plain `Element x y z` lines. If the geometry is inside a Python script, paste its atom lines into a text file first, and add `--unit Bohr` if the script used Bohr. If the original files no longer exist, −709.111974 Eh is unrecoverable and the project should start from the documented geometry.
2. **Next session:** explain the historical number from that file. Then open one small pull request for Milestone 1: README, pinned environment, results format, and a test that checks the H₂ values above to 10⁻⁸ Eh.

## Audit environment and its limits

- Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1 and PySCF 2.14.0, built from the `v2.14.0` git tag (commit c63a953) with xcfun (libxc unavailable; no DFT was used). geomeTRIC from git (commit e31d644, 18 May 2026). Linux x86_64, 2 cores, 8 GB RAM.
- PyPI was unreachable from the audit environment, so Python 3.14 could not be installed and nothing was install-tested on it.
- h5py was replaced by a stub that raises an error on any HDF5 access. Every calculation ran in memory with checkpoint files disabled, so no result depends on it.
- The cysteine structures were built from textbook internal coordinates and then optimized; none were downloaded from a database.
- The results were generated twice: once during the audit, then again from this committed layout. The H₂ output was identical. Every cysteine energy agreed to within 7×10⁻¹² Eh, and the deliberate-error single points to within 4×10⁻¹⁰ Eh, so no quoted number changed.

## Sources

- Repository state: [SchwayX/quantum-insulin](https://github.com/SchwayX/quantum-insulin) (12 Oct 2026)
- PySCF 2.14.0 distribution files: [pypi.org/project/pyscf/2.14.0](https://pypi.org/project/pyscf/2.14.0/#files)
- PyPI project pages: [qiskit](https://pypi.org/project/qiskit/#files), [qiskit-aer](https://pypi.org/project/qiskit-aer/#files), [qiskit-nature](https://pypi.org/project/qiskit-nature/), [openfermion](https://pypi.org/project/openfermion/), [openfermionpyscf](https://pypi.org/project/openfermionpyscf/), [pennylane](https://pypi.org/project/pennylane/), [ffsim](https://pypi.org/project/ffsim/#files), [cirq-core](https://pypi.org/project/cirq-core/)
- Human insulin formula: [Selleckchem datasheet](https://www.selleckchem.com/datasheet/insulin-human-S695502-DataSheet.html); [CymitQuimica, CAS 11061-68-0](https://cymitquimica.com/cas/11061-68-0/)
- WSL2: [.wslconfig defaults](https://learn.microsoft.com/en-us/windows/wsl/wsl-config); [file-system performance](https://learn.microsoft.com/en-us/windows/wsl/filesystems)
- Pauli 1-norm formula: E. Koridon et al., *Phys. Rev. Research* 3, 033127 (2021). The implementation in `scripts/cys_audit.py` matched explicit Jordan–Wigner sums to 10⁻¹² for CAS(2,2)–CAS(10,10).
- Geometry optimizer: L.-P. Wang and C. Song, *J. Chem. Phys.* 144, 214108 (2016) (geomeTRIC).
