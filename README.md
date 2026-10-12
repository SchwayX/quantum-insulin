# quantum-insulin

A small research project on whether, and how, quantum algorithms could contribute to electronic-structure calculations. It starts from benchmark molecules. Molecules relevant to insulin are a long-term motivation, not a near-term target: the first goal is a validated, reproducible pipeline for small systems.

## Status (12 Oct 2026)

**Audit stage. There is no pipeline code yet.** The repository currently holds an independent audit of earlier results that were produced outside it.

| Historical claim | Audit finding |
|---|---|
| H₂/STO-3G Hamiltonian experiment | No code or numbers to check. Reference values were rebuilt and cross-checked four independent ways |
| Cysteine, RHF/STO-3G: 64 electrons, 46 basis functions, 92 spin orbitals | Correct |
| Cysteine E(RHF) ≈ −709.111974 Eh | **Not reproducible.** It is 1.757 Eh above the lowest RHF/STO-3G minimum the audit found (−710.8691575 Eh). The original geometry is missing |

- Audit report (sections A–G): [`audit/2026-10-12/AUDIT.md`](audit/2026-10-12/AUDIT.md)
- Audit scripts, raw results and how to rerun them: [`audit/2026-10-12/README.md`](audit/2026-10-12/README.md)
- Milestones, blockers and the next step: [`docs/roadmap.md`](docs/roadmap.md)

## Rerunning the audit

You need Python 3.10 or newer, NumPy, SciPy and PySCF. The cysteine geometry optimizations also need geomeTRIC (`pip install geometric`). The only tested combination is Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1 and PySCF 2.14.0; nothing has been tested on Python 3.14 yet.

From the repository root:

```bash
python audit/2026-10-12/scripts/h2_audit.py /tmp/rerun/h2_audit.json   # about 10 s
```

The audit README lists every script. Writing to a folder outside the repository, as above, leaves the committed results untouched for comparison.

## What this project does not claim

- No quantum advantage or speedup. Everything so far is a classical simulation of small, textbook systems.
- Nothing about simulating insulin, its biology, or drug discovery. A single gas-phase cysteine calculation does not establish any of that.
