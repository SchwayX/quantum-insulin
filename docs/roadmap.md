# Roadmap and status

_Last updated 12 Oct 2026, after the first audit. Update this file at the end of every work session, so that a new session can resume without repeating the audit._

## Current state

- The repository holds only the [12 Oct 2026 audit](../audit/2026-10-12/AUDIT.md). There is no pipeline code yet.
- **Verified:** the H₂/STO-3G reference energies and qubit Hamiltonian; the cysteine electron, orbital and spin-orbital counts.
- **Not reproducible:** the historical cysteine energy, −709.111974 Eh. It is 1.757 Eh above the lowest RHF/STO-3G minimum the audit found, and the original geometry is missing.

## Milestones

| # | Milestone | Status | Notes |
|---|---|---|---|
| 0 | Audit before implementation | Done, 12 Oct 2026 | [`audit/2026-10-12/`](../audit/2026-10-12/) |
| 1 | Reproducibility and project structure | Not started | Waiting for the historical baseline (next task) |
| 2 | H₂ benchmark | Not started | Reference values and validation checks prototyped in the audit |
| 3 | VQE benchmark | Not started | Needs Milestone 2 and a single quantum SDK |
| 4 | Cysteine calculation audit | Partly done | Counts verified and reference minima found; the historical number stays unexplained until the original geometry is available |
| 5 | Insulin-related feasibility report | Not started | Size and resource numbers in audit sections D and F |

## Blockers

1. The original scripts and cysteine geometry are not in the repository.
2. No quantum SDK has been chosen, and nothing has been install-tested on Python 3.14.

## Next smallest task

**Owner:** push the original files, unchanged, to a `historical-baseline` branch, together with the outputs of `env_report.py` and `cys_check.py`. The commands are in [audit section G](../audit/2026-10-12/AUDIT.md#g-first-concrete-task-capture-the-historical-baseline).

## Decisions for the owner

- **Publishing the originals.** The repository is public. Check the files for keys, tokens and personal information before pushing.
- **License.** There is no license file, so by default others may read the code but not reuse it. Choose one (for example MIT or Apache-2.0) if reuse is intended.
- **Quantum SDK.** Pick one after the Python 3.14 trial (audit section E, item 5).
- **Working location.** Keep the WSL2 copy under `/home/<user>/`, not `/mnt/c/`.

## Session log

| Date | Work | Files | Calculations and tests run | Outcome |
|---|---|---|---|---|
| 2026-10-12 | Audit, then committed the audit record | `audit/2026-10-12/`, `README.md`, `docs/roadmap.md`, `.gitignore` | H₂ audit (218 geometries); cysteine audit (6 optimizations, frequencies, deliberate-error tests, CASCI up to (12,12)); Pauli-count fit; scaling counts. Full regeneration from the committed layout | H₂ verified; cysteine counts verified; historical cysteine energy not reproducible |
