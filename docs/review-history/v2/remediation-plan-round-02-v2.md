# remediation-plan-round-02-v2.md

Produced by Cowork (Claude) from the Codex round-02 review. Approved by Ish
(gate 2) before fixes are implemented.

**Responds to:** Codex round-02 review (findings M-03, M-04).
**Context:** Live PubChem has deprecated `IsomericSMILES`/`CanonicalSMILES`;
stereochemistry is now returned under `SMILES` (verified: `D-ribose` →
`C1[C@H]([C@H](...`), flat connectivity under `ConnectivitySMILES`. The dead keys
return nothing, so any literal read of them yields `""`.
**Ish judgment calls locked:** M-04 → (1a) exclude still-unconverged conformers
from ranking; (2b) if all conformers fail, carry one flagged best-effort seed to
DFT. M-03 dead-key sweep → in scope.

## Decisions on each finding

| ID | Decision | Rationale | Owner | Verify by |
|----|----------|-----------|-------|-----------|
| M-03 | **Accept (+ sweep)** | `score_candidate` reads the dead `IsomericSMILES` key, so ambiguous candidates never get the stereo bonus and the resolver can pick the wrong CID. Root cause (dead keys) affects any other literal reader too. | Claude Code | scoring test + dead-key grep clean |
| M-04 | **Accept** | Unconverged FF optimizations are recorded as optimized and can win a top-3 slot on a meaningless energy, handing DFT a non-minimized seed while logs claim success. Violates provenance + "ran ≠ validated". | Claude Code | convergence tests below |

No findings rejected.

## Ordered fix tasks

### M-03 — route all SMILES reads through the helper; retire the dead keys

1. **Fix the cited line.** In `pipeline/pubchem.py`, `score_candidate()` must read
   the SMILES via `_isomeric_smiles(prop)`, not `prop.get("IsomericSMILES", "")`.
2. **Sweep (in scope).** Replace every remaining literal read of `IsomericSMILES`
   or `CanonicalSMILES` with the helper. This includes the **property-fetch list**
   (originally requesting `IsomericSMILES,CanonicalSMILES`): change it to request
   `SMILES,ConnectivitySMILES`, the current keys.
3. **Confirm helper semantics.** `_isomeric_smiles(prop)` should prefer `SMILES`
   (stereo-bearing) and fall back to the legacy key only if present, returning
   `""` when neither exists (the 2a stereo gate already skips+logs empty).
4. **Tests** (`tests/test_pubchem.py`, offline fixtures):
   - a current-schema record (stereo under `SMILES`) scores with the stereo bonus;
     a lower-CID achiral/unspecified candidate does not win the tie.
   - `_isomeric_smiles` returns the `SMILES` value on current records.

**Verify M-03 by:** `pytest tests/ -q` green with the new scoring test; a repo
grep shows **no** literal `IsomericSMILES`/`CanonicalSMILES` reads outside the
helper's fallback; `python scripts/check_invariants.py` green. Do not mark
resolved on code change alone — run these.

### M-04 — reject/flag unconverged force-field conformers

1. **Capture convergence.** Read the per-conformer `not_converged` flag returned
   by MMFF/UFF optimization. Add a `converged` (bool) column to `conformer_log.csv`
   for every conformer — the log must never imply optimization succeeded when it
   did not (provenance invariant).
2. **Retry once.** For conformers flagged not-converged on the first pass,
   re-run the FF optimization with increased max iterations before judging them
   failed.
3. **(1a) Exclude from ranking.** After the retry, `select_top_n` ranks and
   selects **only converged** conformers. A still-unconverged conformer with an
   unreliable energy can never win a top-3 slot.
4. **(2b) All-fail fallback.** If **no** conformer for a molecule converges after
   retry: carry exactly **one** best-effort geometry (lowest FF energy among the
   unconverged set) into Gaussian — not three. Mark its log row `converged=False`,
   emit a warning, and write an explicit convergence marker into the `.com` title
   line (e.g. a `UNCONVERGED_FF_SEED` tag) so the unminimized start is visible on
   inspection. The recorded energy stays labeled unreliable.
5. **Tests** (`tests/test_conformers.py`, mostly pure/synthetic):
   - filtering: given mixed converged/unconverged conformers, only converged ones
     reach `select_top_n`.
   - retry: a conformer unconverged first-pass but converged after retry is
     included.
   - log: an unconverged conformer is recorded `converged=False`.
   - all-fail branch (2b): all-unconverged input → exactly one row carried,
     `converged=False`, warning emitted, and the `.com` title carries the marker.

**Verify M-04 by:** `pytest tests/ -q` green including the four cases above;
an offline run over a normal molecule (top-3 converged, unchanged behavior) and a
mocked all-unconverged molecule (one flagged seed, marker in the `.com` title).
Confirm the log's `converged` column reflects reality. Not "it ran" — check the
outputs.

## Regression checks
- Normal (all-converged) molecules still yield the same top-3 as round 01;
  the new `converged` column is additive.
- The M-03 fetch-key change doesn't break `download_sdfs` or the stereo gate
  (stereo still arrives via `SMILES`).

## Findings deliberately rejected
- None.

---
### Implementation evidence (fill in as fixes land — do NOT mark resolved on code change alone; run the verify step)
| ID | Commit | Verification run | Result |
|----|--------|------------------|--------|
| M-03 | `5f635bd` | `pytest tests/ -q` → 68 passed (incl. `TestScoreCandidateCurrentSchema`); dead-key grep — only remaining `.get("IsomericSMILES")` reads are the helper's own fallback (`pubchem.py:127`) and the documented molecule-table-column read (`conformers.py`, our DataFrame, not a PubChem dict); no `CanonicalSMILES` reads; `check_invariants.py` → passed | **pass** |
| M-04 | `016c197` | `pytest tests/ -q` → 79 passed (incl. the four convergence cases); `check_invariants.py` → passed; offline run — ribose → top-3 all `converged=True`, unchanged vs round-01, no marker; mocked all-unconverged → exactly one `converged=False` seed, warning emitted, `UNCONVERGED_FF_SEED` in the `.com` title, Link1 intact | **pass** |

**M-03 note.** The property-fetch list already requested the current
`SMILES,ConnectivitySMILES` keys (fixed in round-01 B-01); confirmed unchanged.
The only code-level dead-key *read* was `score_candidate`, now routed through
`_isomeric_smiles`. The pipeline's own molecule-table column keeps the name
`IsomericSMILES` (established schema, `MOLECULE_TABLE_COLUMNS`) and holds real
stereo SMILES sourced from the live `SMILES` key — reading it is not a dead-key
read; a clarifying comment marks the distinction.

**M-04 note (deviation).** `generate_conformers` now returns a 4-tuple
(`coords, energies_kcal, method, converged`) vs the round-01 3-tuple — required to
carry per-conformer convergence into the log. Recorded in
`docs/implementation-status-v2.md` §3.
