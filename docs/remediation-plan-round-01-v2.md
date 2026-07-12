# remediation-plan-round-01-v2.md

Produced by Cowork (Claude) from the Codex review on the v2 conformer-search PR
(round 01). Approved by Ish (gate 2) before fixes are implemented.

**Responds to:** Codex round-01 review (findings B-01, M-02).
**Ish judgment calls locked:** M-02 → update the notebook (1a). B-01 stereo →
skip + log molecules with empty/undefined-stereo SMILES (2a).

## Decisions on each finding

| ID | Decision | Rationale | Owner | Verify by |
|----|----------|-----------|-------|-----------|
| B-01 | **Accept** | `IsomericSMILES` is fetched but never carried into resolved rows, so every molecule dead-ends at "missing IsomericSMILES" and no inputs are produced. Pipeline-breaking. | Claude Code | offline schema + end-to-end tests below |
| M-02 | **Accept** | README v2 promises conformer behavior the documented notebook never runs; users would submit v1.1 single-geometry inputs believing the mitigation ran. | Claude Code | local notebook run + offline path smoke test |

No findings rejected.

## Ordered fix tasks

### B-01 — carry IsomericSMILES into resolved rows, skip undefined stereo, test offline

1. **Carry the field.** In `pipeline/pubchem.py` `build_molecule_table()`, add
   `"IsomericSMILES": prop.get("IsomericSMILES", "")` to the resolved-row dict
   (currently `name, pubchem_query, cid, formula, iupac_name, title, status,
   warnings`, ~lines 318-327). Use `IsomericSMILES` — not `CanonicalSMILES`,
   which drops stereochemistry.
2. **Refactor for testability.** Extract resolved-row assembly into a pure helper
   (e.g. `_resolved_row(prop, ...) -> dict`) so the schema can be asserted
   offline without hitting PubChem.
3. **Stereo gate (2a).** In `pipeline/conformers.py`, before embedding, validate
   the SMILES per molecule:
   - empty/missing `IsomericSMILES` → **skip**, log to
     `conformer_search_failed.csv` with reason `"no IsomericSMILES"`.
   - parses but has ≥1 *unspecified* stereo element → **skip**, log with reason
     `"undefined stereochemistry"`. Use RDKit `FindPotentialStereo`
     (or `FindMolChiralCenters(includeUnassigned=True)`), flagging any element
     whose stereo is unspecified. A molecule with **no** stereocenters (adenine,
     water) is NOT skipped — proceed normally.
   - Never let RDKit auto-assign stereo silently; skipping is the safe failure.
4. **Tests** (`tests/test_pubchem.py`, `tests/test_conformers.py`, offline):
   - schema guard: `_resolved_row` output contains `IsomericSMILES`.
   - valid defined-stereo SMILES → conformer XYZ + log rows produced.
   - empty SMILES → skipped + logged `"no IsomericSMILES"`, no XYZ written.
   - undefined-stereo SMILES (e.g. a sugar with unspecified centers) → skipped +
     logged `"undefined stereochemistry"`, no XYZ written.
   - no-stereocenter molecule (adenine) → proceeds, produces conformer(s).
5. **Regression check.** Confirm adding the column doesn't break existing
   consumers of `build_molecule_table` (e.g. `download_sdfs`).

**Verify B-01 by:** `pytest tests/ -q` green including the new cases;
`python scripts/check_invariants.py` green; and an offline end-to-end run of a
2-row table (one defined-stereo, one undefined-stereo) producing conformer XYZ +
Gaussian `.com` for the first and a logged skip for the second. Do not mark
resolved because code changed — run these.

### M-02 — wire the notebook to the conformer path

1. **Rewire the walkthrough** in `notebooks/run_pipeline.ipynb`: replace the
   geometry+gaussian cells (~lines 184, 208, 234) so the default flow is
   `build_molecule_table → search_conformers(...) →
   write_gaussian_coms_from_conformers(...) → write_slurm_scripts`. Remove the
   `download_sdfs → convert_sdfs_to_xyz → write_gaussian_coms` path from the
   default walkthrough (or demote it to a clearly-labeled "v1.1 legacy" appendix
   cell, not the main flow).
2. **Expose config** in the notebook's configuration cell: `n_generate`,
   `top_n`, `rmsd_prune`, `seed`.
3. **Match the narrative** to behavior: describe the top-3 conformer stage and
   the skip-on-undefined-stereo rule, so a skipped molecule isn't a surprise.
4. **README consistency.** With the notebook fixed the v2 claim becomes true;
   still verify README wording matches actual behavior (top 3, skip+log stereo)
   and keeps the three honest limitations (MMFF ranking, gas-phase FF vs
   solution-phase DFT, fixed sampling).
5. **Offline path smoke test.** Add a test that runs the exact code path the
   notebook uses (`search_conformers → write_gaussian_coms_from_conformers`) on a
   hardcoded defined-stereo SMILES, no network, asserting per-conformer `.com`
   files named `{base}_c00/_c01/_c02_F.com` with intact Link1 and ΔE in the
   title. This proves the notebook's code path without executing PubChem in CI.

**Verify M-02 by:** the offline smoke test passing in CI; plus a **manual local
run** of the notebook top-to-bottom on the demo molecules, confirming it calls
the conformer path and emits per-conformer `.com` files (notebook uses live
PubChem, so it stays a local check, consistent with the no-network CI rule).
"It executed" is not sufficient — confirm the outputs are per-conformer files.

## Findings deliberately rejected
- None.

---
### Implementation evidence (fill in as fixes land — do NOT mark resolved on code change alone; run the verify step)
| ID | Commit | Verification run | Result |
|----|--------|------------------|--------|
| B-01 | <hash> | `pytest tests/ -q`; offline 2-row end-to-end | pass/fail |
| M-02 | <hash> | offline smoke test; manual notebook run | pass/fail |
