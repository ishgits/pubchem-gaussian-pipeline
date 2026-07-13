# implementation-status.md

> Canonical status doc — the AGENTS.md §5 merge gate, so it must always reflect
> the real current status; `scripts/check_invariants.py` fails the objective
> floor if this file ever drifts back to the empty template. The full per-round
> working history (commit hashes, verify runs) is archived under
> `docs/review-history/v2/` (`implementation-status-v2.md` and each
> `remediation-plan-round-0N-v2.md`).

Maintained by the implementing agents. Reviewers must verify these claims
against the code, not trust them.

**PR:** #3   **Branch:** `feat/conformer-search-v2`   **Round:** 8
(v2 conformer search + remediation rounds 01–08; released as v2.0.0). Works
against `docs/architecture.md` and `docs/implementation-plan.md`.

## 1. What was implemented

- **Conformer search stage (`pipeline/conformers.py`, new).** RDKit ETKDGv3
  `EmbedMultipleConfs` (RMSD-pruned) then MMFF94 optimize/rank (UFF logged
  fallback when MMFF params are missing); keep the top `TOP_N=3` lowest-energy
  distinct conformers per molecule. `select_top_n` is a pure, RDKit-free ranking
  helper. `search_conformers` is the batch driver: it reads the molecule table's
  `IsomericSMILES`, writes `{base}_c{ii}.xyz` per kept conformer, and appends one
  provenance row per conformer to `conformer_log.csv`. Energies are kcal/mol;
  `rel_energy_kcalmol` is ΔE from the molecule's lowest carried conformer.
- **Stereo/validity gate (round 01, B-01).** `check_conformer_eligibility` runs
  before embedding: empty or unparseable SMILES, or a molecule with undefined
  stereochemistry, is skipped and logged to `conformer_search_failed.csv` rather
  than letting RDKit assign an arbitrary stereoisomer. No-stereocenter molecules
  proceed normally.
- **PubChem SMILES sourcing and scoring (round 01/02/06, B-01/M-03/M-13).** The
  resolved molecule table carries the stereo-bearing SMILES into an
  `IsomericSMILES` column via the `_isomeric_smiles` helper, and
  `score_candidate` reads stereo through the same helper. Both use PubChem's
  current `SMILES` key (see §3 deviation). The stereo bonus recognizes
  tetrahedral `@` and both directional-bond E/Z tokens (`/` and `\`), while
  `ConnectivitySMILES` remains excluded from scoring.
- **FF convergence handling (round 02, M-04/M-05).** `generate_conformers`
  captures each conformer's `not_converged` flag and retries only the failed
  conformers once with more iterations (`_optimize_single_conf`), so
  already-converged conformers keep the first-pass energy that matches their
  written geometry. `select_converged_top_n` ranks only converged conformers;
  if none converge, exactly one lowest-energy best-effort seed is carried with
  `converged=False`, a warning, and an `UNCONVERGED_FF_SEED` marker in the `.com`
  title. `conformer_log.csv` records a `converged` column.
- **Version + commit provenance (round 03, M-06).** `pipeline.__version__` is
  `"2.0.0"` (round-04 MOD-01); `pipeline_provenance()` returns that version plus a best-effort git
  short SHA (`.dirty` suffix on an uncommitted tree, empty string when git is
  absent). `search_conformers` records `pipeline_version` and `pipeline_commit`
  on every `conformer_log.csv` row and appends RDKit, pipeline-version, and
  pipeline-commit tokens to each per-conformer XYZ comment line.
- **Gaussian writer (extended).** `write_gaussian_com` gained optional
  `conformer_id`, `rel_energy_kcalmol`, and `unconverged` parameters:
  per-conformer basenames `{base}_c{ii}_F.com` / `.chk`, the ΔE (kcal/mol) and any
  `UNCONVERGED_FF_SEED` marker in the title. The Link1 opt→freq checkpoint
  contract is untouched. `write_gaussian_coms_from_conformers` writes one `.com`
  per `conformer_log.csv` row.
- **Self-contained conformer COM provenance (round-07 M-14).** Every v2
  conformer-derived `.com` copies `pipeline_version`, `pipeline_commit`, and
  `rdkit_version` from its conformer-log row into a separate `provenance` line in
  the Gaussian title section. A missing commit is explicit as
  `commit=unavailable`; the same fields are retained in `com_write_log.csv`,
  including its header-only zero-job schema. Legacy v1.1 output has no provenance
  line unless those optional arguments are explicitly supplied. Route lines,
  coordinates, charge/multiplicity, checkpoint names, and Link1 behavior are
  unchanged.
- **Required conformer provenance at the Gaussian boundary (round-08 M-16).**
  `write_gaussian_coms_from_conformers` rejects every nonempty external log that
  lacks a nonblank `pipeline_version` or `rdkit_version`, before deleting a stale
  failure log or writing any COM/log output. These source versions describe
  conformer generation and are never inferred from the downstream writer's
  environment. Empty zero-job logs remain valid; an absent/blank source commit
  remains explicit in each COM as `commit=unavailable`.
- **Required provenance for direct conformer COM writes (round-08 M-17).** A
  direct `write_gaussian_com()` call with `conformer_id` now requires nonblank
  `pipeline_version` and `rdkit_version` before output-directory creation or file
  writing. This closes the direct-writer bypass around the protected batch path.
  Legacy calls with `conformer_id=None` retain optional provenance and unchanged
  v1.1 naming/text behavior; a missing conformer source commit remains explicit
  as `commit=unavailable`.
- **Notebook (round 01, M-02).** `notebooks/run_pipeline.ipynb` runs the
  conformer path by default; the v1.1 single-geometry path is a commented-out
  labeled legacy appendix.
- **Canonical status doc + drift guard (round 03, M-07).** This file is populated
  and synced from `implementation-status-v2.md`; `scripts/check_invariants.py`
  gained a guard that fails when this file still contains template markers.
- **Identity/config-validated resume (M-09 + round-04 B-02).** `search_conformers`
  skips a molecule only when all its existing rows match this run's run-level
  config (`seed`, `n_generate`, `top_n`, `rmsd_prune`, `pipeline_version`,
  **`rdkit_version`**) *and*, when available, a matching clean
  **`pipeline_commit`**, plus the requested per-molecule identity (**`cid`**,
  **`smiles`**), and every recorded **`xyz_path` exists and is non-empty**. A
  `.dirty` commit is never reusable because the marker cannot identify the
  working-tree content; if either commit is unavailable, the documented
  pipeline-version fallback applies. Rows from a different config/identity, a
  stale RDKit build, a missing geometry, or a pre-provenance log are dropped and
  regenerated with a warning. By default (round-04 M-02) the log holds exactly
  the molecules requested this run; `append=True` retains unrequested molecules'
  rows only after the same validation.
- **Complete-group resume validation (round-05 M-12).** A requested molecule is
  resumable only when every row has a parseable, agreeing `n_kept`; the group has
  exactly that many rows; `conformer_id` values are the unique contiguous set
  `0..n_kept-1`; and `xyz_path` values are unique, present, and non-empty. A
  truncated, duplicated, or otherwise damaged group is dropped and regenerated
  rather than accepted because its surviving rows happen to be valid.
- **Collision-safe output labels (round-05 B-05).** Before any output directory,
  failure log, or geometry is mutated, `search_conformers` rejects distinct
  molecule labels that collapse to the same `sanitize_basename()` value (including
  punctuation/whitespace and case-only collisions). This prevents one molecule's
  XYZ/COM/checkpoint/SLURM files from overwriting another's.
- **Append-union output identity (round-06 B-07).** `append=True` loads the
  existing conformer log before any output mutation and validates one union of
  current labels plus all unrequested labels that would be retained. The pure
  `validate_unique_output_basenames()` helper rejects a new/current label that
  collides with retained chemistry and also rejects an already-corrupt append
  log. On failure, the prior conformer log, XYZ files, and failure log remain
  unchanged.
- **Append carry-forward integrity (round-08 M-15).** Every unrequested group
  retained by `append=True` must pass the same complete-group/XYZ and current
  config/version checks used for resume, plus internal name/CID/SMILES
  consistency and explicit `pipeline_commit` field presence. Invalid retained
  groups abort append mode before any output mutation: they are never copied,
  silently dropped, partially retained, or regenerated without a current
  molecule-table identity.
- **Physical-line XYZ parsing (round-04 B-01).** `xyz_to_gaussian_coords` reads
  line 1 = atom count, line 2 = comment (may be empty), then exactly that many
  coordinate rows; a count mismatch or malformed row raises `ValueError` instead
  of silently dropping an atom — an empty comment line no longer truncates the
  geometry sent to Gaussian.
- **Run-scoped, submission-independent SLURM (round-04 B-03, M-01, M-03).**
  `write_slurm_scripts` defaults to the current run's `com_write_log.csv` (stale
  `.com` files are never picked up; a `com_dir` glob is the explicit legacy mode)
  and overwrites `.sh` files (reporting `WROTE`/`OVERWROTE`) so a rerun cannot
  leave stale SBATCH directives. Round-05 B-06 additionally removes `.sh` files
  not represented by the current COM set, including all old scripts on a zero-job
  rerun, so the documented submission glob agrees with `slurm_write_log.csv`.
  Each script resolves its `.com` relative to its own location and runs `g16` on
  the basename, so `sbatch` works from any directory.
- **Validated log-driven SLURM inputs (round-08 M-18).** Before creating the
  SLURM directory, pruning stale scripts, or rewriting `slurm_write_log.csv`, the
  default log-driven path now requires every `com_path` to be nonblank and point
  to an existing file. A damaged/stale COM log fails with row-level details
  rather than reporting `WROTE` for a job that cannot start. The explicit legacy
  `com_dir` glob remains unchanged.
- **Complete XYZ software provenance (round-08 M-19).** Generated conformer XYZ
  comments now record the exact RDKit version alongside the force-field method,
  seed, pipeline version, and pipeline commit, making each geometry artifact
  self-describing without its CSV sidecar.
- **Commit-aware conformer reuse (round-08 M-20).** Resume and append reuse now
  require matching clean nonblank pipeline commits when both are available.
  Any `.dirty` commit on the retained row or current run forces regeneration;
  missing commit provenance retains the documented pipeline-version fallback.
- **Stable zero-result schemas (round-05 M-11).** Both the legacy and conformer
  Gaussian batch writers emit header-only, readable COM logs when zero files are
  written (including all-write-failure batches). The SLURM writer likewise emits
  a fixed-schema empty log and completes with zero scripts, so a scientifically
  valid "no eligible jobs" run proceeds through the default notebook stages.
- **Early parameter validation + stale-log hygiene (round-04 MIN-03, MIN-02).**
  `search_conformers` rejects `n_generate<1`, `top_n<1`, `rmsd_prune<0`, duplicate
  molecule labels, colliding sanitized basenames, and empty sanitized filenames
  at entry; it and the Gaussian writers clear a stale v2 conformer or COM
  `*_failed.csv` at stage start.
- **Repo hygiene + release versioning (round-04 B-04, MOD-01).** Generated
  outputs are untracked and gitignored (`git ls-files -ci --exclude-standard`
  empty), enforced by a new `review-readiness.yml` step; `pipeline.__version__ =
  "2.0.0"` and the PubChem User-Agent is `gaussian-input-pipeline/2.0`.
- **Env / CI / README.** `rdkit` in `environment.yml` and the review-readiness CI
  install step; README leads with the RDKit conformer flow, with Open Babel demoted
  to a labeled legacy v1.1 section.

Required checks locally green after round 08 (RDKit 2025.03.3):
`pytest tests/ -q` → **235 passed**; `python scripts/check_invariants.py` →
**passed**; `git diff --check` → **passed**;
`git ls-files -ci --exclude-standard` → empty.

### Round-08 append and provenance boundary hardening

- **M-15 — Resolved:** append-mode carry-forward groups now undergo
  complete-group, XYZ-existence, current-config/version, internal-identity, and
  provenance checks. Invalid retained groups abort append mode before any output
  mutation because unrequested molecules cannot be regenerated safely from the
  current input.
- **M-16 — Resolved:** the conformer-to-Gaussian batch writer rejects nonempty
  logs with missing/blank pipeline or RDKit provenance before creating any COM.
  Missing source commits remain explicit as `commit=unavailable`.
- **M-17 — Resolved:** direct `write_gaussian_com()` calls using
  `conformer_id` now require nonblank pipeline and RDKit source versions before
  output mutation. This closes the direct-writer bypass around the protected
  conformer-log batch path. Legacy calls with `conformer_id=None` remain
  unchanged; missing commits remain explicit as `commit=unavailable`.
- **M-18 — Resolved:** log-driven SLURM generation now rejects blank or missing
  `com_path` entries before output-directory creation, stale-script pruning, or
  SLURM-log rewrites. This prevents fresh scripts and `WROTE` records that point
  to absent Gaussian inputs; the explicit legacy glob path remains unchanged.
- **M-19 — Resolved:** conformer XYZ comments now include the source RDKit
  version together with the existing pipeline version/commit, force-field
  method, and seed provenance.
- **M-20 — Resolved:** retained conformers produced by a different clean commit
  are regenerated, and `.dirty` commit provenance is never reusable. When one
  side has no commit, exact pipeline-version and RDKit-version matching remains
  the conservative available fallback.

## 2. What was NOT implemented (and why)

- Out-of-scope-for-v2 items: rotatable-bond gating, xTB/CREST, energy-window
  logic, Boltzmann/entropy weighting, solvent-aware search, and any change to the
  level of theory or the Link1 contract. Deliberately omitted per architecture-v2.
- Version/commit provenance in the legacy `sdf_download_log.csv` and unrelated
  legacy v1.1 outputs remains deferred. Round-07 M-14 resolves the default v2
  conformer path by stamping provenance into each COM and `com_write_log.csv`.
- Post-optimization RMSD re-pruning is not added; distinctness relies on
  `pruneRmsThresh` at embed time (see §6).
- **Remaining recorded deferrals:** the per-study `runs/` run-directory redesign;
  replacing label-based filenames with a collision-proof internal identifier
  rather than rejecting collisions; and extending version/commit provenance to
  `sdf_download_log.csv` and unrelated legacy outputs. The round-04
  `n_rows == n_kept` reconciliation is no longer deferred: round-05 M-12
  implements it with ID/path integrity checks; v2 COM provenance is likewise no
  longer deferred after round-07 M-14, and round-08 M-20 makes commit provenance
  part of resume validation.
- **Physical cleanup of stale XYZ/COM files is not implemented.** Reduced reruns
  make the current CSV logs and pruned `slurm_scripts/` authoritative, so old
  files in `conformer_xyz/` or `gaussian_inputs/` cannot enter the documented
  submission path. They may remain on disk; README documents this behavior and
  recommends fresh output directories when studies require physical separation.
- **Old or damaged append inputs require explicit repair/regeneration.** Round-08
  M-15 deliberately refuses to infer or backfill retained-group identity and
  provenance. Include the affected molecule in a current run so it can be
  regenerated, repair the log/XYZ set, or use `append=False`.

## 3. Deviations from architecture / plan

- **PubChem SMILES property rename (B-01/M-03).** PubChem renamed its SMILES
  properties in 2025: the stereo-bearing SMILES now arrives under the `SMILES`
  key and the flat one under `ConnectivitySMILES` (the legacy `IsomericSMILES` /
  `CanonicalSMILES` keys return nothing). Verified against live PubChem
  (L-alanine, D-glucose). The code reads stereo via `_isomeric_smiles` (prefer
  `SMILES`, fall back to legacy `IsomericSMILES`, never the flat key). The
  pipeline's own molecule-table column keeps the name `IsomericSMILES` but holds
  real stereo SMILES. No scientific assumption changed — stereochemistry is
  preserved, not dropped.
- **`generate_conformers` returns a 4-tuple** `(coords_list, energies_kcal,
  method, converged)` rather than the 2-tuple in the plan. `method` records the
  force field actually used (MMFF94 or UFF); `converged` records per-conformer FF
  convergence. Both are provenance needed to honor the "record software
  versions" and "ran is not validated" invariants. API-shape only.
- **`conformer_log.csv` column set is a superset of the plan.** Columns:
  `name, cid, smiles, conformer_id, rel_energy_kcalmol, xyz_path, rdkit_version,
  pipeline_version, pipeline_commit, seed, n_generate, top_n, method, n_generated,
  n_kept, rmsd_prune, converged`. `converged` (M-04), `pipeline_version` /
  `pipeline_commit` (M-06), and the requested `n_generate` / `top_n` (M-09) are
  additive; existing columns and the top-3 selection for a given config are
  unchanged. The XYZ comment format is likewise extended with `rdkit=`, `pver=`,
  and `pcommit=` tokens (M-06/M-19). `n_generate` and `top_n` record the
  *requested* search knobs (distinct from the result columns `n_generated` /
  `n_kept`) so a resumed run can validate that stale rows match the current
  config (M-09); clean commit identity is also checked when available, while
  dirty commits are never reused (M-20).
- **All-fail best-effort seed (M-04 decision 2b).** When no conformer converges,
  one real (not fabricated) FF geometry is still handed to DFT, explicitly flagged
  `converged=False`, warned at runtime, and tagged `UNCONVERGED_FF_SEED`. An
  intentional, visible exception to "no placeholder science"; the FF energy is
  labeled unreliable and DFT is expected to refine the geometry.
- **XYZ output directory** defaults to `conformer_xyz/` (overridable), mirroring
  the existing `pubchem_xyz/` convention. Not a scientific change.

No scientific invariant (AGENTS.md §2) was altered: route lines, units,
charge/multiplicity handling, and the Link1 contract are unchanged. Conformer
energies are labeled kcal/mol at every surface (CSV column, XYZ comment, `.com`
title) and never mixed with DFT Hartree values.

## 4. Tests added

- `tests/test_conformers.py` — ranking (`TestSelectTopN`); seeded embedding and
  determinism (`TestGenerateConformers`); batch behavior incl. adenine collapse,
  resume, and the stereo skip cases (`TestSearchConformers`, `TestCheckEligibility`);
  reproducibility (`TestReproducibility`); the offline notebook code path
  (`TestNotebookPathOffline`); convergence selection and retry-merge
  (`TestSelectConvergedTopN`, `TestFinalizeConvergence`); retry alignment so
  recorded energy matches written geometry (`TestRetryAlignment`); batch
  convergence incl. the flagged best-effort seed (`TestConvergenceBatch`);
  version/commit provenance in the log and XYZ header (`TestProvenanceLogging`);
  and identity/config-validated resume — matching config resumes; a changed
  seed/`top_n`/`cid`/`smiles`/RDKit-version, a deleted XYZ, or a pre-provenance log
  regenerates stale rows with a warning; unrequested molecules are dropped by
  default and retained with `append=True`
  (`TestRowConfigMatches`, `TestRowIdentityAndXyz`, `TestResumePartition`,
  `TestResumeConfigValidationBatch`, `TestPreserveUnrequestedBatch`);
  early parameter validation (`TestParameterValidation`) and stale-`*_failed.csv`
  clearing (`TestStaleFailedCsvCleared`).
- `tests/test_conformers.py` round-05 — sanitized-basename collisions fail before
  writes while distinct basenames succeed (B-05); complete resume groups are
  accepted, while truncation, missing/duplicate IDs, inconsistent `n_kept`, and
  duplicate XYZ paths invalidate the group; a truncated three-conformer batch is
  regenerated to the full set (M-12).
- `tests/test_conformers.py` round-06 — append-mode punctuation/whitespace and
  case-only collisions raise before generation; prior log, XYZ, and failure-log
  bytes remain unchanged; distinct basenames append successfully; and an already
  corrupt retained-label set is rejected before mutation (B-07).
- `tests/test_conformers.py` round-08 — pure retained-group checks cover complete
  groups, config/version equality, internal CID/SMILES consistency, and explicit
  commit-field presence. Batch tests corrupt each resume config/version field,
  truncate a group, remove an XYZ, remove provenance columns, and mix CID/SMILES
  identities; every case asserts the conformer log, surviving XYZ files, failure
  log, and full output file set remain byte-for-byte unchanged (M-15). A valid
  retained group still appends successfully. Follow-up coverage verifies RDKit
  in every generated XYZ comment, different clean commits regenerating, dirty
  commits never resuming or carrying forward, and missing commits retaining the
  version-based fallback (M-19/M-20).
- `tests/test_gaussian.py` — physical-line XYZ parsing: empty comment keeps all
  atoms, count-mismatch (either direction) and malformed/non-integer rows raise,
  trailing blank tolerated (`TestXyzParsingByPhysicalLine`).
- `tests/test_slurm.py` — script resolves its `.com` from a sibling directory
  (`TestSlurmScriptResolvesInput`); log-driven default vs legacy glob
  (`TestWriteSlurmScriptsLogDriven`); overwrite-on-rerun
  (`TestWriteSlurmScriptsOverwrite`); smaller reruns prune prior scripts and a
  zero-job rerun removes all `.sh` files while keeping the fixed log schema
  (`TestWriteSlurmScriptsCurrentRunCleanup`).
- `tests/test_slurm.py` round-08 follow-up — blank, whitespace-only, NaN/empty,
  and nonexistent logged COM paths fail before directory/log creation. A mixed
  valid/missing log preserves prior scripts and `slurm_write_log.csv`
  byte-for-byte and writes no partial new job (M-18).
- `tests/test_pubchem.py` — SMILES key handling (`TestIsomericSmiles`), resolved-row
  schema (`TestResolvedRow`), and current-schema scoring with the stereo bonus
  (`TestScoreCandidateCurrentSchema`), including `/` and backslash-only E/Z
  markers, selection over a lower-CID stereo-free candidate, and explicit
  exclusion of `ConnectivitySMILES` from scoring (M-13).
- `tests/test_gaussian.py` — per-conformer filenames, ΔE in title, Link1 intact
  (`TestWriteGaussianComConformer`, `TestWriteGaussianComsFromConformers`);
  header-only legacy/conformer COM logs, all-write-failure behavior, and the
  all-ineligible conformer → Gaussian → SLURM zero-job path
  (`TestEmptyComLogSchemas`). Round-07 tests verify direct and batch provenance,
  explicit unavailable commits, expanded COM-log schemas, unchanged route/
  charge/checkpoint/Link1 fields, and unchanged legacy output (M-14).
- `tests/test_gaussian.py` round-08 — nonempty logs missing either/both required
  source-version columns or containing blank/NaN values fail with row details
  before any output mutation; missing commits still emit
  `commit=unavailable`; empty logs and current valid logs still succeed (M-16).
- `tests/test_gaussian.py` round-08 follow-up — direct conformer-specific calls
  missing either/both pipeline and RDKit versions, or containing blank,
  whitespace-only, or NaN versions, fail before creating the output directory.
  Valid direct calls retain filename, ΔE, checkpoint, route, charge/multiplicity,
  and Link1 behavior while adding the required provenance line; legacy calls
  without `conformer_id` remain unchanged (M-17).
- `tests/test_utils.py` — the offline provenance helper: git absent / non-zero /
  timeout give empty string, clean tree gives the SHA, dirty tree appends
  `.dirty` (`TestGitShortSha`, `TestPipelineProvenance`). No test asserts a
  concrete SHA.
- `tests/test_check_invariants.py` — the status-doc drift guard fires on a
  template and passes on a populated file (`TestStatusDocDriftGuard`); AST-based
  guards verify M-14 provenance threading/title tokens, M-16 required-version
  batch validation before mutation, M-17 conditional direct-writer validation
  before mutation, M-15 complete/config/identity/provenance carry-forward
  validation plus fail-before-mutation ordering, M-19 RDKit tokens in generated
  XYZ comments, and M-20 commit-aware/dirty-rejecting resume configuration.

RDKit-dependent tests use `pytest.importorskip("rdkit")` so a bare environment
still runs the pure tests; CI installs rdkit so they execute there.

## 5. Known limitations

- MMFF ranking is unreliable for intramolecular-H-bonding species (sugars,
  nucleosides); carrying top 3 mitigates this, and an xTB rerank is future work.
- FF ranking is gas-phase while the DFT default is IEFPCM water, so the FF-lowest
  conformer may not be the solution-phase minimum.
- Fixed `N_GENERATE=20` may under-sample very flexible molecules.
- These are FF starting geometries, not optimized minima — DFT makes the final
  call among the carried candidates.
- Reduced reruns may leave prior pipeline-generated XYZ/COM files on disk. They
  are not referenced by the current logs or pruned SLURM directory and therefore
  are not submitted by the documented workflow; use fresh output directories per
  study if those artifacts must be physically isolated.
- **Reproducibility of `pipeline_commit`.** It pins code identity only when the
  tree is clean. A `.dirty` suffix means uncommitted edits produced the output,
  so that output is not fully reproducible from the commit alone; such rows are
  never reused by resume or append. An empty `pipeline_commit` (no git) falls
  back to `pipeline_version`, which is only as precise as manual version bumping.
  A recorded commit is therefore not by itself a guarantee of exact code.
- **Nonempty pre-round-08 conformer logs without source-version provenance cannot
  be converted to v2 COM files.** They must be regenerated or repaired from
  trustworthy records; the Gaussian stage intentionally does not infer historical
  RDKit/pipeline versions from its current environment.

## 6. Questions requiring scientific judgment  ← Ish reads this FIRST

1. **PubChem `SMILES` as the stereochemistry source of record.** The current
   `SMILES` property is the stereo-bearing one we consume; `ConnectivitySMILES`
   is stereo-free and deliberately unused (verified for L-alanine and D-glucose).
   Confirm this is acceptable. If PubChem ever emits a `SMILES` without stereo for
   a molecule that has stereocenters, the eligibility gate correctly skips it as
   "undefined stereochemistry" rather than embedding an arbitrary isomer.
2. **Distinctness relies on `pruneRmsThresh=0.5 Å` at embed time.** After MMFF
   optimization two kept conformers could in principle relax toward each other;
   v2 does not re-prune post-optimization. Acceptable, or add a post-optimization
   RMSD check (new scope)?
3. **UFF-fallback energies vs MMFF94 energies are on different scales.** Both are
   labeled kcal/mol and `method` is recorded per row, but ΔE from a UFF molecule
   is not comparable to MMFF94 ΔE from another molecule. Confirm per-molecule ΔE
   (never cross-molecule) is the only intended comparison.
4. **Charge/multiplicity for conformers.** The conformer stage builds neutral
   molecules from `IsomericSMILES` (RDKit default protonation). Non-default
   charge/multiplicity must still be passed to the Gaussian writer separately, as
   in v1.1. Confirm this is expected.

## Provenance

- pipeline version: `2.0.0` (`pipeline.__version__`; recorded per row in
  `conformer_log.csv` `pipeline_version`, and best-effort git commit in
  `pipeline_commit`). PubChem User-Agent `gaussian-input-pipeline/2.0`. Branch
  `feat/conformer-search-v2`; base v1.1 (Zenodo 10.5281/zenodo.18894724).
- RDKit version used for the round-08 local test run: 2025.03.3; round-07 used
  2025.09.3. Runtime conformer generation records its actual version per row in
  `conformer_log.csv` `rdkit_version`.
- Conformer stage config: `N_GENERATE=20`, `TOP_N=3`, `RMSD_PRUNE=0.5 Å`,
  `SEED=42`, ranking MMFF94 (UFF logged fallback), energies kcal/mol.
- Open Babel: not used on the conformer path (RDKit consumes `IsomericSMILES`
  directly). The v1.1 SDF-to-XYZ path still uses Open Babel where invoked.
- Gaussian route lines, charge, multiplicity, nproc: unchanged from v1.1; still
  passed in from user config, not constructed in `pipeline/`.
