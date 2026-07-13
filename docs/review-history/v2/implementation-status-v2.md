# implementation-status-v2.md

Maintained by the implementing agent (Claude Code). The reviewer (Codex) must
verify these claims against the code, not trust them.

**PR:** #3   **Branch:** `feat/conformer-search-v2`   **Round:** 1 (+ round-01, 02 & 03 remediation)
Works against `docs/architecture-v2.md` and `docs/implementation-plan-v2.md`.

This is the **working** status file. It is synced into the canonical
`docs/implementation-status.md` (the AGENTS.md §5 merge gate) each remediation
round; `scripts/check_invariants.py` fails if the canonical file drifts back to
the template.

## 0a. Round-04 remediation (PR #3 v2-release review — 12 + 3 findings)

Implements `docs/remediation-plan-round-04-v2.md` in three commits. **Resolved
only after the listed verify step passed**, not on code change alone.

### Commit 1 — correctness hardening (B-01, B-02, B-03, M-01, M-02, M-03, MIN-03)

Verified with `pytest tests/ -q` → **142 passed** and
`python scripts/check_invariants.py` → passed (rdkit 2025.09.3, offline).

| ID | What changed | Verify (tests) |
|----|--------------|----------------|
| B-01 | `xyz_to_gaussian_coords` parses by **physical line** (L1 = count `N`, L2 = comment which may be empty, next `N` = coords); count mismatch or malformed row raises `ValueError` instead of silently dropping atoms. Only trailing blank lines tolerated. | `tests/test_gaussian.py::TestXyzParsingByPhysicalLine` (empty comment keeps all atoms; count≠rows both directions raise; malformed/non-integer raise; trailing blank tolerated) |
| B-02 | Resume now invalidates on per-molecule **`cid`/`smiles`** identity (matched against the requested molecule) and run-level **`rdkit_version`**, and requires each recorded `xyz_path` to exist and be non-empty. `pipeline_commit`-in-key / row-count reconcile / duplicate-label keying explicitly deferred (call 1b). | `tests/test_conformers.py::TestRowIdentityAndXyz`, `::TestRowConfigMatches::test_rdkit_version_mismatch`, `::TestResumePartition`, `::TestResumeConfigValidationBatch` (changed cid/smiles/rdkit/deleted-xyz each regenerate; unchanged skips) |
| B-03 | SLURM template resolves the `.com` **relative to the script's own location** (`os.path.relpath(com_path, slurm_dir)`), `cd`s into the input's directory and runs `g16` on the basename — submission is directory-independent, custom dir names preserved (no hardcoded `../gaussian_inputs`). | `tests/test_slurm.py::TestSlurmScriptResolvesInput` |
| M-01 | `write_slurm_scripts` **defaults** to consuming the current run's `com_write_log.csv`; `com_dir` glob kept as an explicit non-default mode; notebook default switched. | `tests/test_slurm.py::TestWriteSlurmScriptsLogDriven` (3 logged + 1 stale `.com` → 3 scripts; legacy glob still reachable) |
| M-02 | `_resume_partition` gains `preserve_unrequested: bool = False` (default drops unrequested molecules from the new log); `search_conformers(append=True)` retains carry-forward (call 2a). | `tests/test_conformers.py::TestPreserveUnrequestedBatch`, `::TestResumePartition::test_*_unrequested` |
| M-03 | `write_slurm_scripts` **overwrites** `.sh` by default (dropped `SKIPPED_EXISTS`); status column reports `WROTE`/`OVERWROTE` (call 4a). | `tests/test_slurm.py::TestWriteSlurmScriptsOverwrite` |
| MIN-03 | `search_conformers` validates at entry: `n_generate<1`, `top_n<1`, `rmsd_prune<0`, duplicate molecule labels, empty sanitized filename → `ValueError`. | `tests/test_conformers.py::TestParameterValidation` |

No scientific invariant (§2) changed: route lines, charge/mult, units, and the
Link1 opt→freq contract are untouched. B-01/B-02 make previously-silent geometry
corruption/reuse loud (`ValueError` / regeneration), strengthening the
reproducibility and "no silent corruption" invariants rather than altering them.

### Commit 2 — repository cleanup (B-04, MIN-02)

Verified with `git ls-files -ci --exclude-standard` → empty and
`pytest tests/ -q` → **143 passed**.

| ID | What changed | Verify |
|----|--------------|--------|
| B-04 | `.gitignore` gains the genuinely-missing patterns (`.claude/settings.local.json`, `conformer_xyz/`, `pubchem_sdf/`, `pubchem_xyz/`, `molecules_pubchem_*.csv`) and broadens `com_write_failed.csv`→`*_failed.csv`; no existing pattern duplicated. `git rm -r --cached` untracked the already-tracked generated dirs/files (`pipeline/__pycache__`, `notebooks/.ipynb_checkpoints`, `.pubchem_cache`, `gaussian_inputs`, `slurm_scripts`, `pubchem_sdf`, `pubchem_xyz`, and the tracked `notebooks/*_log.csv` + `molecules_pubchem_*.csv`) — local copies kept. New CI step in `review-readiness.yml` fails the floor if any ignored file is tracked. | `git ls-files -ci --exclude-standard` prints nothing; a fresh clone carries no `*.com`/`*.sdf`/`*.xyz` under `notebooks/` |
| MIN-02 | `search_conformers` and both `write_gaussian_coms*` clear a stale `*_failed.csv` at stage start; rewritten only if the run actually fails. | `tests/test_conformers.py::TestStaleFailedCsvCleared` |

### Commit 3 — docs + versioning (MOD-01..MOD-04, MIN-01; MOD-05 PR body)

| ID | What changed |
|----|--------------|
| MOD-01 | `pipeline.__version__ = "2.0.0"`; PubChem UA → `gaussian-input-pipeline/2.0 (research use)`. Deliberately changes the `pipeline_version` stamped into provenance logs (call 5a). |
| MOD-02 | `docs/architecture.md`, `docs/implementation-plan.md`, `docs/implementation-status.md` are now the authoritative v2 files (no template placeholders); the `-v2` architecture/plan/status drafts and every `remediation-plan-round-0N-v2.md` moved under `docs/review-history/v2/`. `check_invariants.py` drift-guard message updated to the canonical file. |
| MOD-03 | README rewritten around the RDKit conformer flow as primary (names → CID+stereo SMILES → ETKDGv3 → MMFF94/UFF rank → top-N XYZ → Gaussian opt→freq → run-scoped SLURM); RDKit added to requirements, `conformers.py` in the tree; Open Babel demoted to a labeled Legacy v1.1 section. |
| MOD-04 | `WORKFLOW.md` rewritten for the native Codex review (no `codex-review.yml`/`review.md`/API-key path); `AGENTS.md` §1 makes the v2 conformer pathway primary. No `OPENAI_API_KEY` path introduced. |
| MIN-01 | Pinned/locked env shipped for v2.0.0 (`environment-lock.yml`); `numpy` removed from explicit deps after confirming it is unused outside `pipeline/`. |
| MOD-05 | PR #3 body rewritten last to match the final diff (test count, rewired notebook, changelog; stray `EOF`/`)` removed). |

This is the final working snapshot for the v2.0.0 release; archived under
`docs/review-history/v2/`. The canonical `docs/implementation-status.md` is the
live merge-gate doc going forward.

## 0c. Round-03 remediation (Codex findings M-06, M-07)

Both findings from `docs/remediation-plan-round-03-v2.md` are addressed as
separate commits. **Resolved only after the listed verify step passed** — not on
code change alone.

| ID | Commit | Verification run | Result |
|----|--------|------------------|--------|
| M-06 | `83cc915` | `pytest tests/ -q` → 91 passed (incl. `TestProvenanceLogging`, `TestGitShortSha`, `TestPipelineProvenance`); `check_invariants.py` → passed; offline ribose run — `pipeline_version=0.2.0` + `pipeline_commit` on every `conformer_log.csv` row, XYZ header carries `pver=/pcommit=`, no-git branch → `""` | **Resolved** |
| M-07 | `b98581f` | canonical `implementation-status.md` populated (grep: no `<` placeholders, no empty bullets); drift guard fails on the template (12 violations) / passes on the populated file; `pytest tests/ -q` → 95 passed (incl. `TestStatusDocDriftGuard`); `check_invariants.py` → passed | **Resolved** |

- **M-06** — `pipeline.__version__ = "0.2.0"` + a best-effort, offline-safe
  `pipeline_provenance` helper; `conformer_log.csv` gains `pipeline_version` and
  `pipeline_commit` columns and each XYZ comment gains `pver=`/`pcommit=` tokens.
  Conformer path only (com/sdf logs out of scope). Additive; ΔE and top-3
  selection unchanged. Reproducibility caveat recorded in §3/§5.
- **M-07** — canonical `docs/implementation-status.md` populated and synced from
  this file; `scripts/check_invariants.py` drift guard makes canonical staleness
  an objective floor failure; `WORKFLOW.md` records the sync step.

## 0b. Round-02 remediation (Codex findings M-03, M-04)

Both findings from `docs/remediation-plan-round-02-v2.md` are addressed as
separate commits. **Resolved only after the listed verify step passed** — not on
code change alone.

| ID | Commit | Verification run | Result |
|----|--------|------------------|--------|
| M-03 | `5f635bd` | `pytest tests/ -q` → 68 passed (incl. `TestScoreCandidateCurrentSchema`); dead-key grep clean (only helper fallback + documented table-column read; no `CanonicalSMILES` reads); `check_invariants.py` → passed | **Resolved** |
| M-04 | `016c197` | `pytest tests/ -q` → 79 passed (incl. the four convergence cases); `check_invariants.py` → passed; offline run — ribose top-3 all `converged=True` (unchanged vs round-01); mocked all-unconverged → one `converged=False` seed + warning + `UNCONVERGED_FF_SEED` in the `.com` title | **Resolved** |
| M-05 | `c5a6974` | `pytest tests/ -q` → 81 passed (incl. `TestRetryAlignment`); `check_invariants.py` → passed; offline alignment run — forced single-conformer retry on ribose → every recorded energy equals the MMFF energy recomputed from the written coordinates (all 8 conformers, <1e-3 kcal/mol) | **Resolved** |

- **M-03** — `score_candidate` now reads the stereo-bearing SMILES via
  `_isomeric_smiles(prop)` instead of the dead `IsomericSMILES` key, so ambiguous
  candidates earn the stereo bonus and the resolver stops picking the wrong CID.
  Repo swept: the property-fetch list already requested `SMILES,ConnectivitySMILES`
  (round-01); the only other dead-key read (`score_candidate`) is fixed. The
  molecule-table column keeps the name `IsomericSMILES` but holds real stereo
  SMILES (from the live `SMILES` key); reading it is not a dead-key read.
- **M-04** — FF convergence is now first-class: `generate_conformers` captures the
  per-conformer `not_converged` flag, retries unconverged conformers once with
  more iterations, and `conformer_log.csv` records a `converged` column.
  `select_converged_top_n` ranks **only converged** conformers (1a); if none
  converge, exactly one lowest-energy best-effort seed is carried with
  `converged=False`, a warning, and an `UNCONVERGED_FF_SEED` marker in the `.com`
  title (2b). The unreliable FF energy stays labeled unreliable.
- **M-05** — correctness follow-up to M-04's retry: the retry now re-optimizes
  **only the conformers that failed** the first pass (`_optimize_single_conf`),
  so already-converged conformers keep the first-pass energy that matches their
  (untouched) geometry. Previously a whole-ensemble retry could move converged
  geometries while retaining their first-pass energies, letting the ranked
  `rel_energy_kcalmol` describe a different geometry than the written XYZ/`.com`.

## 0. Round-01 remediation (Codex findings B-01, M-02)

Both findings from `docs/remediation-plan-round-01-v2.md` are addressed as
separate commits. **Resolved only after the listed verify step passed** — not on
code change alone.

| ID | Commit | Verification run | Result |
|----|--------|------------------|--------|
| B-01 | `6831e30` | `pytest tests/ -q` → 66 passed; `check_invariants.py` → passed; offline 2-row end-to-end (defined-stereo → XYZ + `.com`; undefined-stereo → logged skip, no XYZ) | **Resolved** |
| M-02 | `21dd80d` | offline `TestNotebookPathOffline` in `pytest tests/ -q` → 66 passed; manual notebook run top-to-bottom on demo molecules → per-conformer `_c00_F.com` via conformer path | **Resolved** |

- **B-01** — the resolved molecule table now carries the stereo-bearing SMILES
  into an `IsomericSMILES` column (pure `_resolved_row` helper), and
  `search_conformers` runs `check_conformer_eligibility` before embedding:
  empty/unparseable SMILES or **undefined stereochemistry** → skip + log to
  `conformer_search_failed.csv` (never silent RDKit stereo auto-assignment);
  no-stereocenter molecules proceed. See §3 for the PubChem key-rename deviation
  that this fix had to absorb to actually produce inputs.
- **M-02** — `notebooks/run_pipeline.ipynb` now runs the conformer path by
  default (`build_molecule_table → search_conformers →
  write_gaussian_coms_from_conformers → write_slurm_scripts`); the v1.1
  single-geometry path is a commented-out labeled legacy appendix. README updated
  to match; an offline test exercises the notebook's exact code path in CI.

## 1. What was implemented

- **`pipeline/conformers.py` (new).**
  - `select_top_n(energies_kcal, n) -> list[int]` — pure ranking helper, lowest
    energy first, deterministic stable tie-break, no RDKit needed. (Plan Task 1.)
  - `generate_conformers(smiles, n_generate, rmsd_prune, seed)` — RDKit ETKDGv3
    `EmbedMultipleConfs(numConfs=n_generate, params)` with `randomSeed=seed` and
    `pruneRmsThresh=rmsd_prune`; MMFF94 optimize+score via
    `MMFFOptimizeMoleculeConfs`; UFF fallback with a **logged** warning when
    `MMFFHasAllMoleculeParams` is false. Returns `(coords_list, energies_kcal,
    method)`. Energies in kcal/mol. (Plan Task 2.)
  - `search_conformers(molecule_table, ...)` — batch driver. Reads
    `IsomericSMILES` per row, generates/ranks, keeps top `TOP_N=3` distinct,
    writes `{base}_c{ii}.xyz`, appends one provenance row per kept conformer to
    `conformer_log.csv`. Failures (including missing SMILES) → 
    `conformer_search_failed.csv`. Resume-safe: molecules already in the log are
    skipped. `rel_energy_kcalmol` is ΔE from that molecule's lowest conformer
    (min = 0.0). (Plan Task 3.)
- **`pipeline/gaussian.py` (extended).**
  - `write_gaussian_com` gained optional `conformer_id` and `rel_energy_kcalmol`.
    When `conformer_id` is set, the basename becomes `{base}_c{ii}` →
    `{base}_c{ii}_F.com` / `{base}_c{ii}_F.chk`, and the title records the id and
    `dE=<val> kcal/mol`. When `None`, v1.1 single-geometry naming is unchanged.
    The Link1 opt→freq checkpoint contract is untouched.
  - `write_gaussian_coms_from_conformers(conformer_log_csv, ...)` — new batch
    consumer that writes one `.com` per conformer-log row. (Plan Task 4.)
- **`pipeline/__init__.py`** — exports `generate_conformers`, `search_conformers`,
  `select_top_n`, `write_gaussian_coms_from_conformers`; docstring pipeline steps
  updated. The v1.1 SDF→XYZ path (`download_sdfs`, `convert_sdfs_to_xyz`) is
  retained.
- **`environment.yml`** and **`.github/workflows/review-readiness.yml`** — added
  `rdkit` (offline, no cluster). (Plan Task 6.)
- **README** — "Conformer searching (v2)" section replaces the old "not included"
  caveat; states the three retained limitations honestly. (Plan Task 7.)
- **Tests** — `tests/test_conformers.py` (new) and conformer cases added to
  `tests/test_gaussian.py`. See §4.
- **Round-03 (M-06)** — `pipeline.__version__ = "0.2.0"`; a `pipeline_provenance`
  helper (`utils.py`) returns `(version, best-effort git commit)`, offline-safe
  (`""` when git is absent, `.dirty` on an uncommitted tree). `search_conformers`
  records `pipeline_version` + `pipeline_commit` on every `conformer_log.csv` row
  and appends `pver=`/`pcommit=` to each XYZ comment line. Conformer path only.
- **Round-03 (M-07)** — canonical `docs/implementation-status.md` is populated
  and synced from this file; `scripts/check_invariants.py` gained a drift guard
  that fails when the canonical file still contains template markers; `WORKFLOW.md`
  records the sync step.
- **M-09** — config-validated resume. `search_conformers` records the requested
  `n_generate` and `top_n` per row and, on rerun, skips a molecule only when all
  its rows match this run's config (`seed`, `n_generate`, `top_n`, `rmsd_prune`,
  `pipeline_version`). Config/version drift, or a pre-provenance log, drops the
  stale rows and regenerates them with a warning; molecules not in the current
  table are left untouched. Prevents downstream Gaussian inputs from being built
  on outdated conformers after a settings change.

Required checks locally green after round-03 remediation (incl. M-06/M-07/M-09):
`pytest tests/ -q` → **110 passed**; `python scripts/check_invariants.py` →
**passed**.

## 2. What was NOT implemented (and why)

- Everything under the plan's "Explicitly out of scope": rotatable-bond gating,
  xTB/CREST, energy-window logic, Boltzmann/entropy weighting, solvent-aware
  search, changes to level of theory or the Link1 contract. Deliberately omitted
  per architecture-v2 "Out of scope for v2".
- ~~The notebook was not rewired.~~ **Resolved in round-01 remediation (M-02,
  commit `21dd80d`):** `notebooks/run_pipeline.ipynb` now runs the conformer path
  by default; the v1.1 path is a labeled legacy appendix.

## 3. Deviations from architecture-v2.md / plan

- **`generate_conformers` returns a 4-tuple, not the 2-tuple in the plan.** The
  plan (Task 2) writes `-> (coords_list, energies_kcal)`. It now returns
  `(coords_list, energies_kcal, method, converged)`: `method` (`"MMFF94"`/`"UFF"`)
  was added in round-01 for FF provenance, and `converged` (per-conformer bool)
  was added in **round-02 M-04** because the provenance / "ran ≠ validated"
  invariants require recording whether each FF optimization actually converged.
  API-shape deviation only; no scientific assumption changed. Documented in the
  function docstring.
- **XYZ output directory naming.** Architecture-v2 says "write XYZ per conformer"
  and the plan says `{base}_c{ii}.xyz` without pinning a directory. I default to
  `conformer_xyz/` (parameter `xyz_dir`, overridable), mirroring the existing
  `pubchem_xyz/` convention. Not a scientific change.
- **`conformer_log.csv` column set.** Implemented columns:
  `name, cid, smiles, conformer_id, rel_energy_kcalmol, xyz_path, rdkit_version,
  pipeline_version, pipeline_commit, seed, n_generate, top_n, method, n_generated,
  n_kept, rmsd_prune, converged`. A superset of the plan's required provenance
  columns: `cid` / `smiles` / `n_kept` for traceability
  name→CID→SMILES→conformer, `converged` for M-04, `pipeline_version` /
  `pipeline_commit` for M-06, and the requested `n_generate` / `top_n` for M-09
  (the config a resumed run validates against; distinct from the result columns
  `n_generated` / `n_kept`). All additive; no required column omitted.

### Round-02 remediation deviations

- **M-04 — `generate_conformers` returns a 4-tuple.** The convergence retry
  re-optimizes **only the conformers that failed** the first pass, one at a time
  (`_optimize_single_conf`), so already-converged conformers are never disturbed
  and `_finalize_convergence` keeps their first-pass energies consistent with
  their first-pass geometry (see M-05 below). All-converged ensembles never retry
  and are byte-identical to the pre-M-04 behavior. Implementation detail, not a
  scientific change.
- **M-05 — retry alignment (correctness fix, supersedes the earlier
  whole-ensemble retry).** The initial M-04 retry re-optimized the whole molecule
  in place, which could move already-converged geometries while their first-pass
  energies were retained — leaving ranked `rel_energy_kcalmol` describing a
  different geometry than the written XYZ/`.com`. Now only failed conformers are
  retried, guaranteeing recorded energy ↔ written geometry alignment. Not a
  scientific-assumption change; it removes a geometry/energy mismatch.
- **M-04 — all-fail fallback carries an unminimized geometry (decision 2b).** When
  no conformer converges, one best-effort geometry is still handed to DFT. This is
  an intentional, **flagged** exception to "no placeholder science": it is a real
  FF geometry (not fabricated), explicitly marked `converged=False` in the log,
  warned at runtime, and tagged `UNCONVERGED_FF_SEED` in the `.com` title. The DFT
  optimization is expected to refine it; the FF energy is labeled unreliable.

### Round-03 remediation deviations

- **M-06 — additive provenance columns + extended XYZ comment format.**
  `conformer_log.csv` gains `pipeline_version` and `pipeline_commit`; each
  per-conformer XYZ comment line gains `pver=<version> pcommit=<commit>` tokens.
  Both are additive: existing columns, ordering, ΔE values, and the top-3
  selection are unchanged versus round 02. Log/geometry-comment only — no route
  line, unit, charge/multiplicity, or Link1 edit.
- **M-06 — reproducibility caveat.** `pipeline_commit` pins code identity **only
  when the tree is clean**. A `.dirty` suffix means uncommitted edits produced the
  output, so that output is *not* fully reproducible from the commit alone — the
  flag makes this visible rather than hiding it. An empty `pipeline_commit`
  (no git) falls back to `pipeline_version`, which is only as precise as manual
  bumping. A recorded commit is therefore **not** a guarantee of exact code.
- **M-09 — config-validated resume (behavioral change to resume-safety).** Resume
  previously skipped any molecule already in the log. It now skips only when the
  recorded config (`seed`, `n_generate`, `top_n`, `rmsd_prune`, `pipeline_version`)
  matches this run; otherwise the stale rows are dropped and the molecule is
  regenerated (with a warning). This means a rerun with changed settings now does
  real work where it previously no-oped — the point of the fix, so downstream
  inputs are never built on outdated conformers. `pipeline_version` (not the
  per-commit `pipeline_commit`) is the code-change signal used, since commit is
  empty without git; a `.dirty` edit that changes generation logic without a
  version bump is therefore not auto-detected (bump `__version__` to force it).

### Round-01 remediation deviations

- **B-01 — PubChem SMILES property rename (deviation from the plan's literal
  step).** The remediation plan (B-01 step 1) said to add
  `prop.get("IsomericSMILES", "")`. Implementing exactly that still produced
  **zero** conformers: PubChem renamed its SMILES properties in 2025
  (`IsomericSMILES`→`SMILES`, `CanonicalSMILES`→`ConnectivitySMILES`), so
  `prop.get("IsomericSMILES")` is now always empty and every molecule dead-ended
  at `"no IsomericSMILES"`. **Verified against live PubChem** (CID 5950 L-alanine):
  the stereo-bearing SMILES arrives under `"SMILES"` (`C[C@@H](N)C(=O)O`), while
  `"ConnectivitySMILES"` drops stereo (`CC(N)C(=O)O`). The fix reads the
  stereo-bearing SMILES via `_isomeric_smiles()` — prefer `"SMILES"`, fall back to
  legacy `"IsomericSMILES"`, **never** `ConnectivitySMILES`/`CanonicalSMILES`.
  `get_props_by_cids` now requests the current property names. This is faithful to
  B-01's intent (get the stereo SMILES into resolved rows) and changes **no**
  scientific assumption — stereochemistry is still preserved, not dropped. This is
  the item the reviewer should scrutinize; see §6 Q5.
- **B-01 — stereo skip gate.** `check_conformer_eligibility` skips (and logs)
  molecules with undefined stereochemistry rather than letting RDKit assign a
  stereoisomer. This is the plan's locked decision 2a, not a deviation, but it is
  a deliberate **behavioral** choice: such molecules produce no inputs by design.

No scientific invariant (§2 of AGENTS.md) was altered: route lines, units,
charge/multiplicity handling, and the Link1 contract are unchanged. Conformer
energies are labeled kcal/mol at every surface (CSV column name, XYZ comment,
`.com` title) and never mixed with DFT Hartree values.

## 4. Tests added

- `tests/test_conformers.py`
  - `TestSelectTopN::*` — ordering, top-N truncation, n≥len, n≤0, empty input,
    deterministic tie-break, single conformer. (No RDKit.)
  - `TestGenerateConformers::test_butane_seeded_deterministic` — n-butane yields
    ≥1 conformer, `method == "MMFF94"`, atom tuples are `(sym, x, y, z)` floats,
    and the seeded lowest-energy index + full energy list reproduce exactly.
  - `TestGenerateConformers::test_invalid_smiles_raises` — bad SMILES → ValueError.
  - `TestSearchConformers::test_two_molecule_table` — adenine collapses to 1 row,
    ribose yields 1–3 rows, provenance columns populated, per-molecule min ΔE = 0
    and all ΔE ≥ 0, XYZ files exist.
  - `TestSearchConformers::test_missing_smiles_logged_no_isomericsmiles` — empty
    SMILES → skipped + logged `"no IsomericSMILES"`, no XYZ written. **(B-01)**
  - `TestSearchConformers::test_undefined_stereo_is_skipped_and_logged` —
    undefined-stereo sugar → logged `"undefined stereochemistry"`, no XYZ, while a
    no-stereo molecule in the same table proceeds. **(B-01)**
  - `TestSearchConformers::test_resume_skips_completed` — rerun appends nothing.
  - `TestCheckEligibility::*` — empty/None/whitespace (pure), no-stereocenter and
    defined-stereo eligible, undefined-stereo and unparseable skipped. **(B-01)**
  - `TestReproducibility::test_same_seed_same_selected_conformers` — same seed →
    identical `conformer_id` + `rel_energy_kcalmol` (validation ≠ "it ran").
  - `TestNotebookPathOffline::test_conformer_path_produces_per_conformer_coms` —
    the exact notebook code path → per-conformer `_c{ii}_F.com`, Link1 + ΔE.
    **(M-02)**
  - `TestSelectConvergedTopN::*` — only converged conformers eligible; top-N caps;
    all-unconverged → one best-effort seed + `all_failed=True`; empty input.
    (Pure.) **(M-04 1a/2b)**
  - `TestFinalizeConvergence::*` — first-pass converged kept; retry-converged
    included; still-unconverged-after-retry flagged False; mixed uses retry only
    for the failed conformer. (Pure.) **(M-04 retry)**
  - `TestRetryAlignment::*` — the retry re-optimizes only the failed conformer id
    (converged conformers keep first-pass energy/geometry); no retry call when all
    converge. Guards recorded-energy ↔ written-geometry alignment. **(M-05)**
  - `TestConvergenceBatch::*` (RDKit stubbed): mixed batch keeps only converged
    (logged `converged=True`); all-unconverged logs one `converged=False` seed;
    all-unconverged carries exactly one row, emits a warning, and the `.com` title
    carries `UNCONVERGED_FF_SEED` (Link1 intact). **(M-04)**
  - `TestGenerateConformers::test_butane_seeded_deterministic` also now asserts the
    4th return value `converged` is an all-True bool list. **(M-04)**
- `tests/test_pubchem.py` **(B-01, M-03)**
  - `TestIsomericSmiles::*` — prefers current `"SMILES"` key, falls back to legacy
    `"IsomericSMILES"`, never returns the stereo-free `ConnectivitySMILES`. **(B-01)**
  - `TestResolvedRow::*` — resolved rows carry stereo SMILES in the
    `IsomericSMILES` column; unresolved rows carry an empty column; schema matches
    `MOLECULE_TABLE_COLUMNS`; downstream-consumer columns present. **(B-01)**
  - `TestScoreCandidateCurrentSchema::*` — a current-schema record (stereo under
    `"SMILES"`, no legacy key) earns the stereo bonus, and a chiral higher-CID
    candidate outscores a lower-CID achiral one; both fail under the old dead-key
    read. **(M-03)**
- `tests/test_conformers.py` (continued)
  - `TestProvenanceLogging::*` — `conformer_log.csv` has `pipeline_version` ==
    `pipeline.__version__` on every row, `pipeline_commit` present and a `str`
    (never a concrete SHA asserted), and the XYZ comment carries `pver=`/`pcommit=`
    tokens. **(M-06)**
  - `TestRowConfigMatches::*` / `TestResumePartition::*` — pure config-match and
    partition logic: seed/`n_generate`/`top_n`/`rmsd_prune`/`pipeline_version`
    mismatch (and a missing column) count as drift; all rows of a molecule must
    match; unrequested molecules are carried forward. **(M-09)**
  - `TestResumeConfigValidationBatch::*` — batch (RDKit stubbed): matching config
    resumes without regenerating; changed seed / `top_n` and a pre-provenance log
    regenerate stale rows with a warning; a molecule absent from the current table
    is preserved. **(M-09)**
- `tests/test_utils.py` **(M-06)**
  - `TestGitShortSha::*` — `""` when git is absent / non-zero / times out; the SHA
    on a clean tree; SHA + `.dirty` on a mocked dirty tree. Offline (subprocess
    mocked); no concrete SHA asserted.
  - `TestPipelineProvenance::*` — returns `(pipeline.__version__, commit)` and
    falls back to `""` commit when git fails.
- `tests/test_check_invariants.py` **(M-07)**
  - `TestStatusDocDriftGuard::*` — the guard fires on a template file (angle-bracket
    placeholders + empty bullets) and passes on a populated one.
- `tests/test_gaussian.py`
  - `TestWriteGaussianComConformer::*` — conformer filename/chk (`ribose_c00_F`),
    ΔE in title (`dE=1.2345 kcal/mol`), Link1 preserved, and `conformer_id=None`
    keeps v1.1 naming.
  - `TestWriteGaussianComsFromConformers::test_three_conformer_rows_three_coms` —
    3 log rows → 3 correctly named `.com` files, each with intact Link1 and its ΔE.

RDKit-dependent tests use `pytest.importorskip("rdkit")` so a bare environment
still runs the pure tests; CI installs rdkit so they execute there.

## 5. Known limitations (retained honestly, not solved in v2)

- MMFF ranking is unreliable for intramolecular-H-bonding species (sugars,
  nucleosides). Carrying top 3 mitigates; xTB rerank is future work.
- Force-field ranking is gas-phase while the DFT default is IEFPCM water; the
  FF-lowest conformer may not be the solution-phase minimum.
- Fixed `N_GENERATE=20` may under-sample very flexible molecules.
- These are FF starting geometries, not optimized minima — DFT makes the final
  call among carried candidates.
- `pipeline_commit` reproducibility (M-06): a recorded commit pins code only when
  clean; `.dirty` means not reproducible from the commit alone; empty means the
  fallback to the manually bumped `pipeline_version`. See §3 round-03 deviations.

## 6. Questions requiring scientific judgment  ← Ish reads this FIRST

1. ~~Should the notebook be wired to the conformer stage?~~ **Resolved by M-02
   (Ish decision 1a):** the notebook now runs the conformer path by default.
2. **PubChem SMILES property rename (B-01) — please sanity-check the chemistry.**
   PubChem's current `"SMILES"` property is the stereo-bearing SMILES we now
   consume; `"ConnectivitySMILES"` is stereo-free and deliberately unused. I
   verified this against live PubChem for L-alanine and D-glucose (both show `@`
   in `"SMILES"`, none in `"ConnectivitySMILES"`). Confirm you are comfortable
   treating PubChem `"SMILES"` as the stereochemistry source of record. If PubChem
   ever emits a `"SMILES"` without stereo for a molecule that *has* stereocenters,
   the eligibility gate will (correctly) skip it as "undefined stereochemistry"
   rather than embed an arbitrary isomer.
3. **`TOP_N=3` distinctness relies solely on `pruneRmsThresh=0.5 Å` at embed
   time.** After MMFF optimization, two kept conformers could in principle relax
   toward each other; v2 does not re-prune post-optimization (matches "no separate
   duplicate-removal code" in architecture-v2). Acceptable, or do you want a
   post-optimization RMSD check? (Would be new scope.)
4. **UFF fallback energies vs MMFF94 energies are on different scales.** Both are
   labeled kcal/mol and `method` is recorded per row, but ΔE values from a UFF
   molecule are not comparable to MMFF94 ΔE from another molecule. Confirm that
   per-molecule ΔE (never cross-molecule) is the only intended comparison.
5. **Charge/multiplicity for conformers.** The conformer stage builds neutral
   molecules from `IsomericSMILES` (RDKit default protonation). Non-default
   charge/multiplicity still must be passed to the Gaussian writer separately, as
   in v1.1. No per-molecule charge handling was added. Confirm this is expected.

## Provenance

- pipeline version: `pipeline.__version__ = "0.2.0"` (M-06), recorded per-row in
  `conformer_log.csv` `pipeline_version`, with a best-effort git short SHA in
  `pipeline_commit` (`.dirty` on an uncommitted tree, `""` when git is absent).
  Branch `feat/conformer-search-v2`; base v1.1, Zenodo 10.5281/zenodo.18894724.
- RDKit version used for local test run: 2025.09.3 (recorded per-row in
  `conformer_log.csv` `rdkit_version` at runtime).
- Conformer stage config: `N_GENERATE=20`, `TOP_N=3`, `RMSD_PRUNE=0.5 Å`,
  `SEED=42`, ranking MMFF94 (UFF logged fallback), energies kcal/mol.
- Open Babel: not used on the conformer path (RDKit consumes `IsomericSMILES`
  directly). The v1.1 SDF→XYZ path still uses Open Babel where invoked.
- Gaussian route lines / charge / mult / nproc: unchanged from v1.1; still passed
  in from user config, not constructed in `pipeline/`.
