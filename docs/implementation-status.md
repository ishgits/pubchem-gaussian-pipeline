# implementation-status.md

> **Synced each remediation round** from `docs/implementation-status-v2.md` (the
> working file). AGENTS.md §5 names *this* file as the merge gate, so it must
> always reflect the real current status — `scripts/check_invariants.py` fails
> the objective floor if this file ever drifts back to the empty template.

Maintained by the implementing agent (Claude Code). The reviewer (Codex) must
verify these claims against the code, not trust them.

**PR:** #3   **Branch:** `feat/conformer-search-v2`   **Round:** 3
(v2 conformer search + remediation rounds 01–03). Works against
`docs/architecture-v2.md` and `docs/implementation-plan-v2.md`. The per-round
verification evidence (commit hashes, verify runs) lives in
`implementation-status-v2.md` and each `docs/remediation-plan-round-0N-v2.md`.

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
- **PubChem SMILES sourcing (round 01/02, B-01/M-03).** The resolved molecule
  table carries the stereo-bearing SMILES into an `IsomericSMILES` column via the
  `_isomeric_smiles` helper, and `score_candidate` reads stereo through the same
  helper. Both use PubChem's current `SMILES` key (see §3 deviation).
- **FF convergence handling (round 02, M-04/M-05).** `generate_conformers`
  captures each conformer's `not_converged` flag and retries only the failed
  conformers once with more iterations (`_optimize_single_conf`), so
  already-converged conformers keep the first-pass energy that matches their
  written geometry. `select_converged_top_n` ranks only converged conformers;
  if none converge, exactly one lowest-energy best-effort seed is carried with
  `converged=False`, a warning, and an `UNCONVERGED_FF_SEED` marker in the `.com`
  title. `conformer_log.csv` records a `converged` column.
- **Version + commit provenance (round 03, M-06).** `pipeline.__version__` is
  `"0.2.0"`; `pipeline_provenance()` returns that version plus a best-effort git
  short SHA (`.dirty` suffix on an uncommitted tree, empty string when git is
  absent). `search_conformers` records `pipeline_version` and `pipeline_commit`
  on every `conformer_log.csv` row and appends `pver=`/`pcommit=` tokens to each
  per-conformer XYZ comment line.
- **Gaussian writer (extended).** `write_gaussian_com` gained optional
  `conformer_id`, `rel_energy_kcalmol`, and `unconverged` parameters:
  per-conformer basenames `{base}_c{ii}_F.com` / `.chk`, the ΔE (kcal/mol) and any
  `UNCONVERGED_FF_SEED` marker in the title. The Link1 opt→freq checkpoint
  contract is untouched. `write_gaussian_coms_from_conformers` writes one `.com`
  per `conformer_log.csv` row.
- **Notebook (round 01, M-02).** `notebooks/run_pipeline.ipynb` runs the
  conformer path by default; the v1.1 single-geometry path is a commented-out
  labeled legacy appendix.
- **Canonical status doc + drift guard (round 03, M-07).** This file is populated
  and synced from `implementation-status-v2.md`; `scripts/check_invariants.py`
  gained a guard that fails when this file still contains template markers.
- **Config-validated resume (M-09).** `search_conformers` records the requested
  `n_generate` and `top_n` per row and, on rerun, skips a molecule only when all
  its existing rows match this run's config (`seed`, `n_generate`, `top_n`,
  `rmsd_prune`, `pipeline_version`). Rows from a different config, or a
  pre-provenance log, are treated as stale, dropped, and regenerated with a
  warning, so downstream Gaussian inputs are never built on outdated conformers.
  Molecules present in the log but not in the current table are left untouched.
- **Env / CI / README.** `rdkit` added to `environment.yml` and the
  review-readiness CI install step; README documents the conformer stage, the
  skip-on-undefined-stereo rule, and the retained limitations.

Required checks locally green after round 03 + M-09: `pytest tests/ -q` →
**110 passed**; `python scripts/check_invariants.py` → **passed**.

## 2. What was NOT implemented (and why)

- Out-of-scope-for-v2 items: rotatable-bond gating, xTB/CREST, energy-window
  logic, Boltzmann/entropy weighting, solvent-aware search, and any change to the
  level of theory or the Link1 contract. Deliberately omitted per architecture-v2.
- Version/commit provenance in `com_write_log.csv` and `sdf_download_log.csv`
  (and the `.com` title line) is deferred: round 03 (M-06) is scoped to the
  conformer path only. AGENTS.md §3 applies to every generated output, so this is
  tracked as a next-round candidate, not a silent omission.
- Post-optimization RMSD re-pruning is not added; distinctness relies on
  `pruneRmsThresh` at embed time (see §6).

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
  unchanged. The XYZ comment format is likewise extended with `pver=`/`pcommit=`
  tokens (M-06). `n_generate` and `top_n` record the *requested* search knobs
  (distinct from the result columns `n_generated` / `n_kept`) so a resumed run can
  validate that stale rows match the current config (M-09).
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
  and config-validated resume — matching config resumes, a changed
  seed/`top_n`/pre-provenance log regenerates stale rows with a warning, and an
  unrequested molecule is preserved (`TestRowConfigMatches`, `TestResumePartition`,
  `TestResumeConfigValidationBatch`).
- `tests/test_pubchem.py` — SMILES key handling (`TestIsomericSmiles`), resolved-row
  schema (`TestResolvedRow`), and current-schema scoring with the stereo bonus
  (`TestScoreCandidateCurrentSchema`).
- `tests/test_gaussian.py` — per-conformer filenames, ΔE in title, Link1 intact
  (`TestWriteGaussianComConformer`, `TestWriteGaussianComsFromConformers`).
- `tests/test_utils.py` — the offline provenance helper: git absent / non-zero /
  timeout give empty string, clean tree gives the SHA, dirty tree appends
  `.dirty` (`TestGitShortSha`, `TestPipelineProvenance`). No test asserts a
  concrete SHA.
- `tests/test_check_invariants.py` — the status-doc drift guard fires on a
  template and passes on a populated file (`TestStatusDocDriftGuard`).

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
- **Reproducibility of `pipeline_commit`.** It pins code identity only when the
  tree is clean. A `.dirty` suffix means uncommitted edits produced the output,
  so that output is not fully reproducible from the commit alone. An empty
  `pipeline_commit` (no git) falls back to `pipeline_version`, which is only as
  precise as manual version bumping. A recorded commit is therefore not a
  guarantee of exact code.

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

- pipeline version: `0.2.0` (`pipeline.__version__`; recorded per row in
  `conformer_log.csv` `pipeline_version`, and best-effort git commit in
  `pipeline_commit`). Branch `feat/conformer-search-v2`; base v1.1
  (Zenodo 10.5281/zenodo.18894724).
- RDKit version used for local test runs: 2025.09.3 (recorded per row in
  `conformer_log.csv` `rdkit_version` at runtime).
- Conformer stage config: `N_GENERATE=20`, `TOP_N=3`, `RMSD_PRUNE=0.5 Å`,
  `SEED=42`, ranking MMFF94 (UFF logged fallback), energies kcal/mol.
- Open Babel: not used on the conformer path (RDKit consumes `IsomericSMILES`
  directly). The v1.1 SDF-to-XYZ path still uses Open Babel where invoked.
- Gaussian route lines, charge, multiplicity, nproc: unchanged from v1.1; still
  passed in from user config, not constructed in `pipeline/`.
