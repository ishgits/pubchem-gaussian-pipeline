# implementation-plan.md

> Canonical v2 plan. Approved by Ish before implementation; works against
> `architecture.md` (v2). The original draft and the per-round remediation plans
> are archived under `docs/review-history/v2/`.

## Objective

Add a simple RDKit conformer-search stage: generate an ensemble, MMFF-rank it,
carry the top 3 distinct conformers per molecule into Gaussian input generation,
with recorded provenance. No gating, no xTB, no energy-window logic. Ship it as
the released v2.0.0 default alongside run-scoped, submission-independent SLURM
scripts and identity/config-aware resume.

## Locked defaults

`N_GENERATE=20`, `TOP_N=3`, `RMSD_PRUNE=0.5` Å, `SEED=42`, ranking=MMFF94
(UFF logged fallback). Energy units kcal/mol throughout the conformer stage.

## Tasks (ordered, each independently verifiable)

1. **`pipeline/conformers.py` — pure ranking helper.** `select_top_n(energies_kcal,
   n)` returns lowest-energy-first indices. Acceptance: unit test on a synthetic
   list; no RDKit needed.
2. **RDKit embed + rank core.** `generate_conformers(smiles, n_generate,
   rmsd_prune, seed)`: ETKDGv3 `EmbedMultipleConfs(pruneRmsThresh, randomSeed)`,
   MMFF94 optimize + score, UFF fallback with a logged warning. Acceptance: a
   seeded n-butane test asserts ≥1 conformer and a deterministic lowest-energy
   index.
3. **Batch `search_conformers(molecule_table, ...)`.** Per row: read
   `IsomericSMILES`, generate, rank, keep top 3 distinct, write `{base}_c{ii}.xyz`,
   append provenance rows to `conformer_log.csv`; failures →
   `conformer_search_failed.csv`. Acceptance: a 2-molecule run (adenine, ribose)
   — adenine collapses to 1 row, ribose ≤3 rows, provenance populated, ΔE kcal/mol.
4. **Extend `gaussian.py` to consume conformers.** Write one `.com` per row as
   `{base}_c{ii}_F.com`, preserve the Link1 section, write conformer id + ΔE into
   the title. Acceptance: 3 rows → 3 distinct `.com` files, intact Link1.
5. **Reproducibility check.** Validation = same seed reproduces the same selected
   conformers (not "embedding succeeded"). Acceptance: rerun-with-same-seed test.
6. **Environment + CI.** Add `rdkit` to `environment.yml` and the CI install step.
   Acceptance: `review-readiness` stays green.
7. **Release hardening (v2.0.0).** XYZ parsing is by physical line (never drops an
   atom); resume is identity/config/RDKit-version aware and XYZ-existence checked;
   SLURM scripts resolve their `.com` relative to their own location and default
   to the current run's `com_write_log.csv`, overwriting stale scripts; parameters
   are validated at entry; generated outputs are untracked and gitignored.
   Acceptance: the round-04 verify gate (see `implementation-status.md`).
8. **Docs + status + versioning.** `pipeline.__version__ = "2.0.0"`, PubChem UA
   `gaussian-input-pipeline/2.0`; README leads with the RDKit flow (Open Babel
   demoted to a labeled legacy section); `implementation-status.md` current.

## Explicitly out of scope

Rotatable-bond gating; xTB/CREST; energy-window logic; Boltzmann/entropy
weighting; solvent-aware search; running or parsing Gaussian; changing the level
of theory or the Link1 contract; the `runs/<study>/` directory redesign; a fuller
resume key (`pipeline_commit` in-key, row-count reconcile, duplicate-label
re-keying); version/commit provenance in `com_write_log.csv` /
`sdf_download_log.csv`.

## Files expected to change

- `pipeline/conformers.py` (new), `pipeline/gaussian.py`, `pipeline/slurm.py`,
  `pipeline/pubchem.py`, `pipeline/__init__.py`.
- `environment.yml`, `.github/workflows/review-readiness.yml`, `.gitignore`.
- `tests/test_conformers.py` (new), `tests/test_gaussian.py`, `tests/test_slurm.py`.
- `notebooks/run_pipeline.ipynb`, `README.md`, `WORKFLOW.md`, `AGENTS.md`,
  `docs/architecture.md`, `docs/implementation-plan.md`,
  `docs/implementation-status.md`.

## Escalate, don't decide

If MMFF params are missing for a target species, record whether UFF-fallback was
used or the molecule was skipped — never silently substitute chemistry.
