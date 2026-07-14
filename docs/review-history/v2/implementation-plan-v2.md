# implementation-plan.md

> Human gate: approved by Ish before implementation. Works against architecture.md (v2).
> Scope is intentionally small. Do not expand it without recording a deviation.

## Objective

Add a simple RDKit conformer-search stage: generate an ensemble, MMFF-rank it,
carry the top 3 distinct conformers per molecule into Gaussian input generation,
with recorded provenance. No gating, no xTB, no energy-window logic.

## Locked defaults (confirm at approval)

`N_GENERATE=20`, `TOP_N=3`, `RMSD_PRUNE=0.5` Å, `SEED=42`, ranking=MMFF94
(UFF logged fallback). Energy units kcal/mol throughout the conformer stage.

## Tasks (ordered, each independently verifiable)

1. **`pipeline/conformers.py` — pure ranking helper.**
   - `select_top_n(energies_kcal, n) -> ordered_indices` (lowest energy first).
   - Acceptance: unit test on a synthetic energy list; no RDKit needed.

2. **RDKit embed + rank core.**
   - `generate_conformers(smiles, n_generate, rmsd_prune, seed) -> (coords_list, energies_kcal)`:
     ETKDGv3 `EmbedMultipleConfs(pruneRmsThresh=rmsd_prune, randomSeed=seed)`,
     MMFF94 optimize + score; UFF fallback with a logged warning if MMFF params
     are missing.
   - Acceptance: one seeded integration test on n-butane asserts ≥1 conformer and
     a deterministic lowest-energy index under the fixed seed.

3. **Batch `search_conformers(molecule_table, ...)`.**
   - For each row: read `IsomericSMILES`, generate, rank, keep top 3 distinct,
     write `{base}_c{ii}.xyz`, append to `conformer_log.csv` with provenance
     columns (rdkit version, seed, method, n_generated, n_kept; per conformer:
     id, rel_energy_kcalmol, xyz_path). Failures → `conformer_search_failed.csv`.
     Rerun skips molecules already in the log.
   - Acceptance: run on a 2-molecule table (adenine, ribose); adenine collapses
     to 1 row, ribose yields ≤3 rows; provenance columns populated; ΔE in kcal/mol.

4. **Extend `gaussian.py` to consume conformers.**
   - Read `conformer_log.csv`; write one `.com` per row as `{base}_c{ii}_F.com`;
     preserve the Link1 section; write conformer id + ΔE (kcal/mol) into the title.
   - Acceptance: 3 conformer rows → 3 distinct `.com` files, correct names, intact
     Link1, ΔE recorded in each title.

5. **Reproducibility check.**
   - Do NOT treat "embedding succeeded" as validation. Validation = same seed
     reproduces the same selected conformers.
   - Acceptance: rerun-with-same-seed test yields identical selected conformers.

6. **Environment + CI.**
   - Add `rdkit` to `environment.yml` and to the CI install step in
     `.github/workflows/review-readiness.yml` (offline, no cluster — fine per the
     no-network test rule).
   - Acceptance: `review-readiness` stays green.

7. **Docs + status.**
   - README "Important Caveats": conformer search now included; state the three
     retained limitations honestly (MMFF unreliable for H-bonded species;
     gas-phase FF vs solution-phase DFT; fixed sampling). Keep it short.
   - Update `docs/architecture.md` change log; write `docs/implementation-status.md`
     with the §6 scientific-judgment items.

## Files expected to change

- `pipeline/conformers.py` (new)
- `pipeline/gaussian.py` (consume conformer log; conformer filenames + title)
- `pipeline/__init__.py` (exports)
- `environment.yml`, `.github/workflows/review-readiness.yml` (add rdkit)
- `tests/test_conformers.py` (new), `tests/test_gaussian.py` (conformer cases)
- `README.md`, `docs/architecture.md`, `docs/implementation-status.md`

## Explicitly out of scope (do not add)

Rotatable-bond gating; xTB/CREST; energy windows; Boltzmann weighting;
solvent-aware search; changing level of theory or the Link1 contract.

## Escalate, don't decide

If MMFF params are missing for a target species, record whether UFF-fallback was
used or the molecule was skipped — do not silently substitute chemistry.
