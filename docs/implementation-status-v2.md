# implementation-status-v2.md

Maintained by the implementing agent (Claude Code). The reviewer (Codex) must
verify these claims against the code, not trust them.

**PR:** #<pending>   **Branch:** `feat/conformer-search-v2`   **Round:** 1
Works against `docs/architecture-v2.md` and `docs/implementation-plan-v2.md`.

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

Required checks locally green: `pytest tests/ -q` → **48 passed**;
`python scripts/check_invariants.py` → **passed**.

## 2. What was NOT implemented (and why)

- Everything under the plan's "Explicitly out of scope": rotatable-bond gating,
  xTB/CREST, energy-window logic, Boltzmann/entropy weighting, solvent-aware
  search, changes to level of theory or the Link1 contract. Deliberately omitted
  per architecture-v2 "Out of scope for v2".
- The notebook (`notebooks/run_pipeline.ipynb`) was **not** rewired to call the
  conformer stage. The plan's "Files expected to change" list does not include
  the notebook, and the module API + README are the deliverable surface. The new
  functions are importable and documented; wiring the notebook end-to-end is left
  as a follow-up to avoid scope creep. Flagged as a question in §6.

## 3. Deviations from architecture-v2.md / plan

- **`generate_conformers` returns a 3-tuple, not the 2-tuple in the plan.** The
  plan (Task 2) writes `-> (coords_list, energies_kcal)`. I added `method`
  (`"MMFF94"`/`"UFF"`) as a third element because provenance (architecture-v2
  "Units & provenance") requires recording *which* force field actually ran, and
  that decision is made inside this function. No scientific assumption changed;
  this is an API-shape deviation only. Documented in the function docstring.
- **XYZ output directory naming.** Architecture-v2 says "write XYZ per conformer"
  and the plan says `{base}_c{ii}.xyz` without pinning a directory. I default to
  `conformer_xyz/` (parameter `xyz_dir`, overridable), mirroring the existing
  `pubchem_xyz/` convention. Not a scientific change.
- **`conformer_log.csv` column set.** Implemented columns:
  `name, cid, smiles, conformer_id, rel_energy_kcalmol, xyz_path, rdkit_version,
  seed, method, n_generated, n_kept, rmsd_prune`. This is a superset of the plan's
  required provenance columns (it adds `cid`, `smiles`, `n_kept` for traceability
  name→CID→SMILES→conformer). No required column omitted.

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
  - `TestSearchConformers::test_missing_smiles_is_logged_not_skipped` — missing
    SMILES is recorded in the failed CSV, not silently dropped.
  - `TestSearchConformers::test_resume_skips_completed` — rerun appends nothing.
  - `TestReproducibility::test_same_seed_same_selected_conformers` — same seed →
    identical `conformer_id` + `rel_energy_kcalmol` (validation ≠ "it ran").
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

## 6. Questions requiring scientific judgment  ← Ish reads this FIRST

1. **Should the notebook be wired to the conformer stage in this PR, or in a
   follow-up?** I kept the notebook unchanged (not in the plan's file list) and
   delivered the module API + README. Confirm the intended integration surface.
2. **`TOP_N=3` distinctness relies solely on `pruneRmsThresh=0.5 Å` at embed
   time.** After MMFF optimization, two kept conformers could in principle relax
   toward each other; v2 does not re-prune post-optimization (matches "no separate
   duplicate-removal code" in architecture-v2). Acceptable, or do you want a
   post-optimization RMSD check? (Would be new scope.)
3. **UFF fallback energies vs MMFF94 energies are on different scales.** Both are
   labeled kcal/mol and `method` is recorded per row, but ΔE values from a UFF
   molecule are not comparable to MMFF94 ΔE from another molecule. Confirm that
   per-molecule ΔE (never cross-molecule) is the only intended comparison.
4. **Charge/multiplicity for conformers.** The conformer stage builds neutral
   molecules from `IsomericSMILES` (RDKit default protonation). Non-default
   charge/multiplicity still must be passed to the Gaussian writer separately, as
   in v1.1. No per-molecule charge handling was added. Confirm this is expected.

## Provenance

- pipeline version: v2.0 (branch `feat/conformer-search-v2`; base v1.1,
  Zenodo 10.5281/zenodo.18894724).
- RDKit version used for local test run: 2025.09.3 (recorded per-row in
  `conformer_log.csv` `rdkit_version` at runtime).
- Conformer stage config: `N_GENERATE=20`, `TOP_N=3`, `RMSD_PRUNE=0.5 Å`,
  `SEED=42`, ranking MMFF94 (UFF logged fallback), energies kcal/mol.
- Open Babel: not used on the conformer path (RDKit consumes `IsomericSMILES`
  directly). The v1.1 SDF→XYZ path still uses Open Babel where invoked.
- Gaussian route lines / charge / mult / nproc: unchanged from v1.1; still passed
  in from user config, not constructed in `pipeline/`.
