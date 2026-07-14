# architecture.md

> Human gate: approved by Ish before implementation. v2 revision.
> Scope deliberately minimal: one simple conformer-search stage, top-3 to DFT.

## What v2 adds

A conformer-search stage between molecule resolution and Gaussian input
generation. For each molecule the pipeline generates an RDKit conformer
ensemble, MMFF-optimizes and ranks it, and carries the **top 3 lowest-energy,
distinct** conformers forward as DFT starting geometries — each with recorded
provenance (seed, method, relative energy).

## Why

`geometry.py` currently hands DFT a single force-field conformer. For flexible
molecules that conformer is often not the global minimum, so Gaussian can
optimize into a local minimum and yield wrong relative energies — corrupting the
comparative thermodynamics this pipeline feeds. Carrying the top 3 lets the DFT
step, not the force field, make the final call among the best candidates.

## Deliberately simple: no gating

There is **no rotatable-bond branching**. Every molecule goes through the same
path. RMSD pruning + "keep top 3 distinct" makes rigid molecules collapse to a
single conformer on their own (adenine's conformers are ~identical → pruned to
1 → one DFT job). So v1.1 behavior for rigid molecules is preserved emergently,
with no special-case code.

## Data flow (delta in **bold**)

```
names → [pubchem] CID + IsomericSMILES  (SMILES already retrieved today)
      → **[conformers] RDKit ETKDGv3 embed (RMSD-pruned) → MMFF94 optimize/rank**
        **→ keep top 3 distinct by energy → write XYZ per conformer**
      → [gaussian] XYZ→.com per conformer (Link1 opt→freq preserved)
      → [slurm]    .com→.sh
```

Input is the `IsomericSMILES` already present in `build_molecule_table` output
(stereochemistry preserved; no Open Babel needed on this path).

## New module: `pipeline/conformers.py`

- `search_conformers(molecule_table, n_generate=20, top_n=3, rmsd_prune=0.5, seed=42)`
  → writes `conformer_log.csv` (one row **per kept conformer**) and the XYZ files.
- Small pure helper(s) — e.g. `select_top_n(energies_kcal, n)` — split out so the
  ranking logic is unit-testable without RDKit (matches the repo's testing style).
- Failures logged to `conformer_search_failed.csv` (mirrors `com_write_failed.csv`).

Ensemble generation uses RDKit `EmbedMultipleConfs` with `pruneRmsThresh`, so
duplicate removal is a built-in parameter, not separate code. MMFF94 optimizes
and scores; UFF is a logged fallback only if MMFF params are unavailable.

## Integration with `gaussian.py`

Batch writer consumes `conformer_log.csv` (multiple rows per molecule) and writes
one `.com` per conformer as `{base}_c{ii}_F.com` (e.g. `ribose_c00_F.com`),
extending `sanitize_basename`. The Link1 opt→freq checkpoint contract is
unchanged. The conformer id and its relative energy are written into the `.com`
title line for traceability.

## Units & provenance (invariant-critical, not cut for simplicity)

- MMFF energies and ΔE-from-minimum are **kcal/mol**, labeled as such, never
  mixed with DFT Hartree values downstream.
- `conformer_log.csv` records: RDKit version, random seed, method (MMFF94/UFF),
  n_generated, n_kept, rmsd_prune; and per conformer: id, ΔE (kcal/mol). A rerun
  with the same seed reproduces the ensemble (reproducibility invariant).

## Known limitations (state honestly in README; not solved in v2)

- MMFF ranking is unreliable for intramolecular-H-bonding species (sugars,
  nucleosides); a semi-empirical rerank (xTB) is a possible future addition.
- MMFF ranking is gas-phase; the DFT default is IEFPCM water, so the FF-lowest
  conformer may not be the solution-phase minimum. Carrying top 3 mitigates but
  does not eliminate this.
- Fixed `n_generate=20` may under-sample very flexible molecules.

## Out of scope for v2

Rotatable-bond gating; xTB/CREST (required or optional); energy-window logic;
Boltzmann/entropy weighting; solvent-aware search; running/parsing Gaussian;
changing the level of theory or the Link1 contract.

## Change log

- v1.1 — PubChem → Open Babel single geometry → Gaussian opt+freq
  (Zenodo 10.5281/zenodo.18894724).
- v2.0 — adds simple RDKit conformer search; top 3 distinct conformers per
  molecule to DFT; no gating.
