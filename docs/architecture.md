# architecture.md

> Canonical architecture (v2.0). Approved by Ish before implementation. The
> earlier v1.1 architecture and the v2 draft/review history are archived under
> `docs/review-history/v2/`.

## What the system is

An automated pipeline that turns a list of molecule names into ready-to-submit
Gaussian quantum-chemistry jobs, removing manual structure building. It produces
*inputs* to Gaussian (`.com`) plus SLURM submission scripts (`.sh`); it does not
run Gaussian.

## Why it exists

Building Gaussian inputs by hand (name → structure → conformers → coordinates →
formatted `.com` → SLURM script) is slow and error-prone at the scale of
prebiotic-chemistry surveys. The pipeline makes input generation reproducible and
traceable, and — new in v2 — hands DFT a *ranked ensemble* of conformers rather
than a single force-field guess.

## What v2 adds

A conformer-search stage between molecule resolution and Gaussian input
generation. For each molecule the pipeline embeds an RDKit conformer ensemble,
MMFF94-optimizes and ranks it, and carries the **top 3 lowest-energy, distinct**
conformers forward as DFT starting geometries — each with recorded provenance
(RDKit version, seed, method, relative energy).

**Why:** v1.1 handed DFT a single Open Babel conformer. For flexible molecules
that conformer is often not the global minimum, so Gaussian can optimize into a
local minimum and yield wrong relative energies — corrupting the comparative
thermodynamics this pipeline feeds. Carrying the top 3 lets the DFT step, not the
force field, make the final call.

## Deliberately simple: no gating

There is **no rotatable-bond branching**. Every molecule takes the same path.
RMSD pruning at embed time plus "keep top 3 distinct" makes rigid molecules
collapse to a single conformer on their own (adenine's conformers are
near-identical → pruned to 1 → one DFT job), so v1.1 behavior for rigid molecules
is preserved emergently, with no special-case code.

## Data flow

```
names → [pubchem]     CID + IsomericSMILES (stereo-bearing) + properties
      → [conformers]  RDKit ETKDGv3 embed (RMSD-pruned) → MMFF94 optimize/rank
                      → keep top 3 distinct by energy → one XYZ per conformer
      → [gaussian]    XYZ → .com per conformer (Link1 opt→freq, checkpoint-linked)
      → [slurm]       .com → .sh per conformer (run-scoped, submission-independent)
```

Input to the conformer stage is the `IsomericSMILES` already present in
`build_molecule_table` output (stereochemistry preserved; no Open Babel needed on
this path). Every stage writes a log CSV so any output traces back to its input
(molecule name → CID → SMILES → conformer XYZ → `.com` → `.sh`).

## Modules

- `pipeline/pubchem.py` — name→CID scoring/fallback, property + stereo-SMILES
  retrieval, SDF download, caching, rate limiting.
- `pipeline/conformers.py` — RDKit ETKDGv3 embed + MMFF94/UFF rank; `select_top_n`
  is a pure, RDKit-free ranking helper; `search_conformers` is the batch driver
  that writes `conformer_log.csv` (one row per kept conformer) and the XYZ files.
- `pipeline/gaussian.py` — XYZ→.com; owns the route lines and the Link1 contract;
  `write_gaussian_coms_from_conformers` writes one `.com` per conformer row.
- `pipeline/slurm.py` — .com→.sh from an editable cluster template; scripts
  resolve their `.com` relative to their own location so submission works from any
  directory.
- `pipeline/geometry.py` — legacy v1.1 SDF→XYZ via Open Babel (single-geometry
  path, retained but not the default).
- `pipeline/utils.py` — sanitization, dir handling, CID normalization, provenance.

## Reproducibility & resume

- A rerun with the same seed reproduces the ensemble. `conformer_log.csv` records
  RDKit version, seed, method (MMFF94/UFF), `n_generated`, `n_kept`, `rmsd_prune`,
  `pipeline_version`, `pipeline_commit`, and per conformer the id + ΔE (kcal/mol).
- Resume is identity- and config-aware: a molecule is skipped only when its
  recorded run-level config (seed, `n_generate`, `top_n`, `rmsd_prune`, pipeline
  and RDKit version) **and** per-molecule identity (CID, SMILES) match this run,
  and every recorded XYZ still exists. Any drift regenerates the molecule rather
  than reusing a stale geometry. By default the log holds exactly the molecules
  requested this run (`append=True` retains carry-forward).
- SLURM scripts default to the current run's `com_write_log.csv`, so stale `.com`
  files on disk are never turned into jobs; `.sh` files are overwritten so a
  rerun cannot leave stale SBATCH directives behind.

## Scientific assumptions (invariants live in AGENTS.md §2)

- Default level of theory: B3LYP/6-311++G(2df,2p), IEFPCM water, 298 K.
- Conformer energies are MMFF94 (or UFF fallback) in **kcal/mol**, never mixed
  with DFT Hartree values.
- RDKit/MMFF geometries are DFT *starting* points, not minima; the DFT step makes
  the final call among the carried conformers.

## Out of scope

Rotatable-bond gating; xTB/CREST; energy-window logic; Boltzmann/entropy
weighting; solvent-aware search; running or parsing Gaussian; changing the level
of theory or the Link1 contract.

## Change log

- v1.1 — PubChem → Open Babel single geometry → Gaussian opt+freq
  (Zenodo 10.5281/zenodo.18894724).
- v2.0 — adds RDKit conformer search (top 3 distinct conformers per molecule to
  DFT, no gating); run-scoped, submission-independent SLURM scripts;
  identity/config-aware resume. `pipeline.__version__ = "2.0.0"`.
