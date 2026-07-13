[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18894724.svg)](https://doi.org/10.5281/zenodo.18894724)



# PubChem → Gaussian Input File Pipeline

Automated generation of **Gaussian input files** (`.com`) and **SLURM submission scripts** (`.sh`) from molecule names — no manual structure building required.

Give it a list of molecule names, and the pipeline resolves each to a PubChem record, runs an **RDKit conformer search**, and writes ready-to-submit Gaussian opt→freq inputs (one per conformer) plus cluster submission scripts.

## Pipeline Overview (v2 — the default)

```
names → PubChem              → RDKit conformer search           → Gaussian .com        → SLURM .sh
──────────────────────────    ───────────────────────────────    ─────────────────────  ─────────────
"Adenine"  CID + stereo        ETKDGv3 embed (RMSD-pruned) →      opt+freq per conformer  run-scoped,
"Glycine"  IsomericSMILES      MMFF94/UFF optimize + rank →       (Link1 pattern)         submission-
"Water"    + properties        keep top 3 distinct → XYZ each     {base}_c00_F.com …      independent
```

**Step 1 — Resolve.** Molecule names → PubChem CIDs with a scoring heuristic that picks the best candidate record (handles ambiguous names, stereochemistry, multiple CIDs) and retrieves the stereo-bearing `IsomericSMILES` + properties.

**Step 2 — Conformer search** (`pipeline/conformers.py`). For each molecule, embed an RDKit ETKDGv3 ensemble (`N_GENERATE=20`, RMSD-pruned), MMFF94-optimize and rank it (UFF is a logged fallback only when MMFF params are unavailable), and carry the **top 3 lowest-energy, distinct** conformers forward — one XYZ each, with provenance (RDKit version, seed, method, ΔE in kcal/mol) in `conformer_log.csv`. Rigid molecules collapse to a single conformer on their own.

**Step 3 — Gaussian inputs.** Write one `.com` per conformer as `{base}_c{ii}_F.com` using the **Link1 pattern**: the optimization job writes a checkpoint, the frequency job reads it via `Geom=AllChk Guess=Read`. The conformer id and its ΔE (kcal/mol) go in the title line.

**Step 4 — SLURM scripts.** One `.sh` per conformer. Scripts default to the current run's `com_write_log.csv` (stale `.com` files on disk are never picked up) and resolve their input **relative to their own location**, so `sbatch slurm_scripts/*.sh` works from any directory.

Every step writes a log CSV, so any output traces back to its input. Reruns are resume-safe: a molecule is regenerated (never silently reused) if its search knobs, pipeline/RDKit version, or CID/SMILES identity changed, or if its recorded XYZ is missing.

## Quick Start

### 1. Clone & install dependencies

```bash
git clone https://github.com/ishgits/pubchem-gaussian-pipeline.git
cd pubchem-gaussian-pipeline
conda env create -f environment.yml
conda activate gaussian-pipeline
# For the exact tested v2.0.0 stack: python -m pip install -r requirements-lock.txt
```

### 2. Edit the notebook

Open `notebooks/run_pipeline.ipynb` and:

1. **Configuration cell** — set your Gaussian method/basis set, conformer knobs (`N_GENERATE`, `TOP_N`, `RMSD_PRUNE`, `SEED`), SLURM account, nproc, memory, and walltime.
2. **Molecule cell** — replace the demo molecules with your own list.
3. Run all cells top to bottom.

### 3. Submit to your cluster

```bash
# Transfer the gaussian_inputs/ and slurm_scripts/ directories to your HPC system, then:
for f in slurm_scripts/*.sh; do sbatch "$f"; done
```

## Repository Structure

```
pubchem-gaussian-pipeline/
├── pipeline/                   # Reusable Python modules
│   ├── pubchem.py              #   Name resolution, scoring, SMILES/SDF retrieval
│   ├── conformers.py           #   RDKit ETKDGv3 embed + MMFF94/UFF rank (v2 core)
│   ├── gaussian.py             #   XYZ → .com (Link1 opt+freq)
│   ├── slurm.py                #   .com → .sh (run-scoped, submission-independent)
│   ├── geometry.py             #   Legacy v1.1 SDF → XYZ (Open Babel)
│   └── utils.py                #   Shared helpers, provenance
├── notebooks/
│   └── run_pipeline.ipynb      # Main walkthrough notebook (conformer path default)
├── examples/
│   └── nucleobases_nucleosides.py  # Real-world molecule list example
├── tests/                      # Unit tests (pytest, offline)
├── docs/                       # architecture.md, implementation-plan.md, status
├── environment.yml             # Conda environment
├── requirements-lock.txt       # Pinned v2.0.0 stack
├── LICENSE
└── README.md
```

## Customizing for Your Project

### Molecules

The notebook's molecule cell accepts any name PubChem can resolve. For ambiguous names, use the **alias** dictionary to map your label to a specific PubChem query; for tricky names, add **fallback queries** (tried in order). See `examples/nucleobases_nucleosides.py`.

### Conformer search

Tune `N_GENERATE` (ensemble size), `TOP_N` (conformers carried to DFT), `RMSD_PRUNE` (Å, duplicate threshold), and `SEED` (fixed → reproducible ensemble) in the configuration cell. Molecules whose `IsomericSMILES` has **undefined stereochemistry** are **skipped and logged** to `conformer_search_failed.csv` (reason `"undefined stereochemistry"`) rather than letting RDKit guess a stereoisomer.

### Level of theory

Edit `ROUTE_OPT` and `ROUTE_FREQ` in the configuration cell. The defaults are:

```
# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read
```

Change the functional, basis set, solvation model, or temperature to match your needs. To run opt and freq as separate jobs, set `link1=False` in the Gaussian step.

### SLURM template

The default template in `pipeline/slurm.py` assumes `module load gaussian16` and `g16`. Edit the template string or pass a custom one for a different module system or Gaussian version. `.sh` files are overwritten on rerun, so changing your account/resources always re-stamps the scripts.

### Charge & multiplicity

Pass `charge=` and `multiplicity=` to the Gaussian step. Defaults are `0 1` (neutral singlet). For per-molecule charges, define a mapping and loop manually — the batch function assumes uniform charge/multiplicity.

## Important Caveats

Conformer energies are **MMFF94/UFF in kcal/mol** and are never mixed with the DFT Hartree energies computed downstream. These are force-field **starting** geometries, not optimized minima — the DFT step makes the final call among the carried conformers. Three limitations are retained honestly and **not** solved in v2:

- **MMFF ranking is unreliable for intramolecular-H-bonding species** (sugars, nucleosides): the force field over-stabilizes internal H-bonds, so the FF ranking may not match the DFT ranking. Carrying the top 3 mitigates this; a semi-empirical rerank (e.g. xTB) is possible future work.
- **The force field is gas-phase; the DFT default is IEFPCM water.** The FF-lowest conformer may not be the solution-phase minimum.
- **Sampling is fixed at `N_GENERATE=20`**, which may under-sample very flexible molecules.

When no conformer converges under the force field even after a retry, exactly one lowest-energy best-effort geometry is carried, flagged `converged=False`, warned at runtime, and tagged `UNCONVERGED_FF_SEED` in the `.com` title — its FF energy is labeled unreliable.

### PubChem rate limits

The pipeline includes retry logic with exponential backoff and respects PubChem's rate limits (`~5 requests/second`). The on-disk cache (`.pubchem_cache/`) prevents redundant API calls on reruns.

## Legacy v1.1 workflow (Open Babel single geometry)

Before v2 the pipeline fed DFT a **single** Open Babel geometry per molecule (PubChem 3D SDF → Open Babel `--gen3d --minimize` → one `.com`). This path is **superseded** by the conformer search above but remains available via `pipeline/geometry.py` (`download_sdfs` / `convert_sdfs_to_xyz` / `write_gaussian_coms`) and the commented-out appendix at the bottom of `notebooks/run_pipeline.ipynb`. It requires **Open Babel** and hands DFT only one starting geometry, so for flexible molecules it may miss the global-minimum conformation. Use it only if you deliberately want the single-geometry behavior.

## Running Tests

```bash
cd pubchem-gaussian-pipeline
pytest tests/ -q
```

Tests cover the pure functions and the offline conformer/Gaussian/SLURM paths; they require no network access, no cluster, and no Gaussian. RDKit-dependent tests `importorskip` so a bare environment still runs the pure tests.

## Requirements

- **Python** ≥ 3.11
- **pandas**, **requests**
- **RDKit** — the v2 conformer search (`conda install -c conda-forge rdkit`)
- **Jupyter** (to run the notebook)
- **Open Babel** — only for the optional legacy v1.1 single-geometry path
- **Gaussian 16** (on your HPC cluster, not needed locally)

See `requirements-lock.txt` for the exact pinned versions the v2.0.0 suite was tested against.

## Citation

If you use this workflow or any part of the dataset in your research, please cite:

Ishaan Madan. (2026). ishgits/pubchem-gaussian-pipeline: v1.1 – Initial Release (v1.1). Zenodo. https://doi.org/10.5281/zenodo.18894724

## License

MIT — see [LICENSE](LICENSE).
