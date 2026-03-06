[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18894724.svg)](https://doi.org/10.5281/zenodo.18894724)



# PubChem → Gaussian Input File Pipeline

Automated generation of **Gaussian input files** (`.com`) and **SLURM submission scripts** (`.sh`) from molecule names — no manual structure building required.

Give it a list of molecule names, and the pipeline handles PubChem lookup, 3D structure retrieval, coordinate conversion, and input file formatting so you can go straight to submitting quantum chemistry jobs on your cluster.

## Pipeline Overview

```
Molecule names          PubChem PUG-REST          Open Babel          Gaussian .com + SLURM .sh
─────────────────  ──►  ──────────────────  ──►  ──────────  ──►  ──────────────────────────────
"Adenine"               CID → 3D SDF              SDF → XYZ         opt+freq input files
"Glycine"               (with scoring &            (--gen3d           (Link1 pattern)
"Water"                  fallback queries)          --minimize)       + SLURM scripts
```

**Step 1** — Resolve molecule names to PubChem CIDs with a scoring heuristic that picks the best candidate record (handles ambiguous names, stereochemistry, multiple CIDs).

**Step 2** — Download 3D SDF files from PubChem (falls back to 2D if 3D unavailable).

**Step 3** — Convert SDF → XYZ using Open Babel with forced 3D generation and force-field minimization.

**Step 4** — Write Gaussian `.com` input files with the **Link1 pattern**: optimization job writes a checkpoint, frequency job reads it via `Geom=AllChk Guess=Read`.

**Step 5** — Generate one SLURM `.sh` submission script per `.com` file.

Every step writes a log CSV. If the pipeline is interrupted, rerunning skips completed work automatically.

## Quick Start

### 1. Clone & install dependencies

```bash
git clone https://github.com/YOUR_USERNAME/pubchem-gaussian-pipeline.git
cd pubchem-gaussian-pipeline
conda env create -f environment.yml
conda activate gaussian-pipeline
```

### 2. Edit the notebook

Open `notebooks/run_pipeline.ipynb` and:

1. **Configuration cell** — set your Gaussian method/basis set, SLURM account, nproc, memory, and walltime.
2. **Molecule cell** — replace the demo molecules with your own list.
3. Run all cells top to bottom.

### 3. Submit to your cluster

```bash
# Transfer outputs to your HPC system, then:
for f in slurm_scripts/*.sh; do sbatch "$f"; done
```

## Repository Structure

```
pubchem-gaussian-pipeline/
├── pipeline/                   # Reusable Python modules
│   ├── pubchem.py              #   Name resolution, scoring, SDF download
│   ├── geometry.py             #   SDF → XYZ (Open Babel)
│   ├── gaussian.py             #   XYZ → .com (Link1 opt+freq)
│   ├── slurm.py                #   .com → .sh
│   └── utils.py                #   Shared helpers
├── notebooks/
│   └── run_pipeline.ipynb      # Main walkthrough notebook
├── examples/
│   └── nucleobases_nucleosides.py  # Real-world molecule list example
├── tests/                      # Unit tests (pytest)
├── environment.yml             # Conda environment
├── LICENSE
└── README.md
```

## Customizing for Your Project

### Molecules

The notebook's molecule cell accepts any name PubChem can resolve. For molecules where the common name is ambiguous, use the **alias** dictionary to map your label to a specific PubChem query. For really tricky names, add **fallback queries** — a list of alternative names tried in order. See `examples/nucleobases_nucleosides.py` for a detailed example.

### Level of theory

Edit `ROUTE_OPT` and `ROUTE_FREQ` in the configuration cell. The defaults are:

```
# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read
```

Change the functional, basis set, solvation model, or temperature to match your needs. If you don't need the Link1 opt→freq pattern (e.g., you want to run them as separate jobs), set `link1=False` in the Gaussian step.

### SLURM template

The default template in `pipeline/slurm.py` assumes `module load gaussian16` and `g16`. Edit the template string or pass a custom one if your cluster uses a different module system, Gaussian version, or queue structure.

### Charge & multiplicity

Pass `charge=` and `multiplicity=` to the Gaussian writing step. The defaults are `0` and `1` (neutral singlet). If different molecules in your list need different charges, you'll want to define a mapping and loop manually — the batch function assumes uniform charge/multiplicity.

## Important Caveats

### Starting geometries are approximate

PubChem 3D structures and Open Babel's `--gen3d --minimize` output provide reasonable starting points, but they are **not** DFT-optimized. Gaussian's geometry optimization will refine them, but:

- If optimization doesn't converge, the starting geometry may be in a bad region of the potential energy surface. Inspect the `.xyz` file and adjust manually if needed.
- For flexible molecules with many rotatable bonds, the PubChem geometry may not correspond to the **global minimum** conformation.

### Conformer searching is not included

This pipeline does **not** perform a systematic conformer search (e.g., via RDKit ETKDG + MMFF/UFF). For small, rigid molecules (most nucleobases, amino acids, small heterocycles), the PubChem 3D geometry is typically adequate as a starting point for DFT optimization. For larger, flexible molecules (long-chain lipids, peptides, sugars with multiple ring conformations), you should consider running a conformer search separately before feeding geometries into this pipeline. Tools like [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules), [CREST](https://crest-lab.github.io/crest-docs/), or [Confab](https://open-babel.readthedocs.io/en/latest/3DStructureGen/multipleconformers.html) can help identify the lowest-energy conformer prior to DFT.

### PubChem rate limits

The pipeline includes retry logic with exponential backoff and respects PubChem's rate limits (`~5 requests/second`). If you're resolving hundreds of molecules, the caching system (`.pubchem_cache/`) prevents redundant API calls on reruns.

## Running Tests

```bash
cd pubchem-gaussian-pipeline
pytest tests/ -v
```

Tests cover the pure functions (coordinate formatting, file generation, scoring heuristic) and require no network access or external software.

## Requirements

- **Python** ≥ 3.10
- **pandas**, **requests**, **numpy**
- **Open Babel** (for SDF → XYZ conversion)
- **Jupyter** (to run the notebook)
- **Gaussian 16** (on your HPC cluster, not needed locally)

## Citation

If you use this workflow or any part of the dataset in your research, please cite:

Ishaan Madan. (2026). ishgits/pubchem-gaussian-pipeline: v1.1 – Initial Release (v1.1). Zenodo. https://doi.org/10.5281/zenodo.18894724

## License

MIT — see [LICENSE](LICENSE).
