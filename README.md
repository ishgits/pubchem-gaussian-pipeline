[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18894724.svg)](https://doi.org/10.5281/zenodo.18894724)



# PubChem → Gaussian Input File Pipeline

Automated generation of **Gaussian input files** (`.com`) and **SLURM submission scripts** (`.sh`) from molecule names — no manual structure building required.

Give it a list of molecule names, and the pipeline resolves each to a PubChem record, runs an **RDKit conformer search**, and writes ready-to-submit Gaussian opt→freq inputs (one per conformer) plus cluster submission scripts.

> **Release-candidate status:** the manifest-centric v2 implementation is
> complete and mechanically verified on the release branch. Final merge remains
> gated on the holistic review, green CI, and Ish's approval recorded in
> `docs/implementation-status.md`.

## Pipeline Overview (v2 — the default)

```
names → PubChem              → RDKit conformer search           → Gaussian .com        → SLURM .sh
──────────────────────────    ───────────────────────────────    ─────────────────────  ─────────────
"Adenine"  CID + stereo        ETKDGv3 embed (RMSD-pruned) →      opt+freq per conformer  run-scoped,
"Glycine"  IsomericSMILES      MMFF94/UFF optimize + rank →       (Link1 pattern)         submission-
"Water"    + properties        keep top 3 distinct → XYZ each     {base}_c00_F.com …      independent
```

**Step 1 — Resolve.** Molecule names → PubChem CIDs with a scoring heuristic that picks the best candidate record (handles ambiguous names, stereochemistry, multiple CIDs) and retrieves the stereo-bearing `IsomericSMILES` + properties.

**Step 2 — Conformer search** (`pipeline/conformers.py`). For each molecule, embed an RDKit ETKDGv3 ensemble (`N_GENERATE=20`, RMSD-pruned), MMFF94-optimize and rank it (UFF is a logged fallback only when MMFF params are unavailable), and carry the **top 3 lowest-energy, distinct** conformers forward. Rigid molecules collapse to a single conformer on their own.

**Step 3 — Gaussian inputs.** Write one `.com` per conformer as `{base}_c{ii}_F.com` using the **Link1 pattern**: the optimization job writes a checkpoint, and the frequency job reads it via `Geom=AllChk Guess=Read`.

**Step 4 — SLURM scripts.** Write one unique `.sh` per valid COM. The default path is manifest/log driven, validates every source before mutation, and rejects path collisions rather than silently overwriting jobs.

**Step 5 — Run manifest.** `run_manifest.json` is the authoritative provenance record. It contains the complete conformer, Gaussian, and SLURM configuration plus molecule identity, artifact lineage, relative paths, and file hashes. XYZ, COM, and SH files carry stable identifiers that link back to the manifest.

The supported archive and transfer unit is the complete run package — in v2.1 the whole immutable `runs/<study>/<run_id>/` folder. An isolated artifact remains identifiable through its IDs, but it is not promised to be independently reproducible without the matching manifest. **There is no resume/append: every run is fresh and immutable; to add molecules, start a new run.**


## Provenance and reproducibility contract

v2.1 uses a manifest-centric model defined in
[`docs/release-contract-v2.1.md`](docs/release-contract-v2.1.md) (a deliberate
revision of the frozen [`docs/release-contract-v2.0.md`](docs/release-contract-v2.0.md)):

- `run_manifest.json` stores the complete configuration and artifact hashes;
- each run writes one immutable `runs/<study>/<run_id>/` folder, with COM+SH
  co-located in `gaussian_jobs/`;
- XYZ and COM bodies carry only inline science plus a single `artifact_id`
  back-pointer (everything else is a manifest lookup);
- SLURM scripts identify the exact source COM basename and hash and run
  `g16 <base>_F.com` from their own directory;
- every accepted source record maps to one unique destination path;
- collisions, missing or zero-byte sources, and manifest/hash disagreement fail
  before mutation;
- **resume/append is removed** — re-running against a populated run folder raises
  rather than reusing; start a new run instead;
- undefined-stereo molecules take the **provisional** path (one loudly-flagged
  arbitrated structure, `dE=NA`), never silently guessed.

The complete run package (the `runs/<study>/<run_id>/` folder), not an isolated
COM or XYZ, is the supported archival unit.

## Quick Start

### 1. Clone & install dependencies

```bash
git clone https://github.com/ishgits/pubchem-gaussian-pipeline.git
cd pubchem-gaussian-pipeline
conda env create -f environment.yml
conda activate gaussian-pipeline
# For the pinned v2.0.0 release-target stack: python -m pip install -r requirements-lock.txt
```

### 2. Edit the notebook

Open `notebooks/run_pipeline.ipynb` and:

1. **Configuration cell** — set your Gaussian method/basis set, conformer knobs (`N_GENERATE`, `TOP_N`, `RMSD_PRUNE`, `SEED`), SLURM account, nproc, memory, and walltime.
2. **Molecule cell** — replace the demo molecules with your own list.
3. Set the `STUDY` label and run all cells top to bottom. Each run creates a
   fresh, immutable `runs/<study>/<run_id>/` folder; the notebook writes
   `run_manifest.json` after PubChem resolution and before XYZ, COM, or SH
   artifacts, then verifies all recorded hashes at the end.

### 3. Submit to your cluster

```bash
# Transfer/archive the complete run package — the whole runs/<study>/<run_id>/
# folder (run_manifest.json, the stage-log CSVs, conformer_xyz/, and
# gaussian_jobs/ with COM+SH co-located). Submit each .sh FROM the directory
# that holds its .com (they ship together; the script does not hunt for input):
cd runs/<study>/<run_id>/gaussian_jobs
for f in *.sh; do sbatch "$f"; done
```

`run_manifest.json` is authoritative. The current CSV logs are operational
indexes and must agree with it. Keep the package together: isolated artifacts
contain lookup identifiers, but they require the matching manifest for the
complete configuration and supported reproducibility.

There is **no resume/append** (v2.1). Every run is fresh and immutable; to add
molecules, start a new run in a new `runs/<study>/<run_id>/` folder. Pointing the
pipeline at an already-populated run folder raises rather than reusing or
appending — two clean single-purpose runs are better provenance than a resume
path. This deletes the v2.0 stale-conformer-resurrection finding (B-04) by
construction.

## Repository Structure

```
pubchem-gaussian-pipeline/
├── pipeline/                   # Reusable Python modules
│   ├── pubchem.py              #   Name resolution, scoring, SMILES/SDF retrieval
│   ├── conformers.py           #   RDKit ETKDGv3 embed + MMFF94/UFF rank (v2 core)
│   ├── gaussian.py             #   XYZ → .com (Link1 opt+freq)
│   ├── slurm.py                #   .com → .sh (co-located; submit from the .com's dir)
│   ├── manifest.py             #   Canonical config hash, IDs, lineage, file hashes
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

Tune `N_GENERATE` (ensemble size), `TOP_N` (conformers carried to DFT), `RMSD_PRUNE` (Å, duplicate threshold), and `SEED` (fixed → reproducible ensemble) in the configuration cell.

**Undefined stereochemistry (v2.1 provisional path).** Molecules whose `IsomericSMILES` leaves stereochemistry unspecified are **no longer skipped**. RDKit embeds **one** provisional structure from the PubChem SMILES with an **arbitrary** choice at the undefined centre(s), applies a light MMFF/UFF cleanup, and records it loudly: `provenance_status=provisional_undefined_stereo`, `undefined_centers`, `pubchem_smiles`, and the post-embed `arbitrated_smiles` in the manifest and `conformer_log.csv`; `dE=NA` and a `PROVISIONAL: stereo arbitrated at …` marker in the XYZ/COM; and a loud console warning. This arbitrated structure is an **unvalidated DFT starting geometry, not the compound's real configuration** — one arbitrary pick among 2^k for k undefined centres. It must never be treated as the defined stereoisomer; stereoisomer enumeration is out of scope for v2.1.

### Level of theory

Edit `ROUTE_OPT` and `ROUTE_FREQ` in the configuration cell. The defaults are:

```
# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read
```

Change the functional, basis set, solvation model, or temperature only as an
explicit scientific choice and record it in the manifest. The supported v2.0
path preserves the Link1 opt→freq contract; decoupling the jobs is a scientific
architecture deviation and must be documented before use.

### SLURM template

The default template in `pipeline/slurm.py` assumes `module load gaussian16` and
`g16`. Edit the template string or pass a custom one for a different module
system or Gaussian version. The exact template hash and resources are part of
the immutable run configuration, so changing them requires a new manifest/run.

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

> **Deprecated compatibility path:** this workflow is explicitly exempt from the
> strict v2 manifest and artifact-linkage guarantees. Do not treat its outputs as
> satisfying the v2 provenance contract.

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

See `requirements-lock.txt` for the pinned v2.0.0 release-target stack and its recorded verification status.

## Citation

If you use this workflow or any part of the dataset in your research, please cite:

Ishaan Madan. (2026). ishgits/pubchem-gaussian-pipeline: v1.1 – Initial Release (v1.1). Zenodo. https://doi.org/10.5281/zenodo.18894724

## License

MIT — see [LICENSE](LICENSE).
