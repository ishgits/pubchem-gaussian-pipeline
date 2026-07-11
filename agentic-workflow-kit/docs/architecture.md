# architecture.md

> Human gate: this document must be approved by Ish before implementation.

## What the system is

An automated pipeline that turns a list of molecule names into ready-to-submit
Gaussian quantum-chemistry jobs, removing manual structure building.

## Why it exists

Building Gaussian inputs by hand (name → structure → coordinates → formatted
`.com` → SLURM script) is slow and error-prone at the scale of prebiotic
chemistry surveys. The pipeline makes input generation reproducible and traceable.

## Data flow

```
names → [pubchem] CID+props → SDF (3D, 2D fallback)
      → [geometry] Open Babel SDF→XYZ (--gen3d --minimize)
      → [gaussian] XYZ→.com  (Link1 opt→freq, checkpoint-linked)
      → [slurm]    .com→.sh   (per-job submission scripts)
```

Every stage writes a log CSV so any output traces back to its input.

## Modules

- `pipeline/pubchem.py` — name→CID scoring/fallback, SDF download, caching, rate limiting.
- `pipeline/geometry.py` — SDF→XYZ via Open Babel.
- `pipeline/gaussian.py` — XYZ→.com; owns the route lines and the Link1 contract.
- `pipeline/slurm.py` — .com→.sh from an editable cluster template.
- `pipeline/utils.py` — sanitization, dir handling, CID normalization.

## Scientific assumptions (invariants live in AGENTS.md §2)

- Default level of theory: B3LYP/6-311++G(2df,2p), IEFPCM water, 298 K.
- PubChem+Open Babel geometries are DFT *starting* points, not minima.
- No conformer search is performed; adequate only for small rigid molecules.

## Out of scope

Running Gaussian; parsing Gaussian output; conformer searching; energy analysis.

## Change log

- v1.1 — initial released architecture (Zenodo 10.5281/zenodo.18894724).
- v2.x — <describe the change this workflow round introduces>
