# architecture.md

> Canonical v2.0 architecture, revised and frozen by Ish on 2026-07-13. The
> exact normative release boundary is `docs/release-contract-v2.0.md`. Earlier
> drafts and remediation history remain under `docs/review-history/v2/`.

## System purpose

The pipeline converts molecule names into ranked conformer starting geometries,
Gaussian Link1 opt/freq inputs, and SLURM scripts. It does not run Gaussian.

The primary v2 flow is:

```text
names
  -> PubChem CID + stereo-bearing SMILES
  -> RDKit ETKDGv3 ensemble
  -> MMFF94 optimization/ranking, UFF recorded fallback
  -> top-N carried conformers
  -> conformer XYZ artifacts
  -> Gaussian COM artifacts
  -> SLURM SH artifacts
```

Rigid molecules naturally collapse to fewer conformers through embed-time RMSD
pruning. There is no rotatable-bond gate.

## Provenance architecture

v2.0 uses an authoritative **manifest-centric** model.

```text
                         run_manifest.json
                         /       |        \
                conformer rows  COM rows  SLURM rows
                    |              |          |
                   XYZ ----------> COM ------> SH
```

`run_manifest.json` records the complete run configuration, resolved molecular
identity, conformer generation and selection data, Gaussian settings, SLURM
settings, artifact lineage, relative paths, and SHA-256 file hashes.

Every v2 artifact carries stable linkage:

- `run_id` identifies the execution;
- `config_hash` identifies the canonical complete configuration;
- `artifact_id` identifies one artifact record.

Stage CSVs remain operational indexes and must agree with the manifest. They are
not the sole provenance authority.

## Artifact contract

### XYZ and COM

XYZ and COM are scientific artifacts. They carry the minimal self-identifying
metadata defined in `docs/release-contract-v2.0.md`, including manifest linkage,
conformer identity, and source software identity. They do not duplicate every
conformer-search knob; the complete configuration is in the manifest.

### SLURM scripts

SH files are operational artifacts. Each script identifies its source COM path
and hash plus run/artifact identity. SBATCH resources and the execution command
remain visible in the script itself.

### Supported package

The supported transfer/archive unit is the complete current run package:
manifest, stage logs, XYZ directory, COM directory, and SLURM directory. An
isolated artifact is attributable but not promised to be independently
reproducible without the manifest.

## Modules

- `pipeline/pubchem.py` — name resolution, CID selection, stereo-SMILES and
  property retrieval, caching, and legacy SDF download.
- `pipeline/conformers.py` — eligibility, ETKDGv3 generation, MMFF94/UFF
  optimization, ranking, convergence handling, selection, XYZ writing, and
  conformer-stage records.
- `pipeline/gaussian.py` — physical-line XYZ parsing, Gaussian COM generation,
  routes, charge/multiplicity, checkpoints, and Link1 behavior.
- `pipeline/slurm.py` — validated COM-to-SH mapping and cluster template output.
- `pipeline/utils.py` — shared normalization, hashing, provenance, identifiers,
  and filesystem helpers.
- a manifest module or equivalent shared implementation — canonical
  configuration serialization, IDs, hashes, lineage, and final manifest write.

## Identity and one-to-one mapping

Each stage must preserve one-to-one lineage. Every accepted source record maps to
one unique destination record and path. Distinct inputs must never collapse to a
single filename or script.

v2.0 fails before mutation on:

- sanitized molecule-label collisions;
- duplicate source paths or artifact IDs;
- duplicate destination paths;
- two COMs mapping to one SH path;
- blank, missing, or zero-byte required inputs;
- manifest/path/hash disagreement.

Automatic filename disambiguation is deferred to v2.1.

## Reproducibility, resume, and append

A conformer group may be reused only when:

- molecular identity matches;
- search configuration matches;
- the group is complete and internally consistent;
- referenced files exist, are nonempty, and match manifest hashes;
- pipeline and RDKit versions match;
- the recorded and current pipeline commits are the same clean, nonblank commit.

Dirty trees, source archives without git metadata, and installed packages without
commit identity disable reuse and force regeneration. Version-only fallback is
not permitted under the frozen v2.0 contract.

Append mode applies the same integrity and provenance checks to retained groups.
An invalid unrequested group aborts before mutation because it cannot be safely
regenerated without a current molecule-table identity.

## Scientific choices

- PubChem's stereo-bearing SMILES is authoritative; undefined stereo is skipped.
- ETKDG embed-time RMSD pruning defines distinctness in v2.0; no post-MMFF
  re-pruning is performed.
- MMFF94 is preferred; UFF is a recorded fallback. Force-field ΔE comparisons are
  within one molecule and one method only.
- RDKit geometries are starting points, not minima.
- Charge and multiplicity are explicit Gaussian inputs; conformer generation
  does not infer multiplicity.
- Default Gaussian chemistry and Link1 semantics remain unchanged.

## Directory scope

The `runs/<study>/<run_id>/` redesign is deferred to v2.1. v2.0 retains flat
output directories but treats the current manifest and its referenced artifact
set as authoritative. Physical stale-file cleanup is not required if stale files
cannot enter the manifest-driven submission path.

## Legacy boundary

The Open Babel v1.1 single-geometry path is deprecated and explicitly exempt
from the strict v2 manifest contract. It must be labeled legacy and must not be
presented as having v2 provenance guarantees.

## Out of scope for v2.0

Rotatable-bond gating, xTB/CREST, energy windows, Boltzmann or entropy weighting,
solvent-aware conformer search, post-optimization RMSD re-pruning, running or
parsing Gaussian, auto-disambiguated internal filenames, and per-study run
directories.
