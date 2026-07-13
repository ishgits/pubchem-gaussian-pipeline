# v2.0 release contract

> **Status:** frozen and approved by Ish on 2026-07-13. This document defines
> the exact provenance, artifact, reuse, collision, scientific-judgment, and
> review boundaries for the v2.0.0 release candidate. Reviewers must test the
> implementation against this contract; they must not silently broaden it.

Normative terms such as **MUST**, **MUST NOT**, **SHOULD**, and **MAY** are used
in their ordinary requirements sense.

## 1. Provenance model

v2.0 uses a **manifest-centric** provenance model.

- `run_manifest.json` is the authoritative, complete record for one pipeline
  execution.
- Stage CSVs remain useful operational indexes, but they are subordinate to the
  manifest and MUST agree with it.
- Individual XYZ and COM files are scientific artifacts. They MUST carry enough
  stable identity to locate their exact manifest record, but they do **not**
  need to duplicate every search knob.
- SLURM scripts are operational artifacts. Their resource directives remain
  visible in the script, and they MUST identify the exact source COM record.

Absence of duplicated full search configuration inside every XYZ or COM is not a
contract violation when the required manifest linkage below is present and
valid.

## 2. Supported archive and transfer unit

The supported archive/transfer unit is the **complete run package**:

```text
run_manifest.json
conformer_log.csv
com_write_log.csv
slurm_write_log.csv
conformer_xyz/
gaussian_inputs/
slurm_scripts/
```

An isolated XYZ, COM, or SH remains attributable through its identifiers, but it
is not promised to be independently reproducible without the matching manifest.
Documentation MUST tell users to transfer and archive the package together.

The `runs/<study>/<run_id>/` directory redesign is deferred to v2.1. v2.0 may
retain flat output directories, provided the current package is unambiguous and
all destination collisions fail before mutation.

## 3. Stable identifiers and hashes

Every v2 run MUST define:

- `manifest_schema`: the manifest schema version;
- `run_id`: an opaque, immutable, collision-resistant identifier generated once
  for the execution;
- `config_hash`: SHA-256 over a canonical JSON representation of the complete
  scientific and operational configuration, excluding timestamps, absolute
  machine-specific paths, and output file hashes;
- `artifact_id`: an immutable, unique identifier for each XYZ, COM, and SH
  record within the run.

The canonical configuration used for `config_hash` MUST include:

- requested molecule labels and resolved CID/stereo SMILES identity;
- pipeline version and clean source commit;
- RDKit version;
- conformer method policy, seed, `n_generate`, `top_n`, and `rmsd_prune`;
- Gaussian routes, charge, multiplicity, nproc, and Link1 choice;
- SLURM account/resources and template identity.

The manifest MUST store SHA-256 file hashes for generated XYZ, COM, and SH
artifacts after writing them. File hashes are not included in `config_hash`.

## 4. Authoritative manifest contents

`run_manifest.json` MUST record at least:

### Run-level identity and configuration

```text
manifest_schema
run_id
config_hash
pipeline_version
pipeline_commit
rdkit_version
conformer configuration
Gaussian configuration
SLURM configuration
```

### Molecule and conformer records

```text
molecule_name
CID
IsomericSMILES
molecule identity hash
conformer_id
method
seed
n_generate
n_generated
top_n
n_kept
rmsd_prune
relative_energy_kcalmol
converged
```

### Artifact records

```text
artifact_id
artifact kind
parent artifact_id where applicable
relative path
SHA-256 file hash
```

The manifest MUST preserve the lineage:

```text
molecule identity -> conformer record -> XYZ -> COM -> SLURM script
```

## 5. Per-artifact metadata contract

### Conformer XYZ scientific artifact

The XYZ comment MUST carry:

```text
run_id
artifact_id
config_hash
conformer_id
relative_energy_kcalmol
method
pipeline_version
rdkit_version
```

The full ensemble and selection configuration belongs in the manifest.

### Gaussian COM scientific artifact

The Gaussian title section MUST carry:

```text
run_id
artifact_id
config_hash
conformer_id
relative_energy_kcalmol
pipeline_version
rdkit_version
```

The COM already directly records routes, coordinates, charge/multiplicity,
processor count, checkpoint directives, and Link1 behavior. The full conformer
search configuration belongs in the manifest.

### SLURM SH operational artifact

The script header MUST carry:

```text
run_id
artifact_id
source COM relative path
source COM SHA-256
```

SBATCH resources and execution commands are already visible in the script.

## 6. Output mapping and collision behavior

v2.0 MUST fail loudly before mutation when any of the following occurs:

- distinct molecule labels sanitize to one output basename;
- duplicate source artifact paths are supplied;
- two source artifacts map to one destination path;
- two COM paths map to one SLURM script path;
- a required source artifact is missing, blank, or zero bytes;
- a manifest or stage log contains duplicate `artifact_id` values;
- an artifact path or hash disagrees with the manifest.

v2.0 MUST NOT auto-disambiguate filenames. Collision-proof internal filenames
are deferred to the v2.1 run-directory redesign.

## 7. Resume and append policy

Resume or append reuse is allowed only when all existing integrity, identity,
configuration, file-existence, and complete-group checks pass **and** both the
recorded and current pipeline commits are the same clean, nonblank commit.

- A dirty tree MUST disable reuse.
- A missing git commit MUST disable reuse.
- A source ZIP or installed package without git metadata MUST regenerate rather
  than use a version-only fallback.
- The manifest and artifact hashes MUST agree before reuse.

This strict policy favors honest regeneration over unverifiable cache reuse.

## 8. Legacy v1.1 boundary

The Open Babel single-geometry pathway is deprecated compatibility functionality.
It is explicitly exempt from the strict v2 manifest and per-artifact metadata
contract unless a future release deliberately upgrades it.

Documentation MUST label it legacy and MUST NOT imply that its outputs satisfy
the v2 provenance guarantees.

## 9. Frozen scientific judgments

The following decisions are approved for v2.0:

1. PubChem's stereo-bearing `SMILES` is the source of record. Molecules with
   undefined stereochemistry are skipped and logged rather than allowing RDKit
   to guess.
2. Distinctness relies on ETKDG embed-time RMSD pruning. Post-optimization RMSD
   re-pruning is deferred.
3. UFF is permitted only as a clearly recorded fallback. Relative force-field
   energies are interpreted only within one molecule and one method, never
   across molecules or between MMFF94 and UFF.
4. RDKit builds the structure represented by the supplied SMILES/protonation
   state. Gaussian charge and multiplicity are supplied explicitly; the
   conformer stage does not infer electronic multiplicity.

## 10. Review boundary and stop rule

A release-blocking finding must demonstrate one of the following:

- wrong chemistry or violation of a scientific invariant;
- silent data loss, overwrite, or one-to-many/many-to-one mapping corruption;
- molecular identity mismatch;
- violation of this exact provenance contract;
- a supported workflow producing an invalid or non-runnable job;
- regression in required tests or mechanical invariants.

A proposed improvement that expands this frozen contract is recorded for v2.1
unless Ish explicitly approves the expansion. It is not automatically a v2.0
Blocker or Major.

After implementation aligns with this contract:

1. run one holistic base-to-head audit;
2. fix actual contract-level Blocker and Major findings;
3. run one final re-review;
4. defer non-contract recommendations to v2.1;
5. merge when CI is green, no contract-level Blocker or Major remains, and Ish
   gives final approval.
