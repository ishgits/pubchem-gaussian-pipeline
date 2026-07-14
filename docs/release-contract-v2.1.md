# v2.1 release contract

> **Status:** approved by Ish as a deliberate revision of the frozen
> `release-contract-v2.0.md`. This document defines the exact provenance,
> artifact, collision, scientific-judgment, and review boundaries for the v2.1
> release candidate. It is derived from the v2.0 contract by applying the
> "Contract deltas" section of `docs/architecture-v2.1.md`. Reviewers must test
> the implementation against this contract; they must not silently broaden it.
>
> **What changed from v2.0.** Per-artifact metadata is reduced to one stable
> back-pointer plus inline-useful science (§5); flat output directories are
> replaced by one immutable `runs/<study>/<run_id>/` folder with COM+SH
> co-located in `gaussian_jobs/` (§2); resume/append reuse is removed entirely
> (§7); and undefined-stereo molecules now produce one loudly-flagged
> provisional starting structure instead of being skipped (§9). Level of theory,
> route lines, charge/multiplicity, and the Link1 opt→freq checkpoint contract
> are unchanged.

Normative terms such as **MUST**, **MUST NOT**, **SHOULD**, and **MAY** are used
in their ordinary requirements sense.

## 1. Provenance model

v2.1 uses a **manifest-centric** provenance model, unchanged from v2.0.

- `run_manifest.json` is the authoritative, complete record for one pipeline
  execution.
- Stage CSVs remain useful operational indexes, but they are subordinate to the
  manifest and MUST agree with it.
- Individual XYZ and COM files are scientific artifacts. They MUST carry enough
  stable identity to locate their exact manifest record, but they do **not**
  duplicate the search configuration. In v2.1 that stable identity is a single
  `artifact_id` (see §5).
- SLURM scripts are operational artifacts. Their resource directives remain
  visible in the script, and they MUST identify the exact source COM record.

Absence of duplicated full search configuration inside every XYZ or COM is not a
contract violation when the required manifest linkage below is present and
valid. In v2.1 the per-artifact matrix is deliberately reduced to the minimum in
§5; a reviewer MUST NOT report the absence of the removed fields as a finding.

## 2. Supported archive and transfer unit

The supported archive/transfer unit is the **complete run package**, which in
v2.1 is literally one immutable per-run directory:

```text
runs/<study>/<run_id>/
  run_manifest.json
  conformer_log.csv
  com_write_log.csv
  slurm_write_log.csv
  conformer_search_failed.csv        # only if failures
  conformer_xyz/
    <base>_c00.xyz …
  gaussian_jobs/                      # COM + SH co-located
    <base>_c00_F.com
    <base>_c00_F.sh
```

- `<study>` is a user-supplied label in the config cell (e.g. `nucleobases`);
  the `run_id` remains the opaque, immutable UUID generated once per run.
- The complete run package is exactly the `runs/<study>/<run_id>/` folder,
  which is cleaner to archive and transfer than the former flat layout.
- COM and SH for a conformer MUST be written into the **same** directory
  (`gaussian_jobs/`) so the operational pair travels together.
- An isolated XYZ, COM, or SH remains attributable through its identifiers, but
  it is not promised to be independently reproducible without the matching
  manifest. Documentation MUST tell users to transfer and archive the package
  together, and MUST document that SLURM scripts are submitted from the
  directory that holds their `.com` (see §5).
- Adding molecules later means starting a **new run** in a new `run_id/` folder.
  Two clean single-purpose runs are the supported way to grow a study; there is
  no resume path (see §7).

The `create_run_manifest` `FileExistsError` guard remains as a safety net. It
never fires in normal use because every run writes a fresh directory; it is not
weakened to permit overwriting an existing run folder.

## 3. Stable identifiers and hashes

Every v2.1 run MUST define:

- `manifest_schema`: the manifest schema version;
- `run_id`: an opaque, immutable, collision-resistant identifier generated once
  for the execution;
- `config_hash`: SHA-256 over a canonical JSON representation of the complete
  scientific and operational configuration, excluding timestamps, absolute
  machine-specific paths, and output file hashes;
- `artifact_id`: an immutable, unique identifier for each XYZ, COM, and SH
  record within the run. In v2.1 the `artifact_id` is the **single** reverse
  link recorded in artifact bodies; from it the manifest yields `run_id`,
  `config_hash`, versions, and the full search configuration.

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

For a molecule taken through the **provisional undefined-stereo** path (§9), the
molecule record MUST additionally carry:

```text
provenance_status = provisional_undefined_stereo
undefined_centers
pubchem_smiles
arbitrated_smiles
```

Molecules on the normal path carry `provenance_status = normal` and do not carry
the provisional-only fields.

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

## 5. Per-artifact metadata contract (reduced from v2.0)

The manifest keeps the full provenance matrix. Artifact bodies keep only inline
science plus one `artifact_id`. This is the deliberate reduction of the v2.0
matrix; the removed fields are **manifest-only** and their absence from artifact
bodies is not a finding.

### Conformer XYZ scientific artifact

The XYZ comment (line 2) MUST carry exactly:

```text
dE=<kcal/mol> method=<MMFF94|UFF> artifact_id=xyz-…
```

- `dE` stays in kcal/mol and labeled; `method` stays because a chemist opening
  the XYZ needs to know MMFF vs UFF. These are inline-useful science.
- For a provisional undefined-stereo structure (§9), the comment carries
  `dE=NA` (a single structure has no ensemble reference — a `dE=0.000` MUST NOT
  be fabricated) and a visible `PROVISIONAL: stereo arbitrated at <atoms>`
  marker.
- Removed (manifest-only): `run_id`, `config_hash`, `pipeline_version`,
  `rdkit_version`, `conformer_id`, `relative_energy_kcalmol` as a separate
  field.

### Gaussian COM scientific artifact

The Gaussian title section MUST be a single line carrying exactly:

```text
<name> PCM 298 K <basis> dE=<kcal/mol> artifact_id=com-…
```

- For a provisional undefined-stereo structure (§9), the title carries `dE=NA`
  and a visible `PROVISIONAL: stereo arbitrated at <atoms>` marker.
- Removed (manifest-only): the entire second `provenance … commit=…` line, and
  the former first-line linkage fields `run_id`, `config_hash`, `conformer_id`,
  `relative_energy_kcalmol`, `pipeline_version`, `rdkit_version`.
- The COM still directly records routes, coordinates, charge/multiplicity,
  processor count, checkpoint directives, and Link1 behavior. Those are the
  science of the job and are unchanged.

### SLURM SH operational artifact

The script header MUST carry exactly:

```text
# artifact_id=sh-…
# source_com=<basename> sha256=<…>
```

- The `source_com` sha256 stays — it is the one operationally load-bearing
  header (it lets the job confirm it is running the intended COM).
- Removed (manifest-only): `run_id` as a separate header field, and the former
  `source_com_relative_path` (the COM is now co-located, referenced by
  basename).
- The SH body drops all path-resolution machinery (no `SCRIPT_DIR`, `cd`, or
  `../`). Because COM and SH are co-located in `gaussian_jobs/` (§2), the script
  assumes the `.com` is in the current working directory and ends with
  `ml gaussian16` then `g16 <base>_F.com`.

  The bundled default uses Purdue RCAC/Lmod syntax; a custom template MAY use
  the module-loading command required by another site, provided Gaussian 16 is
  loaded before `g16 <base>_F.com`.

**Operational contract (MUST be documented in the notebook next-steps cell):**
on the HPC, submit the `.sh` from the directory that contains its `.com`. The
COM+SH pair travels together in `gaussian_jobs/`; the script does not hunt for
its input. This is a deliberate trade — submit-from-anywhere robustness is
dropped in exchange for a clean, obvious script.

## 6. Output mapping and collision behavior

v2.1 MUST fail loudly before mutation when any of the following occurs, now
evaluated **within a single run folder**:

- distinct molecule labels sanitize to one output basename;
- duplicate source artifact paths are supplied;
- two source artifacts map to one destination path;
- two COM paths map to one SLURM script path;
- a required source artifact is missing, blank, or zero bytes;
- a manifest or stage log contains duplicate `artifact_id` values;
- an artifact path or hash disagrees with the manifest.

v2.1 MUST NOT auto-disambiguate filenames. The one-to-one source→destination
mapping (molecule identity → conformer → XYZ → COM → SH) is preserved.

## 7. Resume and append policy (removed)

Resume and append reuse are **removed** in v2.1. Every run is fresh and
immutable; to add molecules, start a new run.

- `search_conformers` MUST NOT expose an `append` parameter or any resume/reuse
  path.
- Pointing a run at an existing populated run folder MUST raise (the immutable
  `run_manifest.json` `FileExistsError` guard) rather than reusing, appending
  to, or repairing prior conformers.
- The v2.0 dirty/no-git reuse language is retired with the feature it governed.
  There is nothing to reuse, so there is no clean-commit precondition to check.

This is a deliberate scope reduction: it deletes a large, correct-but-fragile
body of code whose only purpose was to avoid recomputing conformer searches
(seconds of ETKDG against hours of downstream DFT). The v2.0 carried finding
B-04 (stale conformers resurrected on a dirty/blank-commit resume) is removed by
construction — the code path that carried it no longer exists.

## 8. Legacy v1.1 boundary

Unchanged from v2.0. The Open Babel single-geometry pathway is deprecated
compatibility functionality. It is explicitly exempt from the strict v2 manifest
and per-artifact metadata contract unless a future release deliberately upgrades
it. Documentation MUST label it legacy and MUST NOT imply that its outputs
satisfy the v2 provenance guarantees.

## 9. Frozen scientific judgments

The following decisions are approved for v2.1:

1. **PubChem's stereo-bearing `SMILES` is the source of record. Amended:**
   molecules with undefined stereochemistry are **no longer skipped**. They are
   taken through the provisional path (judgment #5) — provisionally embedded and
   loudly flagged, never silently guessed and never presented as the defined
   stereoisomer.
2. Distinctness relies on ETKDG embed-time RMSD pruning. Post-optimization RMSD
   re-pruning is deferred.
3. UFF is permitted only as a clearly recorded fallback. Relative force-field
   energies are interpreted only within one molecule and one method, never
   across molecules or between MMFF94 and UFF.
4. RDKit builds the structure represented by the supplied SMILES/protonation
   state. Gaussian charge and multiplicity are supplied explicitly; the
   conformer stage does not infer electronic multiplicity.
5. **Provisional structure for undefined stereo (new).** When a molecule's
   PubChem `IsomericSMILES` leaves one or more stereocentres unspecified, the
   pipeline produces **one** provisional starting structure instead of skipping:
   - the undefined stereocentre(s) are detected via `FindPotentialStereo` and
     the affected atom(s) are recorded as `undefined_centers`;
   - **one** conformer is embedded from the PubChem `IsomericSMILES` (single
     ETKDG embed, fixed seed) with a **light MMFF/UFF cleanup** — explicitly not
     an ensemble conformer search;
   - RDKit necessarily fixes the undefined centre(s) to an arbitrary
     configuration at embed time. The post-embed isomeric SMILES is read back
     and recorded as `arbitrated_smiles`, distinct from the underspecified
     `pubchem_smiles`;
   - the molecule is marked PROVISIONAL everywhere it surfaces:
     `provenance_status = provisional_undefined_stereo` in the manifest molecule
     record and the `conformer_log.csv` `provenance_status` column; a
     `PROVISIONAL: stereo arbitrated at <atoms>` marker and `dE=NA` in the XYZ
     comment and COM title; and a loud console warning;
   - **honesty guardrails (invariant-critical):** the arbitrated structure is an
     unvalidated DFT starting geometry, not the compound's real configuration,
     and is one arbitrary pick among 2^k for k undefined centres. It MUST never
     be presented as the defined stereoisomer. With k > 1 the warning states
     "k undefined centres → 2^k isomers possible; one arbitrated"; exactly one
     structure is still emitted (stereoisomer enumeration is deferred).
   - provisional molecules flow through the Gaussian/SLURM stages unchanged,
     carrying the marker via the one `provenance_status` column — no
     special-case downstream branch.

If MMFF parameters are missing for a provisional molecule, the record MUST state
whether UFF fallback was used or the cleanup was skipped; chemistry is never
silently substituted.

## 10. Review boundary and stop rule

A release-blocking finding must demonstrate one of the following:

- wrong chemistry or violation of a scientific invariant;
- silent data loss, overwrite, or one-to-many/many-to-one mapping corruption;
- molecular identity mismatch;
- violation of this exact provenance contract;
- a supported workflow producing an invalid or non-runnable job;
- regression in required tests or mechanical invariants.

A proposed improvement that expands this contract is recorded for a future
release unless Ish explicitly approves the expansion. It is not automatically a
Blocker or Major. In particular, the absence of the v2.0 per-artifact fields
that §5 deliberately removed MUST NOT be reported as a finding.

After implementation aligns with this contract:

1. run one holistic base-to-head audit;
2. fix actual contract-level Blocker and Major findings;
3. run one final re-review;
4. defer non-contract recommendations;
5. merge when CI is green, no contract-level Blocker or Major remains, and Ish
   gives final approval.
