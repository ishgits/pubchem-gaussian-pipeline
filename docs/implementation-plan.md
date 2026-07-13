# implementation-plan.md

> Frozen v2.0 release-completion plan approved by Ish on 2026-07-13. It
> supersedes the earlier open-ended provenance interpretation while preserving
> the implemented conformer science. Historical plans remain under
> `docs/review-history/v2/`.

## Objective

Complete v2.0.0 against the manifest-centric contract in
`docs/release-contract-v2.0.md`, close the remaining one-to-one SLURM mapping
gap, remove unverifiable no-git cache reuse, and run one holistic final audit.

The goal is not to place every conformer-search knob inside every artifact. The
goal is a complete authoritative manifest plus stable artifact linkage.

## Locked scientific defaults

```text
N_GENERATE=20
TOP_N=3
RMSD_PRUNE=0.5 Å
SEED=42
MMFF94 preferred
UFF recorded fallback
force-field energy unit=kcal/mol
```

Approved judgments: PubChem stereo SMILES is authoritative; undefined stereo is
skipped; no post-optimization RMSD re-pruning; force-field ΔE is compared only
within one molecule and one method; Gaussian charge/multiplicity remains
explicit.

## Ordered release tasks

### 1. Canonical manifest utilities

Implement `run_manifest.json` with:

- schema version, immutable `run_id`, and deterministic `config_hash`;
- complete run, molecule, conformer, Gaussian, and SLURM configuration;
- artifact lineage, relative paths, and SHA-256 file hashes;
- canonical JSON serialization that excludes timestamps, absolute
  machine-specific paths, and output file hashes from `config_hash`.

Acceptance:

- identical canonical configuration yields the same `config_hash`;
- scientifically relevant configuration changes alter it;
- ordering-only differences do not alter it;
- duplicate IDs or records are rejected.

### 2. Conformer-stage manifest linkage

Generate stable conformer and XYZ artifact IDs. Write the required XYZ metadata:

```text
run_id artifact_id config_hash conformer_id relative_energy_kcalmol
method pipeline_version rdkit_version
```

Record the complete search configuration and XYZ SHA-256 in the manifest.
Maintain `conformer_log.csv` as a manifest-consistent operational index.

Acceptance:

- every conformer row has one unique XYZ artifact;
- manifest and CSV identity/config values agree;
- XYZ file hash verifies;
- the full search configuration is present in the manifest, not necessarily
  duplicated in the XYZ comment.

### 3. Gaussian-stage manifest linkage

Require and forward run/artifact/config identity for every conformer-derived COM.
Write the required COM title metadata:

```text
run_id artifact_id config_hash conformer_id relative_energy_kcalmol
pipeline_version rdkit_version
```

Record parent XYZ artifact ID, full Gaussian configuration, relative COM path,
and COM SHA-256 in the manifest and `com_write_log.csv`.

Acceptance:

- a COM cannot be generated without valid manifest linkage;
- route lines, charge/multiplicity, coordinates, checkpoint names, and Link1 are
  unchanged;
- absence of duplicated `n_generate`, `top_n`, or `rmsd_prune` in the COM is
  allowed because the manifest contains them.

### 4. SLURM one-to-one mapping and linkage

Before mutation, validate every log-driven COM input:

- nonblank, existing regular file, and size greater than zero;
- unique normalized source path;
- unique destination SH path;
- no two COM basenames collapse to one script;
- COM hash agrees with the manifest.

Write `run_id`, SH `artifact_id`, source COM relative path, and source COM hash in
each script header. Record the SH hash in the manifest.

Acceptance:

- two different `same.com` paths fail before mutation;
- duplicate source paths fail;
- zero-byte COM fails;
- valid N COM inputs create N unique SH files and N unique log/manifest records;
- failing validation preserves existing scripts and logs byte-for-byte.

### 5. Strict reuse policy

Remove version-only resume and append fallback. Reuse requires the same clean,
nonblank pipeline commit on both the retained record and current run, plus all
existing identity/config/group/file/hash checks.

Acceptance:

- matching clean commit may reuse;
- changed commit regenerates;
- dirty commit regenerates;
- missing commit regenerates;
- source ZIP or no-git execution never reuses cached conformers.

### 6. Notebook and user documentation

Update the notebook and README so users:

- generate the manifest before artifacts and finalize it with file hashes;
- transfer/archive the complete run package;
- understand that isolated artifacts require the matching manifest;
- understand that legacy v1.1 is exempt from v2 guarantees;
- understand that no-git execution disables reuse.

### 7. Mechanical invariants and tests

Add offline tests and invariant guards for:

- manifest schema, canonical hash, IDs, lineage, and file hashes;
- required XYZ/COM/SH linkage fields;
- one-to-one mappings and destination collisions;
- zero-byte inputs;
- failure before mutation;
- no-clean-commit regeneration;
- zero-job manifest behavior;
- clean-archive reproducibility.

The invariant checker should enforce the frozen matrix, not an ever-expanding
requirement to duplicate all manifest fields into every artifact.

### 8. Status, finding disposition, and final review

Update `docs/implementation-status.md` with verified evidence.

- M-21 and M-22 are **rejected as written** because they require duplicating the
  full search configuration inside each COM/XYZ, contrary to the approved
  manifest-centric contract. Their underlying traceability concern is accepted
  and resolved only when tasks 1–3 are complete.
- The SLURM basename collision and zero-byte behavior are accepted and must be
  fixed under task 4.
- Recommendations outside the frozen contract are recorded for v2.1.

Run one holistic base-to-head audit and one final re-review under the stop rule.

## Required verification

```bash
pytest tests/ -q
python scripts/check_invariants.py
git diff --check
test -z "$(git ls-files -ci --exclude-standard)"
```

Repeat tests and invariants from a clean `git archive`. Before tagging, test the
actual pinned dependency stack or change the lock-file wording so it does not
claim an unperformed verification.

## Explicitly deferred to v2.1

`runs/<study>/<run_id>/`, automatic collision-proof filenames, xTB/CREST,
energy-window selection, post-optimization RMSD re-pruning, solvent-aware search,
Gaussian execution/parsing, and upgrading the legacy Open Babel path to the v2
manifest contract.
