# implementation-status.md

> Canonical merge-gate status. Reviewers must verify every claim against code
> and tests. Historical remediation evidence remains under
> `docs/review-history/v2/`.

**PR:** #3
**Branch:** `feat/conformer-search-v2`
**Current phase:** v2.0 final remediation
**Merge status:** **local remediation complete; requires push, green remote CI, one current-head re-review, and Ish's final approval**

This status is governed by `docs/architecture.md`,
`docs/implementation-plan.md`, and `docs/release-contract-v2.0.md`.

## 1. Ordered release-task completion

### Task 1 — canonical manifest utilities: complete

`pipeline/manifest.py` implements schema `2.0`, UUID run identity,
order-independent canonical JSON, deterministic SHA-256 `config_hash`, stable
conformer/artifact IDs, molecule identity hashes, relative package paths,
atomic manifest writes, duplicate-record rejection, exact lineage validation,
and SHA-256 verification. Timestamps, machine-specific absolute paths,
run/artifact IDs, and output hashes are excluded from `config_hash`; the SLURM
template identity remains included because it is operational configuration.

### Task 2 — conformer linkage: complete

Every retained conformer has one stable conformer record and one unique XYZ
artifact. `conformer_log.csv`, the XYZ comment, and the manifest agree on
`run_id`, `artifact_id`, `config_hash`, conformer identity, relative energy,
method, pipeline version, and RDKit version. The manifest records the complete
search configuration and verified XYZ file hash. Resume/append validation checks
the manifest record and current file bytes, not only CSV values.

### Task 3 — Gaussian linkage: complete

Conformer-derived COM generation requires a real manifest-linked XYZ parent.
Both batch and direct conformer-specific writes validate source identity,
versions, configuration, path, hash, conformer record, energy, and convergence
before output mutation. The COM title carries the approved minimal linkage
matrix. The manifest and `com_write_log.csv` record parent XYZ artifact ID,
relative COM path, and COM hash. Route lines, coordinates, charge/multiplicity,
checkpoint names, processor directives, and Link1 opt→freq construction remain
unchanged.

### Task 4 — SLURM one-to-one mapping: complete

Log-driven generation validates all COM inputs before mutation: nonblank,
regular, nonzero files; unique normalized source paths; unique script
destinations; no basename collapse; and exact manifest path/hash agreement.
Each script header records run/artifact identity plus source COM relative path
and hash. The manifest and log record each SH hash. Collision, duplicate,
zero-byte, missing, and hash-failure tests verify byte-for-byte preservation of
prior scripts and logs.

### Task 5 — strict reuse: complete

Version-only fallback was removed. Reuse requires the same clean, nonblank
pipeline commit plus all prior identity/config/group/file/hash checks. Changed,
dirty, or missing commits regenerate. Source archives and no-git installations
therefore never reuse cached conformers.

### Task 6 — notebook and user documentation: complete

The notebook creates `run_manifest.json` after PubChem identity resolution and
before XYZ/COM/SH artifacts, passes the manifest through all v2 stages, and
finalizes it by verifying every artifact hash. README/notebook guidance defines
the complete run package as the supported archive/transfer unit, warns that
isolated artifacts require their manifest, labels v1.1 exempt/deprecated, and
documents no-git regeneration.

### Task 7 — mechanical invariants and tests: complete

Offline coverage now includes schema/canonical hashing, ordering stability,
scientific config drift, duplicate records/IDs, stable lineage, file hashes,
XYZ/COM/SH linkage fields, direct COM linkage rejection, source/destination
collisions, duplicate normalized paths, zero-byte COMs, failure-before-mutation,
strict clean-commit reuse, and zero-job manifests. The invariant checker enforces
the frozen per-artifact matrix without requiring full search-knob duplication.

### Task 8 — finding disposition and review: complete locally

One holistic contract audit was performed after tasks 1–7. It identified and
resolved two new Major traceability issues before the final re-review:

- **M-24 — Resolved.** Evidence before remediation: the direct conformer COM
  writer accepted nonblank identifiers without proving they belonged to a real
  manifest parent. Consequence: an apparently linked COM could be orphaned or
  falsely attributed. Required remediation: validate full manifest identity,
  parent path/hash/record/config, then record the COM artifact. Implemented in
  `pipeline/gaussian.py:100-219` and `pipeline/gaussian.py:287-455`;
  `tests/test_gaussian.py:318-339` rejects a tampered COM artifact ID before
  directory creation.
- **M-25 — Resolved.** Evidence before remediation: relative artifact paths could
  contain `..` and leave the run package. Consequence: the advertised archive
  could omit a recorded scientific artifact. Required remediation: constrain
  every artifact to the manifest package root and reject path traversal.
  Implemented in `pipeline/manifest.py:397-425` and
  `pipeline/manifest.py:488-506`, with regression coverage in
  `tests/test_manifest.py:147-170`.

The subsequent current-head review identified three additional contract-level
gaps. They are resolved locally in this remediation pass:

- **B-08 — Resolved locally.** Manifest creation and validation now reject a
  Link1 frequency route unless it contains both `Geom=AllChk` and `Guess=Read`
  (case-insensitive, ordinary whitespace allowed). This prevents the generated
  title/charge/coordinate-free frequency section from becoming an invalid job.
- **M-26 — Resolved locally.** Before any Gaussian output or manifest mutation,
  `conformer_log.csv` must contain exactly the manifest's XYZ artifact-ID set.
  Valid-looking truncated and empty subset logs are rejected, preserving prior
  COM files, logs, failure logs, and manifest lineage byte-for-byte.
- **M-27 — Resolved locally.** Before script pruning or writing,
  `com_write_log.csv` must contain exactly the manifest's COM artifact-ID set.
  Valid-looking truncated and empty subset logs are rejected, preserving prior
  scripts, SLURM logs, and manifest lineage byte-for-byte.

The next exact-head review identified two more stage-ordering violations. They
are also resolved locally:

- **B-09 — Resolved locally.** The documented root-level `run_manifest.json` is
  now ignored as a generated run artifact. Creating the immutable manifest in a
  clean checkout no longer makes the source tree appear dirty before conformer
  generation recomputes pipeline provenance. A real temporary Git-repository
  regression test exercises the exact manifest-then-provenance sequence.
- **M-28 — Resolved locally.** Manifest creation now applies the same nonempty,
  unique sanitized-basename rule as downstream conformer generation before the
  manifest is written. Distinct labels such as `Water` and `water` therefore
  fail without leaving an immutable manifest for an invalid one-to-one mapping.

- **B-10 — Resolved locally.** `search_conformers()` now requires the runtime
  molecule identities to cover the immutable manifest exactly when
  `append=False`. With `append=True`, current rows plus fully validated retained
  groups must account for every manifest molecule. Missing configured molecules
  fail before lineage, logs, failure records, or XYZ outputs are mutated.
- **M-29 — Resolved locally.** The Gaussian stage now requires a nonblank,
  parseable convergence flag for every manifest-linked conformer-log row and
  compares it exactly with the manifest conformer record during preflight.
  Missing, malformed, or altered convergence metadata fails before prior COM/SH
  lineage, COM files, logs, or failure records are changed.
- **M-30 — Resolved locally.** Each v2 stage now validates that its output root
  and authoritative write log stay inside the manifest run package *before* the
  first mutation, reusing the existing `relative_artifact_path` boundary check
  (no manifest-schema or path-helper change). `search_conformers` checks
  `xyz_dir` and `conformer_log.csv` before any conformer-lineage removal,
  directory creation, failure-log deletion, XYZ write, or log rewrite;
  `write_gaussian_coms_from_conformers` checks `outdir` right after manifest
  configuration and the `com_write_log.csv` after source/artifact-set preflight,
  before any COM/SH lineage removal or write; `write_slurm_scripts` checks
  `slurm_dir` and `slurm_write_log.csv` in manifest-driven mode before directory
  creation, SH-lineage removal, stale-script pruning, or script writing (the
  zero-job path is covered because it still creates `slurm_dir`). The legacy
  Open Babel v1.1 and explicit `com_dir` paths remain exempt and unchanged. Six
  new parameterized failure-atomicity regression cases (two per stage) start
  from an already-valid run, monkeypatch the writer to fail if reached, and
  prove that an outside destination raises the package-boundary `ValueError`
  while leaving the manifest, prior artifacts, and all logs byte-for-byte
  unchanged and creating no outside file. One pre-existing MIN-02 test
  (`test_failed_csv_cleared_on_clean_rerun`) placed its shared outputs outside a
  subdirectory manifest package; its fixture was corrected to root both
  manifests at the shared package directory, preserving the test's intent.

The targeted reviewed-head remediation is also complete:

- **B-11 — Resolved.** PubChem JSON caches now use a readable prefix plus the
  SHA-256 of the complete canonical GET request, store a schema-2 request/response
  envelope, and treat legacy, malformed, or request-mismatched entries as misses.
- **M-31 — Resolved.** Artifact resolution rechecks the final real path against
  the manifest package root, rejecting direct and parent-directory symlink escapes
  before verification or finalization while retaining internal symlink support.
- **M-32 — Resolved.** One shared strict boolean parser now governs conformer
  resume, direct and batch Gaussian preflight, manifest recording, and manifest
  validation; blank, missing, NaN, infinite, malformed, and ambiguous values
  cannot pass by Python truthiness.
- **M-33 — Resolved.** Every manifest conformer requires the complete frozen
  record schema, strict types and finite values, stable identity, agreement with
  authoritative search configuration, and exactly one matching XYZ artifact;
  shared, orphaned, missing, wrong-kind, or mismatched XYZ lineage is rejected.
- **MOD-06 / MIN-07 / MIN-08 — Resolved.** Molecule rows record the PubChem query
  that actually succeeded, DOI targets are consistent, and the notebook states
  the supported Python >=3.11 floor.

The holistic re-review additionally closed **M-34**, a frozen-contract manifest
consistency gap: duplicated per-conformer search knobs could previously disagree
with `configuration.conformer` without invalidating the manifest. Validation now
rejects that internal provenance drift. Regression coverage exercises every item
above on both normal and failure-atomic paths.

**M-35 — Resolved.** Conformer ensembles are now staged and published atomically
at molecule-group scope. A complete candidate manifest is validated before one
canonical manifest write; ordinary file-placement or manifest-write exceptions
restore prior XYZ bytes and manifest state. Every nonempty group must have exact
cardinality, contiguous IDs, shared search metadata, the approved convergence
shape, and exact XYZ lineage. This is deliberately not described as a journaled
crash transaction across multiple filesystem operations.

The final holistic audit also closed **M-36 — Resolved**, an additional
frozen-contract publication finding discovered during implementation: crossed or
duplicate staged source/destination mappings could otherwise mutate a protected
final path, and replacing a group with fewer conformers could leave obsolete XYZ
artifacts in the run package. The group publisher now rejects those mappings
before mutation and transactionally retires obsolete prior XYZ files while
preserving the current policy for unreferenced downstream COM/SH files.

**M-37 — Resolved.** The final re-review found that publishing through an
already-supported internal XYZ symlink could replace the link while recording
its resolved target, leaving the authoritative hash attached to old bytes. The
publisher now moves staged bytes to the same resolved in-package target used by
manifest path validation, preserving the internal symlink and verified lineage.

**M-38 — Resolved.** The post-M-37 review found that molecule-group publication
was atomic only across XYZ files and `run_manifest.json`; the stage rewrote
`conformer_log.csv` separately after all group attempts. A caught staged-write
or publication failure could therefore retain the prior complete ensemble in
the manifest while publishing an empty/subset CSV, and a final CSV-write failure
could leave newly committed manifest/XYZ state without its required subordinate
index. `record_conformer_group()` now accepts a staged complete conformer log,
requires its artifact-ID set to equal the complete candidate manifest XYZ set,
and publishes the CSV through the same backup/placement/rollback path as the XYZ
ensemble and manifest. `search_conformers()` reconstructs canonical log rows
from manifest authority, repairs a stale/missing subordinate index atomically
when the existing XYZ package verifies, stages the complete candidate log for
each group, and returns only the committed on-disk log state. Ordinary staged
CSV write, CSV replacement, XYZ placement, and manifest-write failures preserve
the previous complete manifest/XYZ/log package. Successful smaller-group
replacement removes obsolete XYZ files and leaves exact manifest/log XYZ
artifact-ID equality.

No known local frozen-contract Blocker or Major remains. This exact patched head
has received a final local holistic re-review; remote CI and Ish's human merge
decision remain required.

## 2. Prior finding disposition

- **M-21 — Rejected as written; underlying concern Resolved.** Full conformer
  search configuration is authoritative in the manifest; COM carries the exact
  approved linkage matrix. Requiring all knobs in every COM contradicts the
  frozen manifest-centric contract.
- **M-22 — Rejected as written; underlying concern Resolved.** Full conformer
  search configuration is authoritative in the manifest; XYZ carries the exact
  approved linkage matrix. Requiring all knobs in every XYZ contradicts the
  frozen contract.
- **M-23 — Resolved.** Same-basename COM paths, duplicate normalized sources,
  duplicate destinations, zero-byte inputs, and manifest/hash disagreement all
  fail before SLURM mutation.
- **MIN-06 — Resolved.** The exact pinned stack was installed and exercised on
  2026-07-13 rather than described as an unperformed target.

## 3. Verification evidence

Pinned release-target environment used locally:

```text
Python 3.12.13
pandas 3.0.3
requests 2.34.2
RDKit 2025.09.3
pytest 9.1.1
```

Current local results:

```text
pytest tests/ -q: 412 passed
python scripts/check_invariants.py: passed
Python compilation: passed
notebook JSON validation: passed
clean `git archive`: 412 passed; invariant checks and compilation passed
```

The M-38 remediation and clean-archive verification were run on CPython 3.13.5
with pandas 2.2.3, requests 2.32.5, RDKit 2025.09.4, and pytest 9.0.2. This was
not an independent rerun of the pinned Python 3.12 release-target environment.

The M-35 through M-37 verification used the pinned dependency versions listed above on
CPython 3.12.10.

The M-30 re-verification above (284 passed, six new package-boundary cases)
was run with the pinned dependency versions (pandas 3.0.3, requests 2.34.2,
RDKit 2025.09.3, pytest 9.1.1) on CPython 3.13.3 rather than the 3.12.13
release-target interpreter recorded above.

The publication checkout also passed `git diff --check` and the
ignored-but-tracked hygiene guard. Remote CI and Ish's final human merge gate
remain; successful execution is not treated as scientific validation.

## 4. Architecture deviations and scientific invariants

No deviation from the frozen v2.0 architecture or release contract remains.
There was no silent change to route lines, units, charge/multiplicity, coordinate
format, checkpoint naming, or the Link1 contract. Starting geometries remain
explicitly approximate force-field DFT starts; MMFF94/UFF energies remain
kcal/mol and are never mixed with DFT Hartree values.

## 5. Known limitations and v2.1 recommendations

Deferred exactly as approved: `runs/{study}/{run_id}/`, collision-proof internal
filenames, xTB/CREST, energy windows, post-optimization RMSD re-pruning,
solvent-aware conformer search, Gaussian execution/parsing, and upgrading legacy
Open Babel output to the v2 manifest contract. Flat v2.0 package directories
therefore require a fresh directory per independent study.

## 6. Questions requiring scientific judgment

None. The frozen judgments in `docs/release-contract-v2.0.md` remain unchanged.

## 7. Human gate

Ish's final merge approval is still required even after green local/remote
checks. Successful execution alone is not treated as scientific validation.
