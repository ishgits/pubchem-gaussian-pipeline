# M-38 remediation — atomic conformer package publication

## Status

**Implemented and locally verified.**

This remediation closes the review finding that conformer-group publication was
atomic across the XYZ ensemble and `run_manifest.json`, but not across the
required subordinate index `conformer_log.csv`.

## Root cause

`record_conformer_group()` transactionally published staged XYZ files and the
candidate manifest. `search_conformers()` accumulated log rows in memory and
rewrote `conformer_log.csv` only after all molecule attempts completed.

That split transaction boundary allowed two invalid package outcomes:

1. A failed replacement could restore the old manifest and XYZ ensemble while
   the final batch-level CSV rewrite omitted the restored group.
2. A final `DataFrame.to_csv()` failure could leave a newly committed
   manifest/XYZ ensemble with an absent or stale conformer log.

The manifest remained authoritative, but the package no longer satisfied the
frozen requirement that its subordinate stage index agree exactly with it.

## Implemented design

### 1. Complete candidate log is staged per conformer group

Before publication, `search_conformers()` now:

- derives all existing conformer-log rows from current manifest authority;
- removes the molecule being replaced;
- adds the complete newly staged group;
- writes the complete candidate CSV to a same-package temporary file.

The candidate is not a partial per-molecule append. It represents the exact
complete conformer index that must exist after the group commits.

### 2. Manifest/CSV artifact-set equality is validated before mutation

`record_conformer_group()` now accepts:

- `conformer_log_path`; and
- `staged_conformer_log_path`.

When supplied, both are required. The publisher reads the staged CSV and rejects
it before mutation unless:

- it contains `artifact_id`;
- artifact IDs are unique; and
- its artifact-ID set exactly equals the candidate manifest's complete XYZ
  artifact-ID set.

### 3. XYZ, conformer log, and manifest share one rollback path

The staged conformer log is now included in the same placement sequence as the
staged XYZ files. Existing destination bytes are backed up before replacement.
If XYZ placement, CSV replacement, or manifest writing fails, the publisher:

- removes newly placed files;
- restores prior XYZ bytes;
- restores the prior conformer log; and
- leaves the prior manifest unchanged.

This remains an ordinary exception-safe filesystem transaction, not a journaled
crash-consistency protocol.

### 4. The subordinate index is reconstructible from manifest authority

`search_conformers()` now has a canonical manifest-to-log reconstruction helper.
It derives all overlapping fields from the manifest, including:

- run/config/artifact identity;
- molecule and conformer identity;
- search metadata;
- XYZ paths and manifest hashes; and
- pipeline/RDKit provenance.

When all currently published XYZ files verify, a missing or stale log is repaired
atomically before regeneration begins. Resume decisions remain based on the
input CSV, so a damaged log still triggers the intended regeneration attempt;
the reconstructed complete log becomes the rollback state if that attempt fails.

### 5. Returned results represent committed state only

The stage no longer performs an unconditional batch-level `to_csv()` after the
loop. It reads and returns the committed on-disk conformer log. Failed attempted
publications therefore cannot leak hypothetical or partial in-memory rows into
the reported result.

### 6. Failed regeneration preserves the previous complete ensemble

Eligibility or conformer-generation failure no longer independently removes
prior manifest lineage before a replacement can commit. For a previously valid
published ensemble, failure preserves the old complete manifest/XYZ/log state.
For a molecule with no prior ensemble, no conformer records are exposed.

## Files changed

- `pipeline/manifest.py`
  - expanded `record_conformer_group()` transaction to include the staged
    conformer log;
  - added exact candidate-manifest/CSV XYZ artifact-set validation;
  - extended rollback messaging from XYZ-only to complete package rollback.

- `pipeline/conformers.py`
  - added canonical manifest-derived conformer-log reconstruction;
  - added atomic CSV staging/replacement helpers;
  - stages a complete candidate CSV for every group publication;
  - removed the separate final batch-level CSV rewrite;
  - returns committed package state.

- `tests/test_conformers.py`
  - updated staged XYZ failure expectation to preserve the previous complete
    group rather than return an empty log;
  - added regression coverage for CSV staging failure, CSV replacement failure,
    XYZ placement failure, manifest-write failure, and successful smaller-group
    replacement with exact manifest/log artifact-ID equality.

## Required state invariant

After every completed or caught ordinary failure, the package must be in one of
these states:

```text
A. prior complete manifest group + prior XYZ files + matching prior log rows
B. new complete manifest group + new XYZ files + matching new log rows
C. no manifest group + no corresponding log rows
```

The implementation rejects or rolls back any attempted state where a manifest
XYZ artifact is omitted from the complete conformer index.

## Verification

Environment used for this remediation:

```text
Python 3.13.5
pandas 2.2.3
requests 2.32.5
RDKit 2025.09.4
pytest 9.0.2
```

Results:

```text
pytest -q: 412 passed
python scripts/check_invariants.py: passed
Python compilation: passed
notebook JSON validation: passed
git diff --check: passed
ignored-but-tracked hygiene: passed
clean git archive: 412 passed; invariants and compilation passed
```

The pinned Python 3.12 release environment was not independently rerun in this
remediation session. Remote CI and Ish's final merge decision remain required.
