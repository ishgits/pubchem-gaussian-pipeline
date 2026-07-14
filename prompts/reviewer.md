# Reviewer brief — Codex

Perform one holistic base-to-head audit against:

- `AGENTS.md`;
- `docs/release-contract-v2.0.md`;
- `docs/architecture.md`;
- `docs/implementation-plan.md`.

Inspect code and tests directly. Do not trust `docs/implementation-status.md`.

## Frozen v2.0 provenance interpretation

The architecture is manifest-centric:

- `run_manifest.json` holds the complete configuration and artifact hashes;
- XYZ and COM carry the exact minimal linkage fields in the release contract;
- SH identifies its source COM and hash;
- individual artifacts do not need every manifest field duplicated inside them.

Do not report missing duplicated search knobs in XYZ or COM as a finding when the
required per-artifact matrix and valid complete manifest linkage exist.

## Holistic audit checklist

Review every supported public entry point and stage boundary for:

- scientific-invariant preservation;
- molecular identity;
- complete conformer groups;
- manifest schema, canonical config hash, artifact IDs, lineage, and file hashes;
- one-to-one source-record to destination-path mapping;
- duplicate, colliding, blank, missing, zero-byte, stale, and damaged inputs;
- validation before mutation and byte-preserving failure;
- resume/append behavior under matching, changed, dirty, and unavailable commits;
- zero-job behavior;
- clean-archive reproducibility;
- documentation and actual behavior agreement.

## Finding eligibility

A release-blocking finding must demonstrate wrong chemistry, silent data loss or
overwrite, identity corruption, a frozen-contract violation, an invalid
supported job, or required-check regression.

Label recommendations that expand the frozen contract as **Future contract
proposal**. They are not automatically v2.0 Blockers or Majors.

## Output

For each finding provide:

- stable ID;
- severity: Blocker, Major, Moderate, Minor, or Verified Strength;
- file and line/function location;
- evidence;
- consequence;
- required remediation;
- verification;
- contract mapping: exact violated section, or `Future contract proposal`.

On re-review classify every prior finding as Resolved, Partially resolved,
Unresolved, Rejected with justification, or Regressed. Verify rejected findings
against the frozen contract rather than repeating the previous requested design.

After the final re-review, list non-contract recommendations separately for
v2.1 and state whether any frozen-contract Blocker or Major remains.
