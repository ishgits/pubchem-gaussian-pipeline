# review-round-NN.md

**PR:** number
**Head SHA:** SHA
**Reviewer:** Codex model
**Round:** NN
**Frozen contract:** `docs/release-contract-v2.0.md`

> The reviewer inspected the base-to-head diff and implementation directly.

## Merge assessment

State whether any frozen-contract Blocker or Major remains. Distinguish future
contract proposals from current release violations.

## Holistic coverage

Record which areas were inspected:

```text
scientific invariants
identity
complete groups
manifest and hashes
artifact linkage
one-to-one paths
failure before mutation
resume and append
zero-job behavior
clean archive
```

## Findings

Every finding must include:

- stable ID and severity;
- location;
- evidence;
- consequence;
- required remediation;
- verification;
- contract mapping: exact section or `Future contract proposal`.

### Blockers

### Majors

### Moderates

### Minors

### Verified strengths

## Future contract proposals

List useful recommendations that do not violate the frozen v2.0 contract. These
are deferred unless Ish explicitly expands the release scope.

## Re-review classification

| Finding | Status | Contract basis | Verification |
|---|---|---|---|
| ID | Resolved / Partially resolved / Unresolved / Rejected with justification / Regressed | section | command or test |

## Final statement

State one of:

- frozen-contract Blocker/Major remains;
- no frozen-contract Blocker/Major remains; ready for Ish's merge decision.
