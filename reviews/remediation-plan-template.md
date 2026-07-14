# remediation-plan-round-NN.md

**Responds to:** review round NN
**Frozen contract:** `docs/release-contract-v2.0.md`

## Finding decisions

| ID | Contract violation or expansion | Accept / Reject requested design | Approved response | Verify by |
|---|---|---|---|---|
| ID | contract section or future proposal | decision | task or deferral | test/check |

Rejecting a requested design is legitimate when it conflicts with the frozen
architecture. The plan must still address any accepted underlying risk through
the approved architecture.

## Ordered accepted tasks

1. Finding ID — concrete change and adjacent-boundary regression tests.

## Rejected requested designs

For each rejection, cite the frozen contract and explain how the approved design
addresses the underlying concern or why it is deferred.

## Deferred v2.1 proposals

List recommendations that expand the frozen contract.

## Implementation evidence

| ID | Commit | Verification | Result |
|---|---|---|---|
| ID | SHA | command/test | pass/fail |

Do not mark a finding resolved until verification passes.
