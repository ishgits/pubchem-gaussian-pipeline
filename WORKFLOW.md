# Agentic review workflow

This repository uses a contract-first handoff between Ish, the architect or
remediation planner, the implementer, CI, and an independent reviewer.

## Governing artifacts

- `AGENTS.md` — universal operating and review rules.
- `docs/release-contract-v2.0.md` — frozen v2.0 provenance, artifact, collision,
  reuse, scientific-judgment, and review boundary.
- `docs/architecture.md` — system structure and data flow.
- `docs/implementation-plan.md` — ordered release tasks and acceptance tests.
- `docs/implementation-status.md` — honest current merge-gate status.
- `reviews/*-template.md` and `prompts/*.md` — stable handoff formats.

## Workflow

```text
Ish approves architecture and release contract
        |
implementer changes code on a branch
        |
objective floor: tests + invariants + diff + repo hygiene
        |
one holistic base-to-head Codex review
        |
Ish/remediation planner accepts contract violations or rejects contract expansions
        |
implementer fixes accepted findings and records verification
        |
one final re-review
        |
Ish makes the merge decision
```

## Architecture freeze

Once Ish freezes a release contract:

- reviewers assess conformance to that contract;
- reviewers may recommend broader improvements, but must label them proposed
  contract expansions;
- proposed expansions are deferred to the next release unless Ish explicitly
  accepts them;
- remediation plans may legitimately reject a finding whose requested design
  conflicts with the frozen architecture, while still accepting and solving the
  underlying risk through the approved architecture.

For v2.0, full conformer configuration belongs in `run_manifest.json`; XYZ and
COM files carry stable manifest linkage rather than every duplicated knob.

## Objective floor

Before review:

```bash
pytest tests/ -q
python scripts/check_invariants.py
git diff --check
test -z "$(git ls-files -ci --exclude-standard)"
```

Before release, repeat the checks from a clean `git archive` and accurately
record the dependency environment used.

## Review scope

The holistic review covers every supported public entry point and handoff:

```text
identity
complete groups
manifest lineage and hashes
unique source and destination paths
failure before mutation
resume and append
zero-job behavior
scientific invariants
clean-archive reproducibility
```

It includes malformed, blank, missing, zero-byte, duplicate, colliding, stale,
damaged, dirty-git, and no-git cases where relevant.

## Native Codex review

Comment `@codex review` on the PR. Codex reads `AGENTS.md` and the frozen release
contract. Every finding must state whether it is:

- an actual frozen-contract violation; or
- a proposed expansion for a future release.

Re-review classifies prior findings as Resolved, Partially resolved, Unresolved,
Rejected with justification, or Regressed.

## Stop rule

The v2.0 loop ends after:

1. the implementation aligns with the frozen contract;
2. one holistic review is completed;
3. actual contract-level Blocker/Major findings are fixed;
4. one final re-review finds no contract-level Blocker/Major regression;
5. non-contract recommendations are recorded for v2.1;
6. Ish approves the merge.

Small documentation-only PRs may use the lightweight path: green floor and one
human review, unless they change a governing contract.
