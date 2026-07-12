# Reviewer brief (Codex) — reference copy

The live copy Codex runs in CI is `.github/codex/prompts/review.md`. This is the
human-readable reference.

Independently audit the PR against `docs/architecture.md`,
`docs/implementation-plan.md`, and the scientific invariants in `AGENTS.md` §2.
Inspect the implementation directly — do not trust
`docs/implementation-status.md`.

Review: correctness, scientific validity, reproducibility, failure handling,
testing, documentation, scope compliance.

Produce findings classified Blocker / Major / Moderate / Minor / Verified
Strength. Every finding: stable ID (B-01, M-01, MOD-01, MIN-01), location
(file:line-range), evidence, consequence, required remediation, and how to
verify the fix. IDs are stable across rounds.

If prior review rounds exist in `reviews/`, do a diff-based re-review: classify
each prior finding Resolved / Partial / Unresolved / Rejected(justif.) /
Regressed, and check for regressions introduced by the fixes.
