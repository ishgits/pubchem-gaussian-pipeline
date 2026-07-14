# Implementer brief

Implement `docs/implementation-plan.md` against `docs/architecture.md` and the
frozen `docs/release-contract-v2.0.md`.

- Read `AGENTS.md` first.
- Work on a branch, never directly on `main`.
- Implement tasks in plan order.
- Do not broaden the frozen contract without Ish's explicit approval.
- Validate all complete inputs before mutating outputs or logs.
- Preserve one-to-one source-record to destination-path mapping.
- Add tests for normal, missing, blank, zero-byte, duplicate, colliding, damaged,
  dirty-git, and no-git cases where the task touches those boundaries.
- Maintain `docs/implementation-status.md` honestly.
- Run:
  ```bash
  pytest tests/ -q
  python scripts/check_invariants.py
  git diff --check
  test -z "$(git ls-files -ci --exclude-standard)"
  ```
- Repeat the relevant checks from a clean `git archive` before requesting final
  review.
- Open or update the PR. Do not merge; Ish is the human gate.
