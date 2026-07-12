# Implementer brief (Claude Code / Codex)

Implement `docs/implementation-plan.md` against `docs/architecture.md`.

- Read `AGENTS.md` first and obey its scientific invariants (§2) and development
  rules (§3).
- Work on a new branch named `round-<NN>/<short-slug>`.
- Implement the plan's tasks in order.
- Run the required checks (`pytest tests/ -q` and
  `python scripts/check_invariants.py`) and make them green before opening the PR.
- Create/maintain `docs/implementation-status.md` recording: what was and wasn't
  implemented, deviations from architecture, tests added, known limitations, and
  questions requiring scientific judgment.
- Do NOT alter scientific assumptions without recording the deviation in the
  status doc.
- Open a pull request. Do not merge. Stop.
