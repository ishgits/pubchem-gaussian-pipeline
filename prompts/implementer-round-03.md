# Implementer brief — round 03 (Claude Code)

Implement `docs/remediation-plan-round-03-v2.md` (Codex findings **M-06, M-07**).
That plan is the source of truth for this round; this brief only sets the rules.

## Rules
- Read `AGENTS.md` first. Obey the scientific invariants (§2) and development
  rules (§3). This is a docs + additive-provenance change — **no route line,
  unit, charge/multiplicity, or Link1 edits**.
- **Branch:** continue on the open PR branch `feat/conformer-search-v2`
  (PR #3), as with rounds 01–02. Do **not** cut a new branch. One commit per
  finding (`M-06: …`, `M-07: …`), plus a docs commit if needed.
- Implement the plan's tasks in order. The three locked judgment calls are
  already decided in the plan header — do not re-litigate them:
  1. **M-06 value** = both `pipeline.__version__` (`"0.2.0"`) **and** best-effort
     git short-SHA + `.dirty` marker, empty string when git is absent/errors.
  2. **M-06 scope** = conformer path only (`conformer_log.csv` + XYZ comment).
     Do **not** touch `com_write_log.csv` or `sdf_download_log.csv`.
  3. **M-07** = keep `implementation-status-v2.md` as the working file; sync it
     into the canonical `docs/implementation-status.md`; add the
     `check_invariants.py` drift guard (now **in scope**).
- The provenance helper must be import-safe and offline-safe: wrap the git
  subprocess in `try/except`, return `""` on any failure, never raise, never
  require network or a `.git` dir. Tests run offline (§4) and must not assert a
  concrete SHA.

## Required checks (must be green before the PR)
- `pytest tests/ -q` — including the new M-06 provenance tests and the M-07
  drift-guard unit test.
- `python scripts/check_invariants.py` — including the new drift guard. Sanity-
  check it both ways: reverting `implementation-status.md` to the template makes
  it **fail**; the populated file passes.

## Status doc (AGENTS.md §5)
- Update **both** `docs/implementation-status.md` (canonical gate) and
  `docs/implementation-status-v2.md` (working file) for this round: what was /
  wasn't implemented, deviations (new log columns; extended XYZ comment format),
  tests added, known limitations, scientific-judgment questions, and the
  Provenance block (now naming pipeline version + commit).
- Record the reproducibility caveat from the plan: `pipeline_commit` pins code
  only when clean; `.dirty` = not reproducible from the commit alone; empty =
  fell back to the manually-bumped `__version__`. Do not imply a recorded SHA
  guarantees exact code.
- Fill the round-03 **Implementation evidence** table in the remediation plan as
  fixes land. Do **not** mark a finding resolved on code change alone — run the
  listed verify step first.

## Definition of done (AGENTS.md §5)
Required checks green; both status docs updated; deviations recorded, not silent;
M-06 and M-07 each Resolved with evidence.

Open the pull request. Do **not** merge. Stop. (Codex reviews; Ish holds the
merge gate.)
