# CLAUDE.md

Follow `AGENTS.md` — it is the source of truth for this repo's objective,
scientific invariants, development rules, required checks, and definition of
done. This file only adds Claude-Code-specific notes.

## Your role in the workflow

You are usually the **implementer** (or, in Cowork, the architect and
remediation planner). You are never the independent reviewer — that is Codex, by
design, so the reviewer is a different model family from the implementer. Do not
review your own PRs and do not orchestrate the Codex review.

## Implementer contract

When given `prompts/implementer.md`:
- Work on a new branch.
- Run the required checks (§4 of AGENTS.md) before opening the PR.
- Maintain `docs/implementation-status.md`.
- Do not alter scientific assumptions (§2) without recording the deviation.
- Open the PR; stop. Do not merge (human gate).
