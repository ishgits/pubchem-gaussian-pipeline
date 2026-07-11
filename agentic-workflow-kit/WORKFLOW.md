# Agentic review workflow — how this repo runs

A handoff substrate for going back and forth between Claude (implementer /
architect / remediation planner) and Codex (independent reviewer) over GitHub
pull requests, instead of ZIPs and chat-trapped context.

## The loop

```
        ┌─────────────── HUMAN GATE 1 ───────────────┐
        │  Ish approves architecture.md + plan        │
Cowork  ▼                                             │
(Claude) architecture.md + implementation-plan.md ────┘
        │
Claude  ▼  prompts/implementer.md
Code     branch → implement → run floor → status doc → open PR
        │
CI      ▼  review-readiness.yml  (tests + invariants)   ← objective floor
        │        green?  ──no──► fix; not review-ready
        │        yes
Codex   ▼  codex-review.yml  (prompt: full taxonomy, stable IDs)
         posts review comment + artifact  → review-round-NN
        │
        ├──────────────── HUMAN GATE 2 ───────────────┐
        │  Ish reviews substantive findings           │
Cowork  ▼  prompts/remediation-planner.md             │
(Claude) remediation-plan-round-NN.md  (accept/reject)┘
        │
Claude  ▼  implement accepted fixes; record commit+verification per finding
Code     push → CI floor reruns → Codex diff-based re-review (round NN+1)
        │        classifies each prior finding
        │
        └──────────────── HUMAN GATE 3 ───────────────┐
                          Ish makes final merge call   ┘
```

## What each file is

- `AGENTS.md` — source of truth: objective, scientific invariants, dev rules,
  required checks, definition of done, Codex review guidelines.
- `CLAUDE.md` — thin pointer + Claude-Code implementer role.
- `docs/architecture.md` / `implementation-plan.md` — the two gate-1 artifacts.
- `docs/implementation-status.md` — implementer's self-report (reviewer verifies,
  never trusts). §6 "questions requiring scientific judgment" is what Ish reads
  first.
- `reviews/*-template.md` — the review and remediation artifact formats.
- `prompts/*.md` — the versioned brief each agent gets, so handoffs are identical
  every time (kills "no universal handoff format" + "context trapped in chats").
- `scripts/check_invariants.py` — the mechanical floor.
- `.github/workflows/review-readiness.yml` — runs the floor on every PR.
- `.github/workflows/codex-review.yml` — runs Codex review after the floor is green.
- `.github/codex/prompts/review.md` — the live reviewer prompt Codex runs.

## Two ways to run the Codex review

1. **Native (5-min, no secrets):** turn on Code review for this repo in Codex
   settings, then comment `@codex review` on a PR. Codex reads `AGENTS.md` review
   guidelines. Fast, but posts P0/P1 inline comments only — taxonomy collapses,
   no committed artifact. Good for proving the loop.
2. **Structured (this kit's `codex-review.yml`):** Codex runs the custom prompt,
   emits the full Blocker→Verified-Strength taxonomy with stable IDs + code
   anchors, posts a review comment, and uploads the review as an artifact.
   Requires the `OPENAI_API_KEY` repo secret. This is what your round-2
   diff-based re-review depends on.

Both read the same `AGENTS.md`, so no wasted setup.

## A lightweight path (don't over-ceremony small PRs)

The full 2-round loop is for PRs touching scientific logic (route lines,
geometry, the Link1 contract, energetics). For docstrings, plotting, or README
tweaks: floor green + one human glance is enough. Reserve the ceremony for
load-bearing changes.
