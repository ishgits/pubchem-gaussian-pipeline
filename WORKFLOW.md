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
CI      ▼  review-readiness.yml  (tests + invariants + repo hygiene)  ← floor
        │        green?  ──no──► fix; not review-ready
        │        yes
Codex   ▼  native `@codex review` on the PR (reads AGENTS.md §6)
         posts findings → review-round-NN
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
- `docs/architecture.md` / `implementation-plan.md` — the two gate-1 artifacts
  (canonical v2). Earlier drafts and per-round history live under
  `docs/review-history/v2/`.
- `docs/implementation-status.md` — implementer's self-report and the AGENTS.md
  §5 **merge gate** (reviewer verifies, never trusts). §6 "questions requiring
  scientific judgment" is what Ish reads first. It must never sit as the empty
  template — `scripts/check_invariants.py` fails the floor if this file still
  holds template placeholders. The full per-round evidence (commit hashes, verify
  runs) is archived in `docs/review-history/v2/`.
- `reviews/*-template.md` — the review and remediation artifact formats.
- `prompts/*.md` — the versioned brief each agent gets, so handoffs are identical
  every time (kills "no universal handoff format" + "context trapped in chats").
- `scripts/check_invariants.py` — the mechanical floor.
- `.github/workflows/review-readiness.yml` — runs the floor on every PR
  (tests + scientific-invariant checks + repo-hygiene: no tracked generated files).

## Running the Codex review (native, no secrets)

Turn on Code review for this repo in Codex settings, then comment `@codex review`
on the PR. Codex reads the review guidelines in `AGENTS.md` §6 (full
Blocker→Verified-Strength taxonomy, stable IDs, code anchors) and posts its
findings on the PR. No API key, no extra workflow, and no third-party secret is
required — the review runs entirely through Codex's native GitHub integration
against the same `AGENTS.md`, so there is no separate setup to keep in sync.

On a re-review Codex does a diff-based pass and classifies each prior finding
Resolved / Partial / Unresolved / Rejected(justif.) / Regressed.

## A lightweight path (don't over-ceremony small PRs)

The full 2-round loop is for PRs touching scientific logic (route lines,
geometry, the Link1 contract, energetics). For docstrings, plotting, or README
tweaks: floor green + one human glance is enough. Reserve the ceremony for
load-bearing changes.
