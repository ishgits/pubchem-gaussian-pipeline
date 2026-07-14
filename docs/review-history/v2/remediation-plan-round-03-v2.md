# remediation-plan-round-03-v2.md

Produced by Cowork (Claude) from the Codex round-03 review. Approved by Ish
(gate 2) before fixes are implemented.

**Responds to:** Codex round-03 review (findings M-06, M-07).
**Branch:** remediation commits land on the open PR branch
`feat/conformer-search-v2` (PR #3), as with rounds 01–02 — not a new branch.
**Context:** The repo currently has **no notion of its own version**. No
`pipeline.__version__`, no runtime git-SHA capture. So `conformer_log.csv` rows
from two different pipeline revisions with the same RDKit/seed/method are
byte-identical, and the canonical §5 gate file was never populated.

**Ish judgment calls locked:**
- **M-06 value** → record *both* a static `pipeline.__version__` **and** a
  best-effort git short-SHA with a `.dirty` marker (empty string when git is
  unavailable). Rationale: `__version__` alone doesn't distinguish two un-bumped
  revisions (the exact failure Codex cites); a bare SHA breaks §4 offline tests
  and installed/HPC copies with no `.git`. Both together are offline-safe *and*
  per-commit unique.
- **M-06 scope** → **conformer path only** this round (`conformer_log.csv` + the
  per-conformer XYZ comment line). `com_write_log.csv` and `sdf_download_log.csv`
  are explicitly **out of scope** — tracked as a next-round candidate, not fixed
  here.
- **M-07** → **keep** `implementation-status-v2.md` as the working file; **sync**
  its current content into the canonical `docs/implementation-status.md` (the
  file AGENTS.md §5 names) so the merge gate is never stale.

## Decisions on each finding

| ID | Decision | Rationale | Owner | Verify by |
|----|----------|-----------|-------|-----------|
| M-06 | **Accept** | `conformer_log.csv` is the machine-readable provenance tying generated XYZ/`.com` back to the code that produced them. Without a code-identity field, two runs from different revisions are indistinguishable — violates the §3 "record relevant software versions" invariant. | Claude Code | provenance tests below |
| M-07 | **Accept** | AGENTS.md §5 hard-codes `docs/implementation-status.md` as the merge gate; it is still the empty template while all real status lives in the `-v2` file. Reviewers/maintainers following the mandated gate miss the PubChem SMILES rename and every v2 deviation. | Claude Code | canonical file populated + no template placeholders |

No findings rejected.

## Ordered fix tasks

### M-06 — record pipeline version + commit in conformer provenance

1. **Define the version.** Add `__version__` to `pipeline/__init__.py`
   (start at `"0.2.0"` for the v2 conformer track — Ish may pick a different
   number; the *value* is a judgment call, the *mechanism* is not).
2. **Add an import-safe, offline-safe provenance helper** (in `conformers.py`
   or `utils.py`) that returns `(version, commit)`:
   - `version` = `pipeline.__version__`.
   - `commit` = best-effort `git rev-parse --short HEAD`, with `.dirty`
     appended when `git status --porcelain` is non-empty; **`""`** if git
     errors, is absent, or the code isn't in a repo. Wrap the whole subprocess
     path in `try/except` returning `""` — it must **never** raise and must not
     require network or a `.git` dir (AGENTS.md §4: tests run offline).
3. **Record it, conformer path only.** Add `pipeline_version` and
   `pipeline_commit` to `_LOG_COLUMNS` and to every row appended in
   `search_conformers` (near the existing `rdkit_version` at
   `conformers.py:523`). Also append `pver=<version> pcommit=<commit>` to the
   per-conformer **XYZ comment line** (same conformer output; ties the geometry
   file, not just the CSV, back to the code). Capture the `(version, commit)`
   pair **once per run** (like `rdkit_ver`), not per row.
4. **Do NOT touch** `gaussian.py` / `pubchem.py` logs this round (scope locked).
5. **Tests** (`tests/test_conformers.py`, offline):
   - `conformer_log.csv` has `pipeline_version` == `pipeline.__version__` on
     every row.
   - `pipeline_commit` column is present and is a `str` (may be empty) — do
     **not** assert a specific SHA (non-deterministic across environments).
   - the provenance helper returns `("<version>", "")` and does not raise when
     git is mocked absent/failing; returns a `.dirty`-suffixed commit when
     `git status --porcelain` is mocked non-empty.
   - the XYZ comment line contains the `pver=`/`pcommit=` tokens.

**Verify M-06 by:** `pytest tests/ -q` green including the cases above;
`python scripts/check_invariants.py` green; an offline run over a normal
molecule shows `pipeline_version` populated on every `conformer_log.csv` row and
the token in each XYZ header. Not "it ran" — open the CSV and an XYZ file and
read the columns/comment. Confirm the no-git branch yields `""` without error
(e.g. run the helper from a temp dir outside any repo).

### M-07 — populate and sync the canonical implementation-status.md

1. **Sync content.** Fill every section of `docs/implementation-status.md` (the
   §5 gate file) from the current `docs/implementation-status-v2.md`: what was /
   wasn't implemented, deviations (incl. the PubChem `SMILES` rename and the
   `generate_conformers` 4-tuple deviation), tests added, known limitations,
   the §6 scientific-judgment questions, and the Provenance block (now including
   `pipeline version` + commit from M-06). No `<...>` placeholders or empty
   bullets may remain.
2. **Record the M-06 change** in both status files (canonical + `-v2`) as this
   round's implementation, with its deviation note (new log columns; XYZ comment
   format extended) and the reproducibility caveat below.
3. **Make the sync durable, not a one-off.** Add a one-line header to
   `docs/implementation-status.md` stating it is synced each remediation round
   from `implementation-status-v2.md` (the working file), and add the sync step
   to `WORKFLOW.md` / the implementer checklist so future rounds don't
   re-introduce drift.
4. **Drift guard (in scope).** Add a lightweight check to
   `scripts/check_invariants.py` that fails when
   `docs/implementation-status.md` still contains template markers (`<...>`
   angle-bracket placeholders, `**PR:** #<n>`, empty `- ` / `<none | ...>`
   bullets). This converts "canonical went stale" from a reviewer catch into an
   objective floor failure (AGENTS.md §4), so it can never regress silently to
   the template again. Add a unit test asserting the guard fires on a template
   file and passes on a populated one.

**Verify M-07 by:** `docs/implementation-status.md` reads as the real current
status (grep shows **no** `<` placeholder tokens, no empty bullets); it agrees
with `-v2` on the PubChem rename and the v2 deviations; the Provenance block
names the pipeline version/commit. The new drift guard is exercised: reverting
the canonical file to the template makes `python scripts/check_invariants.py`
**fail** (confirm, then restore); with the populated file it is green, and
`pytest tests/ -q` (incl. the guard's unit test) passes.

## Reproducibility caveat (record in the status doc)

`pipeline_commit` pins code identity **only when the tree is clean**. A
`.dirty` suffix means uncommitted edits produced the output, so that output is
*not* fully reproducible from the commit alone — the flag makes this visible
rather than hiding it. Empty `pipeline_commit` (no git) falls back to
`pipeline_version`, which is only as precise as manual bumping. State this
explicitly so no one reads a recorded commit as a guarantee of exact code.

## Regression checks
- Normal molecules: the two new columns are **additive** — existing
  `conformer_log.csv` columns, ordering, ΔE values, and the top-3 selection are
  unchanged vs round 02.
- Offline test suite still passes with **no** git available (helper returns
  `""`); no test asserts a concrete SHA.
- The Link1 opt→freq contract and route lines are untouched (docs/log-only +
  additive columns; no route-line edits).

## Findings deliberately rejected
- None.

## Out of scope (next-round candidates)
- Pipeline version/commit in `com_write_log.csv` and `sdf_download_log.csv`
  (and the `.com` title line) — AGENTS.md §3 applies to *every* generated
  output, so expect Codex to raise these next unless deferred deliberately.

---
### Implementation evidence (fill in as fixes land — do NOT mark resolved on code change alone; run the verify step)
| ID | Commit | Verification run | Result |
|----|--------|------------------|--------|
| M-06 | `83cc915` | `pytest tests/ -q` → 91 passed (incl. `TestProvenanceLogging`, `TestGitShortSha`, `TestPipelineProvenance`); `check_invariants.py` → passed; offline run over ribose — `conformer_log.csv` `pipeline_version=0.2.0` and `pipeline_commit` populated on every row, XYZ header carries `pver=0.2.0 pcommit=…`, and `git_short_sha` from a non-repo temp dir returns `""` without error | **Resolved** |
| M-07 | `b98581f` | canonical `implementation-status.md` populated — grep shows **no** `<` placeholder tokens, no empty bullets; agrees with `-v2` on the PubChem rename and v2 deviations; Provenance block names version + commit; drift guard **fails on the template** (12 violations) and **passes on the populated file**; `pytest tests/ -q` → 95 passed (incl. `TestStatusDocDriftGuard`); `check_invariants.py` → passed | **Resolved** |
