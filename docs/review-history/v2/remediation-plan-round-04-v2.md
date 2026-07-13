# remediation-plan-round-04-v2.md

Produced by Cowork (Claude) from the PR #3 v2 release review. Approved by Ish
(gate 2) before fixes are implemented.

**Responds to:** PR #3 v2 release review (12 numbered findings + 3 soft items).
**Branch:** remediation commits land on the open PR branch
`feat/conformer-search-v2` (PR #3), as with rounds 01–03 — not a new branch.
**Base commit at triage:** `4d17713`.
**Context:** The conformer science is release-quality; the remaining defects let
the pipeline silently ship a truncated molecule to Gaussian, reuse a stale
geometry for a changed molecule, submit legacy v1 jobs, and emit SLURM scripts
that cannot find their inputs. Docs still advertise the retired Open Babel flow.

## Ish judgment calls locked (gate 2)

- **1b (B-02 resume key):** resume must invalidate on **CID, SMILES, and RDKit
  version** in addition to the current search knobs, and must confirm the
  referenced XYZ exists/non-empty. Rationale: ETKDGv3+MMFF geometry is
  RDKit-version-dependent, so reusing a cached geometry from a different RDKit
  silently breaks the reproducibility that `pipeline_version` claims. **Deferred
  as scope creep:** `pipeline_commit` in the key, `n_rows == n_kept` reconcile,
  and duplicate-label / key-by-identity — logged for a later round.
- **2a (M-02):** `preserve_unrequested` default → **False**; the conformer log
  represents the molecules requested this run. `append=True` stays available,
  non-default.
- **3a (scope):** the `runs/<study>/` directory redesign is **deferred** to a
  future architecture round. Not implemented here.
- **4a (M-03):** **overwrite** `.sh` files by default (cheap regeneration,
  emergent simplicity). No config-stamp-and-compare.
- **5a (MOD-01):** `__version__ = "2.0.0"`, release tag `v2.0.0`, PubChem
  User-Agent → `gaussian-input-pipeline/2.0`. This deliberately changes the
  `pipeline_version` string written into provenance logs going forward; recorded
  here so the change is not silent.

## Scope-expansion flags (not silently widened)

1. `runs/<study>/` run-directory redesign — **deferred** (call 3a).
2. Fuller resume key (`pipeline_commit`, row-count reconcile, duplicate-label
   rejection, identity-keying) — **deferred** (call 1b).
3. Extending version/commit provenance to `com_write_log.csv` /
   `sdf_download_log.csv` — still out of scope (carried from round 03).

Anything beyond the tasks below is out of scope for round 04.

## Decisions on each finding

| ID | Review # | Sev | Decision | Verify by |
|----|----------|-----|----------|-----------|
| B-01 | 7  | Blocker | Accept | XYZ parse tests below |
| B-02 | 4  | Blocker | Accept (call 1b) | resume-identity tests below |
| B-03 | 2  | Blocker | Accept | SLURM path test below |
| B-04 | 1  | Blocker | Accept | `git ls-files -ci` clean + fresh-clone count |
| M-01 | 3  | Major | Accept | COM-log discovery test |
| M-02 | 5  | Major | Accept (call 2a) | preserve-unrequested test |
| M-03 | 6  | Major | Accept (call 4a) | overwrite test |
| MOD-01 | 11 | Moderate | Accept (call 5a) | version + UA grep |
| MOD-02 | 8  | Moderate | Accept | canonical docs current, no placeholders |
| MOD-03 | 9  | Moderate | Accept | README primary flow = RDKit |
| MOD-04 | 10 | Moderate | Accept | WORKFLOW.md/AGENTS.md match reality |
| MOD-05 | 12 | Moderate | Accept | PR body matches final diff |
| MIN-01 | soft | Minor | Accept (+verify numpy) | pinned env documented |
| MIN-02 | soft | Minor | Accept | stale `*_failed.csv` cleared |
| MIN-03 | soft | Minor | Accept | invalid-param tests |

No findings rejected. All twelve + three were reproduced against the code at
`4d17713`.

---

## Blockers

### B-01 — XYZ blank-comment parser silently drops atoms
- **Evidence:** `pipeline/gaussian.py:32-35` — `lines = [ln.strip() for ln in
  f.readlines() if ln.strip()]` filters blank lines *before* `body = lines[2:]`.
  A valid XYZ with an empty comment line loses that line, so `lines[2:]` skips
  the atom count **and** the first atom.
- **Consequence:** Gaussian receives a molecule missing ≥1 atom, with no error.
  Silent scientific corruption of the input geometry — the most severe defect in
  the review.
- **Remediation:** Parse by physical line, not filtered line. Read line 1 as the
  atom count `N`; line 2 is the comment (may be empty); parse the next lines as
  coordinates. Validate that exactly `N` coordinate rows were read; raise
  `ValueError` on mismatch or malformed rows rather than dropping lines.
- **Verify by:** New tests in `tests/test_gaussian.py`: (a) empty-comment XYZ →
  all atoms present; (b) declared-count ≠ actual rows → raises; (c) existing
  well-formed cases still pass. `pytest tests/test_gaussian.py -q` green.

### B-02 — Resume reuses a stale geometry for a changed molecule
- **Evidence:** `pipeline/conformers.py:399` `_RESUME_CONFIG_FIELDS =
  ("seed","n_generate","top_n","rmsd_prune","pipeline_version")`;
  `_row_config_matches` (`:400-423`) compares only those. `cid`, `smiles`, and
  `rdkit_version` are recorded in `_LOG_COLUMNS` (`:376-390`) but never checked;
  `xyz_path` existence is never checked.
- **Consequence:** Same molecule `name` + same search knobs but a different
  CID/SMILES (a corrected structure) is treated as complete → the old geometry
  is silently reused. Chemistry-level wrong result.
- **Remediation (call 1b):**
  - Compare **`cid`** and **`smiles`** per molecule. These are per-molecule, not
    run-level: match each recorded row against the *requested* row's cid/smiles
    from `molecule_table`, not against a single run-level `config` dict. Update
    `_resume_partition` to pass the requested molecule's identity into the match.
  - Add **`rdkit_version`** to the invalidation set (run-level, like
    `pipeline_version`).
  - Confirm every recorded `xyz_path` **exists and is non-empty**; if any is
    missing, treat the molecule as stale → regenerate.
  - **Do not** add `pipeline_commit` to the key, reconcile `n_rows == n_kept`, or
    reject/re-key duplicate labels this round (deferred).
- **Verify by:** New `tests/test_conformers.py` cases — changed CID, changed
  SMILES, changed RDKit version, deleted/empty XYZ each force regeneration;
  unchanged identity + knobs still skips as complete.

### B-03 — SLURM scripts can't locate their `.com` files
- **Evidence:** `pipeline/slurm.py:34` template body is bare `g16 {jobname}.com`;
  submission from the parent dir (`for f in slurm_scripts/*.sh; do sbatch "$f"`)
  makes Gaussian look for `./{jobname}.com`, but inputs live in
  `gaussian_inputs/`.
- **Consequence:** "Ready-to-submit" scripts fail immediately on the cluster.
- **Remediation:** Pass the real `.com` path into the template and make the
  script submission-directory-independent. Compute the path with
  `os.path.relpath(com_path, slurm_dir)` (preserves custom dir names — do **not**
  hardcode `../gaussian_inputs`). In the script, resolve relative to the script's
  own location, `cd` there, and run `g16` on the basename.
- **Verify by:** New `tests/test_slurm.py` case creating sibling
  `gaussian_inputs/` + `slurm_scripts/` and asserting the generated script
  resolves the real input file path.

### B-04 — Tracked generated outputs contaminate a fresh clone
- **Evidence:** `git ls-files` tracks `pipeline/__pycache__/*.pyc`,
  `notebooks/.ipynb_checkpoints/`, `notebooks/.pubchem_cache/`,
  `notebooks/gaussian_inputs/*.com` (incl. legacy `adenine_F.com`),
  `notebooks/pubchem_sdf/`, `notebooks/pubchem_xyz/`, and the log/resolution
  CSVs. `.gitignore` names some of these but they are already tracked, so
  ignoring is a no-op.
- **Consequence:** A fresh clone ships legacy v1 `.com` files; the v2 notebook
  then reports 6 SLURM scripts when only 3 current jobs exist (review's
  fresh-clone simulation). Coupled with M-01 this means **legacy jobs get
  submitted**.
- **Remediation:**
  - `.gitignore` **already** ignores `__pycache__/`, `*.pyc`,
    `.ipynb_checkpoints/`, `.pubchem_cache/`, `gaussian_inputs/`,
    `slurm_scripts/`, `*_log.csv`, `com_write_failed.csv` — these are ineffective
    only because the files are already tracked. Do **not** duplicate them. Add
    the genuinely-missing patterns: `.claude/settings.local.json`,
    `conformer_xyz/`, `pubchem_sdf/`, `pubchem_xyz/`, `molecules_pubchem_*.csv`,
    and broaden `com_write_failed.csv` → `*_failed.csv`.
  - `git rm -r --cached` the tracked generated dirs/files (keep local copies):
    `pipeline/__pycache__`, `notebooks/.ipynb_checkpoints`,
    `notebooks/.pubchem_cache`, `notebooks/gaussian_inputs`,
    `notebooks/slurm_scripts`, `notebooks/pubchem_sdf`, `notebooks/pubchem_xyz`,
    and the tracked `notebooks/*_log.csv`, `molecules_pubchem_*.csv` files.
  - Add a CI hygiene check to `review-readiness.yml`:
    `test -z "$(git ls-files -ci --exclude-standard)"`.
- **Verify by:** `git ls-files -ci --exclude-standard` prints nothing; CI step
  present and passing; fresh clone contains no `*.com`/`*.sdf`/`*.xyz` under
  `notebooks/`.

---

## Major

### M-01 — SLURM stage globs the whole output dir instead of the current run
- **Evidence:** `pipeline/slurm.py:97` `com_files = sorted(glob.glob(
  os.path.join(com_dir, "*.com")))` — every `.com` on disk becomes a job.
- **Consequence:** Stale `.com` files from a previous molecule list are picked up
  and submitted even after the B-04 cleanup.
- **Remediation:** Default the SLURM stage to consume the current run's
  `com_write_log.csv` (interface e.g. `write_slurm_scripts(com_log_csv=...,
  slurm_dir=...)`). Keep the legacy `com_dir` glob as an explicit non-default
  mode. Notebook default switches to the log-driven path.
- **Verify by:** Test — a `com_write_log.csv` listing 3 jobs plus an extra stale
  `.com` on disk yields exactly 3 scripts; legacy glob mode still reachable.

### M-02 — Unrequested molecules preserved by default
- **Evidence:** `_resume_partition` (`conformers.py:445-446`) carries rows for
  molecules `name not in requested_names` into `kept_rows`.
- **Consequence:** Changing the molecule list never yields a clean current run;
  prior molecules accumulate in the log indefinitely and flow downstream.
- **Remediation (call 2a):** Add `preserve_unrequested: bool = False` (default
  drop) to the resume path; unrequested molecules are excluded from the new log.
  `append=True` opt-in retains the carry-forward behavior. **No** `runs/<study>/`
  redesign (call 3a).
- **Verify by:** Test — run {water, glycine}, then run {adenine}; default log
  contains only adenine; with `append=True` it contains all three.

### M-03 — Stale SLURM settings survive a rerun
- **Evidence:** `slurm.py` marks existing scripts `SKIPPED_EXISTS`
  (`:106`); a second call with a new `--account`/resources leaves the old script.
- **Consequence:** Jobs submit with stale SBATCH directives (wrong account,
  wrong walltime/resources).
- **Remediation (call 4a):** Overwrite `.sh` files by default (regeneration is
  cheap). Drop the `SKIPPED_EXISTS` default path; keep the log status column but
  report `WROTE`/`OVERWROTE`.
- **Verify by:** Test — write with `account="old"`, rewrite with `account="new"`;
  file contains `--account=new`, log reports an overwrite.

---

## Moderate — documentation, versioning, release hygiene

### MOD-01 — Version contradiction + stale PubChem UA
- **Evidence:** `pipeline/__init__.py:25` `__version__ = "0.2.0"`;
  `pipeline/pubchem.py:27` UA `gaussian-input-pipeline/1.0`.
- **Consequence:** A `v2.0.0` release writes `pipeline_version=0.2.0` into
  scientific provenance logs; UA misidentifies the client version.
- **Remediation (call 5a):** `__version__ = "2.0.0"`; UA →
  `gaussian-input-pipeline/2.0 (research use)`; release tag `v2.0.0`.
- **Verify by:** `grep __version__` = `2.0.0`; UA grep = `/2.0`; a fresh run
  records `pipeline_version=2.0.0`.

### MOD-02 — Canonical architecture/plan out of sync
- **Evidence:** `docs/architecture.md` still describes PubChem SDF → Open Babel,
  no conformer search, conformer search "out of scope";
  `docs/implementation-plan.md` is the untouched template placeholder while the
  real content lives in the `-v2` files.
- **Consequence:** The earlier Codex architecture finding is marked resolved on
  GitHub but not in the repo; the mandated merge-gate docs are wrong.
- **Remediation:** Make `docs/architecture.md`, `docs/implementation-plan.md`,
  `docs/implementation-status.md` the authoritative current (v2) versions. Move
  the `-v2` plan/status and remediation records under
  `docs/review-history/v2/` (preferred) or delete the duplicates. The canonical
  plan must hold the real v2 plan, not placeholders.
- **Verify by:** No template-placeholder strings in canonical files; canonical
  architecture describes the RDKit conformer pathway.

### MOD-03 — README leads with the retired Open Babel flow
- **Evidence:** README opens on the Open Babel pipeline as primary, omits
  `conformers.py` from the tree, lists Open Babel as required but omits RDKit,
  cites only v1.1, and doesn't explain stale-run behavior.
- **Consequence:** Public docs contradict the implemented default workflow.
- **Remediation:** Primary diagram = names → CID + stereo SMILES → RDKit ETKDGv3
  ensemble → MMFF94/UFF rank → top-N XYZ → Gaussian opt→freq → run-scoped SLURM.
  RDKit in requirements; `conformers.py` in the tree; Open Babel path demoted to
  a labeled **Legacy v1.1 workflow** section.
- **Verify by:** README primary flow is the RDKit path; RDKit listed; legacy
  section clearly labeled.

### MOD-04 — Stale agent-workflow docs
- **Evidence:** `WORKFLOW.md` references `.github/workflows/codex-review.yml`,
  `.github/codex/prompts/review.md`, and an API-key review path — all removed in
  PR #2. `AGENTS.md` still centers the v1 Open Babel architecture.
- **Consequence:** Contributors follow a review process and file layout that no
  longer exist.
- **Remediation:** Rewrite `WORKFLOW.md` to describe the **native Codex review**
  (no API key). Update `AGENTS.md` so the v2 conformer pathway is the primary
  architecture. Do **not** introduce any `OPENAI_API_KEY` / API-key path.
- **Verify by:** No references to removed files; native-review process described;
  no API-key language anywhere.

### MOD-05 — PR description stale
- **Evidence:** PR body says 48 tests (branch has 110), claims the notebook was
  not rewired (it is), ends with stray `EOF`/`)`.
- **Consequence:** Reviewer/merge context is inaccurate.
- **Remediation:** Rewrite the PR body **last**, after code+docs land, to match
  the final diff (test count, rewired notebook, changelog). Remove stray text.
- **Verify by:** PR body matches final diff; no stray tokens.

---

## Minor — recommended, not blockers

### MIN-01 — Dependency reproducibility
- **Evidence:** `environment.yml` is unpinned; RDKit version affects embedding.
  `numpy` has **no direct use in `pipeline/`** (grep clean at `4d17713`) — verify
  notebooks/tests before removing.
- **Remediation:** Ship a pinned/locked environment for v2.0.0 (env export,
  `conda-lock`, or documented exact versions). Remove `numpy` from explicit deps
  only after confirming it's unused outside `pipeline/`.
- **Verify by:** Pinned env artifact present/documented; if numpy removed, full
  suite + notebook still run.

### MIN-02 — Stale failure logs
- **Evidence:** A prior run's `*_failed.csv` survives a later successful run.
- **Remediation:** Clear/overwrite stale `*_failed.csv` at the start of the
  relevant stage so the notebook doesn't surface old failures.
- **Verify by:** Test — a failing run then a clean run leaves no stale
  `*_failed.csv` content.

### MIN-03 — Early parameter validation
- **Evidence:** No guard on `n_generate < 1`, `top_n < 1`, `rmsd_prune < 0`,
  duplicate molecule labels, empty sanitized filename.
- **Remediation:** Validate at entry; raise `ValueError` with a clear message.
- **Verify by:** Tests asserting each invalid input raises.

---

## Recommended implementation sequence

**Commit 1 — correctness hardening (blockers/majors):** B-01, B-02, B-03, M-01,
M-02, M-03, plus MIN-03.
**Commit 2 — repository cleanup:** B-04 (untrack + `.gitignore` + CI hygiene
check), MIN-02.
**Commit 3 — docs + versioning:** MOD-01..MOD-04, MIN-01. MOD-05 (PR body) last.

## Final verification gate (must run — code changing ≠ finding resolved)

```bash
pytest tests/ -q                                   # expect all green
python scripts/check_invariants.py                 # expect passed
test -z "$(git ls-files -ci --exclude-standard)"   # ignored-but-tracked → empty
grep -n '__version__' pipeline/__init__.py         # 2.0.0
grep -n 'User-Agent' pipeline/pubchem.py           # /2.0
```

Then run the notebook from a genuinely clean clone and confirm: exactly the
requested molecules in `conformer_log.csv` → same in `com_write_log.csv` → same
in `slurm_write_log.csv`; every SLURM script resolves its `.com`; changing
molecule list / alias / CID / SMILES / account / resources cannot reuse stale
output. Finally trigger a fresh `@codex review` against the final head. Only a
clean re-review clears gate 3 (merge + tag `v2.0.0`).
