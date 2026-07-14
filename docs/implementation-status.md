# implementation-status.md

> Canonical merge-gate status. Reviewers must verify every claim against code
> and tests. Historical v2.0 remediation evidence remains under
> `docs/review-history/v2/`.

**PR:** (opened for v2.1 on branch `feat/pipeline-v2.2`)
**Branch:** `feat/pipeline-v2.2`
**Current phase:** targeted PR #4 remediation implemented; awaiting verification,
push, and Ish's final approval
**Merge status:** local checks green (pytest, check_invariants, git diff --check,
ignored-tracked-files check); requires push, green remote CI, an independent
review, and Ish's final approval

This status is governed by `docs/architecture-v2.1.md`,
`docs/implementation-plan-v2.1.md`, and the normative
`docs/release-contract-v2.1.md` (a deliberate, Ish-approved revision of the
frozen `docs/release-contract-v2.0.md`).

## 1. Ordered task completion (v2.1)

### Task 1 — release contract + AGENTS.md: complete

`docs/release-contract-v2.1.md` is written from the v2.0 contract with the
architecture "Contract deltas" applied: run-folder layout (section 2), reduced
per-artifact matrix (section 5), resume/append removal (section 7), and the
amended judgment 1 plus new judgment 5 (section 9). `AGENTS.md` section 3 now
references v2.1 and drops the resume/dirty-git reuse bullet.

### Task 2 — trimmed per-artifact metadata: complete

- XYZ line 2 carries `dE=... method=... artifact_id=xyz-...` only (plus the
  unconverged/provisional markers where applicable). The manifest-only fields
  (run_id, config_hash, conformer_id, pipeline_version, rdkit_version, the
  separate relative_energy field) are removed from the comment.
- The Gaussian COM title is a single line: name, level-of-theory suffix, `dE`,
  `artifact_id`. The entire second `provenance ...` line is deleted.
- The SLURM header carries `artifact_id` and `source_com=BASENAME sha256=...`
  only; run_id and the source-COM relative path are manifest-only.
- `scripts/check_invariants.py` is updated to the reduced matrix and now asserts
  the removed tokens are absent from each artifact body.

### Task 3 — SH simplification: complete

`DEFAULT_TEMPLATE` drops the SCRIPT_DIR/COM_PATH/cd path-resolution block and
ends with `ml gaussian16` then `g16 JOBNAME.com`. The script assumes
the COM is in the current working directory. `write_slurm_script` no longer
threads a COM relative path.

### Task 4 — per-run directory layout: complete

The pipeline stage functions already resolve every output path relative to the
manifest package root, so the run-folder layout is realized by rooting the
manifest at `runs/STUDY/RUN_ID/run_manifest.json` and placing `conformer_xyz/`
and `gaussian_jobs/` beside it. The notebook adds a `STUDY` field, generates the
immutable `RUN_ID` up front, and writes COM and SH into the same
`gaussian_jobs/` directory. `create_run_manifest`'s `FileExistsError` guard and
all section 6 collision checks (distinct labels to one basename, one-to-one
COM to SH mapping) are preserved and now evaluated within the run folder.

### Task 5 — resume/append removed (B-04 deleted by construction): complete

`search_conformers` no longer accepts an `append` parameter or runs any resume
path. Deleted: `_resume_partition`, `_row_config_matches`, `_commit_key`,
`_row_manifest_matches`, `_row_identity_matches`, `_row_xyz_present`,
`_integer_key`, `_resume_group_is_complete`, `_group_identity_is_consistent`,
`_carry_forward_group_is_valid`, `_RESUME_CONFIG_FIELDS`, the canonical-log
repair helpers (`_write_conformer_log_atomically`, `_conformer_log_matches_rows`),
and every resume/append/stale/dirty-git/no-git reuse test. `search_conformers`
now raises when pointed at an already-populated run folder (manifest already
holds conformer records or artifacts). The invariant checker's append-integrity
section is replaced by a no-resume guard that fails if a resume helper or the
`append` parameter is reintroduced or if the populated-run guard is removed.

**B-04 disposition:** deferred-then-removed, not lost. The v2.0 finding lived
entirely on the resume path (a failed regeneration into a dirty/blank-commit
manifest could resurrect stale conformers). With no resume path that code is
unreachable and has been deleted. `tests/test_conformers.py` proves it:
`search_conformers` exposes no `append` parameter and a populated run folder
raises rather than reusing.

### Task 6 — provisional undefined-stereo structure (judgment 5): complete

`search_conformers` no longer skips undefined stereo. It detects the undefined
centre(s) via `FindPotentialStereo`, embeds one structure from the PubChem
IsomericSMILES (`generate_provisional_conformer`: single ETKDG embed, fixed
seed, light MMFF/UFF cleanup), reads back the post-embed `arbitrated_smiles`,
and emits one `BASE_c00.xyz` marked PROVISIONAL with `dE=NA`. The manifest
molecule record and `conformer_log.csv` carry `provenance_status`,
`undefined_centers`, `pubchem_smiles`, and `arbitrated_smiles`; the manifest
validator requires the arbitrated SMILES to differ from the PubChem SMILES. The
COM title and XYZ comment carry the `PROVISIONAL: stereo arbitrated at ...`
marker; downstream Gaussian and SLURM consume the row unchanged via the one
`provenance_status` column (no special-case branch). A loud console warning is
emitted, and for k undefined centres greater than one it states that 2^k
isomers are possible with one arbitrated. Verified on D-Ribose (CID 10975657):
one PROVISIONAL XYZ, `provenance_status=provisional_undefined_stereo` in log and
manifest, `undefined_centers` recorded, `arbitrated_smiles != pubchem_smiles`,
`dE=NA` in artifacts, and a COM+SH pair generated carrying the marker.

Honesty guardrail: the arbitrated structure is an unvalidated DFT starting
geometry, one arbitrary pick among 2^k, and is never presented as the defined
stereoisomer — stated in code, the contract, the notebook, and the README. If
MMFF params are missing, UFF is used for the cleanup and the fallback is logged
(never a silent chemistry substitution).

### Task 7 — notebook, README, status doc: complete

The notebook adds the `STUDY` config field, roots all outputs under
`runs/STUDY/RUN_ID/`, writes COM and SH together in `gaussian_jobs/`, documents
the HPC submission contract (submit the `.sh` from the directory holding its
`.com`), and describes the provisional behavior honestly. The README updates the
provenance/reproducibility section, the run-folder archive/submission steps, the
resume-removal note (with B-04 disposition), and the provisional undefined-stereo
caveat.

### PR #4 remediation — complete

- The conformer stage rejects existing conformer logs, failure logs, or a
  nonempty `conformer_xyz/` before it creates staging state or publishes output;
  untracked final XYZ destinations hard-fail before publication mutation.
- Manifest-driven Gaussian COM and SLURM defaults both use `gaussian_jobs/`, so
  the default SLURM command finds its source COM beside the script.
- Gaussian preflight requires conformer-log provisional provenance fields to
  match the authoritative manifest molecule exactly before a COM/log/manifest
  mutation; provisional COMs therefore retain `dE=NA` and the visible marker.
- Manifest validation requires an undefined-stereo provisional molecule to carry
  exactly one conformer, including for hand-edited or damaged manifests.
- Manifest-driven SLURM generation rejects COM/SH directory separation before
  mutation; direct v2 COM writes derive provisional metadata from the manifest;
  and all three conformer-stage logs/directories are package-contained.
- The bundled default SLURM template intentionally uses Purdue RCAC/Lmod syntax,
  `ml gaussian16`; other clusters may provide a site-specific custom template.

## 2. Scientific invariants (unchanged)

Level of theory, route lines, charge/multiplicity, coordinate/energy units, and
the Link1 opt→freq checkpoint contract (`Geom=AllChk Guess=Read`) are unchanged.
No section 2 invariant was touched; no deviation was required.

## 3. Required checks

Locally green with the pinned stack (pandas 3.0.3, rdkit 2025.09.3,
pytest 9.1.1, requests 2.34.2) on CPython 3.12:

- `pytest tests/ -q` — all tests pass
- `python scripts/check_invariants.py` — passes
- `git diff --check` — clean
- `test -z "$(git ls-files -ci --exclude-standard)"` — no ignored tracked files

The clean-`git archive` repeat of the same checks is recorded at PR time.

## 4. Pending / deferred / rejected

- **Pending:** push, green remote CI, independent Codex review, Ish's final
  merge approval.
- **Deferred (explicitly out of scope for v2.1):** stereoisomer enumeration
  (both anomers), xTB/CREST rerank, solvent-aware search, energy-window/Boltzmann
  logic, running or parsing Gaussian, and submit-from-anywhere SH robustness
  (deliberately dropped for the clean co-located script).
- **Rejected:** none.
