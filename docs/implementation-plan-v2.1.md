# implementation-plan-v2.1.md

> Human gate: Ish approves this before implementation. Works against
> `docs/architecture-v2.1.md`. This plan deliberately revises the frozen v2.0
> contract; the revision is the point, not a violation. Do not expand scope
> beyond what is listed here without recording a deviation.

## Objective

Ship v2.1: minimal per-artifact metadata, per-run `runs/<study>/<run_id>/`
directories with co-located COM+SH, removal of resume/append (which deletes
B-04), and a provisional starting structure for undefined-stereo molecules. No
change to level of theory, route lines, charge/multiplicity, or the Link1
opt→freq contract.

## First: write the contract, then code to it

1. **`docs/release-contract-v2.1.md`** — copy `release-contract-v2.0.md` and apply
   the "Contract deltas" section of `architecture-v2.1.md` (§2 layout, §5 reduced
   matrix, §7 resume removed, §9 judgments #1 amended + #5 added). Update
   `AGENTS.md` §3 to reference v2.1 and drop the resume/dirty-git reuse bullet in
   §3. This file is normative; everything below is tested against it.

## Tasks (ordered, each independently verifiable)

2. **Trim per-artifact metadata.**
   - XYZ line 2 → `dE=<kcal/mol> method=<MMFF94|UFF> artifact_id=xyz-…`
     (`pipeline/conformers.py` `_write_xyz` comment builder).
   - COM title → single line `<name> PCM 298 K <basis> dE=<kcal/mol>
     artifact_id=com-…`; delete the second `provenance …` line
     (`pipeline/gaussian.py` title construction).
   - SH header → `# artifact_id=sh-…` and `# source_com=<basename> sha256=<…>`
     only (`pipeline/slurm.py`).
   - Verify-by: tests assert the exact new strings and assert removed fields
     (run_id, config_hash, pipeline_version, rdkit_version, the provenance line)
     are absent. Update `scripts/check_invariants.py` to the reduced matrix.

3. **Simplify the SH script.**
   - Remove the `SCRIPT_DIR`/`COM_PATH` resolution block and its comment; the
     script ends with `ml gaussian16` then `g16 <base>_F.com`, assuming
     the COM is in the cwd.
   - Verify-by: generated SH contains no `SCRIPT_DIR`/`cd`/`../`; a test asserts
     the two-line tail and that `source_com`/`sha256` header is intact.

4. **Per-run directory layout.**
   - Add a `study` config field and route all outputs under
     `runs/<study>/<run_id>/` with `conformer_xyz/` and `gaussian_jobs/`
     (COM+SH co-located) as in `architecture-v2.1.md` Change 3.
   - Keep the `create_run_manifest` `FileExistsError` guard.
   - Preserve all §6 collision checks, now evaluated within the run folder
     (distinct labels → one basename still fails before mutation; COM→SH mapping
     stays one-to-one).
   - Verify-by: a run creates the full tree; COM and SH for a conformer sit in the
     same `gaussian_jobs/` dir; a second run creates a new `run_id/` and leaves the
     first byte-for-byte unchanged.

5. **Remove resume/append (deletes B-04).**
   - Delete the `append` parameter and all append validation from
     `search_conformers`; delete `_resume_partition`, `_row_config_matches`,
     `_commit_key`, `_RESUME_CONFIG_KEYS`, the canonical-log repair block, and the
     `_manifest_conformer_log_rows(verify_xyz=True)` reuse path.
   - Remove the corresponding resume/append/stale/dirty-git/no-git *reuse* tests
     and the invariant checker's reuse assertions.
   - Verify-by (B-04, must run): a test asserts `search_conformers` exposes no
     `append` path and that pointing a run at an existing populated run folder
     raises rather than reusing stale conformers. This is the resolution of B-04 —
     the path is gone, verified by test, not merely refactored.
   - Scope flag: this is an intentional feature removal and a large test deletion.
     Record it in `implementation-status.md` as deferred-then-removed, not lost.

6. **Provisional undefined-stereo structure (judgment #5).**
   - In `conformers.py`, replace the skip-only branch: when
     `check_conformer_eligibility` finds undefined stereo, record the undefined
     atom(s), embed one conformer from the PubChem `IsomericSMILES`, apply a light
     MMFF/UFF cleanup, read back the post-embed `arbitrated_smiles`, and emit one
     `<base>_c00.xyz` marked PROVISIONAL with `dE=NA`.
   - Add `provenance_status` to `conformer_log.csv` (`normal` |
     `provisional_undefined_stereo`) and to the manifest molecule record along with
     `undefined_centers`, `pubchem_smiles`, `arbitrated_smiles`.
   - COM title and XYZ comment carry `PROVISIONAL: stereo arbitrated at <atoms>`;
     downstream gaussian/slurm consume the row unchanged, carrying the marker.
   - Loud console warning per `architecture-v2.1.md` Change 5; k>1 centres →
     "2^k isomers possible, one arbitrated."
   - Honesty guardrail in code + docs: arbitrated structure is an unvalidated DFT
     start, never the compound's real configuration.
   - Verify-by: D-Ribose (CID 10975657) test — one PROVISIONAL XYZ,
     `provenance_status=provisional_undefined_stereo` in log+manifest,
     `undefined_centers` recorded, `arbitrated_smiles != pubchem_smiles`, `dE=NA`,
     and a COM+SH pair generated with the marker.

7. **Notebook + docs.**
   - Config cell: add `study`. Next-steps markdown: state the HPC submission
     contract — submit the `.sh` from the directory holding its `.com`
     (they ship together in `gaussian_jobs/`); document the provisional-structure
     behavior and its PROVISIONAL labeling honestly.
   - README caveats: add the provisional undefined-stereo behavior; note resume was
     removed in favor of one-run-per-folder.
   - Write `docs/implementation-status.md` distinguishing implemented / pending /
     rejected / deferred, including the resume removal and B-04 disposition.

## Required checks (AGENTS.md §5 — must be green before PR)

```bash
pytest tests/ -q
python scripts/check_invariants.py
git diff --check
test -z "$(git ls-files -ci --exclude-standard)"
```

Re-run the same checks from a clean `git archive` before requesting final review;
record whether the pinned dependency stack was used.

## Files expected to change

- `pipeline/conformers.py` (metadata trim; remove resume subsystem; provisional
  branch; `provenance_status`)
- `pipeline/gaussian.py` (single-line title; provisional marker; `dE=NA` path)
- `pipeline/slurm.py` (minimal SH header + body)
- `pipeline/manifest.py` (run-dir paths; molecule record fields; remove reuse
  helpers)
- `scripts/check_invariants.py` (reduced matrix; drop reuse assertions)
- `notebooks/run_pipeline.ipynb` (`study` config + next-steps markdown)
- `tests/…` (rewrite artifact-format + provisional tests; remove resume tests)
- `docs/release-contract-v2.1.md` (new), `docs/architecture-v2.1.md`,
  `docs/implementation-plan-v2.1.md`, `docs/implementation-status.md`
- `AGENTS.md` (§3 reference to v2.1), `README.md`

## Explicitly out of scope (do not add)

Stereoisomer enumeration; xTB/CREST; solvent-aware search; energy-window/Boltzmann
logic; running/parsing Gaussian; any route/basis/solvent/charge/multiplicity or
Link1 change; submit-from-anywhere SH robustness.

## Escalate, don't decide

If MMFF params are missing for a provisional molecule, record whether UFF-fallback
was used or the cleanup was skipped — do not silently substitute chemistry. If any
task appears to require touching a scientific invariant (§2), stop and record a
deviation rather than proceeding.
