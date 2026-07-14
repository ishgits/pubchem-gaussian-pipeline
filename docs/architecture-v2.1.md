# architecture-v2.1.md

> Human gate: Ish approves this before implementation. Supersedes the frozen
> `release-contract-v2.0.md` where noted. v2.1 is a deliberate revision of the
> v2.0 provenance/artifact contract, prompted by user testing of the merged v2.0
> pipeline (PR #3).
> Scope: cleaner artifacts, per-run directories, no resume machinery, and a
> provisional path for undefined-stereo molecules. No change to level of theory,
> route lines, charge/multiplicity, or the Link1 opt→freq contract.

## Why v2.1 (what testing exposed)

The v2.0 pipeline works and the config cells are good, but real use surfaced four
things, three of which trace back to the v2.0 contract being over-built for the
actual goal — *reduce manual structure prep so we can screen hundreds of
molecules*:

1. **Artifacts are cluttered.** XYZ/COM/SH bodies carry the full provenance
   matrix mandated by contract §5. In practice this is noise the manifest already
   records; it hurts the UX with no scientific gain.
2. **`run_manifest.json` "won't overwrite"** blocks re-runs in a shared working
   directory. The immutability guard is correct; the flat-directory layout is the
   problem.
3. **D-Ribose (and every sugar/anomer) silently drops out.** Correct per v2.0
   judgment #1 (skip undefined stereo), but at scale it means whole compound
   classes yield nothing to run QC on.
4. **B-04 (carried, unresolved):** on a resume into a dirty/blank-commit manifest,
   a failed regeneration leaves stale conformers presented as valid, defeating the
   §7 "dirty disables reuse" policy. This lives entirely on the resume/append
   path.

## Design principle

The manifest is already the authoritative, complete provenance record (v2.0
§1/§4). Everything duplicated into an artifact body is redundant. So: **an
artifact carries one stable back-pointer into the manifest plus only the
information a human or tool reading that specific file actually uses.** One key,
not seven. Everything else is a manifest lookup away.

Two structural simplifications fall out of that principle and the "one run = one
folder" model, and together they delete the most fragile code in the repo:

- Per-run output directories remove the overwrite friction *and* co-locate COM+SH.
- Dropping resume/append removes `_resume_partition`, the manifest-driven
  canonical-log repair, and stale detection — which is exactly where B-04 lives.
  B-04 is not fixed; the code that carries it is deleted.

## Change 1 — Minimal per-artifact metadata (revises contract §5)

The manifest keeps carrying the full matrix. Artifact bodies keep only inline
science + one `artifact_id`.

| Artifact | Comment/title now carries | Removed (manifest-only) |
|---|---|---|
| **XYZ** line 2 | `dE=<kcal/mol> method=<MMFF94\|UFF> artifact_id=xyz-…` | run_id, config_hash, pipeline_version, rdkit_version |
| **COM** title | one line: `<name> PCM 298 K <basis> dE=<kcal/mol> artifact_id=com-…` | the entire second `provenance … commit=…` line; config_hash, run_id, versions |
| **SH** header | `# artifact_id=sh-…`, `# source_com=<basename> sha256=<…>` | run_id as a separate field |

Notes:

- `dE` stays kcal/mol and labeled; `method` stays because a chemist opening the
  XYZ needs to know MMFF vs UFF. These are inline-useful science, not provenance.
- The SH `source_com` sha256 stays — it is the one operationally load-bearing
  header (lets the job confirm it is running the intended COM).
- `artifact_id` is the single reverse link; from it the manifest yields run_id,
  config_hash, versions, and full search configuration.

## Change 2 — SH simplification (revises contract §5 SH block)

COM and SH are placed in the **same directory** (see Change 3). The SH therefore
drops all path-resolution machinery — the `SCRIPT_DIR`/`COM_PATH` block and its
explanatory comment (the old "B-03" logic) are removed entirely. The script
assumes the `.com` is in the current working directory:

```bash
#!/bin/bash
# artifact_id=sh-…
# source_com=<base>_F.com sha256=<…>
#SBATCH --account=…
#SBATCH --job-name=<base>_F
# … resources …
ml gaussian16
g16 <base>_F.com
```

The bundled default uses Purdue RCAC/Lmod syntax (`ml gaussian16`). A different
cluster may provide a custom template with its site-specific module command.

**Operational contract (must be documented in the notebook next-steps cell):** on
the HPC, submit the `.sh` from the directory that contains its `.com`. The COM+SH
pair travels together; the script does not hunt for its input. This is a
deliberate trade — we drop submit-from-anywhere robustness in exchange for a
clean, obvious script, and the per-run `gaussian_jobs/` directory makes "they're
together" the default.

## Change 3 — Per-run directories (activates deferred v2.0 §2 redesign)

Replace flat `conformer_xyz/ gaussian_inputs/ slurm_scripts/` with one immutable
directory per run:

```text
runs/<study>/<run_id>/
  run_manifest.json
  conformer_log.csv
  com_write_log.csv
  slurm_write_log.csv
  conformer_search_failed.csv        # only if failures
  conformer_xyz/
    <base>_c00.xyz …
  gaussian_jobs/                      # COM + SH co-located
    <base>_c00_F.com
    <base>_c00_F.sh
```

- `<study>` is a user-supplied label in the config cell (e.g. `nucleobases`); the
  `run_id` remains the opaque immutable UUID.
- `create_run_manifest`'s `FileExistsError` guard stays as a safety net but never
  fires in normal use, because every run writes a fresh directory. This resolves
  the "won't overwrite" friction structurally rather than by weakening
  immutability.
- The complete run package (v2.0 §2) is now literally the run folder — cleaner to
  archive and transfer.
- Adding molecules later = a new run in a new `run_id/` folder. Two clean
  single-purpose runs are better provenance than a resume path.

## Change 4 — Drop resume/append; delete B-04's home (removes contract §7)

Every run is fresh and immutable. Remove the entire reuse subsystem:

- `search_conformers(..., append=...)` — the `append` parameter and all append
  validation.
- `_resume_partition`, `_row_config_matches`, `_commit_key`, the
  `_RESUME_CONFIG_KEYS` machinery.
- The manifest-driven canonical-log repair block (current conformers.py
  ~1150–1159) and `_manifest_conformer_log_rows(verify_xyz=True)` reuse path.
- Every stale/dirty-git/no-git *reuse* test and the invariant checker's reuse
  assertions.

**B-04 disposition:** removed by construction, not patched. With no resume path,
the stale-conformer resurrection is unreachable. The verify-by (below) proves the
path is gone, not merely that a branch changed.

This is an explicit scope reduction: it deletes a large, correct-but-fragile body
of code whose only purpose was to avoid recomputing conformer searches, which are
seconds of ETKDG against hours of downstream DFT.

## Change 5 — Provisional structure for undefined stereo (new judgment #5, amends #1)

Today undefined-stereo molecules (D-Ribose: undefined anomeric C in the PubChem
pyranose SMILES) are skipped with no output. v2.1 instead produces **one
provisional starting structure**, loudly labeled:

- Detect the undefined stereocentre(s) via `FindPotentialStereo` (as now) and
  record which atom(s) are unspecified.
- Embed **one** conformer from the PubChem `IsomericSMILES` (single ETKDG embed,
  fixed seed) and apply a **light MMFF/UFF cleanup** — explicitly *not* an
  ensemble conformer search.
- RDKit necessarily fixes the undefined centre(s) to an arbitrary configuration
  at embed time. Read back the post-embed isomeric SMILES and record it as the
  *arbitrated* structure, distinct from the underspecified PubChem SMILES.
- Mark the molecule PROVISIONAL everywhere it surfaces:
  - manifest molecule record: `provenance_status = provisional_undefined_stereo`,
    plus `undefined_centers`, `pubchem_smiles`, `arbitrated_smiles`;
  - `conformer_log.csv`: a `provenance_status` column (`normal` vs
    `provisional_undefined_stereo`);
  - XYZ comment and COM title: a visible `PROVISIONAL: stereo arbitrated at
    <atoms>` marker, and `dE=NA` (a single structure has no ensemble reference —
    we do not fabricate `dE=0.000`);
  - a loud console warning: conformer search skipped for undefined stereo; single
    structure embedded from PubChem SMILES with an arbitrary choice at the
    undefined centre(s).

**Honesty guardrails (invariant-critical):**

- The arbitrated structure is an unvalidated DFT starting geometry, **not** the
  compound's real configuration, and is one arbitrary pick among 2^k for k
  undefined centres. It must never be presented as the defined stereoisomer.
- With k > 1 the warning states "k undefined centres → 2^k isomers possible; one
  arbitrated." We still emit exactly one structure; enumeration of stereoisomers
  is deferred (not in v2.1).
- Provisional molecules flow through gaussian/slurm unchanged, carrying the
  marker — emergent behavior via one status column, no special-case downstream
  branch.

Judgment #1 is amended accordingly: undefined stereo is no longer skipped; it is
provisionally embedded and loudly flagged, never silently guessed.

## Contract deltas (for the versioned `release-contract-v2.1.md`)

- **§2** — flat directories replaced by `runs/<study>/<run_id>/`; the run folder
  is the archive/transfer unit; COM+SH co-located in `gaussian_jobs/`.
- **§5** — per-artifact matrix reduced to the table in Change 1; SH loses path
  machinery per Change 2.
- **§7** — resume/append removed. Replaced by: *every run is fresh and immutable;
  to add molecules, start a new run.* The dirty/no-git language is retired with
  the feature it governed.
- **§9** — judgment #1 amended (provisional embed, not skip); new judgment #5
  (provisional undefined-stereo structure, per Change 5).
- **§6 collision rules** — unchanged in intent; now evaluated within a single
  run folder. Distinct labels sanitizing to one basename still fail loudly before
  mutation.

## Verify-by (each change must be proven, not assumed)

- **Change 1/2:** a test asserts the exact new XYZ line-2 / COM title / SH header
  strings and asserts the removed fields are absent; `check_invariants.py` updated
  to the reduced matrix.
- **Change 3:** a run writes `runs/<study>/<run_id>/…` with COM+SH in
  `gaussian_jobs/`; re-running writes a new `run_id/` and never touches the prior.
- **Change 4 (B-04):** a test asserts `search_conformers` has no `append`
  parameter and no resume branch, and that pointing a run at an existing
  populated run folder raises rather than reusing. The old dirty-git-resume
  regression scenario is removed with the feature.
- **Change 5:** a test on D-Ribose (CID 10975657) asserts one PROVISIONAL XYZ is
  produced, `provenance_status=provisional_undefined_stereo` in log+manifest,
  `undefined_centers` recorded, `arbitrated_smiles` differs from `pubchem_smiles`,
  `dE=NA` in artifacts, and a COM+SH pair is generated carrying the marker.

## Explicitly out of scope for v2.1

Stereoisomer enumeration (both anomers); xTB/CREST rerank; solvent-aware search;
energy-window/Boltzmann logic; running or parsing Gaussian; any change to route
lines, basis set, solvent, charge/multiplicity, or the Link1 opt→freq checkpoint
contract; submit-from-anywhere SH robustness (deliberately dropped).

## Change log

- v2.0 — simple RDKit conformer search; top-3 distinct conformers to DFT;
  manifest-centric provenance with full per-artifact metadata; resume/append with
  strict git-clean reuse; flat output directories.
- v2.1 — minimal per-artifact metadata; per-run `runs/<study>/<run_id>/`
  directories with co-located COM+SH; resume/append removed (deletes B-04);
  provisional starting structure for undefined-stereo molecules.
