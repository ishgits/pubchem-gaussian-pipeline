# Architecture

## System purpose

The pipeline converts molecule names into ranked conformer starting geometries,
Gaussian Link1 opt/freq inputs, and SLURM scripts. It does not run Gaussian.

```text
names
  -> PubChem CID + stereo-bearing SMILES
  -> RDKit ETKDGv3 ensemble
  -> MMFF94 optimization/ranking, UFF recorded fallback
  -> top-N carried conformers
  -> conformer XYZ artifacts
  -> Gaussian COM artifacts
  -> SLURM SH artifacts
```

Rigid molecules naturally collapse to fewer conformers through embed-time RMSD
pruning. There is no rotatable-bond gate.

The Open Babel single-geometry pathway (`pipeline/geometry.py`) is deprecated
v1.1 compatibility functionality. It is not the default and is exempt from the
provenance model below.

## Design principle: the manifest is authoritative

`run_manifest.json` is the complete provenance record for a run: the full
conformer, Gaussian, and SLURM configuration, resolved molecular identity,
artifact lineage, relative paths, and SHA-256 file hashes.

Because the manifest carries everything, artifact bodies carry only what a
person or tool reading that specific file actually uses, plus one stable
back-pointer:

| Artifact | Comment/title carries |
|---|---|
| **XYZ** line 2 | `dE=<kcal/mol> method=<MMFF94\|UFF> artifact_id=xyz-…` |
| **COM** title | `<name> PCM 298 K <basis> dE=<kcal/mol> artifact_id=com-…` |
| **SH** header | `# artifact_id=sh-…`, `# source_com=<basename> sha256=<…>` |

`dE` and `method` stay in the XYZ because they are inline-useful science. The
SH `source_com` hash stays because it lets a job confirm it is running the
intended COM. Everything else — run ID, config hash, software versions, full
search configuration — is a manifest lookup via `artifact_id`. See
[`provenance.md`](provenance.md) for the full model.

## Per-run directories

Every run writes one immutable directory:

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

- `<study>` is a user-supplied label; `<run_id>` is an opaque immutable UUID.
- COM and SH are co-located in `gaussian_jobs/`, so a SLURM script needs no
  path-resolution machinery: submit the `.sh` from the directory that holds
  its `.com`, and it runs `g16 <base>_F.com` from there.
- The complete run folder, not an isolated artifact, is the supported archive
  and transfer unit.

## No resume or append

Every run is fresh and immutable. There is no resume/append path: pointing the
pipeline at an already-populated run folder raises rather than reusing it, and
adding molecules means starting a new run. Recomputing a conformer search
costs seconds of ETKDG against hours of downstream DFT, so the simplicity is
worth more than the reuse.

## Provisional structures for undefined stereochemistry

Molecules whose PubChem SMILES leaves stereochemistry unspecified are not
skipped. The pipeline embeds **one** provisional structure (a single ETKDG
embed with a light MMFF/UFF cleanup, not an ensemble search), in which RDKit
necessarily fixes the undefined centre(s) to an arbitrary configuration.

The provisional structure is flagged everywhere it surfaces:
`provenance_status=provisional_undefined_stereo`, `undefined_centers`,
`pubchem_smiles`, and the post-embed `arbitrated_smiles` in the manifest and
`conformer_log.csv`; a `PROVISIONAL: stereo arbitrated at …` marker and
`dE=NA` in the XYZ/COM; and a console warning.

The arbitrated structure is an unvalidated DFT starting geometry — one
arbitrary pick among 2^k for k undefined centres — and must never be treated
as the defined stereoisomer. Stereoisomer enumeration is out of scope.

## Deliberately out of scope

xTB/CREST reranking; solvent-aware conformer search; energy-window/Boltzmann
selection; running or parsing Gaussian; submit-from-anywhere SLURM scripts
(dropped in favour of the clean co-located script).

## Version history

- **v1.1** — single Open Babel geometry per molecule (now the legacy path).
- **v2.0** — RDKit conformer search; top-3 distinct conformers to DFT;
  manifest-centric provenance with full per-artifact metadata; resume/append;
  flat output directories.
- **v2.1** — minimal per-artifact metadata (one `artifact_id` back-pointer);
  per-run `runs/<study>/<run_id>/` directories with co-located COM+SH;
  resume/append removed; provisional starting structure for undefined-stereo
  molecules.
