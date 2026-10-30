# Provenance and reproducibility

This page states the guarantees a run from this pipeline gives you, in plain
terms. It is distilled from the v2.1 design; [`architecture.md`](architecture.md)
explains how the pieces fit together.

## The manifest is the record

`run_manifest.json` is the authoritative, complete record of one run:

- the complete conformer, Gaussian, and SLURM configuration, hashed as a
  canonical `config_hash` (identical configuration → identical hash;
  scientifically relevant changes alter it; ordering-only differences do not);
- each molecule's resolved identity (PubChem CID and stereo-bearing SMILES);
- every artifact's lineage — molecule → conformer → XYZ → COM → SH — with
  relative paths and SHA-256 file hashes.

The stage CSVs (`conformer_log.csv`, `com_write_log.csv`, `slurm_write_log.csv`)
are operational indexes. They must agree with the manifest; where they differ,
the manifest wins.

## Artifacts carry one back-pointer

XYZ, COM, and SH files carry only inline science plus a single `artifact_id`
that links them to their manifest record (see the table in
[`architecture.md`](architecture.md)). They deliberately do not duplicate the
full search configuration — that lives once, in the manifest. An isolated
artifact is identifiable through its ID, but it is not independently
reproducible without the matching manifest.

## The archive unit is the run folder

Archive and transfer the complete `runs/<study>/<run_id>/` folder. That folder
— manifest, stage logs, `conformer_xyz/`, and `gaussian_jobs/` with COM+SH
co-located — is the smallest unit the reproducibility guarantee covers.

## Failures happen before mutation

A run fails loudly, before writing or overwriting anything, when:

- distinct molecule labels sanitize to the same output basename;
- two source artifacts map to one destination path, or two COM files map to
  one SLURM script;
- a required source artifact is missing, blank, or zero bytes;
- duplicate `artifact_id` values appear;
- an artifact's path or hash disagrees with the manifest.

Filenames are never auto-disambiguated to work around a collision.

## Runs are immutable; there is no resume

Every run is fresh. Pointing the pipeline at an already-populated run folder
raises rather than reusing, appending to, or repairing it. To add molecules,
start a new run — two clean single-purpose runs are better provenance than a
resume path.

## Undefined stereochemistry is flagged, never guessed

A molecule whose PubChem SMILES leaves stereochemistry unspecified produces
one provisional starting structure: the undefined centre(s) are recorded as
`undefined_centers`, the post-embed arbitrated SMILES is recorded alongside
the original PubChem SMILES, and the XYZ/COM carry a
`PROVISIONAL: stereo arbitrated at …` marker with `dE=NA`. That structure is
an unvalidated starting geometry — one arbitrary pick among 2^k for k
undefined centres — not the compound's real configuration.

## Legacy path

The v1.1 Open Babel single-geometry path (`pipeline/geometry.py`) predates
this model and is exempt from it. Do not treat its outputs as carrying the
guarantees above.
