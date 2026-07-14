"""
Gaussian input file (.com) generation from XYZ coordinates.

Supports the Link1 pattern for combined opt → freq calculations in a
single submission (optimization writes checkpoint, frequency job reads it
via Geom=AllChk Guess=Read).
"""

from __future__ import annotations

import glob
import os

import pandas as pd

from .conformers import UNCONVERGED_FF_SEED
from .manifest import (
    assert_stage_configuration,
    find_artifact,
    find_conformer_record,
    record_child_artifact,
    relative_artifact_path,
    remove_artifacts_by_kind,
    require_exact_artifact_id_set,
    sha256_file,
    stable_record_id,
)
from .utils import ensure_dir, parse_strict_bool, sanitize_basename


# Fixed schemas keep a scientifically valid zero-job run machine-readable for
# the next stage (M-11). The legacy writer has no conformer identifier; the v2
# writer preserves it for traceability.
_LEGACY_COM_LOG_COLUMNS = ["name", "xyz_path", "com_path"]
_CONFORMER_COM_LOG_COLUMNS = [
    "run_id",
    "artifact_id",
    "config_hash",
    "name",
    "conformer_id",
    "conformer_record_id",
    "xyz_artifact_id",
    "xyz_path",
    "com_path",
    "com_sha256",
    "pipeline_version",
    "pipeline_commit",
    "rdkit_version",
]
_REQUIRED_CONFORMER_PROVENANCE_COLUMNS = (
    "run_id",
    "artifact_id",
    "config_hash",
    "pipeline_version",
    "rdkit_version",
    "xyz_sha256",
)


def _optional_text(value) -> str | None:
    """Normalize an optional scalar/CSV field to non-empty text or ``None``."""
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text or None




def _parse_converged_flag(value, *, row_index: int) -> bool:
    """Parse one manifest-linked convergence flag without permissive fallback."""
    return parse_strict_bool(
        value, field_name=f"Conformer row {row_index} converged"
    )

def _validate_required_conformer_provenance(conf_log: pd.DataFrame) -> None:
    """Reject nonempty v2 logs missing source-version provenance (M-16)."""
    if conf_log.empty:
        return

    missing_columns = [
        column
        for column in _REQUIRED_CONFORMER_PROVENANCE_COLUMNS
        if column not in conf_log.columns
    ]
    if missing_columns:
        raise ValueError(
            "Nonempty conformer log is missing required provenance column(s): "
            + ", ".join(missing_columns)
        )

    problems = []
    for column in _REQUIRED_CONFORMER_PROVENANCE_COLUMNS:
        bad_rows = [
            int(index)
            for index, value in conf_log[column].items()
            if _optional_text(value) is None
        ]
        if bad_rows:
            problems.append(f"{column} missing at row(s) {bad_rows}")

    if problems:
        raise ValueError(
            "Nonempty conformer log has missing required provenance: "
            + "; ".join(problems)
        )


def _validate_direct_conformer_provenance(
    conformer_id: int | None,
    pipeline_version,
    rdkit_version,
    run_id,
    artifact_id,
    config_hash,
    manifest_path,
    parent_artifact_id,
    conformer_record_id,
) -> None:
    """Require source versions for direct conformer-specific COM writes (M-17)."""
    if conformer_id is None:
        return

    missing = []
    if _optional_text(pipeline_version) is None:
        missing.append("pipeline_version")
    if _optional_text(rdkit_version) is None:
        missing.append("rdkit_version")
    if _optional_text(run_id) is None:
        missing.append("run_id")
    if _optional_text(artifact_id) is None:
        missing.append("artifact_id")
    if _optional_text(config_hash) is None:
        missing.append("config_hash")
    if _optional_text(manifest_path) is None:
        missing.append("manifest_path")
    if _optional_text(parent_artifact_id) is None:
        missing.append("parent_artifact_id")
    if _optional_text(conformer_record_id) is None:
        missing.append("conformer_record_id")
    if missing:
        raise ValueError(
            "Conformer-specific Gaussian inputs require nonblank provenance: "
            + ", ".join(missing)
        )


def xyz_to_gaussian_coords(xyz_path: str) -> str:
    """
    Read an XYZ file and return the coordinate block formatted for a
    Gaussian input file.

    XYZ format expected::

        <atom_count>
        <comment line>
        Element  x  y  z
        ...

    Parsing is by **physical line**, never by "non-blank line" (B-01): line 1 is
    the atom count ``N``, line 2 is the comment (which may legitimately be empty),
    and the next ``N`` lines are coordinates. Filtering blank lines first would
    drop an empty comment line and silently shift the atom count / first atom out
    of the geometry — corrupting the molecule sent to Gaussian. Any count mismatch
    or malformed coordinate row raises ``ValueError`` rather than dropping atoms.
    """
    with open(xyz_path, "r") as f:
        raw_lines = f.read().splitlines()

    if len(raw_lines) < 2:
        raise ValueError(
            f"XYZ file {xyz_path!r} is too short: expected an atom-count line, a "
            f"comment line, then coordinates (got {len(raw_lines)} line(s))."
        )

    count_str = raw_lines[0].strip()
    try:
        n_atoms = int(count_str)
    except ValueError:
        raise ValueError(
            f"XYZ file {xyz_path!r} line 1 is not an integer atom count: "
            f"{count_str!r}."
        )
    if n_atoms < 1:
        raise ValueError(
            f"XYZ file {xyz_path!r} declares a non-positive atom count: {n_atoms}."
        )

    # Line 2 is the comment (may be empty); everything after it is coordinates.
    # Only purely-trailing blank lines are tolerated (a trailing newline is
    # normal); any other count mismatch raises rather than dropping/padding atoms.
    coord_lines = raw_lines[2:]
    while coord_lines and coord_lines[-1].strip() == "":
        coord_lines.pop()
    if len(coord_lines) != n_atoms:
        raise ValueError(
            f"XYZ file {xyz_path!r} declares {n_atoms} atom(s) but {len(coord_lines)} "
            f"coordinate row(s) are present (declared count ≠ actual rows)."
        )

    out_lines = []
    for i, ln in enumerate(coord_lines, 1):
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(
                f"XYZ file {xyz_path!r} coordinate row {i} is malformed "
                f"(need 'Element x y z'): {ln!r}."
            )
        sym = parts[0]
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            raise ValueError(
                f"XYZ file {xyz_path!r} coordinate row {i} has non-numeric "
                f"coordinates: {ln!r}."
            )
        out_lines.append(f"{sym:<2} {x:>16.8f} {y:>12.8f} {z:>12.8f}")
    return "\n".join(out_lines)


def write_gaussian_com(
    name: str,
    xyz_path: str,
    outdir: str,
    route_opt: str,
    route_freq: str,
    title_suffix: str = "",
    charge: int = 0,
    multiplicity: int = 1,
    nproc: int = 16,
    link1: bool = True,
    conformer_id: int | None = None,
    rel_energy_kcalmol: float | None = None,
    unconverged: bool = False,
    pipeline_version: str | None = None,
    pipeline_commit: str | None = None,
    rdkit_version: str | None = None,
    run_id: str | None = None,
    artifact_id: str | None = None,
    config_hash: str | None = None,
    manifest_path: str | None = None,
    parent_artifact_id: str | None = None,
    conformer_record_id: str | None = None,
    provenance_status: str = "normal",
    undefined_centers: str | None = None,
) -> str:
    """
    Write a Gaussian .com input file from an XYZ file.

    Parameters
    ----------
    name : str
        Molecule label (used for filenames and title line).
    xyz_path : str
        Path to the .xyz coordinate file.
    outdir : str
        Directory to write the .com file into.
    route_opt : str
        Gaussian route line for the optimization job
        (e.g., ``"# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)"``).
    route_freq : str
        Gaussian route line for the frequency job
        (e.g., ``"# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) Geom=AllChk Guess=Read"``).
    title_suffix : str
        Appended to the title line (e.g., ``"PCM 298 K 6-311++G(2df,2p)"``).
    charge : int
        Molecular charge.
    multiplicity : int
        Spin multiplicity.
    nproc : int
        Number of processors (%nprocshared).
    link1 : bool
        If True, append a --Link1-- section for the frequency job that reads
        geometry from the checkpoint file.
    conformer_id : int, optional
        Conformer index (v2 conformer stage). When given, the basename becomes
        ``{base}_c{ii}`` so each conformer gets its own ``.com``/``.chk`` pair
        (e.g. ``ribose_c00_F.com``), and the id is recorded in the title line for
        traceability. Nonblank ``pipeline_version`` and ``rdkit_version`` are
        required for this v2 path (M-17). When ``None`` the v1.1 single-geometry
        naming and optional-provenance behavior are preserved.
    rel_energy_kcalmol : float, optional
        Force-field ΔE (kcal/mol) of this conformer relative to the molecule's
        lowest-energy conformer. Recorded in the title line for traceability.
        Explicitly labeled kcal/mol so it is never mixed with DFT Hartree values.
    unconverged : bool
        If True, the starting geometry came from an FF optimization that did NOT
        converge (M-04 decision 2b best-effort seed). An ``UNCONVERGED_FF_SEED``
        marker is written into the title line so the unminimized start — and its
        unreliable FF energy — are visible on inspection.
    pipeline_version : str, optional
        Pipeline version that produced the conformer. Required when
        ``conformer_id`` is supplied; otherwise optional. When supplied, recorded
        on a separate ``provenance`` line in the Gaussian title section (M-14/M-17).
    pipeline_commit : str, optional
        Source commit that produced the conformer. If missing while another
        provenance field is supplied, the title records ``commit=unavailable``.
    rdkit_version : str, optional
        RDKit version that generated/ranked the starting geometry. Required when
        ``conformer_id`` is supplied; otherwise optional. When supplied, recorded
        on the title-section provenance line.

    Returns
    -------
    str
        Path to the written .com file.
    """
    # M-17: `conformer_id` selects the v2 scientific-output path even for direct
    # calls that bypass the validated batch writer. Require its source versions
    # before creating a directory or writing any file. Legacy v1.1 calls keep
    # optional provenance exactly as before.
    unconverged_value = unconverged
    if conformer_id is not None:
        unconverged_value = parse_strict_bool(
            unconverged, field_name="Direct COM unconverged"
        )

    _validate_direct_conformer_provenance(
        conformer_id,
        pipeline_version,
        rdkit_version,
        run_id,
        artifact_id,
        config_hash,
        manifest_path,
        parent_artifact_id,
        conformer_record_id,
    )
    pipeline_version = _optional_text(pipeline_version)
    pipeline_commit = _optional_text(pipeline_commit)
    rdkit_version = _optional_text(rdkit_version)
    run_id = _optional_text(run_id)
    artifact_id = _optional_text(artifact_id)
    config_hash = _optional_text(config_hash)
    manifest_path = _optional_text(manifest_path)
    parent_artifact_id = _optional_text(parent_artifact_id)
    conformer_record_id = _optional_text(conformer_record_id)

    if conformer_id is not None:
        manifest = assert_stage_configuration(
            manifest_path,
            "gaussian",
            {
                "route_opt": route_opt,
                "route_freq": route_freq,
                "title_suffix": title_suffix,
                "charge": charge,
                "multiplicity": multiplicity,
                "nproc": nproc,
                "link1": link1,
            },
        )
        if run_id != manifest["run_id"] or config_hash != manifest["config_hash"]:
            raise ValueError("Direct COM identity disagrees with run manifest.")
        if pipeline_version != manifest["pipeline_version"] or rdkit_version != manifest["rdkit_version"]:
            raise ValueError("Direct COM software versions disagree with run manifest.")
        if (pipeline_commit or "") != manifest["pipeline_commit"]:
            raise ValueError("Direct COM source commit disagrees with run manifest.")
        parent = find_artifact(manifest, parent_artifact_id)
        if parent["kind"] != "xyz" or parent.get("conformer_record_id") != conformer_record_id:
            raise ValueError("Direct COM parent lineage is invalid.")
        if relative_artifact_path(xyz_path, manifest_path) != parent["relative_path"]:
            raise ValueError("Direct COM XYZ path disagrees with run manifest.")
        if sha256_file(xyz_path) != parent["sha256"]:
            raise ValueError("Direct COM XYZ hash disagrees with run manifest.")
        expected_artifact_id = stable_record_id(
            manifest["run_id"], "com", parent_artifact_id
        )
        if artifact_id != expected_artifact_id:
            raise ValueError("Direct COM artifact ID is not stable for its lineage.")
        molecule_record, conformer_record = find_conformer_record(
            manifest, conformer_record_id
        )
        if str(name) != molecule_record["molecule_name"]:
            raise ValueError("Direct COM molecule name disagrees with run manifest.")
        if int(conformer_id) != conformer_record["conformer_id"]:
            raise ValueError("Direct COM conformer ID disagrees with run manifest.")
        if rel_energy_kcalmol is None or abs(
            float(rel_energy_kcalmol)
            - float(conformer_record["relative_energy_kcalmol"])
        ) > 1e-9:
            raise ValueError("Direct COM relative energy disagrees with run manifest.")
        if unconverged_value == conformer_record["converged"]:
            raise ValueError("Direct COM convergence marker disagrees with run manifest.")

    base = sanitize_basename(name)
    if conformer_id is not None:
        # Extend the basename per conformer (architecture v2): {base}_c{ii}.
        base = f"{base}_c{conformer_id:02d}"
    chk_name = f"{base}_F.chk"
    com_path = os.path.join(outdir, f"{base}_F.com")

    if conformer_id is not None:
        # B-01: preflight the COM destination against the run package before any
        # filesystem mutation. relative_artifact_path() raises ValueError on an
        # out-of-package path, so the direct v2 entry point fails before
        # ensure_dir()/open() can leave an external artifact behind.
        relative_artifact_path(com_path, manifest_path)

    ensure_dir(outdir)

    coords = xyz_to_gaussian_coords(xyz_path)

    # v2.1 per-artifact metadata (contract §5, architecture Change 1): the COM
    # title is a SINGLE line carrying only inline science plus the one
    # artifact_id back-pointer. run_id/config_hash/conformer_id/versions and the
    # former second `provenance …` line are manifest-only now. This must never
    # alter route lines, checkpoint directives, charge/multiplicity, coordinates,
    # or the Link1 frequency section.
    is_provisional = provenance_status == "provisional_undefined_stereo"
    title = f"{base} {title_suffix}".strip()
    if is_provisional:
        # A single arbitrated structure has no ensemble reference; never fabricate
        # dE=0.000 (contract §9 honesty guardrail).
        title = f"{title} dE=NA".strip()
    elif rel_energy_kcalmol is not None:
        # Units labeled explicitly; FF energy, never a DFT Hartree value.
        title = f"{title} dE={rel_energy_kcalmol:.4f} kcal/mol".strip()
    if artifact_id:
        title = f"{title} artifact_id={artifact_id}".strip()
    if unconverged_value:
        # Make the unconverged FF start explicit on the input itself (M-04 2b).
        title = f"{title} {UNCONVERGED_FF_SEED}".strip()
    if is_provisional:
        centers = undefined_centers if undefined_centers else "unspecified center(s)"
        title = f"{title} PROVISIONAL: stereo arbitrated at {centers}".strip()
    title_block = title

    text = (
        f"%nprocshared={nproc}\n"
        f"%chk={chk_name}\n"
        f"{route_opt}\n\n"
        f"{title_block}\n\n"
        f"{charge} {multiplicity}\n"
        f"{coords}\n\n"
    )

    if link1:
        text += (
            f"--Link1--\n"
            f"%nprocshared={nproc}\n"
            f"%chk={chk_name}\n"
            f"{route_freq}\n\n"
        )

    with open(com_path, "w") as f:
        f.write(text)

    if conformer_id is not None:
        record_child_artifact(
            manifest_path,
            kind="com",
            artifact_id=artifact_id,
            parent_artifact_id=parent_artifact_id,
            conformer_record_id=conformer_record_id,
            path=com_path,
        )

    return com_path


def write_gaussian_coms(
    xyz_log_csv: str,
    outdir: str = "gaussian_inputs",
    log_csv: str = "com_write_log.csv",
    **kwargs,
) -> pd.DataFrame:
    """
    Batch-write Gaussian .com files for every XYZ in *xyz_log_csv*.

    All keyword arguments are forwarded to :func:`write_gaussian_com`.
    """
    # Clear any stale failure log from a prior run (MIN-02); rewritten below only
    # if this run actually has failures.
    if os.path.exists("com_write_failed.csv"):
        os.remove("com_write_failed.csv")

    xyz_log = pd.read_csv(xyz_log_csv)
    written = []
    failed = []

    for _, row in xyz_log.iterrows():
        name = row["name"]
        xyz_path = row["xyz_path"]
        try:
            com_path = write_gaussian_com(name, xyz_path, outdir=outdir, **kwargs)
            written.append({"name": name, "xyz_path": xyz_path, "com_path": com_path})
        except Exception as e:
            failed.append({"name": name, "xyz_path": xyz_path, "error": repr(e)})

    out_df = pd.DataFrame(written, columns=_LEGACY_COM_LOG_COLUMNS)
    out_df.to_csv(log_csv, index=False)

    if failed:
        fail_df = pd.DataFrame(failed)
        fail_df.to_csv("com_write_failed.csv", index=False)
        print(f"WARNING: {len(failed)} .com writes failed — see com_write_failed.csv")
    else:
        print("All Gaussian .com files written successfully.")

    print(f"Wrote: {log_csv}")
    return out_df


def write_gaussian_coms_from_conformers(
    conformer_log_csv: str,
    outdir: str = "gaussian_inputs",
    log_csv: str = "com_write_log.csv",
    manifest_path: str = "run_manifest.json",
    **kwargs,
) -> pd.DataFrame:
    """
    Batch-write one Gaussian ``.com`` per conformer from a ``conformer_log.csv``
    (the v2 conformer stage output; multiple rows per molecule).

    Each row must carry ``name``, ``xyz_path``, and ``conformer_id``; the ΔE
    column (``rel_energy_kcalmol``) is recorded in the title line when present.
    A ``converged`` column (M-04), when present and False, tags the title with
    ``UNCONVERGED_FF_SEED`` so an unminimized best-effort start is visible.
    Every nonempty log row must contain nonblank ``pipeline_version`` and
    ``rdkit_version`` values (M-16); they describe conformer generation and are
    never inferred from the current Gaussian-writer environment.
    ``pipeline_version``, ``pipeline_commit``, and ``rdkit_version`` are copied
    from each row into the COM title section and the COM write log (M-14); a
    missing commit is recorded in the COM as ``commit=unavailable``.
    Files are written as ``{base}_c{ii}_F.com``. The Link1 opt→freq checkpoint
    contract is unchanged — every keyword argument is forwarded to
    :func:`write_gaussian_com`, exactly as :func:`write_gaussian_coms` does.
    """
    conf_log = pd.read_csv(conformer_log_csv)
    # M-16: validate source-version provenance before deleting a stale failure
    # log, creating an output directory, or writing any COM/log file. These
    # versions belong to conformer generation and cannot be backfilled safely
    # from the environment running this downstream stage.
    _validate_required_conformer_provenance(conf_log)

    gaussian_config = {
        "route_opt": kwargs.get("route_opt"),
        "route_freq": kwargs.get("route_freq"),
        "title_suffix": kwargs.get("title_suffix", ""),
        "charge": kwargs.get("charge", 0),
        "multiplicity": kwargs.get("multiplicity", 1),
        "nproc": kwargs.get("nproc", 16),
        "link1": kwargs.get("link1", True),
    }
    manifest = assert_stage_configuration(
        manifest_path, "gaussian", gaussian_config
    )
    # M-30: the COM output root must stay inside the manifest package. Validate it
    # before any lineage removal, directory creation, or COM write so an outdir
    # outside the package fails before prior COM/SH lineage is removed.
    relative_artifact_path(outdir, manifest_path)

    # Validate all source linkage, hashes, and destination mappings before the
    # first directory/log/manifest mutation.  A conformer-derived COM has no
    # supported v2 meaning without its exact XYZ manifest parent.
    prepared = []
    seen_rows = set()
    seen_destinations = set()
    for index, row in conf_log.iterrows():
        if str(row["run_id"]) != manifest["run_id"] or str(row["config_hash"]) != manifest["config_hash"]:
            raise ValueError(f"Conformer row {int(index)} disagrees with manifest identity.")
        xyz_artifact_id = str(row["artifact_id"])
        if xyz_artifact_id in seen_rows:
            raise ValueError(f"Duplicate XYZ artifact record in conformer log: {xyz_artifact_id!r}.")
        seen_rows.add(xyz_artifact_id)
        xyz_artifact = find_artifact(manifest, xyz_artifact_id)
        if xyz_artifact["kind"] != "xyz":
            raise ValueError(f"Conformer row {int(index)} does not reference an XYZ artifact.")
        molecule_record, conformer_record = find_conformer_record(
            manifest, xyz_artifact["conformer_record_id"]
        )
        if str(row["name"]) != molecule_record["molecule_name"]:
            raise ValueError(f"Conformer row {int(index)} molecule name disagrees with manifest.")
        row_cid = None if pd.isna(row.get("cid")) else int(float(row.get("cid")))
        row_smiles = "" if pd.isna(row.get("smiles")) else str(row.get("smiles"))
        if row_cid != molecule_record["CID"] or row_smiles != molecule_record["IsomericSMILES"]:
            raise ValueError(f"Conformer row {int(index)} molecule identity disagrees with manifest.")
        if int(row["conformer_id"]) != conformer_record["conformer_id"]:
            raise ValueError(f"Conformer row {int(index)} ID disagrees with manifest.")
        converged = _parse_converged_flag(row.get("converged"), row_index=int(index))
        if converged is not conformer_record["converged"]:
            raise ValueError(
                f"Conformer row {int(index)} convergence disagrees with manifest."
            )
        rel_energy = row.get("rel_energy_kcalmol")
        if pd.isna(rel_energy) or abs(
            float(rel_energy) - float(conformer_record["relative_energy_kcalmol"])
        ) > 1e-9:
            raise ValueError(f"Conformer row {int(index)} energy disagrees with manifest.")
        for field in ("pipeline_version", "pipeline_commit", "rdkit_version"):
            row_value = _optional_text(row.get(field)) or ""
            if row_value != str(manifest[field]):
                raise ValueError(
                    f"Conformer row {int(index)} {field} disagrees with manifest."
                )
        xyz_path = str(row["xyz_path"])
        if relative_artifact_path(xyz_path, manifest_path) != xyz_artifact["relative_path"]:
            raise ValueError(f"Conformer row {int(index)} XYZ path disagrees with manifest.")
        if str(row["xyz_sha256"]) != xyz_artifact["sha256"] or sha256_file(xyz_path) != xyz_artifact["sha256"]:
            raise ValueError(f"Conformer row {int(index)} XYZ hash disagrees with manifest.")
        conformer_id = int(row["conformer_id"])
        base = f"{sanitize_basename(str(row['name']))}_c{conformer_id:02d}"
        com_path = os.path.join(outdir, f"{base}_F.com")
        normalized_destination = os.path.normcase(os.path.abspath(com_path))
        if normalized_destination in seen_destinations:
            raise ValueError(f"Duplicate Gaussian destination path: {com_path!r}.")
        seen_destinations.add(normalized_destination)
        com_artifact_id = stable_record_id(
            manifest["run_id"], "com", xyz_artifact_id
        )
        prepared.append((row, xyz_artifact, com_path, com_artifact_id, converged))

    # A stage CSV is a subordinate index, not an independent source of truth.
    # After validating each present row, reject a valid-looking subset or extra
    # row before removing prior COM/SH lineage or touching any output.  Empty
    # zero-job logs remain valid when the manifest contains no XYZ artifacts.
    observed_xyz_ids = (
        conf_log["artifact_id"].tolist()
        if "artifact_id" in conf_log.columns
        else []
    )
    require_exact_artifact_id_set(
        manifest,
        "xyz",
        observed_xyz_ids,
        source_label="conformer_log.csv",
    )

    # M-30: the authoritative COM write log must also stay inside the package.
    # Validate it after all source/artifact-set preflight but before the first
    # mutation, so a log_csv outside the package fails before any failure-log
    # deletion, COM/SH lineage removal, COM write, or log rewrite.
    relative_artifact_path(log_csv, manifest_path)

    # Clear any stale failure log from a prior run (MIN-02); rewritten below only
    # if this run actually has failures.
    if os.path.exists("com_write_failed.csv"):
        os.remove("com_write_failed.csv")

    # Regenerating COMs replaces the downstream layers in manifest order.
    remove_artifacts_by_kind(manifest_path, "sh")
    remove_artifacts_by_kind(manifest_path, "com")

    written = []
    failed = []

    for row, xyz_artifact, expected_com_path, com_artifact_id, converged in prepared:
        name = row["name"]
        xyz_path = row["xyz_path"]
        conformer_id = int(row["conformer_id"])
        rel_e = row.get("rel_energy_kcalmol")
        rel_e = None if pd.isna(rel_e) else float(rel_e)
        unconverged = not converged
        pipeline_version = _optional_text(row.get("pipeline_version"))
        pipeline_commit = _optional_text(row.get("pipeline_commit"))
        rdkit_version = _optional_text(row.get("rdkit_version"))
        run_id = str(row["run_id"])
        config_hash = str(row["config_hash"])
        # v2.1: the provisional undefined-stereo marker rides through on the
        # conformer-log row's provenance_status column — no special-case branch.
        provenance_status = _optional_text(row.get("provenance_status")) or "normal"
        undefined_centers = _optional_text(row.get("undefined_centers"))
        try:
            com_path = write_gaussian_com(
                name,
                xyz_path,
                outdir=outdir,
                conformer_id=conformer_id,
                rel_energy_kcalmol=rel_e,
                unconverged=unconverged,
                pipeline_version=pipeline_version,
                pipeline_commit=pipeline_commit,
                rdkit_version=rdkit_version,
                run_id=run_id,
                artifact_id=com_artifact_id,
                config_hash=config_hash,
                manifest_path=manifest_path,
                parent_artifact_id=xyz_artifact["artifact_id"],
                conformer_record_id=xyz_artifact["conformer_record_id"],
                provenance_status=provenance_status,
                undefined_centers=undefined_centers,
                **kwargs,
            )
            if os.path.normcase(os.path.abspath(com_path)) != os.path.normcase(os.path.abspath(expected_com_path)):
                raise ValueError("Gaussian writer destination changed after validation.")
            com_digest = sha256_file(com_path)
            written.append({
                "run_id": run_id,
                "artifact_id": com_artifact_id,
                "config_hash": config_hash,
                "name": name,
                "conformer_id": conformer_id,
                "conformer_record_id": xyz_artifact["conformer_record_id"],
                "xyz_artifact_id": xyz_artifact["artifact_id"],
                "xyz_path": xyz_path,
                "com_path": com_path,
                "com_sha256": com_digest,
                "pipeline_version": pipeline_version,
                "pipeline_commit": pipeline_commit or "",
                "rdkit_version": rdkit_version,
            })
        except Exception as e:
            failed.append({
                "name": name,
                "conformer_id": conformer_id,
                "xyz_path": xyz_path,
                "error": repr(e),
            })

    out_df = pd.DataFrame(written, columns=_CONFORMER_COM_LOG_COLUMNS)
    out_df.to_csv(log_csv, index=False)

    if failed:
        fail_df = pd.DataFrame(failed)
        fail_df.to_csv("com_write_failed.csv", index=False)
        print(f"WARNING: {len(failed)} .com writes failed — see com_write_failed.csv")
    else:
        print("All Gaussian .com files written successfully.")

    print(f"Wrote: {log_csv}")
    return out_df
