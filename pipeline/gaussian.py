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
from .utils import ensure_dir, sanitize_basename


# Fixed schemas keep a scientifically valid zero-job run machine-readable for
# the next stage (M-11). The legacy writer has no conformer identifier; the v2
# writer preserves it for traceability.
_LEGACY_COM_LOG_COLUMNS = ["name", "xyz_path", "com_path"]
_CONFORMER_COM_LOG_COLUMNS = [
    "name",
    "conformer_id",
    "xyz_path",
    "com_path",
    "pipeline_version",
    "pipeline_commit",
    "rdkit_version",
]


def _optional_text(value) -> str | None:
    """Normalize an optional scalar/CSV field to non-empty text or ``None``."""
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text or None


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
        traceability. When ``None`` the v1.1 single-geometry naming is preserved.
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
        Pipeline version that produced the conformer. When supplied, recorded on
        a separate ``provenance`` line in the Gaussian title section (M-14).
    pipeline_commit : str, optional
        Source commit that produced the conformer. If missing while another
        provenance field is supplied, the title records ``commit=unavailable``.
    rdkit_version : str, optional
        RDKit version that generated/ranked the starting geometry. When supplied,
        recorded on the title-section provenance line.

    Returns
    -------
    str
        Path to the written .com file.
    """
    ensure_dir(outdir)
    base = sanitize_basename(name)
    if conformer_id is not None:
        # Extend the basename per conformer (architecture v2): {base}_c{ii}.
        base = f"{base}_c{conformer_id:02d}"
    chk_name = f"{base}_F.chk"
    com_path = os.path.join(outdir, f"{base}_F.com")

    coords = xyz_to_gaussian_coords(xyz_path)
    title = f"{base} {title_suffix}".strip()
    if rel_energy_kcalmol is not None:
        # Units labeled explicitly; FF energy, never a DFT Hartree value.
        title = f"{title} dE={rel_energy_kcalmol:.4f} kcal/mol".strip()
    if unconverged:
        # Make the unconverged FF start explicit on the input itself (M-04 2b).
        title = f"{title} {UNCONVERGED_FF_SEED}".strip()

    # M-14: keep software provenance self-contained in the Gaussian title
    # section. It must never alter route lines, checkpoint directives,
    # charge/multiplicity, coordinates, or the Link1 frequency section.
    pipeline_version = _optional_text(pipeline_version)
    pipeline_commit = _optional_text(pipeline_commit)
    rdkit_version = _optional_text(rdkit_version)
    provenance_parts = []
    if pipeline_version:
        provenance_parts.append(f"pipeline={pipeline_version}")
    if pipeline_commit:
        provenance_parts.append(f"commit={pipeline_commit}")
    elif pipeline_version or rdkit_version:
        provenance_parts.append("commit=unavailable")
    if rdkit_version:
        provenance_parts.append(f"rdkit={rdkit_version}")
    title_lines = [title]
    if provenance_parts:
        title_lines.append("provenance " + " ".join(provenance_parts))
    title_block = "\n".join(title_lines)

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
    **kwargs,
) -> pd.DataFrame:
    """
    Batch-write one Gaussian ``.com`` per conformer from a ``conformer_log.csv``
    (the v2 conformer stage output; multiple rows per molecule).

    Each row must carry ``name``, ``xyz_path``, and ``conformer_id``; the ΔE
    column (``rel_energy_kcalmol``) is recorded in the title line when present.
    A ``converged`` column (M-04), when present and False, tags the title with
    ``UNCONVERGED_FF_SEED`` so an unminimized best-effort start is visible.
    ``pipeline_version``, ``pipeline_commit``, and ``rdkit_version`` are copied
    from each row into the COM title section and the COM write log (M-14); a
    missing commit is recorded in the COM as ``commit=unavailable``.
    Files are written as ``{base}_c{ii}_F.com``. The Link1 opt→freq checkpoint
    contract is unchanged — every keyword argument is forwarded to
    :func:`write_gaussian_com`, exactly as :func:`write_gaussian_coms` does.
    """
    # Clear any stale failure log from a prior run (MIN-02); rewritten below only
    # if this run actually has failures.
    if os.path.exists("com_write_failed.csv"):
        os.remove("com_write_failed.csv")

    conf_log = pd.read_csv(conformer_log_csv)
    written = []
    failed = []

    for _, row in conf_log.iterrows():
        name = row["name"]
        xyz_path = row["xyz_path"]
        conformer_id = int(row["conformer_id"])
        rel_e = row.get("rel_energy_kcalmol")
        rel_e = None if pd.isna(rel_e) else float(rel_e)
        # Missing/NaN converged column → assume converged (backward compatible).
        # Handle both native-bool and CSV string ("True"/"False") representations.
        conv = row.get("converged")
        if conv is None or (isinstance(conv, float) and pd.isna(conv)):
            unconverged = False
        elif isinstance(conv, str):
            unconverged = conv.strip().lower() in ("false", "0", "no")
        else:
            unconverged = not bool(conv)
        pipeline_version = _optional_text(row.get("pipeline_version"))
        pipeline_commit = _optional_text(row.get("pipeline_commit"))
        rdkit_version = _optional_text(row.get("rdkit_version"))
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
                **kwargs,
            )
            written.append({
                "name": name,
                "conformer_id": conformer_id,
                "xyz_path": xyz_path,
                "com_path": com_path,
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
