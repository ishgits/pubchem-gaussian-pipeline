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

from .utils import ensure_dir, sanitize_basename


def xyz_to_gaussian_coords(xyz_path: str) -> str:
    """
    Read an XYZ file and return the coordinate block formatted for a
    Gaussian input file.

    XYZ format expected::

        <atom_count>
        <comment line>
        Element  x  y  z
        ...
    """
    with open(xyz_path, "r") as f:
        lines = [ln.strip() for ln in f.readlines() if ln.strip()]

    body = lines[2:]  # skip atom count + comment
    out_lines = []
    for ln in body:
        parts = ln.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        x, y, z = map(float, parts[1:4])
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

    Returns
    -------
    str
        Path to the written .com file.
    """
    ensure_dir(outdir)
    base = sanitize_basename(name)
    chk_name = f"{base}_F.chk"
    com_path = os.path.join(outdir, f"{base}_F.com")

    coords = xyz_to_gaussian_coords(xyz_path)
    title = f"{base} {title_suffix}".strip()

    text = (
        f"%nprocshared={nproc}\n"
        f"%chk={chk_name}\n"
        f"{route_opt}\n\n"
        f"{title}\n\n"
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

    out_df = pd.DataFrame(written)
    out_df.to_csv(log_csv, index=False)

    if failed:
        fail_df = pd.DataFrame(failed)
        fail_df.to_csv("com_write_failed.csv", index=False)
        print(f"WARNING: {len(failed)} .com writes failed — see com_write_failed.csv")
    else:
        print("All Gaussian .com files written successfully.")

    print(f"Wrote: {log_csv}")
    return out_df
