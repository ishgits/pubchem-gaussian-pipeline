"""
SDF → XYZ conversion using Open Babel.

Open Babel must be installed and available on PATH (``obabel`` or ``babel``).
Install via conda:  conda install -c conda-forge openbabel
Or on HPC:          module load openbabel
"""

from __future__ import annotations

import os
import shutil
import subprocess

import pandas as pd

from .utils import ensure_dir, sanitize_basename


def _find_obabel() -> str:
    """Locate the Open Babel executable, or raise with install instructions."""
    exe = shutil.which("obabel") or shutil.which("babel")
    if exe is None:
        raise RuntimeError(
            "Open Babel not found on PATH (obabel / babel).\n"
            "Install:  conda install -c conda-forge openbabel\n"
            "Or HPC:   module load openbabel"
        )
    return exe


def sdf_to_xyz(sdf_path: str, xyz_path: str, gen3d: bool = True, minimize: bool = True) -> None:
    """
    Convert an SDF file to XYZ format using Open Babel.

    Parameters
    ----------
    sdf_path : str
        Path to the input .sdf file.
    xyz_path : str
        Path for the output .xyz file.
    gen3d : bool
        If True, force 3D coordinate generation (useful if the SDF is 2D).
    minimize : bool
        If True, run a quick force-field minimization after 3D generation.
    """
    obabel = _find_obabel()
    cmd = [obabel, sdf_path, "-O", xyz_path]
    if gen3d:
        cmd.append("--gen3d")
    if minimize:
        cmd.append("--minimize")
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def convert_sdfs_to_xyz(
    download_log_csv: str,
    xyz_dir: str = "pubchem_xyz",
    log_csv: str = "xyz_convert_log.csv",
    **kwargs,
) -> pd.DataFrame:
    """
    Batch-convert all SDFs listed in *download_log_csv* to XYZ files.

    Resume-safe: skips files that already exist and are non-empty.
    """
    log = pd.read_csv(download_log_csv)
    ensure_dir(xyz_dir)

    conv_ok = []
    conv_fail = []

    for _, row in log.iterrows():
        name = row["name"]
        sdf_path = row["sdf_path"]
        base = sanitize_basename(name)
        xyz_path = os.path.join(xyz_dir, f"{base}.xyz")

        # Resume-safe
        if os.path.exists(xyz_path) and os.path.getsize(xyz_path) > 0:
            conv_ok.append({"name": name, "sdf_path": sdf_path, "xyz_path": xyz_path, "status": "SKIPPED_EXISTS"})
            continue

        try:
            sdf_to_xyz(sdf_path, xyz_path, **kwargs)
            conv_ok.append({"name": name, "sdf_path": sdf_path, "xyz_path": xyz_path, "status": "OK"})
        except Exception as e:
            conv_fail.append({"name": name, "sdf_path": sdf_path, "error": repr(e)})

    ok_df = pd.DataFrame(conv_ok)
    ok_df.to_csv(log_csv, index=False)

    if conv_fail:
        fail_df = pd.DataFrame(conv_fail)
        fail_df.to_csv("xyz_convert_failed.csv", index=False)
        print(f"WARNING: {len(conv_fail)} conversions failed — see xyz_convert_failed.csv")
    else:
        print("All XYZ conversions succeeded.")

    print(f"Wrote: {log_csv}")
    return ok_df
