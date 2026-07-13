"""
SLURM submission script (.sh) generator for Gaussian jobs.

The default template is a starting point — users should customize the
account, module loads, and resource requests for their own cluster.
"""

from __future__ import annotations

import glob
import os

import pandas as pd

from .utils import ensure_dir


_SLURM_LOG_COLUMNS = ["jobname", "com_path", "sh_path", "status"]


# ---------------------------------------------------------------------------
# Default template — EDIT THIS for your cluster
# ---------------------------------------------------------------------------
DEFAULT_TEMPLATE = """\
#!/bin/bash
#SBATCH --account={account}
#SBATCH --job-name={jobname}
#SBATCH --output={jobname}.out
#SBATCH --error={jobname}.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}

module load gaussian16

# Resolve the Gaussian input relative to THIS script's own location, so the job
# runs correctly no matter which directory `sbatch` is invoked from (B-03). The
# .com path is stored relative to the script (e.g. ../gaussian_inputs/x.com); we
# cd into the input's directory and run g16 on the basename so Gaussian's output
# files land beside the input.
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
COM_PATH="$SCRIPT_DIR/{com_relpath}"
cd "$(dirname "$COM_PATH")"
g16 "$(basename "$COM_PATH")"
"""


def _validated_logged_com_paths(com_log: pd.DataFrame) -> list[str]:
    """Return valid logged COM paths or raise before SLURM output mutation (M-18)."""
    if "com_path" not in com_log.columns:
        raise ValueError("COM write log is missing required column: com_path")

    com_paths = []
    problems = []
    for index, value in com_log["com_path"].items():
        if value is None or pd.isna(value) or str(value).strip() == "":
            problems.append(f"blank com_path at row {int(index)}")
            continue
        com_path = str(value)
        if not os.path.isfile(com_path):
            problems.append(f"missing com_path at row {int(index)}: {com_path!r}")
            continue
        com_paths.append(com_path)

    if problems:
        raise ValueError(
            "Cannot write SLURM scripts from invalid COM log entries: "
            + "; ".join(problems)
        )
    return com_paths


def write_slurm_script(
    jobname: str,
    outdir: str,
    com_path: str | None = None,
    template: str = DEFAULT_TEMPLATE,
    account: str = "myaccount",
    cpus: int = 16,
    mem: str = "32G",
    time: str = "24:00:00",
) -> str:
    """
    Write a single SLURM submission script for a Gaussian job.

    Parameters
    ----------
    jobname : str
        Base name of the .com file (without extension).
    outdir : str
        Directory to write the .sh file.
    com_path : str, optional
        Path to the ``.com`` input this job runs. The script references it
        **relative to its own location** (``os.path.relpath(com_path, outdir)``),
        preserving custom directory names, so submitting from any working
        directory still finds the input (B-03). Defaults to
        ``{outdir}/{jobname}.com`` (sibling to the script) for backward
        compatibility.
    template : str
        SLURM template with ``{jobname}``, ``{account}``, ``{cpus}``, ``{mem}``,
        ``{time}``, and ``{com_relpath}`` placeholders.
    account : str
        SLURM account/allocation name.
    cpus, mem, time : resource parameters.

    Returns
    -------
    str
        Path to the written .sh file.
    """
    ensure_dir(outdir)
    sh_path = os.path.join(outdir, f"{jobname}.sh")

    if com_path is None:
        com_path = os.path.join(outdir, f"{jobname}.com")
    com_relpath = os.path.relpath(com_path, outdir)

    text = template.format(
        jobname=jobname,
        account=account,
        cpus=cpus,
        mem=mem,
        time=time,
        com_relpath=com_relpath,
    )

    with open(sh_path, "w") as f:
        f.write(text)

    return sh_path


def write_slurm_scripts(
    com_log_csv: str = "com_write_log.csv",
    slurm_dir: str = "slurm_scripts",
    log_csv: str = "slurm_write_log.csv",
    com_dir: str | None = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Generate one SLURM ``.sh`` script per Gaussian ``.com`` of the current run.

    By default (M-01) the jobs come from the current run's *com_log_csv* (the
    ``com_write_log.csv`` produced by the Gaussian stage; columns
    ``name,xyz_path,com_path``), so only this run's inputs are turned into
    scripts — stale ``.com`` files left on disk from a previous molecule list are
    never picked up. Passing *com_dir* switches to the legacy explicit mode that
    globs every ``*.com`` in that directory instead.

    Before writing, ``.sh`` files absent from the current COM set are removed
    (B-06), so the dedicated output directory agrees with the current-run log and
    a submission glob cannot include stale jobs. Current scripts are
    **overwritten** (M-03): regeneration is cheap and a rerun with new SBATCH
    directives (account, resources) must not leave a stale script behind. The
    log's ``status`` column reports ``WROTE`` for a new file and ``OVERWROTE``
    when a non-empty script was replaced. A zero-job run writes a header-only log
    and leaves no ``.sh`` files (M-11). In log-driven mode, every ``com_path``
    must be nonblank and identify an existing file; an invalid entry aborts before
    creating/pruning scripts or rewriting the SLURM log (M-18).
    """
    if com_dir is not None:
        # Legacy explicit mode: every .com on disk becomes a job.
        com_paths = sorted(glob.glob(os.path.join(com_dir, "*.com")))
    else:
        # Default: consume the current run's com_write_log.csv.
        com_log = pd.read_csv(com_log_csv)
        com_paths = _validated_logged_com_paths(com_log)

    # M-18 validation above intentionally completes before this first mutation.
    # A stale/damaged COM log must not prune good scripts, create a directory, or
    # overwrite the prior SLURM write log with jobs that cannot run.
    ensure_dir(slurm_dir)

    expected_scripts = {
        os.path.normcase(os.path.abspath(os.path.join(
            slurm_dir,
            f"{os.path.splitext(os.path.basename(com_path))[0]}.sh",
        )))
        for com_path in com_paths
    }
    for stale_path in glob.glob(os.path.join(slurm_dir, "*.sh")):
        if os.path.normcase(os.path.abspath(stale_path)) not in expected_scripts:
            os.remove(stale_path)

    rows = []
    for com_path in com_paths:
        jobname = os.path.splitext(os.path.basename(com_path))[0]
        sh_path = os.path.join(slurm_dir, f"{jobname}.sh")

        existed = os.path.exists(sh_path) and os.path.getsize(sh_path) > 0
        write_slurm_script(jobname, slurm_dir, com_path=com_path, **kwargs)
        rows.append({
            "jobname": jobname,
            "com_path": com_path,
            "sh_path": sh_path,
            "status": "OVERWROTE" if existed else "WROTE",
        })

    df = pd.DataFrame(rows, columns=_SLURM_LOG_COLUMNS)
    df.to_csv(log_csv, index=False)
    print(f"Wrote: {log_csv}")
    return df
