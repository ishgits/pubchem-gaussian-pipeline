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
g16 {jobname}.com
"""


def write_slurm_script(
    jobname: str,
    outdir: str,
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
    template : str
        SLURM template with ``{jobname}``, ``{account}``, ``{cpus}``,
        ``{mem}``, and ``{time}`` placeholders.
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

    text = template.format(
        jobname=jobname,
        account=account,
        cpus=cpus,
        mem=mem,
        time=time,
    )

    with open(sh_path, "w") as f:
        f.write(text)

    return sh_path


def write_slurm_scripts(
    com_dir: str = "gaussian_inputs",
    slurm_dir: str = "slurm_scripts",
    log_csv: str = "slurm_write_log.csv",
    **kwargs,
) -> pd.DataFrame:
    """
    Generate one SLURM .sh script per .com file found in *com_dir*.

    Resume-safe: skips .sh files that already exist.
    """
    ensure_dir(slurm_dir)
    com_files = sorted(glob.glob(os.path.join(com_dir, "*.com")))

    rows = []
    for com_path in com_files:
        jobname = os.path.splitext(os.path.basename(com_path))[0]
        sh_path = os.path.join(slurm_dir, f"{jobname}.sh")

        # Resume-safe
        if os.path.exists(sh_path) and os.path.getsize(sh_path) > 0:
            rows.append({"jobname": jobname, "com_path": com_path, "sh_path": sh_path, "status": "SKIPPED_EXISTS"})
            continue

        write_slurm_script(jobname, slurm_dir, **kwargs)
        rows.append({"jobname": jobname, "com_path": com_path, "sh_path": sh_path, "status": "OK"})

    df = pd.DataFrame(rows)
    df.to_csv(log_csv, index=False)
    print(f"Wrote: {log_csv}")
    return df
