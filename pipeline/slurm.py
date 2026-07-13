"""
SLURM submission script (.sh) generator for Gaussian jobs.

The default template is a starting point — users should customize the
account, module loads, and resource requests for their own cluster.
"""

from __future__ import annotations

import glob
import os

import pandas as pd

from .manifest import (
    assert_stage_configuration,
    find_artifact,
    record_child_artifact,
    relative_artifact_path,
    remove_artifacts_by_kind,
    sha256_file,
    slurm_template_identity,
    stable_record_id,
)
from .utils import ensure_dir


_SLURM_LOG_COLUMNS = [
    "run_id", "artifact_id", "config_hash", "jobname", "com_artifact_id",
    "com_path", "com_sha256", "sh_path", "sh_sha256", "status",
]


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


def _validated_logged_com_paths(
    com_log: pd.DataFrame, manifest_path: str, manifest: dict
) -> list[dict]:
    """Return valid logged COM paths or raise before SLURM output mutation (M-18)."""
    required = {
        "run_id", "artifact_id", "config_hash", "conformer_record_id",
        "com_path", "com_sha256",
    }
    missing = sorted(required - set(com_log.columns))
    if missing and not com_log.empty:
        raise ValueError(
            "COM write log is missing required column(s): " + ", ".join(missing)
        )

    prepared = []
    problems = []
    seen_sources = set()
    seen_artifacts = set()
    for index, value in com_log["com_path"].items():
        if value is None or pd.isna(value) or str(value).strip() == "":
            problems.append(f"blank com_path at row {int(index)}")
            continue
        com_path = str(value)
        if not os.path.isfile(com_path):
            problems.append(f"missing com_path at row {int(index)}: {com_path!r}")
            continue
        if os.path.getsize(com_path) == 0:
            problems.append(f"zero-byte com_path at row {int(index)}: {com_path!r}")
            continue
        normalized = os.path.normcase(os.path.abspath(os.path.realpath(com_path)))
        if normalized in seen_sources:
            problems.append(f"duplicate normalized com_path at row {int(index)}: {com_path!r}")
            continue
        seen_sources.add(normalized)
        try:
            if str(com_log.at[index, "run_id"]) != manifest["run_id"] or str(com_log.at[index, "config_hash"]) != manifest["config_hash"]:
                raise ValueError("run/config identity mismatch")
            artifact_id = str(com_log.at[index, "artifact_id"])
            if artifact_id in seen_artifacts:
                raise ValueError("duplicate COM artifact record")
            seen_artifacts.add(artifact_id)
            artifact = find_artifact(manifest, artifact_id)
            if artifact["kind"] != "com":
                raise ValueError("referenced artifact is not COM")
            if artifact.get("conformer_record_id") != str(com_log.at[index, "conformer_record_id"]):
                raise ValueError("conformer lineage mismatch")
            if relative_artifact_path(com_path, manifest_path) != artifact["relative_path"]:
                raise ValueError("COM path mismatch")
            actual_hash = sha256_file(com_path)
            if str(com_log.at[index, "com_sha256"]) != artifact["sha256"] or actual_hash != artifact["sha256"]:
                raise ValueError("COM hash mismatch")
        except (KeyError, OSError, TypeError, ValueError) as exc:
            problems.append(f"manifest mismatch at row {int(index)}: {exc}")
            continue
        prepared.append({"com_path": com_path, "artifact": artifact})

    if problems:
        raise ValueError(
            "Cannot write SLURM scripts from invalid COM log entries: "
            + "; ".join(problems)
        )
    return prepared


def write_slurm_script(
    jobname: str,
    outdir: str,
    com_path: str | None = None,
    template: str = DEFAULT_TEMPLATE,
    account: str = "myaccount",
    cpus: int = 16,
    mem: str = "32G",
    time: str = "24:00:00",
    run_id: str | None = None,
    artifact_id: str | None = None,
    source_com_relative_path: str | None = None,
    source_com_sha256: str | None = None,
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
    linkage = (run_id, artifact_id, source_com_relative_path, source_com_sha256)
    if any(value is not None for value in linkage):
        if any(value is None or not str(value).strip() for value in linkage):
            raise ValueError(
                "Linked SLURM scripts require run_id, artifact_id, source COM "
                "relative path, and source COM SHA-256."
            )
        header = (
            f"# run_id={run_id}\n"
            f"# artifact_id={artifact_id}\n"
            f"# source_com_relative_path={source_com_relative_path}\n"
            f"# source_com_sha256={source_com_sha256}\n"
        )
        if text.startswith("#!") and "\n" in text:
            shebang, remainder = text.split("\n", 1)
            text = f"{shebang}\n{header}{remainder}"
        else:
            text = header + text

    with open(sh_path, "w") as f:
        f.write(text)

    return sh_path


def write_slurm_scripts(
    com_log_csv: str = "com_write_log.csv",
    slurm_dir: str = "slurm_scripts",
    log_csv: str = "slurm_write_log.csv",
    com_dir: str | None = None,
    manifest_path: str = "run_manifest.json",
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
        prepared = [{"com_path": path, "artifact": None} for path in com_paths]
        manifest = None
    else:
        # Default: consume the current run's com_write_log.csv.
        com_log = pd.read_csv(com_log_csv)
        slurm_config = {
            "account": kwargs.get("account", "myaccount"),
            "cpus": kwargs.get("cpus", 16),
            "mem": kwargs.get("mem", "32G"),
            "time": kwargs.get("time", "24:00:00"),
            "template_sha256": slurm_template_identity(
                kwargs.get("template", DEFAULT_TEMPLATE)
            ),
        }
        manifest = assert_stage_configuration(
            manifest_path, "slurm", slurm_config
        )
        prepared = _validated_logged_com_paths(com_log, manifest_path, manifest)
        com_paths = [item["com_path"] for item in prepared]

    # Prove one-to-one source→destination mapping and validate the template for
    # every job before pruning/creating any file or rewriting either log.
    destinations = {}
    basenames = {}
    for item in prepared:
        com_path = item["com_path"]
        jobname = os.path.splitext(os.path.basename(com_path))[0]
        destination = os.path.normcase(
            os.path.abspath(os.path.join(slurm_dir, f"{jobname}.sh"))
        )
        if jobname in basenames:
            raise ValueError(
                f"Two COM inputs collapse to script basename {jobname!r}: "
                f"{basenames[jobname]!r} and {com_path!r}."
            )
        if destination in destinations:
            raise ValueError(
                f"Duplicate SLURM destination {destination!r} for "
                f"{destinations[destination]!r} and {com_path!r}."
            )
        basenames[jobname] = com_path
        destinations[destination] = com_path
        template = kwargs.get("template", DEFAULT_TEMPLATE)
        template.format(
            jobname=jobname,
            account=kwargs.get("account", "myaccount"),
            cpus=kwargs.get("cpus", 16),
            mem=kwargs.get("mem", "32G"),
            time=kwargs.get("time", "24:00:00"),
            com_relpath=os.path.relpath(com_path, slurm_dir),
        )

    # M-18 validation above intentionally completes before this first mutation.
    # A stale/damaged COM log must not prune good scripts, create a directory, or
    # overwrite the prior SLURM write log with jobs that cannot run.
    ensure_dir(slurm_dir)

    if manifest is not None:
        remove_artifacts_by_kind(manifest_path, "sh")

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
    for item in prepared:
        com_path = item["com_path"]
        com_artifact = item["artifact"]
        jobname = os.path.splitext(os.path.basename(com_path))[0]
        sh_path = os.path.join(slurm_dir, f"{jobname}.sh")

        existed = os.path.exists(sh_path) and os.path.getsize(sh_path) > 0
        linked_kwargs = {}
        sh_artifact_id = ""
        if manifest is not None:
            sh_artifact_id = stable_record_id(
                manifest["run_id"], "sh", com_artifact["artifact_id"]
            )
            linked_kwargs = {
                "run_id": manifest["run_id"],
                "artifact_id": sh_artifact_id,
                "source_com_relative_path": com_artifact["relative_path"],
                "source_com_sha256": com_artifact["sha256"],
            }
        write_slurm_script(
            jobname, slurm_dir, com_path=com_path, **linked_kwargs, **kwargs
        )
        sh_digest = ""
        if manifest is not None:
            sh_digest = record_child_artifact(
                manifest_path,
                kind="sh",
                artifact_id=sh_artifact_id,
                parent_artifact_id=com_artifact["artifact_id"],
                conformer_record_id=com_artifact["conformer_record_id"],
                path=sh_path,
            )
        rows.append({
            "run_id": "" if manifest is None else manifest["run_id"],
            "artifact_id": sh_artifact_id,
            "config_hash": "" if manifest is None else manifest["config_hash"],
            "jobname": jobname,
            "com_artifact_id": "" if com_artifact is None else com_artifact["artifact_id"],
            "com_path": com_path,
            "com_sha256": "" if com_artifact is None else com_artifact["sha256"],
            "sh_path": sh_path,
            "sh_sha256": sh_digest,
            "status": "OVERWROTE" if existed else "WROTE",
        })

    df = pd.DataFrame(rows, columns=_SLURM_LOG_COLUMNS)
    df.to_csv(log_csv, index=False)
    print(f"Wrote: {log_csv}")
    return df
