"""Tests for pipeline.slurm"""

import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.slurm import write_slurm_script, write_slurm_scripts


class TestWriteSlurmScript:
    def test_creates_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("adenine_F", tmpdir)
            assert os.path.exists(sh_path)
            assert sh_path.endswith("adenine_F.sh")

    def test_jobname_in_script(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("adenine_F", tmpdir)
            with open(sh_path) as f:
                text = f.read()
            assert "#SBATCH --job-name=adenine_F" in text

    def test_account_placeholder(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("water_F", tmpdir, account="pearce21")
            with open(sh_path) as f:
                text = f.read()
            assert "#SBATCH --account=pearce21" in text

    def test_resource_parameters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script(
                "water_F", tmpdir,
                cpus=32, mem="64G", time="48:00:00",
            )
            with open(sh_path) as f:
                text = f.read()
            assert "--cpus-per-task=32" in text
            assert "--mem=64G" in text
            assert "--time=48:00:00" in text

    def test_g16_runs_the_com_basename(self):
        # B-03: the script cd's to the input's directory and runs g16 on the
        # basename (resolved relative to the script's own location), not a bare
        # `g16 {jobname}.com` that only works from one directory.
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("cytosine_F", tmpdir)
            with open(sh_path) as f:
                text = f.read()
            assert 'g16 "$(basename "$COM_PATH")"' in text
            assert 'COM_PATH="$SCRIPT_DIR/cytosine_F.com"' in text

    def test_custom_template(self):
        custom = "#!/bin/bash\n#SBATCH --job-name={jobname}\necho {jobname}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("test_mol", tmpdir, template=custom)
            with open(sh_path) as f:
                text = f.read()
            assert "echo test_mol" in text


class TestSlurmScriptResolvesInput:
    """B-03: a script in slurm_scripts/ must resolve its .com in a sibling
    gaussian_inputs/ directory, regardless of the submission directory."""

    def test_sibling_dirs_relpath(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_dir = os.path.join(tmpdir, "gaussian_inputs")
            slurm_dir = os.path.join(tmpdir, "slurm_scripts")
            os.makedirs(com_dir)
            os.makedirs(slurm_dir)
            com_path = os.path.join(com_dir, "adenine_F.com")
            with open(com_path, "w") as f:
                f.write("%chk=adenine_F.chk\n")

            sh_path = write_slurm_script("adenine_F", slurm_dir, com_path=com_path)
            with open(sh_path) as f:
                text = f.read()

            # The stored path is relative and preserves the real directory name.
            assert 'COM_PATH="$SCRIPT_DIR/../gaussian_inputs/adenine_F.com"' in text

            # And it actually resolves back to the real input file on disk.
            rel = os.path.join("..", "gaussian_inputs", "adenine_F.com")
            resolved = os.path.normpath(os.path.join(slurm_dir, rel))
            assert resolved == os.path.normpath(com_path)
            assert os.path.exists(resolved)

    def test_custom_com_dir_name_preserved(self):
        # relpath preserves a non-default directory name (not hardcoded).
        with tempfile.TemporaryDirectory() as tmpdir:
            com_dir = os.path.join(tmpdir, "my_inputs")
            slurm_dir = os.path.join(tmpdir, "jobs")
            os.makedirs(com_dir)
            os.makedirs(slurm_dir)
            com_path = os.path.join(com_dir, "water_F.com")
            open(com_path, "w").close()
            sh_path = write_slurm_script("water_F", slurm_dir, com_path=com_path)
            with open(sh_path) as f:
                text = f.read()
            assert "../my_inputs/water_F.com" in text


class TestWriteSlurmScriptsLogDriven:
    """M-01: default consumes the current run's com_write_log.csv, not a glob."""

    def _make_log(self, tmpdir, com_paths):
        rows = [
            {"name": os.path.basename(p), "xyz_path": "x.xyz", "com_path": p}
            for p in com_paths
        ]
        log_csv = os.path.join(tmpdir, "com_write_log.csv")
        pd.DataFrame(rows).to_csv(log_csv, index=False)
        return log_csv

    def test_only_logged_jobs_get_scripts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_dir = os.path.join(tmpdir, "gaussian_inputs")
            slurm_dir = os.path.join(tmpdir, "slurm_scripts")
            os.makedirs(com_dir)
            logged = []
            for name in ("a_F", "b_F", "c_F"):
                p = os.path.join(com_dir, f"{name}.com")
                open(p, "w").close()
                logged.append(p)
            # A stale .com on disk that is NOT in the log must be ignored.
            open(os.path.join(com_dir, "stale_F.com"), "w").close()

            log_csv = self._make_log(tmpdir, logged)
            df = write_slurm_scripts(
                com_log_csv=log_csv,
                slurm_dir=slurm_dir,
                log_csv=os.path.join(tmpdir, "slurm_write_log.csv"),
            )
            assert len(df) == 3
            jobs = sorted(df["jobname"])
            assert jobs == ["a_F", "b_F", "c_F"]
            assert "stale_F" not in jobs
            # Exactly the 3 scripts on disk.
            scripts = sorted(f for f in os.listdir(slurm_dir) if f.endswith(".sh"))
            assert scripts == ["a_F.sh", "b_F.sh", "c_F.sh"]

    def test_legacy_com_dir_glob_still_reachable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_dir = os.path.join(tmpdir, "gaussian_inputs")
            slurm_dir = os.path.join(tmpdir, "slurm_scripts")
            os.makedirs(com_dir)
            for name in ("a_F", "b_F"):
                open(os.path.join(com_dir, f"{name}.com"), "w").close()
            df = write_slurm_scripts(
                com_dir=com_dir,
                slurm_dir=slurm_dir,
                log_csv=os.path.join(tmpdir, "slurm_write_log.csv"),
            )
            assert sorted(df["jobname"]) == ["a_F", "b_F"]


class TestWriteSlurmScriptsOverwrite:
    """M-03: a rerun with new SBATCH directives overwrites the stale script."""

    def test_rewrite_updates_account_and_reports_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_dir = os.path.join(tmpdir, "gaussian_inputs")
            slurm_dir = os.path.join(tmpdir, "slurm_scripts")
            os.makedirs(com_dir)
            com_path = os.path.join(com_dir, "adenine_F.com")
            open(com_path, "w").close()
            log_csv = os.path.join(tmpdir, "com_write_log.csv")
            pd.DataFrame([
                {"name": "adenine", "xyz_path": "x.xyz", "com_path": com_path}
            ]).to_csv(log_csv, index=False)

            slog = os.path.join(tmpdir, "slurm_write_log.csv")
            first = write_slurm_scripts(
                com_log_csv=log_csv, slurm_dir=slurm_dir, log_csv=slog, account="old"
            )
            assert list(first["status"]) == ["WROTE"]

            second = write_slurm_scripts(
                com_log_csv=log_csv, slurm_dir=slurm_dir, log_csv=slog, account="new"
            )
            assert list(second["status"]) == ["OVERWROTE"]
            with open(os.path.join(slurm_dir, "adenine_F.sh")) as f:
                text = f.read()
            assert "#SBATCH --account=new" in text
            assert "#SBATCH --account=old" not in text
