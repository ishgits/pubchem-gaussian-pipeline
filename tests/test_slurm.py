"""Tests for pipeline.slurm"""

import os
import sys
import tempfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.slurm import write_slurm_script, write_slurm_scripts
from manifest_helpers import ensure_manifest, write_linked_com_log

SAMPLE_XYZ = os.path.join(os.path.dirname(__file__), "sample_data", "water.xyz")


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
            from pathlib import Path

            root = Path(tmpdir)
            com_dir = os.path.join(tmpdir, "gaussian_inputs")
            slurm_dir = os.path.join(tmpdir, "slurm_scripts")
            log_csv, manifest_path = write_linked_com_log(
                root,
                [
                    {"name": name, "com_path": root / "gaussian_inputs" / f"{name}_F.com"}
                    for name in ("a", "b", "c")
                ],
                SAMPLE_XYZ,
            )
            # A stale .com on disk that is NOT in the log must be ignored.
            with open(os.path.join(com_dir, "stale_F.com"), "w") as handle:
                handle.write("stale\n")

            df = write_slurm_scripts(
                com_log_csv=log_csv,
                slurm_dir=slurm_dir,
                log_csv=os.path.join(tmpdir, "slurm_write_log.csv"),
                manifest_path=manifest_path,
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

    @pytest.mark.parametrize(
        "bad_path",
        [None, "", "   ", "missing.com"],
    )
    def test_invalid_logged_com_path_fails_before_directory_creation(
        self, tmp_path, bad_path
    ):
        com_log = tmp_path / "com_write_log.csv"
        pd.DataFrame([{
            "run_id": "run-test",
            "artifact_id": "com-test",
            "config_hash": "a" * 64,
            "conformer_record_id": "conformer-test",
            "com_sha256": "b" * 64,
            "name": "bad",
            "xyz_path": "x.xyz",
            "com_path": bad_path,
        }]).to_csv(com_log, index=False)
        slurm_dir = tmp_path / "slurm_scripts"
        slurm_log = tmp_path / "slurm_write_log.csv"
        empty_table = pd.DataFrame(columns=["name", "cid", "IsomericSMILES"])
        manifest_path = ensure_manifest(tmp_path, empty_table)

        with pytest.raises(ValueError, match="invalid COM log entries"):
            write_slurm_scripts(
                com_log_csv=str(com_log),
                slurm_dir=str(slurm_dir),
                log_csv=str(slurm_log),
                manifest_path=manifest_path,
            )

        assert not slurm_dir.exists()
        assert not slurm_log.exists()

    def test_one_missing_logged_com_preserves_all_prior_outputs(self, tmp_path):
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [{"name": "Water", "com_path": tmp_path / "gaussian_inputs" / "water_F.com"}],
            SAMPLE_XYZ,
        )
        rows = pd.read_csv(com_log)
        com_dir = tmp_path / "gaussian_inputs"
        missing_com = com_dir / "missing_F.com"
        missing = rows.iloc[0].copy()
        missing["artifact_id"] = "missing-artifact"
        missing["com_path"] = str(missing_com)
        rows = pd.concat([rows, missing.to_frame().T], ignore_index=True)
        rows.to_csv(com_log, index=False)

        slurm_dir = tmp_path / "slurm_scripts"
        slurm_dir.mkdir()
        prior_script = slurm_dir / "prior_F.sh"
        prior_script.write_bytes(b"prior script\n")
        slurm_log = tmp_path / "slurm_write_log.csv"
        slurm_log.write_bytes(b"prior log\n")

        with pytest.raises(ValueError, match="missing com_path at row 1"):
            write_slurm_scripts(
                com_log_csv=str(com_log),
                slurm_dir=str(slurm_dir),
                log_csv=str(slurm_log),
                manifest_path=manifest_path,
            )

        assert prior_script.read_bytes() == b"prior script\n"
        assert slurm_log.read_bytes() == b"prior log\n"
        assert {path.name for path in slurm_dir.glob("*.sh")} == {"prior_F.sh"}
        assert not (slurm_dir / "water_F.sh").exists()


class TestStrictOneToOneManifestMapping:
    """Frozen v2 mapping failures must happen before any output mutation."""

    @staticmethod
    def _prior_outputs(tmp_path):
        slurm_dir = tmp_path / "slurm_scripts"
        slurm_dir.mkdir()
        prior_script = slurm_dir / "prior.sh"
        prior_script.write_bytes(b"prior script\n")
        slurm_log = tmp_path / "slurm_write_log.csv"
        slurm_log.write_bytes(b"prior log\n")
        return slurm_dir, prior_script, slurm_log

    def test_distinct_same_basename_paths_fail_before_mutation(self, tmp_path):
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [
                {"name": "First", "com_path": tmp_path / "a" / "same.com"},
                {"name": "Second", "com_path": tmp_path / "b" / "same.com"},
            ],
            SAMPLE_XYZ,
        )
        slurm_dir, prior_script, slurm_log = self._prior_outputs(tmp_path)
        with pytest.raises(ValueError, match="collapse to script basename"):
            write_slurm_scripts(
                com_log_csv=com_log,
                slurm_dir=str(slurm_dir),
                log_csv=str(slurm_log),
                manifest_path=manifest_path,
            )
        assert prior_script.read_bytes() == b"prior script\n"
        assert slurm_log.read_bytes() == b"prior log\n"

    def test_duplicate_normalized_source_fails_before_mutation(self, tmp_path):
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [{"name": "Water", "com_path": tmp_path / "inputs" / "water.com"}],
            SAMPLE_XYZ,
        )
        rows = pd.read_csv(com_log)
        rows = pd.concat([rows, rows], ignore_index=True)
        rows.to_csv(com_log, index=False)
        slurm_dir, prior_script, slurm_log = self._prior_outputs(tmp_path)
        with pytest.raises(ValueError, match="duplicate normalized com_path"):
            write_slurm_scripts(
                com_log_csv=com_log,
                slurm_dir=str(slurm_dir),
                log_csv=str(slurm_log),
                manifest_path=manifest_path,
            )
        assert prior_script.read_bytes() == b"prior script\n"
        assert slurm_log.read_bytes() == b"prior log\n"

    def test_zero_byte_com_fails_before_mutation(self, tmp_path):
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [{
                "name": "Water",
                "com_path": tmp_path / "inputs" / "water.com",
                "content": "",
            }],
            SAMPLE_XYZ,
        )
        slurm_dir, prior_script, slurm_log = self._prior_outputs(tmp_path)
        with pytest.raises(ValueError, match="zero-byte com_path"):
            write_slurm_scripts(
                com_log_csv=com_log,
                slurm_dir=str(slurm_dir),
                log_csv=str(slurm_log),
                manifest_path=manifest_path,
            )
        assert prior_script.read_bytes() == b"prior script\n"
        assert slurm_log.read_bytes() == b"prior log\n"

    def test_valid_inputs_have_one_record_and_linked_header_each(self, tmp_path):
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [
                {"name": "Water", "com_path": tmp_path / "inputs" / "water.com"},
                {"name": "Ammonia", "com_path": tmp_path / "inputs" / "ammonia.com"},
            ],
            SAMPLE_XYZ,
        )
        out = write_slurm_scripts(
            com_log_csv=com_log,
            slurm_dir=str(tmp_path / "slurm_scripts"),
            log_csv=str(tmp_path / "slurm_write_log.csv"),
            manifest_path=manifest_path,
        )
        assert len(out) == 2
        assert out["artifact_id"].is_unique
        assert out["sh_path"].is_unique
        for _, row in out.iterrows():
            text = open(row["sh_path"], encoding="utf-8").read()
            assert f"# run_id={row['run_id']}" in text
            assert f"# artifact_id={row['artifact_id']}" in text
            assert f"# source_com_sha256={row['com_sha256']}" in text


class TestWriteSlurmScriptsOverwrite:
    """M-03: a rerun with new SBATCH directives overwrites the stale script."""

    def test_rewrite_updates_account_and_reports_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_dir = os.path.join(tmpdir, "gaussian_inputs")
            slurm_dir = os.path.join(tmpdir, "slurm_scripts")
            os.makedirs(com_dir)
            com_path = os.path.join(com_dir, "adenine_F.com")
            with open(com_path, "w") as handle:
                handle.write("%chk=adenine.chk\n")
            log_csv = os.path.join(tmpdir, "com_write_log.csv")
            pd.DataFrame([
                {"name": "adenine", "xyz_path": "x.xyz", "com_path": com_path}
            ]).to_csv(log_csv, index=False)

            slog = os.path.join(tmpdir, "slurm_write_log.csv")
            first = write_slurm_scripts(
                com_dir=com_dir, slurm_dir=slurm_dir, log_csv=slog, account="old"
            )
            assert list(first["status"]) == ["WROTE"]

            second = write_slurm_scripts(
                com_dir=com_dir, slurm_dir=slurm_dir, log_csv=slog, account="new"
            )
            assert list(second["status"]) == ["OVERWROTE"]
            with open(os.path.join(slurm_dir, "adenine_F.sh")) as f:
                text = f.read()
            assert "#SBATCH --account=new" in text
            assert "#SBATCH --account=old" not in text


class TestWriteSlurmScriptsCurrentRunCleanup:
    """B-06/M-11: disk scripts and log reflect only the current run."""

    def _write_com_log(self, path, com_paths):
        pd.DataFrame(
            [
                {"name": os.path.basename(p), "xyz_path": "x.xyz", "com_path": p}
                for p in com_paths
            ],
            columns=["name", "xyz_path", "com_path"],
        ).to_csv(path, index=False)

    def test_smaller_rerun_prunes_stale_scripts(self, tmp_path):
        slurm_dir = tmp_path / "slurm_scripts"
        slurm_log = tmp_path / "slurm_write_log.csv"
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [{"name": "Water", "com_path": tmp_path / "gaussian_inputs" / "water_F.com"}],
            SAMPLE_XYZ,
        )
        slurm_dir.mkdir()
        (slurm_dir / "stale_F.sh").write_text("#!/bin/bash\n")
        out = write_slurm_scripts(
            com_log_csv=str(com_log),
            slurm_dir=str(slurm_dir),
            log_csv=str(slurm_log),
            manifest_path=manifest_path,
        )

        assert list(out["jobname"]) == ["water_F"]
        assert {p.name for p in slurm_dir.glob("*.sh")} == {"water_F.sh"}
        assert len(list(slurm_dir.glob("*.sh"))) == len(out)

    @pytest.mark.parametrize("retained_rows", [1, 0])
    def test_truncated_com_log_rejected_before_script_or_manifest_mutation(
        self, tmp_path, retained_rows
    ):
        com_log, manifest_path = write_linked_com_log(
            tmp_path,
            [
                {"name": "Water", "com_path": tmp_path / "gaussian_inputs" / "water_F.com"},
                {"name": "Glycine", "com_path": tmp_path / "gaussian_inputs" / "glycine_F.com"},
            ],
            SAMPLE_XYZ,
        )
        slurm_dir = tmp_path / "slurm_scripts"
        slurm_log = tmp_path / "slurm_write_log.csv"
        first = write_slurm_scripts(
            com_log_csv=str(com_log),
            slurm_dir=str(slurm_dir),
            log_csv=str(slurm_log),
            manifest_path=manifest_path,
        )
        manifest_before = open(manifest_path, "rb").read()
        log_before = slurm_log.read_bytes()
        script_bytes = {
            path: path.read_bytes() for path in slurm_dir.glob("*.sh")
        }

        rows = pd.read_csv(com_log)
        rows.iloc[:retained_rows].to_csv(com_log, index=False)

        with pytest.raises(ValueError, match="does not exactly match manifest com artifacts"):
            write_slurm_scripts(
                com_log_csv=str(com_log),
                slurm_dir=str(slurm_dir),
                log_csv=str(slurm_log),
                manifest_path=manifest_path,
            )

        assert open(manifest_path, "rb").read() == manifest_before
        assert slurm_log.read_bytes() == log_before
        assert len(first) == 2
        assert {path: path.read_bytes() for path in slurm_dir.glob("*.sh")} == script_bytes

    def test_zero_job_rerun_prunes_all_and_keeps_log_schema(self, tmp_path):
        slurm_dir = tmp_path / "slurm_scripts"
        slurm_dir.mkdir()
        (slurm_dir / "old_F.sh").write_text("#!/bin/bash\n")
        com_log = tmp_path / "com_write_log.csv"
        self._write_com_log(com_log, [])
        slurm_log = tmp_path / "slurm_write_log.csv"
        empty_table = pd.DataFrame(columns=["name", "cid", "IsomericSMILES"])
        manifest_path = ensure_manifest(tmp_path, empty_table)

        out = write_slurm_scripts(
            com_log_csv=str(com_log),
            slurm_dir=str(slurm_dir),
            log_csv=str(slurm_log),
            manifest_path=manifest_path,
        )

        expected = [
            "run_id", "artifact_id", "config_hash", "jobname",
            "com_artifact_id", "com_path", "com_sha256", "sh_path",
            "sh_sha256", "status",
        ]
        assert out.empty
        assert list(out.columns) == expected
        assert list(pd.read_csv(slurm_log).columns) == expected
        assert list(slurm_dir.glob("*.sh")) == []
