"""Tests for pipeline.slurm"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.slurm import write_slurm_script


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

    def test_g16_command(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("cytosine_F", tmpdir)
            with open(sh_path) as f:
                text = f.read()
            assert "g16 cytosine_F.com" in text

    def test_custom_template(self):
        custom = "#!/bin/bash\n#SBATCH --job-name={jobname}\necho {jobname}\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            sh_path = write_slurm_script("test_mol", tmpdir, template=custom)
            with open(sh_path) as f:
                text = f.read()
            assert "echo test_mol" in text
