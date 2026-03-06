"""Tests for pipeline.gaussian"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.gaussian import xyz_to_gaussian_coords, write_gaussian_com

SAMPLE_XYZ = os.path.join(os.path.dirname(__file__), "sample_data", "water.xyz")


class TestXyzToGaussianCoords:
    def test_three_atoms(self):
        coords = xyz_to_gaussian_coords(SAMPLE_XYZ)
        lines = coords.strip().split("\n")
        assert len(lines) == 3, f"Expected 3 coordinate lines, got {len(lines)}"

    def test_element_symbols(self):
        coords = xyz_to_gaussian_coords(SAMPLE_XYZ)
        symbols = [line.split()[0] for line in coords.strip().split("\n")]
        assert symbols == ["O", "H", "H"]

    def test_coordinate_values_parseable(self):
        coords = xyz_to_gaussian_coords(SAMPLE_XYZ)
        for line in coords.strip().split("\n"):
            parts = line.split()
            assert len(parts) == 4
            # Coordinates should be parseable as floats
            float(parts[1])
            float(parts[2])
            float(parts[3])


class TestWriteGaussianCom:
    def test_creates_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Water",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
            )
            assert os.path.exists(com_path)
            assert com_path.endswith("water_F.com")

    def test_link1_present(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Water",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
                link1=True,
            )
            with open(com_path) as f:
                text = f.read()
            assert "--Link1--" in text

    def test_link1_absent_when_disabled(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Water",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d)",
                link1=False,
            )
            with open(com_path) as f:
                text = f.read()
            assert "--Link1--" not in text

    def test_charge_multiplicity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Water",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d)",
                charge=-1,
                multiplicity=2,
            )
            with open(com_path) as f:
                text = f.read()
            assert "-1 2\n" in text

    def test_chk_filename_matches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Adenine",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
            )
            with open(com_path) as f:
                text = f.read()
            assert "%chk=adenine_F.chk" in text
            # chk should appear twice (opt section + Link1 section)
            assert text.count("%chk=adenine_F.chk") == 2

    def test_nproc_in_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Water",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d)",
                nproc=24,
            )
            with open(com_path) as f:
                text = f.read()
            assert "%nprocshared=24" in text
