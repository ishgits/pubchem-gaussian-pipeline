"""Tests for pipeline.gaussian"""

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.gaussian import (
    xyz_to_gaussian_coords,
    write_gaussian_com,
    write_gaussian_coms,
    write_gaussian_coms_from_conformers,
)

SAMPLE_XYZ = os.path.join(os.path.dirname(__file__), "sample_data", "water.xyz")


def _write_xyz(tmpdir, contents):
    path = os.path.join(tmpdir, "mol.xyz")
    with open(path, "w") as f:
        f.write(contents)
    return path


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


class TestXyzParsingByPhysicalLine:
    """B-01: parse by physical line so an empty comment never drops an atom, and
    a declared/actual count mismatch raises instead of silently truncating."""

    def test_empty_comment_keeps_all_atoms(self):
        # Line 2 (the comment) is legitimately empty. All 3 atoms must survive —
        # the old blank-line filter dropped the count line + first atom here.
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz = _write_xyz(
                tmpdir,
                "3\n\n"
                "O    0.0  0.0  0.117\n"
                "H    0.0  0.757 -0.469\n"
                "H    0.0 -0.757 -0.469\n",
            )
            coords = xyz_to_gaussian_coords(xyz)
            lines = coords.strip().split("\n")
            assert len(lines) == 3
            assert [ln.split()[0] for ln in lines] == ["O", "H", "H"]

    def test_count_greater_than_rows_raises(self):
        # Declares 4 atoms but only 3 coordinate rows are present.
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz = _write_xyz(
                tmpdir,
                "4\nwater\n"
                "O    0.0  0.0  0.117\n"
                "H    0.0  0.757 -0.469\n"
                "H    0.0 -0.757 -0.469\n",
            )
            with pytest.raises(ValueError):
                xyz_to_gaussian_coords(xyz)

    def test_count_less_than_rows_raises(self):
        # Declares 2 atoms but 3 rows follow — a naive slice would silently drop
        # the extra atom; we require an exact match and raise.
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz = _write_xyz(
                tmpdir,
                "2\nwater\n"
                "O    0.0  0.0  0.117\n"
                "H    0.0  0.757 -0.469\n"
                "H    0.0 -0.757 -0.469\n",
            )
            with pytest.raises(ValueError):
                xyz_to_gaussian_coords(xyz)

    def test_trailing_blank_line_tolerated(self):
        # A trailing newline / blank line is normal and must not trip the count.
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz = _write_xyz(
                tmpdir,
                "2\nH2\n"
                "H  0.0 0.0 0.0\n"
                "H  0.0 0.0 0.74\n\n",
            )
            coords = xyz_to_gaussian_coords(xyz)
            assert len(coords.strip().split("\n")) == 2

    def test_malformed_coordinate_row_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz = _write_xyz(
                tmpdir,
                "2\nbad\n"
                "O    0.0  0.0  0.117\n"
                "H    0.0  not_a_number -0.469\n",
            )
            with pytest.raises(ValueError):
                xyz_to_gaussian_coords(xyz)

    def test_noninteger_count_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz = _write_xyz(tmpdir, "notanumber\ncomment\nO 0 0 0\n")
            with pytest.raises(ValueError):
                xyz_to_gaussian_coords(xyz)


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


class TestWriteGaussianComConformer:
    """v2: conformer-aware naming, title ΔE, and Link1 preservation."""

    def test_conformer_filename_and_chk(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Ribose",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
                conformer_id=0,
                rel_energy_kcalmol=0.0,
            )
            assert com_path.endswith("ribose_c00_F.com")
            with open(com_path) as f:
                text = f.read()
            assert "%chk=ribose_c00_F.chk" in text
            assert text.count("%chk=ribose_c00_F.chk") == 2  # opt + Link1
            assert "--Link1--" in text  # Link1 contract preserved

    def test_conformer_delta_energy_in_title(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Ribose",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
                conformer_id=2,
                rel_energy_kcalmol=1.2345,
            )
            with open(com_path) as f:
                text = f.read()
            assert "dE=1.2345 kcal/mol" in text
            assert "ribose_c02" in text  # id in title line

    def test_provenance_is_stamped_only_in_title_section(self):
        route_opt = "# opt b3lyp/6-31g(d)"
        route_freq = "# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read"
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Ribose",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt=route_opt,
                route_freq=route_freq,
                charge=-1,
                multiplicity=2,
                conformer_id=0,
                rel_energy_kcalmol=0.0,
                pipeline_version="2.0.0",
                pipeline_commit="abc1234",
                rdkit_version="2025.09.3",
            )
            with open(com_path) as f:
                text = f.read()

            provenance = (
                "provenance pipeline=2.0.0 commit=abc1234 rdkit=2025.09.3"
            )
            title_and_rest = text.split(f"{route_opt}\n\n", 1)[1]
            title_block, after_title = title_and_rest.split("\n\n", 1)
            assert provenance in title_block
            assert text.count(provenance) == 1
            assert provenance not in after_title
            assert text.count(route_opt) == 1
            assert text.count(route_freq) == 1
            assert "\n-1 2\n" in text
            assert text.count("%chk=ribose_c00_F.chk") == 2
            assert "--Link1--" in text

    def test_missing_commit_is_explicit_when_other_provenance_exists(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Ribose",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
                conformer_id=0,
                pipeline_version="2.0.0",
                pipeline_commit="",
                rdkit_version="2025.09.3",
            )
            with open(com_path) as f:
                text = f.read()
            assert (
                "provenance pipeline=2.0.0 commit=unavailable rdkit=2025.09.3"
                in text
            )

    def test_none_conformer_preserves_v1_naming(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            com_path = write_gaussian_com(
                name="Water",
                xyz_path=SAMPLE_XYZ,
                outdir=tmpdir,
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d)",
            )
            assert com_path.endswith("water_F.com")
            with open(com_path) as f:
                text = f.read()
            assert "provenance " not in text


class TestWriteGaussianComsFromConformers:
    def test_three_conformer_rows_three_coms(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmpdir:
            log_csv = os.path.join(tmpdir, "conformer_log.csv")
            pd.DataFrame([
                {"name": "Ribose", "conformer_id": 0, "rel_energy_kcalmol": 0.0, "xyz_path": SAMPLE_XYZ},
                {"name": "Ribose", "conformer_id": 1, "rel_energy_kcalmol": 0.5, "xyz_path": SAMPLE_XYZ},
                {"name": "Ribose", "conformer_id": 2, "rel_energy_kcalmol": 1.1, "xyz_path": SAMPLE_XYZ},
            ]).to_csv(log_csv, index=False)

            outdir = os.path.join(tmpdir, "gaussian_inputs")
            out = write_gaussian_coms_from_conformers(
                log_csv,
                outdir=outdir,
                log_csv=os.path.join(tmpdir, "com_write_log.csv"),
                route_opt="# opt b3lyp/6-31g(d)",
                route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
            )
            assert len(out) == 3
            names = sorted(os.path.basename(p) for p in out["com_path"])
            assert names == ["ribose_c00_F.com", "ribose_c01_F.com", "ribose_c02_F.com"]

            # Each file has an intact Link1 section and its ΔE in the title.
            deltas = {0: "0.0000", 1: "0.5000", 2: "1.1000"}
            for _, row in out.iterrows():
                with open(row["com_path"]) as f:
                    text = f.read()
                assert "--Link1--" in text
                assert f"dE={deltas[row['conformer_id']]} kcal/mol" in text

    def test_batch_copies_provenance_into_com_and_log(self, tmp_path):
        import pandas as pd

        conformer_log = tmp_path / "conformer_log.csv"
        pd.DataFrame([{
            "name": "Ribose",
            "conformer_id": 0,
            "rel_energy_kcalmol": 0.0,
            "xyz_path": SAMPLE_XYZ,
            "pipeline_version": "2.0.0",
            "pipeline_commit": "abc1234",
            "rdkit_version": "2025.09.3",
        }]).to_csv(conformer_log, index=False)
        com_log = tmp_path / "com_write_log.csv"

        out = write_gaussian_coms_from_conformers(
            str(conformer_log),
            outdir=str(tmp_path / "gaussian_inputs"),
            log_csv=str(com_log),
            route_opt="# opt b3lyp/6-31g(d)",
            route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
        )

        expected_columns = [
            "name", "conformer_id", "xyz_path", "com_path",
            "pipeline_version", "pipeline_commit", "rdkit_version",
        ]
        assert list(out.columns) == expected_columns
        assert out.loc[0, "pipeline_version"] == "2.0.0"
        assert out.loc[0, "pipeline_commit"] == "abc1234"
        assert out.loc[0, "rdkit_version"] == "2025.09.3"
        written_log = pd.read_csv(com_log, dtype=str, keep_default_na=False)
        assert list(written_log.columns) == expected_columns
        assert written_log.loc[0, "pipeline_version"] == "2.0.0"
        assert written_log.loc[0, "pipeline_commit"] == "abc1234"
        assert written_log.loc[0, "rdkit_version"] == "2025.09.3"
        with open(out.loc[0, "com_path"]) as f:
            text = f.read()
        assert (
            "provenance pipeline=2.0.0 commit=abc1234 rdkit=2025.09.3" in text
        )


class TestEmptyComLogSchemas:
    """M-11: zero writes remain valid, readable CSVs for downstream stages."""

    def test_empty_legacy_log_keeps_headers(self, tmp_path, monkeypatch):
        import pandas as pd

        monkeypatch.chdir(tmp_path)
        xyz_log = tmp_path / "xyz_log.csv"
        pd.DataFrame(columns=["name", "xyz_path"]).to_csv(xyz_log, index=False)
        com_log = tmp_path / "com_write_log.csv"

        out = write_gaussian_coms(str(xyz_log), log_csv=str(com_log))

        assert list(out.columns) == ["name", "xyz_path", "com_path"]
        assert out.empty
        assert list(pd.read_csv(com_log).columns) == ["name", "xyz_path", "com_path"]

    def test_empty_conformer_log_keeps_headers(self, tmp_path, monkeypatch):
        import pandas as pd

        monkeypatch.chdir(tmp_path)
        conformer_log = tmp_path / "conformer_log.csv"
        pd.DataFrame(columns=["name", "conformer_id", "xyz_path"]).to_csv(
            conformer_log, index=False
        )
        com_log = tmp_path / "com_write_log.csv"

        out = write_gaussian_coms_from_conformers(
            str(conformer_log), log_csv=str(com_log)
        )

        expected = [
            "name", "conformer_id", "xyz_path", "com_path",
            "pipeline_version", "pipeline_commit", "rdkit_version",
        ]
        assert list(out.columns) == expected
        assert out.empty
        assert list(pd.read_csv(com_log).columns) == expected

    def test_every_conformer_write_failure_keeps_headers(
        self, tmp_path, monkeypatch
    ):
        import pandas as pd
        import pipeline.gaussian as G

        monkeypatch.chdir(tmp_path)
        conformer_log = tmp_path / "conformer_log.csv"
        pd.DataFrame([{
            "name": "Water",
            "conformer_id": 0,
            "xyz_path": "missing.xyz",
        }]).to_csv(conformer_log, index=False)
        monkeypatch.setattr(
            G,
            "write_gaussian_com",
            lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("bad XYZ")),
        )
        com_log = tmp_path / "com_write_log.csv"

        out = G.write_gaussian_coms_from_conformers(
            str(conformer_log), log_csv=str(com_log)
        )

        expected = [
            "name", "conformer_id", "xyz_path", "com_path",
            "pipeline_version", "pipeline_commit", "rdkit_version",
        ]
        assert out.empty
        assert list(pd.read_csv(com_log).columns) == expected
        failures = pd.read_csv(tmp_path / "com_write_failed.csv")
        assert len(failures) == 1

    def test_all_ineligible_molecules_complete_zero_job_pipeline(
        self, tmp_path, monkeypatch
    ):
        import pandas as pd
        import pipeline.conformers as C
        from pipeline.slurm import write_slurm_scripts

        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        monkeypatch.setattr(
            C, "check_conformer_eligibility", lambda smiles: "undefined stereochemistry"
        )
        molecule_table = pd.DataFrame([{
            "name": "Ambiguous",
            "cid": 1,
            "IsomericSMILES": "CC(F)Cl",
        }])
        conformer_log = tmp_path / "conformer_log.csv"
        failed_log = tmp_path / "conformer_search_failed.csv"
        C.search_conformers(
            molecule_table,
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(conformer_log),
            failed_csv=str(failed_log),
        )

        com_log = tmp_path / "com_write_log.csv"
        coms = write_gaussian_coms_from_conformers(
            str(conformer_log), log_csv=str(com_log)
        )
        slurm_dir = tmp_path / "slurm_scripts"
        slurm_dir.mkdir()
        (slurm_dir / "stale_F.sh").write_text("#!/bin/bash\n")
        slurm_log = tmp_path / "slurm_write_log.csv"
        scripts = write_slurm_scripts(
            com_log_csv=str(com_log),
            slurm_dir=str(slurm_dir),
            log_csv=str(slurm_log),
        )

        assert pd.read_csv(conformer_log).empty
        assert len(pd.read_csv(failed_log)) == 1
        assert coms.empty
        assert scripts.empty
        assert list(slurm_dir.glob("*.sh")) == []
        assert list(pd.read_csv(slurm_log).columns) == [
            "jobname", "com_path", "sh_path", "status"
        ]
