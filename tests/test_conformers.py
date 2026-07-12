"""Tests for pipeline.conformers

The pure ranking helper (`select_top_n`) is tested with no RDKit. The embedding
and batch tests require RDKit but are still offline / no-cluster (allowed by the
no-network test rule); they `importorskip` so a bare environment still passes the
pure tests.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.conformers import check_conformer_eligibility, select_top_n


# ---------------------------------------------------------------------------
# Task 1 — pure ranking helper (no RDKit)
# ---------------------------------------------------------------------------

class TestSelectTopN:
    def test_orders_lowest_first(self):
        # energies: index 2 lowest, then 0, then 3, then 1
        energies = [1.0, 5.0, -2.0, 3.0]
        assert select_top_n(energies, 4) == [2, 0, 3, 1]

    def test_keeps_only_top_n(self):
        energies = [1.0, 5.0, -2.0, 3.0]
        assert select_top_n(energies, 3) == [2, 0, 3]

    def test_n_larger_than_available(self):
        energies = [0.5, 0.1]
        assert select_top_n(energies, 5) == [1, 0]

    def test_n_zero_returns_empty(self):
        assert select_top_n([1.0, 2.0], 0) == []

    def test_n_negative_returns_empty(self):
        assert select_top_n([1.0, 2.0], -1) == []

    def test_empty_input(self):
        assert select_top_n([], 3) == []

    def test_tie_break_is_deterministic(self):
        # Equal energies → keep original index order (stable, reproducible).
        energies = [2.0, 2.0, 2.0]
        assert select_top_n(energies, 2) == [0, 1]

    def test_single_conformer(self):
        assert select_top_n([42.0], 3) == [0]


# Flexible sugar with UNspecified stereocenters (no @) — must be skipped.
UNDEFINED_STEREO_SUGAR = "OC1OCC(O)C(O)C1O"


# ---------------------------------------------------------------------------
# B-01 — stereo/validity eligibility gate
# ---------------------------------------------------------------------------

class TestCheckEligibility:
    def test_empty_smiles(self):
        # Pure (no RDKit): empty/None short-circuits before parsing.
        assert check_conformer_eligibility("") == "no IsomericSMILES"

    def test_none_smiles(self):
        assert check_conformer_eligibility(None) == "no IsomericSMILES"

    def test_whitespace_smiles(self):
        assert check_conformer_eligibility("   ") == "no IsomericSMILES"

    def test_no_stereocenter_molecule_eligible(self):
        pytest.importorskip("rdkit")
        # Adenine / water have no stereo elements → eligible (not skipped).
        assert check_conformer_eligibility("c1[nH]cnc2c1ncn2") is None
        assert check_conformer_eligibility("O") is None

    def test_defined_stereo_eligible(self):
        pytest.importorskip("rdkit")
        assert check_conformer_eligibility("C[C@@H](N)C(=O)O") is None

    def test_undefined_stereo_skipped(self):
        pytest.importorskip("rdkit")
        assert check_conformer_eligibility(UNDEFINED_STEREO_SUGAR) == "undefined stereochemistry"

    def test_unparseable_smiles(self):
        pytest.importorskip("rdkit")
        assert check_conformer_eligibility("not-a-smiles(") == "unparseable SMILES"


# ---------------------------------------------------------------------------
# Task 2 — RDKit embed + rank core (offline; requires rdkit)
# ---------------------------------------------------------------------------

class TestGenerateConformers:
    def test_butane_seeded_deterministic(self):
        pytest.importorskip("rdkit")
        from pipeline.conformers import generate_conformers

        coords1, energies1, method1 = generate_conformers("CCCC", seed=42)
        # ≥1 conformer produced
        assert len(coords1) >= 1
        assert len(energies1) == len(coords1)
        assert method1 == "MMFF94"  # butane has full MMFF params
        # atom tuples look like (symbol, x, y, z), coords in Ångström
        sym, x, y, z = coords1[0][0]
        assert isinstance(sym, str)
        assert all(isinstance(v, float) for v in (x, y, z))

        # Deterministic lowest-energy index under the fixed seed.
        coords2, energies2, _ = generate_conformers("CCCC", seed=42)
        assert select_top_n(energies1, 1) == select_top_n(energies2, 1)
        assert energies1 == energies2

    def test_invalid_smiles_raises(self):
        pytest.importorskip("rdkit")
        from pipeline.conformers import generate_conformers

        with pytest.raises(ValueError):
            generate_conformers("this-is-not-smiles(", seed=42)


# ---------------------------------------------------------------------------
# Task 3 — batch driver (offline; requires rdkit + pandas)
# ---------------------------------------------------------------------------

# Rigid, planar aromatic → conformers collapse to ~1 under RMSD pruning.
ADENINE_SMILES = "c1[nH]cnc2c1ncn2"
# Flexible furanose sugar → multiple distinct conformers expected.
RIBOSE_SMILES = "C([C@@H]1[C@H]([C@H]([C@H](O1)O)O)O)O"


class TestSearchConformers:
    def _table(self):
        pd = pytest.importorskip("pandas")
        return pd.DataFrame([
            {"name": "Adenine", "cid": 190, "IsomericSMILES": ADENINE_SMILES},
            {"name": "Ribose", "cid": 5779, "IsomericSMILES": RIBOSE_SMILES},
        ])

    def test_two_molecule_table(self, tmp_path):
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        table = self._table()
        log = search_conformers(
            table,
            xyz_dir=str(tmp_path / "conf_xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "conformer_search_failed.csv"),
        )

        # Adenine (rigid) collapses to a single conformer row.
        adenine_rows = log[log["name"] == "Adenine"]
        assert len(adenine_rows) == 1

        # Ribose (flexible) yields at most TOP_N=3 rows, at least 1.
        ribose_rows = log[log["name"] == "Ribose"]
        assert 1 <= len(ribose_rows) <= 3

        # Provenance columns populated for every kept conformer.
        for col in ["rdkit_version", "seed", "method", "n_generated", "n_kept", "rmsd_prune"]:
            assert log[col].notna().all()

        # Each molecule's lowest conformer has ΔE == 0; all ΔE ≥ 0 (kcal/mol).
        for name in ("Adenine", "Ribose"):
            rows = log[log["name"] == name]
            assert rows["rel_energy_kcalmol"].min() == pytest.approx(0.0, abs=1e-9)
            assert (rows["rel_energy_kcalmol"] >= -1e-9).all()

        # XYZ files actually written.
        for xyz_path in log["xyz_path"]:
            assert os.path.exists(xyz_path)

    def test_missing_smiles_logged_no_isomericsmiles(self, tmp_path):
        pd = pytest.importorskip("pandas")
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        table = pd.DataFrame([
            {"name": "NoSmiles", "cid": None, "IsomericSMILES": None},
        ])
        failed_csv = tmp_path / "conformer_search_failed.csv"
        xyz_dir = tmp_path / "conf_xyz"
        log = search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(failed_csv),
        )
        assert len(log) == 0
        assert failed_csv.exists()
        fail_df = pd.read_csv(failed_csv)
        assert "NoSmiles" in set(fail_df["name"])
        assert set(fail_df["error"]) == {"no IsomericSMILES"}
        # No XYZ files written for a skipped molecule.
        assert not any(os.scandir(xyz_dir)) if os.path.isdir(xyz_dir) else True

    def test_undefined_stereo_is_skipped_and_logged(self, tmp_path):
        pd = pytest.importorskip("pandas")
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        table = pd.DataFrame([
            {"name": "Adenine", "cid": 190, "IsomericSMILES": ADENINE_SMILES},
            {"name": "UndefSugar", "cid": 1, "IsomericSMILES": UNDEFINED_STEREO_SUGAR},
        ])
        failed_csv = tmp_path / "conformer_search_failed.csv"
        xyz_dir = tmp_path / "conf_xyz"
        log = search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(failed_csv),
        )
        # Adenine (no stereo) proceeds; the undefined-stereo sugar is skipped.
        assert set(log["name"]) == {"Adenine"}
        fail_df = pd.read_csv(failed_csv)
        assert dict(zip(fail_df["name"], fail_df["error"])) == {
            "UndefSugar": "undefined stereochemistry"
        }
        # No XYZ written for the skipped molecule.
        written = [os.path.basename(p) for p in log["xyz_path"]]
        assert all(name.startswith("adenine") for name in written)
        assert not any("undefsugar" in f for f in os.listdir(xyz_dir))

    def test_resume_skips_completed(self, tmp_path):
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        log_csv = tmp_path / "conformer_log.csv"
        xyz_dir = tmp_path / "conf_xyz"
        table = self._table()

        first = search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(log_csv),
            failed_csv=str(tmp_path / "fail.csv"),
        )
        # Rerun with the same table: nothing new appended (both already done).
        second = search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(log_csv),
            failed_csv=str(tmp_path / "fail.csv"),
        )
        assert len(second) == len(first)


# ---------------------------------------------------------------------------
# Task 5 — reproducibility (validation ≠ "it ran")
# ---------------------------------------------------------------------------

class TestReproducibility:
    def test_same_seed_same_selected_conformers(self, tmp_path):
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        table_kwargs = dict(
            xyz_dir=str(tmp_path / "xyz_a"),
            failed_csv=str(tmp_path / "fail_a.csv"),
        )
        import pandas as pd
        table = pd.DataFrame([
            {"name": "Ribose", "cid": 5779, "IsomericSMILES": RIBOSE_SMILES},
        ])

        log_a = search_conformers(
            table, log_csv=str(tmp_path / "log_a.csv"), seed=42, **table_kwargs
        )
        log_b = search_conformers(
            table,
            xyz_dir=str(tmp_path / "xyz_b"),
            log_csv=str(tmp_path / "log_b.csv"),
            failed_csv=str(tmp_path / "fail_b.csv"),
            seed=42,
        )

        # Same seed → identical selected conformers (count + ΔE ranking).
        assert list(log_a["conformer_id"]) == list(log_b["conformer_id"])
        assert list(log_a["rel_energy_kcalmol"]) == list(log_b["rel_energy_kcalmol"])


# ---------------------------------------------------------------------------
# M-02 — notebook code path, exercised offline (no PubChem in CI)
# ---------------------------------------------------------------------------

class TestNotebookPathOffline:
    """Runs the exact default-notebook path search_conformers →
    write_gaussian_coms_from_conformers on a hardcoded defined-stereo SMILES,
    proving per-conformer .com files without executing PubChem."""

    def test_conformer_path_produces_per_conformer_coms(self, tmp_path):
        pd = pytest.importorskip("pandas")
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers
        from pipeline.gaussian import write_gaussian_coms_from_conformers

        # Stand-in for the notebook's `df` (build_molecule_table output).
        df = pd.DataFrame([
            {"name": "Ribose", "cid": 5779, "IsomericSMILES": RIBOSE_SMILES},
        ])
        conf_log_csv = tmp_path / "conformer_log.csv"
        conf_log = search_conformers(
            df,
            xyz_dir=str(tmp_path / "conformer_xyz"),
            log_csv=str(conf_log_csv),
            failed_csv=str(tmp_path / "conformer_search_failed.csv"),
            n_generate=20,
            top_n=3,
            rmsd_prune=0.5,
            seed=42,
        )
        n = len(conf_log)
        assert 1 <= n <= 3  # ribose (flexible) → up to TOP_N distinct conformers

        com_log = write_gaussian_coms_from_conformers(
            conformer_log_csv=str(conf_log_csv),
            outdir=str(tmp_path / "gaussian_inputs"),
            log_csv=str(tmp_path / "com_write_log.csv"),
            route_opt="# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)",
            route_freq="# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read",
            title_suffix="PCM 298 K 6-311++G(2df,2p)",
            nproc=16,
        )
        assert len(com_log) == n

        # One .com per conformer, named ribose_c{ii}_F.com contiguously from c00.
        names = sorted(os.path.basename(p) for p in com_log["com_path"])
        assert names == [f"ribose_c{ii:02d}_F.com" for ii in range(n)]

        # Each file: intact Link1 opt→freq contract + ΔE (kcal/mol) in the title.
        for com_path in com_log["com_path"]:
            with open(com_path) as f:
                text = f.read()
            assert "--Link1--" in text
            assert "Geom=AllChk Guess=Read" in text
            assert "kcal/mol" in text
