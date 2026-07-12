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

from pipeline.conformers import select_top_n


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

    def test_missing_smiles_is_logged_not_skipped(self, tmp_path):
        pd = pytest.importorskip("pandas")
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        table = pd.DataFrame([
            {"name": "NoSmiles", "cid": None, "IsomericSMILES": None},
        ])
        failed_csv = tmp_path / "conformer_search_failed.csv"
        log = search_conformers(
            table,
            xyz_dir=str(tmp_path / "conf_xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(failed_csv),
        )
        assert len(log) == 0
        assert failed_csv.exists()
        fail_df = pd.read_csv(failed_csv)
        assert "NoSmiles" in set(fail_df["name"])

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
