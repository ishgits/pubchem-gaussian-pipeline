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

from pipeline.conformers import (
    UNCONVERGED_FF_SEED,
    _finalize_convergence,
    check_conformer_eligibility,
    select_converged_top_n,
    select_top_n,
)


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

        coords1, energies1, method1, converged1 = generate_conformers("CCCC", seed=42)
        # ≥1 conformer produced
        assert len(coords1) >= 1
        assert len(energies1) == len(coords1)
        assert len(converged1) == len(coords1)
        assert method1 == "MMFF94"  # butane has full MMFF params
        assert all(isinstance(c, bool) for c in converged1)
        assert all(converged1)  # butane conformers converge under the FF
        # atom tuples look like (symbol, x, y, z), coords in Ångström
        sym, x, y, z = coords1[0][0]
        assert isinstance(sym, str)
        assert all(isinstance(v, float) for v in (x, y, z))

        # Deterministic lowest-energy index under the fixed seed.
        coords2, energies2, _, _ = generate_conformers("CCCC", seed=42)
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


# ---------------------------------------------------------------------------
# M-04 — FF convergence: reject/flag unconverged conformers
# ---------------------------------------------------------------------------

class TestSelectConvergedTopN:
    """M-04 1a/2b selection logic — pure, no RDKit."""

    def test_only_converged_are_eligible(self):
        # idx1 and idx3 converged; ranked by energy among the converged set.
        keep, all_failed = select_converged_top_n(
            [5.0, 1.0, 2.0, 0.5], [False, True, False, True], 3
        )
        assert all_failed is False
        assert keep == [3, 1]  # 0.5 then 1.0; unconverged idx0/idx2 excluded

    def test_top_n_caps_converged(self):
        keep, all_failed = select_converged_top_n(
            [1.0, 2.0, 3.0], [True, True, True], 2
        )
        assert (keep, all_failed) == ([0, 1], False)

    def test_all_unconverged_returns_one_best_effort(self):
        keep, all_failed = select_converged_top_n(
            [5.0, 3.0, 9.0], [False, False, False], 3
        )
        assert all_failed is True
        assert keep == [1]  # single lowest-energy best-effort seed

    def test_empty_input(self):
        assert select_converged_top_n([], [], 3) == ([], True)


class TestFinalizeConvergence:
    """M-04 retry-merge logic — pure, no RDKit."""

    def test_first_pass_converged_kept(self):
        assert _finalize_convergence([(0, 3.0)]) == [(True, 3.0)]

    def test_retry_converges_is_included(self):
        # Unconverged first pass (nc=1), converged on retry (nc=0) → included.
        assert _finalize_convergence([(1, 5.0)], retry=[(0, 4.2)]) == [(True, 4.2)]

    def test_still_unconverged_after_retry(self):
        assert _finalize_convergence([(1, 5.0)], retry=[(1, 4.9)]) == [(False, 4.9)]

    def test_mixed_uses_retry_only_for_failed(self):
        out = _finalize_convergence([(0, 1.0), (1, 2.0)], retry=[(0, 9.9), (0, 1.5)])
        # idx0 kept first-pass energy (retry ignored); idx1 took retry result.
        assert out == [(True, 1.0), (True, 1.5)]


class TestRetryAlignment:
    """M-05: the retry must re-optimize ONLY the failed conformers, so converged
    conformers keep the first-pass energy that matches their (untouched) geometry
    — ranked energy never describes a different geometry than the written XYZ."""

    def test_retry_touches_only_failed_confs_and_keeps_energies_aligned(self, monkeypatch):
        pytest.importorskip("rdkit")
        import pipeline.conformers as C

        state = {}

        def fake_confs(mol, method, max_iters):
            # First pass: every conformer converged except the last one.
            conf_ids = [c.GetId() for c in mol.GetConformers()]
            state["conf_ids"] = conf_ids
            res = [(0, float(10 + i)) for i in range(len(conf_ids))]
            res[-1] = (1, 999.0)  # last conf failed to converge
            return res

        retried = []

        def fake_single(mol, method, conf_id, max_iters):
            retried.append(conf_id)
            return (0, 42.0)  # converges on retry with a distinctive energy

        monkeypatch.setattr(C, "_optimize_confs", fake_confs)
        monkeypatch.setattr(C, "_optimize_single_conf", fake_single)

        coords, energies, method, converged = C.generate_conformers("CCCC", seed=1)
        n = len(state["conf_ids"])
        assert n >= 1

        # Exactly one retry call, and only for the conformer that failed.
        assert retried == [state["conf_ids"][-1]]
        # The failed conformer carries its retry energy; every already-converged
        # conformer keeps its first-pass energy (never re-optimized → geometry and
        # energy stay in sync).
        assert energies[-1] == 42.0
        for i in range(n - 1):
            assert energies[i] == float(10 + i)
        assert all(converged)
        assert len(coords) == n

    def test_no_retry_when_all_converge(self, monkeypatch):
        pytest.importorskip("rdkit")
        import pipeline.conformers as C

        def fake_confs(mol, method, max_iters):
            return [(0, float(10 + i)) for i in range(mol.GetNumConformers())]

        retried = []
        monkeypatch.setattr(C, "_optimize_confs", fake_confs)
        monkeypatch.setattr(
            C, "_optimize_single_conf",
            lambda *a, **k: retried.append(a) or (0, 0.0),
        )
        _, energies, _, converged = C.generate_conformers("CCCC", seed=1)
        # No conformer failed → single-conf retry never called.
        assert retried == []
        assert all(converged)


class TestConvergenceBatch:
    """M-04 batch behavior — RDKit stubbed out so this is synthetic/offline."""

    def _run(self, tmp_path, monkeypatch, energies, converged, top_n=3):
        import pandas as pd

        import pipeline.conformers as C

        # Stub every RDKit touchpoint so we control convergence deterministically.
        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        coords = [[("C", 0.0, 0.0, float(i))] for i in range(len(energies))]
        monkeypatch.setattr(
            C,
            "generate_conformers",
            lambda smiles, **kw: (coords, list(energies), "MMFF94", list(converged)),
        )
        table = pd.DataFrame([{"name": "Mol", "cid": 1, "IsomericSMILES": "C"}])
        return C.search_conformers(
            table,
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "fail.csv"),
            top_n=top_n,
        )

    def test_mixed_keeps_only_converged_logged_true(self, tmp_path, monkeypatch):
        # 4 confs, 2 converged → only the 2 converged are carried and logged True.
        log = self._run(
            tmp_path, monkeypatch,
            energies=[1.0, 0.5, 2.0, 0.1], converged=[True, False, True, False],
        )
        assert len(log) == 2
        assert all(bool(c) for c in log["converged"])

    def test_unconverged_seed_logged_converged_false(self, tmp_path, monkeypatch):
        # All-fail branch: the single carried seed is recorded converged=False.
        log = self._run(
            tmp_path, monkeypatch,
            energies=[5.0, 3.0, 9.0], converged=[False, False, False],
        )
        assert len(log) == 1
        assert bool(log.iloc[0]["converged"]) is False

    def test_all_unconverged_one_flagged_seed_with_marker(self, tmp_path, monkeypatch, capsys):
        from pipeline.gaussian import write_gaussian_coms_from_conformers

        log = self._run(
            tmp_path, monkeypatch,
            energies=[5.0, 3.0, 9.0], converged=[False, False, False],
        )
        # (2b) exactly one best-effort geometry carried, not three.
        assert len(log) == 1
        row = log.iloc[0]
        assert int(row["conformer_id"]) == 0
        assert bool(row["converged"]) is False

        # A warning was emitted naming the unconverged best-effort seed.
        out = capsys.readouterr().out
        assert "no conformer converged" in out
        assert UNCONVERGED_FF_SEED in out

        # The .com title carries the UNCONVERGED_FF_SEED marker on inspection.
        com_log = write_gaussian_coms_from_conformers(
            conformer_log_csv=str(tmp_path / "conformer_log.csv"),
            outdir=str(tmp_path / "gaussian_inputs"),
            log_csv=str(tmp_path / "com_write_log.csv"),
            route_opt="# opt b3lyp/6-31g(d)",
            route_freq="# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
        )
        assert len(com_log) == 1
        with open(com_log.iloc[0]["com_path"]) as f:
            com_text = f.read()
        assert UNCONVERGED_FF_SEED in com_text
        assert "--Link1--" in com_text  # Link1 contract still intact
