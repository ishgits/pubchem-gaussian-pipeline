"""Tests for pipeline.conformers

The pure ranking helper (`select_top_n`) is tested with no RDKit. The embedding
and batch tests require RDKit but are still offline / no-cluster (allowed by the
no-network test rule); they `importorskip` so a bare environment still passes the
pure tests.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from manifest_helpers import ensure_manifest

from pipeline.conformers import (
    UNCONVERGED_FF_SEED,
    _carry_forward_group_is_valid,
    _finalize_convergence,
    _group_identity_is_consistent,
    _resume_group_is_complete,
    _resume_partition,
    _row_config_matches,
    _row_identity_matches,
    _row_xyz_present,
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
        manifest_path = ensure_manifest(tmp_path, table)
        log = search_conformers(
            table,
            xyz_dir=str(tmp_path / "conf_xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "conformer_search_failed.csv"),
            manifest_path=manifest_path,
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
        manifest_path = ensure_manifest(tmp_path, table)
        log = search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(failed_csv),
            manifest_path=manifest_path,
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
        manifest_path = ensure_manifest(tmp_path, table)
        log = search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(failed_csv),
            manifest_path=manifest_path,
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

    def test_resume_skips_completed(self, tmp_path, monkeypatch):
        pytest.importorskip("rdkit")
        import pipeline.conformers as C

        monkeypatch.setattr(
            C, "pipeline_provenance", lambda: ("2.0.0", "test-clean-commit")
        )

        log_csv = tmp_path / "conformer_log.csv"
        xyz_dir = tmp_path / "conf_xyz"
        table = self._table()
        manifest_path = ensure_manifest(
            tmp_path,
            table,
            pipeline_version="2.0.0",
            pipeline_commit="test-clean-commit",
        )

        first = C.search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(log_csv),
            failed_csv=str(tmp_path / "fail.csv"),
            manifest_path=manifest_path,
        )
        # Rerun with the same table: nothing new appended (both already done).
        second = C.search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(log_csv),
            failed_csv=str(tmp_path / "fail.csv"),
            manifest_path=manifest_path,
        )
        assert len(second) == len(first)


# ---------------------------------------------------------------------------
# Task 5 — reproducibility (validation ≠ "it ran")
# ---------------------------------------------------------------------------

class TestReproducibility:
    def test_same_seed_same_selected_conformers(self, tmp_path):
        pytest.importorskip("rdkit")
        from pipeline.conformers import search_conformers

        import pandas as pd
        table = pd.DataFrame([
            {"name": "Ribose", "cid": 5779, "IsomericSMILES": RIBOSE_SMILES},
        ])

        run_a = tmp_path / "run_a"
        run_b = tmp_path / "run_b"
        run_a.mkdir()
        run_b.mkdir()
        manifest_a = ensure_manifest(run_a, table)
        manifest_b = ensure_manifest(run_b, table)

        log_a = search_conformers(
            table,
            xyz_dir=str(run_a / "xyz"),
            log_csv=str(run_a / "log.csv"),
            failed_csv=str(run_a / "fail.csv"),
            seed=42,
            manifest_path=manifest_a,
        )
        log_b = search_conformers(
            table,
            xyz_dir=str(run_b / "xyz"),
            log_csv=str(run_b / "log.csv"),
            failed_csv=str(run_b / "fail.csv"),
            seed=42,
            manifest_path=manifest_b,
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
        route_opt = "# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)"
        route_freq = "# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read"
        manifest_path = ensure_manifest(
            tmp_path,
            df,
            route_opt=route_opt,
            route_freq=route_freq,
            title_suffix="PCM 298 K 6-311++G(2df,2p)",
        )
        conf_log = search_conformers(
            df,
            xyz_dir=str(tmp_path / "conformer_xyz"),
            log_csv=str(conf_log_csv),
            failed_csv=str(tmp_path / "conformer_search_failed.csv"),
            n_generate=20,
            top_n=3,
            rmsd_prune=0.5,
            seed=42,
            manifest_path=manifest_path,
        )
        n = len(conf_log)
        assert 1 <= n <= 3  # ribose (flexible) → up to TOP_N distinct conformers

        com_log = write_gaussian_coms_from_conformers(
            conformer_log_csv=str(conf_log_csv),
            outdir=str(tmp_path / "gaussian_inputs"),
            log_csv=str(tmp_path / "com_write_log.csv"),
            route_opt=route_opt,
            route_freq=route_freq,
            title_suffix="PCM 298 K 6-311++G(2df,2p)",
            nproc=16,
            manifest_path=manifest_path,
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
        manifest_path = ensure_manifest(
            tmp_path,
            table,
            top_n=top_n,
            rdkit_version="test-rdkit",
        )
        self._manifest_path = manifest_path
        return C.search_conformers(
            table,
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "fail.csv"),
            top_n=top_n,
            manifest_path=manifest_path,
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
            manifest_path=self._manifest_path,
        )
        assert len(com_log) == 1
        with open(com_log.iloc[0]["com_path"]) as f:
            com_text = f.read()
        assert UNCONVERGED_FF_SEED in com_text
        assert "--Link1--" in com_text  # Link1 contract still intact


class TestProvenanceLogging:
    """M-06: pipeline version + commit recorded in conformer_log.csv and XYZ.

    RDKit is stubbed out (synthetic/offline); the real `pipeline_provenance`
    runs, so `pipeline_commit` is whatever git reports here (may be empty) — we
    never assert a concrete SHA."""

    def _run(self, tmp_path, monkeypatch):
        import pandas as pd

        import pipeline.conformers as C

        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        coords = [[("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.96)]]
        monkeypatch.setattr(
            C, "generate_conformers",
            lambda smiles, **kw: (coords, [1.23], "MMFF94", [True]),
        )
        table = pd.DataFrame([{"name": "Water", "cid": 962, "IsomericSMILES": "O"}])
        manifest_path = ensure_manifest(
            tmp_path,
            table,
            rdkit_version="test-rdkit",
        )
        return C.search_conformers(
            table,
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "fail.csv"),
            manifest_path=manifest_path,
        )

    def test_pipeline_version_on_every_row(self, tmp_path, monkeypatch):
        import pipeline

        log = self._run(tmp_path, monkeypatch)
        assert len(log) >= 1
        assert (log["pipeline_version"] == pipeline.__version__).all()

    def test_pipeline_commit_present_and_str(self, tmp_path, monkeypatch):
        # Present and a string (possibly empty) — never assert a concrete SHA.
        log = self._run(tmp_path, monkeypatch)
        assert "pipeline_commit" in log.columns
        for v in log["pipeline_commit"]:
            assert isinstance(v, str)

    def test_xyz_comment_has_provenance_tokens(self, tmp_path, monkeypatch):
        log = self._run(tmp_path, monkeypatch)
        with open(log.iloc[0]["xyz_path"]) as f:
            comment = f.read().splitlines()[1]  # line 2 of an XYZ file is the comment
        for token in (
            "run_id=", "artifact_id=", "config_hash=", "conformer_id=",
            "relative_energy_kcalmol=", "method=", "pipeline_version=",
            "rdkit_version=test-rdkit",
        ):
            assert token in comment


# ---------------------------------------------------------------------------
# M-09 — resume must validate recorded config, not just the molecule name
# ---------------------------------------------------------------------------

_CFG = {"seed": 42, "n_generate": 20, "top_n": 3, "rmsd_prune": 0.5,
        "pipeline_version": "0.2.0", "pipeline_commit": "abc1234",
        "rdkit_version": "2024.03.1"}


class TestRowConfigMatches:
    """Pure config-match predicate (no RDKit)."""

    def _row(self, **override):
        row = dict(_CFG)
        row.update(override)
        return row

    def test_exact_match(self):
        assert _row_config_matches(self._row(), _CFG) is True

    def test_seed_mismatch(self):
        assert _row_config_matches(self._row(seed=7), _CFG) is False

    def test_n_generate_mismatch(self):
        assert _row_config_matches(self._row(n_generate=50), _CFG) is False

    def test_top_n_mismatch(self):
        assert _row_config_matches(self._row(top_n=1), _CFG) is False

    def test_pipeline_version_mismatch(self):
        assert _row_config_matches(self._row(pipeline_version="0.1.0"), _CFG) is False

    def test_rdkit_version_mismatch(self):
        # B-02: ETKDGv3+MMFF geometry is RDKit-version-dependent.
        assert _row_config_matches(self._row(rdkit_version="2020.09.1"), _CFG) is False

    def test_clean_pipeline_commit_mismatch(self):
        assert _row_config_matches(
            self._row(pipeline_commit="def5678"), _CFG
        ) is False

    def test_missing_pipeline_commit_disables_reuse(self):
        row = self._row()
        del row["pipeline_commit"]
        assert _row_config_matches(row, _CFG) is False

    @pytest.mark.parametrize(
        ("row_commit", "config_commit"),
        [
            ("abc1234.dirty", "abc1234"),
            ("abc1234", "abc1234.dirty"),
            ("abc1234.dirty", ""),
            ("", "abc1234.dirty"),
        ],
    )
    def test_dirty_pipeline_commit_never_reuses(self, row_commit, config_commit):
        row = self._row(pipeline_commit=row_commit)
        config = dict(_CFG, pipeline_commit=config_commit)
        assert _row_config_matches(row, config) is False

    def test_rmsd_within_float_tolerance(self):
        assert _row_config_matches(self._row(rmsd_prune=0.5 + 1e-12), _CFG) is True

    def test_rmsd_mismatch(self):
        assert _row_config_matches(self._row(rmsd_prune=0.75), _CFG) is False

    def test_missing_column_is_mismatch(self):
        row = self._row()
        del row["top_n"]  # pre-provenance / older-schema log
        assert _row_config_matches(row, _CFG) is False

    def test_missing_rdkit_version_is_mismatch(self):
        row = self._row()
        del row["rdkit_version"]  # log predates the rdkit_version guard
        assert _row_config_matches(row, _CFG) is False


class TestRowIdentityAndXyz:
    """Pure per-molecule identity + XYZ-existence predicates (B-02, no RDKit)."""

    def test_cid_matches_across_int_and_float_string(self):
        assert _row_identity_matches({"cid": 190, "smiles": "O"}, 190.0, "O") is True
        assert _row_identity_matches({"cid": "190", "smiles": "O"}, 190, "O") is True

    def test_changed_cid_mismatch(self):
        assert _row_identity_matches({"cid": 190, "smiles": "O"}, 191, "O") is False

    def test_changed_smiles_mismatch(self):
        assert _row_identity_matches({"cid": 190, "smiles": "O"}, 190, "[OH2]") is False

    def test_xyz_present_true_for_nonempty(self, tmp_path):
        p = tmp_path / "a.xyz"
        p.write_text("3\n\nO 0 0 0\n")
        assert _row_xyz_present({"xyz_path": str(p)}) is True

    def test_xyz_present_false_for_missing(self, tmp_path):
        assert _row_xyz_present({"xyz_path": str(tmp_path / "nope.xyz")}) is False

    def test_xyz_present_false_for_empty(self, tmp_path):
        p = tmp_path / "empty.xyz"
        p.write_text("")
        assert _row_xyz_present({"xyz_path": str(p)}) is False


class TestResumePartition:
    """Pure partition of an existing log into done / kept / stale (no RDKit)."""

    def _rows_with_xyz(self, tmp_path, *specs):
        """Build rows with a real, non-empty xyz_path per spec (name, **override)."""
        import pandas as pd

        rows = []
        counts = {name: sum(spec_name == name for spec_name, _ in specs) for name, _ in specs}
        next_id = {name: 0 for name, _ in specs}
        for i, (name, override) in enumerate(specs):
            xyz = tmp_path / f"{name}_{i}.xyz"
            xyz.write_text("1\n\nO 0 0 0\n")
            row = dict(_CFG, name=name, cid=1, smiles="O",
                       pipeline_commit="abc1234",
                       conformer_id=next_id[name], n_kept=counts[name],
                       xyz_path=str(xyz))
            next_id[name] += 1
            row.update(override)
            rows.append(row)
        return pd.DataFrame(rows)

    def _requested(self, *names):
        return {n: {"cid": 1, "smiles": "O"} for n in names}

    def test_partition_default_drops_unrequested(self, tmp_path):
        existing = self._rows_with_xyz(
            tmp_path,
            ("A", {}),            # requested, matches
            ("B", {"seed": 7}),   # requested, config drift
            ("C", {}),            # NOT requested
        )
        done, kept, stale, invalid = _resume_partition(
            existing, _CFG, self._requested("A", "B")
        )
        assert done == {"A"}
        assert stale == {"B"}
        assert invalid == {}
        # M-02 default: unrequested C is dropped; only resumed A is kept.
        assert {r["name"] for r in kept} == {"A"}

    def test_preserve_unrequested_keeps_carry_forward(self, tmp_path):
        existing = self._rows_with_xyz(
            tmp_path, ("A", {}), ("C", {}),
        )
        done, kept, stale, invalid = _resume_partition(
            existing, _CFG, self._requested("A"), preserve_unrequested=True
        )
        assert done == {"A"}
        assert invalid == {}
        # append=True: unrequested C carried forward alongside resumed A.
        assert {r["name"] for r in kept} == {"A", "C"}

    def test_changed_identity_is_stale(self, tmp_path):
        # Same name A, but recorded cid differs from requested → stale.
        existing = self._rows_with_xyz(tmp_path, ("A", {"cid": 999}))
        done, kept, stale, invalid = _resume_partition(
            existing, _CFG, self._requested("A")
        )
        assert done == set()
        assert stale == {"A"}
        assert invalid == {}

    def test_missing_xyz_is_stale(self, tmp_path):
        import pandas as pd

        row = dict(_CFG, name="A", cid=1, smiles="O",
                   conformer_id=0, xyz_path=str(tmp_path / "gone.xyz"))
        done, kept, stale, invalid = _resume_partition(
            pd.DataFrame([row]), _CFG, self._requested("A")
        )
        assert stale == {"A"}
        assert done == set()
        assert invalid == {}

    def test_all_rows_of_a_molecule_must_match(self, tmp_path):
        # Two rows for A: one matches, one drifted → A is stale (regenerate all).
        existing = self._rows_with_xyz(
            tmp_path, ("A", {}), ("A", {"top_n": 1}),
        )
        done, kept, stale, invalid = _resume_partition(
            existing, _CFG, self._requested("A")
        )
        assert done == set()
        assert stale == {"A"}
        assert kept == []
        assert invalid == {}

    def test_invalid_unrequested_is_reported_not_kept(self, tmp_path):
        existing = self._rows_with_xyz(
            tmp_path, ("A", {}), ("C", {"seed": 7}),
        )
        done, kept, stale, invalid = _resume_partition(
            existing, _CFG, self._requested("A"), preserve_unrequested=True
        )
        assert done == {"A"}
        assert stale == set()
        assert {row["name"] for row in kept} == {"A"}
        assert set(invalid) == {"C"}


class TestResumeGroupCompleteness:
    """M-12: resume only complete, self-consistent conformer groups."""

    def _rows(self, tmp_path):
        rows = []
        for conformer_id in range(3):
            xyz = tmp_path / f"a_{conformer_id}.xyz"
            xyz.write_text("1\n\nO 0 0 0\n")
            rows.append(dict(
                _CFG,
                name="A",
                cid=1,
                smiles="O",
                conformer_id=conformer_id,
                n_kept=3,
                xyz_path=str(xyz),
            ))
        return rows

    def test_intact_group_is_complete(self, tmp_path):
        assert _resume_group_is_complete(self._rows(tmp_path)) is True

    def test_deleted_last_row_is_incomplete(self, tmp_path):
        rows = self._rows(tmp_path)[:-1]
        assert _resume_group_is_complete(rows) is False

    def test_missing_middle_id_is_incomplete(self, tmp_path):
        rows = self._rows(tmp_path)
        assert _resume_group_is_complete([rows[0], rows[2]]) is False

    def test_duplicate_conformer_row_is_incomplete(self, tmp_path):
        rows = self._rows(tmp_path)
        rows[1] = dict(rows[0])
        assert _resume_group_is_complete(rows) is False

    def test_disagreeing_n_kept_is_incomplete(self, tmp_path):
        rows = self._rows(tmp_path)
        rows[1]["n_kept"] = 2
        assert _resume_group_is_complete(rows) is False

    def test_duplicate_xyz_path_is_incomplete(self, tmp_path):
        rows = self._rows(tmp_path)
        rows[1]["xyz_path"] = rows[0]["xyz_path"]
        assert _resume_group_is_complete(rows) is False


class TestCarryForwardGroupValidation:
    """M-15: retained groups need integrity, identity, config, and provenance."""

    def _rows(self, tmp_path):
        rows = []
        for conformer_id in range(2):
            xyz = tmp_path / f"retained_{conformer_id}.xyz"
            xyz.write_text("1\n\nO 0 0 0\n")
            rows.append(dict(
                _CFG,
                name="Retained",
                cid=1,
                smiles="O",
                pipeline_commit="abc1234",
                conformer_id=conformer_id,
                n_kept=2,
                xyz_path=str(xyz),
            ))
        return rows

    def test_valid_group_passes(self, tmp_path):
        rows = self._rows(tmp_path)
        assert _group_identity_is_consistent(rows) is True
        assert _carry_forward_group_is_valid(rows, _CFG) is True

    @pytest.mark.parametrize("field,value", [("cid", 2), ("smiles", "N")])
    def test_inconsistent_identity_fails(self, tmp_path, field, value):
        rows = self._rows(tmp_path)
        rows[1][field] = value
        assert _group_identity_is_consistent(rows) is False
        assert _carry_forward_group_is_valid(rows, _CFG) is False

    def test_missing_commit_field_fails(self, tmp_path):
        rows = self._rows(tmp_path)
        for row in rows:
            del row["pipeline_commit"]
        assert _carry_forward_group_is_valid(rows, _CFG) is False

    def test_blank_commit_value_disables_reuse(self, tmp_path):
        rows = self._rows(tmp_path)
        for row in rows:
            row["pipeline_commit"] = ""
        assert _carry_forward_group_is_valid(rows, _CFG) is False

    def test_dirty_commit_group_is_not_reusable(self, tmp_path):
        rows = self._rows(tmp_path)
        config = dict(_CFG, pipeline_commit="abc1234.dirty")
        for row in rows:
            row["pipeline_commit"] = "abc1234.dirty"
        assert _carry_forward_group_is_valid(rows, config) is False


class TestResumeConfigValidationBatch:
    """search_conformers end-to-end with RDKit stubbed (synthetic/offline)."""

    def _patch(
        self,
        monkeypatch,
        calls,
        rdkit_version="test-rdkit",
        pipeline_commit="test-clean-commit",
    ):
        import pipeline.conformers as C

        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: rdkit_version)
        monkeypatch.setattr(
            C, "pipeline_provenance", lambda: ("2.0.0", pipeline_commit)
        )
        self._rdkit_version = rdkit_version
        self._pipeline_commit = pipeline_commit

        def gen(smiles, **kw):
            calls.append(smiles)
            return ([[("O", 0.0, 0.0, 0.0)]], [0.0], "MMFF94", [True])

        monkeypatch.setattr(C, "generate_conformers", gen)
        return C

    def _table(self, name="Water", smiles="O", cid=1):
        import pandas as pd

        return pd.DataFrame([{"name": name, "cid": cid, "IsomericSMILES": smiles}])

    def _kw(self, tmp_path, manifest_table=None, **over):
        kw = dict(
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "fail.csv"),
            seed=42, n_generate=20, top_n=3, rmsd_prune=0.5,
        )
        kw.update(over)
        table = self._table() if manifest_table is None else manifest_table
        kw["manifest_path"] = ensure_manifest(
            tmp_path,
            table,
            seed=kw["seed"],
            n_generate=kw["n_generate"],
            top_n=kw["top_n"],
            rmsd_prune=kw["rmsd_prune"],
            pipeline_version="2.0.0",
            pipeline_commit=self._pipeline_commit,
            rdkit_version=self._rdkit_version,
        )
        return kw

    def test_matching_config_skips_regeneration(self, tmp_path, monkeypatch):
        calls = []
        C = self._patch(monkeypatch, calls)
        C.search_conformers(self._table(), **self._kw(tmp_path))
        C.search_conformers(self._table(), **self._kw(tmp_path))  # same config
        assert len(calls) == 1  # second run resumed; molecule NOT regenerated

    @pytest.mark.parametrize(
        "damaged_value", [None, float("nan"), "", "unknown", False]
    )
    def test_damaged_convergence_never_resumes(
        self, tmp_path, monkeypatch, damaged_value
    ):
        calls = []
        C = self._patch(monkeypatch, calls)
        kw = self._kw(tmp_path)
        C.search_conformers(self._table(), **kw)
        rows = pd.read_csv(kw["log_csv"]).astype({"converged": object})
        rows.loc[0, "converged"] = damaged_value
        rows.to_csv(kw["log_csv"], index=False)

        C.search_conformers(self._table(), **kw)
        assert len(calls) == 2

    def test_missing_convergence_column_never_resumes(
        self, tmp_path, monkeypatch
    ):
        calls = []
        C = self._patch(monkeypatch, calls)
        kw = self._kw(tmp_path)
        C.search_conformers(self._table(), **kw)
        rows = pd.read_csv(kw["log_csv"]).drop(columns=["converged"])
        rows.to_csv(kw["log_csv"], index=False)

        C.search_conformers(self._table(), **kw)
        assert len(calls) == 2

    def test_changed_seed_regenerates(self, tmp_path, monkeypatch, capsys):
        calls = []
        C = self._patch(monkeypatch, calls)
        C.search_conformers(self._table(), **self._kw(tmp_path, seed=42))
        log2 = C.search_conformers(self._table(), **self._kw(tmp_path, seed=7))
        assert len(calls) == 2  # regenerated under the new seed
        assert set(log2[log2["name"] == "Water"]["seed"].astype(int)) == {7}
        assert "regenerating" in capsys.readouterr().out

    def test_changed_top_n_regenerates(self, tmp_path, monkeypatch):
        calls = []
        C = self._patch(monkeypatch, calls)
        C.search_conformers(self._table(), **self._kw(tmp_path, top_n=1))
        C.search_conformers(self._table(), **self._kw(tmp_path, top_n=3))
        assert len(calls) == 2  # the TOP_N=1→3 case from the finding

    def test_changed_cid_regenerates(self, tmp_path, monkeypatch):
        # B-02: same name + knobs but a corrected CID must NOT reuse the geometry.
        calls = []
        C = self._patch(monkeypatch, calls)
        first_table = self._table(cid=1)
        second_table = self._table(cid=2)
        C.search_conformers(first_table, **self._kw(tmp_path, manifest_table=first_table))
        log2 = C.search_conformers(second_table, **self._kw(tmp_path, manifest_table=second_table))
        assert len(calls) == 2
        assert set(log2[log2["name"] == "Water"]["cid"].astype(int)) == {2}

    def test_changed_smiles_regenerates(self, tmp_path, monkeypatch):
        # B-02: same name + knobs but a corrected SMILES must regenerate.
        calls = []
        C = self._patch(monkeypatch, calls)
        first_table = self._table(smiles="O")
        second_table = self._table(smiles="[OH2]")
        C.search_conformers(first_table, **self._kw(tmp_path, manifest_table=first_table))
        C.search_conformers(second_table, **self._kw(tmp_path, manifest_table=second_table))
        assert len(calls) == 2

    def test_changed_rdkit_version_regenerates(self, tmp_path, monkeypatch):
        # B-02: a different RDKit build changes ETKDGv3+MMFF geometry → regenerate.
        calls = []
        C = self._patch(monkeypatch, calls, rdkit_version="2024.03.1")
        C.search_conformers(self._table(), **self._kw(tmp_path))
        C = self._patch(monkeypatch, calls, rdkit_version="2020.09.1")
        C.search_conformers(self._table(), **self._kw(tmp_path))
        assert len(calls) == 2

    def test_changed_clean_pipeline_commit_regenerates(self, tmp_path, monkeypatch):
        calls = []
        C = self._patch(monkeypatch, calls, pipeline_commit="abc1234")
        C.search_conformers(self._table(), **self._kw(tmp_path))
        C = self._patch(monkeypatch, calls, pipeline_commit="def5678")
        C.search_conformers(self._table(), **self._kw(tmp_path))
        assert len(calls) == 2

    def test_dirty_pipeline_commit_never_resumes(self, tmp_path, monkeypatch):
        calls = []
        C = self._patch(monkeypatch, calls, pipeline_commit="abc1234.dirty")
        C.search_conformers(self._table(), **self._kw(tmp_path))
        C.search_conformers(self._table(), **self._kw(tmp_path))
        assert len(calls) == 2

    def test_missing_pipeline_commit_regenerates(
        self, tmp_path, monkeypatch
    ):
        calls = []
        C = self._patch(monkeypatch, calls, pipeline_commit="")
        C.search_conformers(self._table(), **self._kw(tmp_path))
        C.search_conformers(self._table(), **self._kw(tmp_path))
        assert len(calls) == 2

    def test_deleted_xyz_regenerates(self, tmp_path, monkeypatch):
        # B-02: identity + config match, but the recorded XYZ is gone → regenerate.
        calls = []
        C = self._patch(monkeypatch, calls)
        log1 = C.search_conformers(self._table(), **self._kw(tmp_path))
        os.remove(log1.iloc[0]["xyz_path"])  # user cleaned the geometry away
        C.search_conformers(self._table(), **self._kw(tmp_path))
        assert len(calls) == 2

    def test_pre_provenance_log_regenerates(self, tmp_path, monkeypatch):
        import pandas as pd

        calls = []
        C = self._patch(monkeypatch, calls)
        log_csv = tmp_path / "conformer_log.csv"
        # Old-schema log: no n_generate/top_n/pipeline_version columns.
        pd.DataFrame([{
            "name": "Water", "cid": 1, "conformer_id": 0,
            "seed": 42, "rmsd_prune": 0.5, "xyz_path": "old.xyz",
        }]).to_csv(log_csv, index=False)
        C.search_conformers(self._table(), **self._kw(tmp_path))
        assert len(calls) == 1  # name present but stale → regenerated

    def test_truncated_three_conformer_log_regenerates(self, tmp_path, monkeypatch):
        import pandas as pd

        calls = []
        C = self._patch(monkeypatch, calls)

        def gen_three(smiles, **kw):
            calls.append(smiles)
            coords = [[("O", float(i), 0.0, 0.0)] for i in range(3)]
            return coords, [0.0, 0.5, 1.0], "MMFF94", [True, True, True]

        monkeypatch.setattr(C, "generate_conformers", gen_three)
        kw = self._kw(tmp_path)
        first = C.search_conformers(self._table(), **kw)
        assert len(first) == 3

        # Simulate a truncated but still parseable log: surviving rows and XYZ
        # files are individually valid, but the group declares n_kept=3.
        first.iloc[[0, 2]].to_csv(kw["log_csv"], index=False)
        second = C.search_conformers(self._table(), **kw)
        assert len(calls) == 2
        assert list(second["conformer_id"]) == [0, 1, 2]
        assert set(second["n_kept"].astype(int)) == {3}

    def test_second_staged_xyz_failure_publishes_no_partial_group(
        self, tmp_path, monkeypatch
    ):
        calls = []
        C = self._patch(monkeypatch, calls)

        def gen_three(smiles, **kw):
            calls.append(smiles)
            coords = [[("O", float(i), 0.0, 0.0)] for i in range(3)]
            return coords, [0.0, 0.5, 1.0], "MMFF94", [True, True, True]

        monkeypatch.setattr(C, "generate_conformers", gen_three)
        kw = self._kw(tmp_path)
        first = C.search_conformers(self._table(), **kw)
        assert len(first) == 3
        manifest_before = open(kw["manifest_path"], "rb").read()
        xyz_before = {
            path: open(path, "rb").read() for path in first["xyz_path"]
        }

        damaged = first.copy()
        damaged.loc[0, "xyz_sha256"] = "0" * 64
        damaged.to_csv(kw["log_csv"], index=False)
        real_write_xyz = C._write_xyz
        staged_writes = 0

        def fail_second(path, coords, comment=""):
            nonlocal staged_writes
            staged_writes += 1
            if staged_writes == 2:
                raise OSError("simulated staged XYZ failure")
            return real_write_xyz(path, coords, comment)

        monkeypatch.setattr(C, "_write_xyz", fail_second)
        second = C.search_conformers(self._table(), **kw)

        assert second.empty
        assert open(kw["manifest_path"], "rb").read() == manifest_before
        for path, expected in xyz_before.items():
            assert open(path, "rb").read() == expected
        assert not any(
            name.startswith(".staging-") for name in os.listdir(kw["xyz_dir"])
        )
        assert pd.read_csv(kw["failed_csv"]).shape[0] == 1


class TestPreserveUnrequestedBatch:
    """M-02 call 2a: the conformer log represents molecules requested this run."""

    def _patch(self, monkeypatch, calls):
        import pipeline.conformers as C

        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        monkeypatch.setattr(
            C, "pipeline_provenance", lambda: ("2.0.0", "test-clean-commit")
        )

        def gen(smiles, **kw):
            calls.append(smiles)
            return ([[("O", 0.0, 0.0, 0.0)]], [0.0], "MMFF94", [True])

        monkeypatch.setattr(C, "generate_conformers", gen)
        return C

    def _kw(self, tmp_path, **over):
        import pandas as pd

        manifest_table = over.pop("manifest_table", None)
        kw = dict(
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "fail.csv"),
            seed=42, n_generate=20, top_n=3, rmsd_prune=0.5,
        )
        kw.update(over)
        if manifest_table is None:
            manifest_table = pd.concat([self._two(), self._one()], ignore_index=True)
        kw["manifest_path"] = ensure_manifest(
            tmp_path,
            manifest_table,
            seed=kw["seed"],
            n_generate=kw["n_generate"],
            top_n=kw["top_n"],
            rmsd_prune=kw["rmsd_prune"],
            pipeline_version="2.0.0",
            pipeline_commit="test-clean-commit",
            rdkit_version="test-rdkit",
        )
        return kw

    def _two(self):
        import pandas as pd

        return pd.DataFrame([
            {"name": "Water", "cid": 1, "IsomericSMILES": "O"},
            {"name": "Glycine", "cid": 2, "IsomericSMILES": "O"},
        ])

    def _one(self):
        import pandas as pd

        return pd.DataFrame([{"name": "Adenine", "cid": 3, "IsomericSMILES": "O"}])

    def test_default_subset_rejected(self, tmp_path, monkeypatch):
        import pandas as pd

        calls = []
        C = self._patch(monkeypatch, calls)
        full = pd.concat([self._two(), self._one()], ignore_index=True)
        kw = self._kw(tmp_path, manifest_table=full)
        C.search_conformers(full, **kw)
        before = (tmp_path / "conformer_log.csv").read_bytes()
        with pytest.raises(ValueError, match="match the run manifest exactly"):
            C.search_conformers(self._one(), **kw)
        assert (tmp_path / "conformer_log.csv").read_bytes() == before

    def test_append_retains_all(self, tmp_path, monkeypatch):
        import pandas as pd

        calls = []
        C = self._patch(monkeypatch, calls)
        full = pd.concat([self._two(), self._one()], ignore_index=True)
        kw = self._kw(tmp_path, manifest_table=full)
        C.search_conformers(full, **kw)
        log2 = C.search_conformers(self._one(), **dict(kw, append=True))
        # append=True may operate on a subset only when valid retained rows account
        # for every other manifest molecule.
        assert set(log2["name"]) == {"Water", "Glycine", "Adenine"}

    @pytest.mark.parametrize(
        "corruption",
        [
            "truncated",
            "missing_xyz",
            "seed",
            "n_generate",
            "top_n",
            "rmsd_prune",
            "pipeline_version",
            "pipeline_commit",
            "rdkit_version",
            "pre_provenance",
            "inconsistent_cid",
            "inconsistent_smiles",
        ],
    )
    def test_invalid_retained_group_fails_before_any_mutation(
        self, tmp_path, monkeypatch, corruption
    ):
        import pandas as pd

        calls = []
        C = self._patch(monkeypatch, calls)
        full = pd.concat([self._two(), self._one()], ignore_index=True)
        kw = self._kw(tmp_path, manifest_table=full)
        C.search_conformers(full, **kw)
        log_path = tmp_path / "conformer_log.csv"
        existing = pd.read_csv(log_path)
        water_index = existing.index[existing["name"] == "Water"][0]

        if corruption == "truncated":
            existing.loc[water_index, "n_kept"] = 2
        elif corruption == "missing_xyz":
            os.remove(existing.loc[water_index, "xyz_path"])
        elif corruption == "pre_provenance":
            existing = existing.drop(
                columns=["pipeline_version", "rdkit_version", "pipeline_commit"]
            )
        elif corruption in {"inconsistent_cid", "inconsistent_smiles"}:
            second = existing.loc[water_index].copy()
            existing.loc[water_index, "n_kept"] = 2
            second["conformer_id"] = 1
            second["n_kept"] = 2
            copied_xyz = tmp_path / "xyz" / "water_c01.xyz"
            copied_xyz.write_text("1\n\nO 1 0 0\n")
            second["xyz_path"] = str(copied_xyz)
            if corruption == "inconsistent_cid":
                second["cid"] = 999
            else:
                second["smiles"] = "N"
            existing = pd.concat([existing, second.to_frame().T], ignore_index=True)
        else:
            mismatches = {
                "seed": 7,
                "n_generate": 50,
                "top_n": 1,
                "rmsd_prune": 0.75,
                "pipeline_version": "1.9.0",
                "pipeline_commit": "other-clean-commit",
                "rdkit_version": "old-rdkit",
            }
            existing.loc[water_index, corruption] = mismatches[corruption]

        if corruption != "missing_xyz":
            existing.to_csv(log_path, index=False)
        else:
            # Preserve the original CSV so it points at the now-missing XYZ.
            assert not os.path.exists(existing.loc[water_index, "xyz_path"])

        failed_path = tmp_path / "fail.csv"
        failed_path.write_bytes(b"prior failure log\n")
        before = {
            path.relative_to(tmp_path): path.read_bytes()
            for path in tmp_path.rglob("*")
            if path.is_file()
        }
        calls_before = len(calls)

        with pytest.raises(ValueError, match="cannot carry forward invalid") as excinfo:
            C.search_conformers(self._one(), **dict(kw, append=True))

        assert "Water" in str(excinfo.value)
        assert len(calls) == calls_before
        after = {
            path.relative_to(tmp_path): path.read_bytes()
            for path in tmp_path.rglob("*")
            if path.is_file()
        }
        assert after == before

    @pytest.mark.parametrize(
        "prior_label,current_label",
        [("A+B", "A B"), ("Water", "water")],
    )
    def test_append_collision_raises_before_any_mutation(
        self, tmp_path, monkeypatch, prior_label, current_label
    ):
        import pandas as pd

        calls = []
        C = self._patch(monkeypatch, calls)
        first = pd.DataFrame([{
            "name": prior_label,
            "cid": 1,
            "IsomericSMILES": "O",
        }])
        second = pd.DataFrame([{
            "name": current_label,
            "cid": 2,
            "IsomericSMILES": "N",
        }])
        with pytest.raises(ValueError, match="both map to output basename"):
            self._kw(
                tmp_path,
                manifest_table=pd.concat([first, second], ignore_index=True),
            )

        # The manifest boundary now rejects the invalid one-to-one mapping
        # before conformer generation or any run-output mutation begins.
        assert calls == []
        assert not list(tmp_path.glob("run_manifest_*.json"))
        assert not (tmp_path / "conformer_log.csv").exists()
        assert not (tmp_path / "xyz").exists()
        assert not (tmp_path / "fail.csv").exists()

    def test_already_corrupt_append_log_is_rejected_without_mutation(
        self, tmp_path
    ):
        import pandas as pd
        import pipeline.conformers as C

        log_path = tmp_path / "conformer_log.csv"
        pd.DataFrame([{"name": "A+B"}, {"name": "A B"}]).to_csv(
            log_path, index=False
        )
        xyz_dir = tmp_path / "xyz"
        xyz_dir.mkdir()
        xyz_path = xyz_dir / "a_b_c00.xyz"
        xyz_path.write_text("prior geometry\n")
        before_log = log_path.read_bytes()
        before_xyz = xyz_path.read_bytes()
        current = pd.DataFrame([{
            "name": "Adenine",
            "cid": 3,
            "IsomericSMILES": "N",
        }])

        with pytest.raises(ValueError, match="both map to output basename"):
            C.search_conformers(
                current,
                xyz_dir=str(xyz_dir),
                log_csv=str(log_path),
                failed_csv=str(tmp_path / "fail.csv"),
                append=True,
            )

        assert log_path.read_bytes() == before_log
        assert xyz_path.read_bytes() == before_xyz
        assert not (tmp_path / "fail.csv").exists()


class TestStaleFailedCsvCleared:
    """MIN-02: a prior run's *_failed.csv must not survive a later clean run."""

    def test_failed_csv_cleared_on_clean_rerun(self, tmp_path, monkeypatch):
        import pandas as pd

        import pipeline.conformers as C

        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        failed_csv = tmp_path / "conformer_search_failed.csv"
        kw = dict(
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(failed_csv),
        )

        # First run fails eligibility → a failure log is written. Both runs share
        # the shared xyz_dir/log_csv/failed_csv paths, so their manifests must be
        # rooted at tmp_path to keep every output inside the run package (M-30);
        # distinct molecule tables give the two manifests distinct digest names.
        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: "no IsomericSMILES")
        bad_table = pd.DataFrame([{"name": "Bad", "cid": 1, "IsomericSMILES": None}])
        bad_kw = dict(kw, manifest_path=ensure_manifest(
            tmp_path, bad_table, rdkit_version="test-rdkit"
        ))
        C.search_conformers(bad_table, **bad_kw)
        assert failed_csv.exists()

        # Second run is clean → the stale failure log is cleared, not left behind.
        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: None)
        monkeypatch.setattr(
            C, "generate_conformers",
            lambda smiles, **k: ([[("O", 0.0, 0.0, 0.0)]], [0.0], "MMFF94", [True]),
        )
        good_table = pd.DataFrame([{"name": "Good", "cid": 2, "IsomericSMILES": "O"}])
        good_kw = dict(kw, manifest_path=ensure_manifest(
            tmp_path, good_table, rdkit_version="test-rdkit"
        ))
        C.search_conformers(good_table, **good_kw)
        assert not failed_csv.exists()


class TestParameterValidation:
    """MIN-03: invalid parameters raise ValueError at entry."""

    def _table(self, rows=None):
        import pandas as pd

        return pd.DataFrame(rows or [{"name": "Water", "cid": 1, "IsomericSMILES": "O"}])

    def _kw(self, tmp_path, **over):
        kw = dict(
            xyz_dir=str(tmp_path / "xyz"),
            log_csv=str(tmp_path / "conformer_log.csv"),
            failed_csv=str(tmp_path / "fail.csv"),
        )
        kw.update(over)
        return kw

    def test_n_generate_below_one_raises(self, tmp_path):
        from pipeline.conformers import search_conformers

        with pytest.raises(ValueError):
            search_conformers(self._table(), **self._kw(tmp_path, n_generate=0))

    def test_top_n_below_one_raises(self, tmp_path):
        from pipeline.conformers import search_conformers

        with pytest.raises(ValueError):
            search_conformers(self._table(), **self._kw(tmp_path, top_n=0))

    def test_negative_rmsd_prune_raises(self, tmp_path):
        from pipeline.conformers import search_conformers

        with pytest.raises(ValueError):
            search_conformers(self._table(), **self._kw(tmp_path, rmsd_prune=-0.1))

    def test_duplicate_labels_raise(self, tmp_path):
        from pipeline.conformers import search_conformers

        table = self._table([
            {"name": "Water", "cid": 1, "IsomericSMILES": "O"},
            {"name": "Water", "cid": 2, "IsomericSMILES": "O"},
        ])
        with pytest.raises(ValueError):
            search_conformers(table, **self._kw(tmp_path))

    def test_empty_sanitized_label_raises(self, tmp_path):
        from pipeline.conformers import search_conformers

        table = self._table([{"name": "!!!", "cid": 1, "IsomericSMILES": "O"}])
        with pytest.raises(ValueError):
            search_conformers(table, **self._kw(tmp_path))

    @pytest.mark.parametrize(
        "first,second",
        [("A+B", "A B"), ("Water", "water")],
    )
    def test_sanitized_basename_collision_raises_before_writes(
        self, tmp_path, first, second
    ):
        from pipeline.conformers import search_conformers

        table = self._table([
            {"name": first, "cid": 1, "IsomericSMILES": "O"},
            {"name": second, "cid": 2, "IsomericSMILES": "N"},
        ])
        kw = self._kw(tmp_path)
        with pytest.raises(ValueError, match="both map to output basename"):
            search_conformers(table, **kw)
        assert not (tmp_path / "xyz").exists()
        assert not (tmp_path / "conformer_log.csv").exists()
        assert not (tmp_path / "fail.csv").exists()

    def test_distinct_sanitized_basenames_succeed(self, tmp_path, monkeypatch):
        import pipeline.conformers as C

        monkeypatch.setattr(C, "check_conformer_eligibility", lambda smiles: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        monkeypatch.setattr(
            C,
            "generate_conformers",
            lambda smiles, **kw: (
                [[("O", 0.0, 0.0, 0.0)]], [0.0], "MMFF94", [True]
            ),
        )
        table = self._table([
            {"name": "Water", "cid": 1, "IsomericSMILES": "O"},
            {"name": "Ammonia", "cid": 2, "IsomericSMILES": "N"},
        ])
        kw = self._kw(tmp_path)
        kw["manifest_path"] = ensure_manifest(
            tmp_path, table, rdkit_version="test-rdkit"
        )
        out = C.search_conformers(table, **kw)
        assert {os.path.basename(path) for path in out["xyz_path"]} == {
            "water_c00.xyz",
            "ammonia_c00.xyz",
        }


class TestManifestMoleculeCoverage:
    """B-10: runtime molecule inputs must account for the immutable manifest."""

    @staticmethod
    def _patch(monkeypatch):
        import pipeline.conformers as C

        monkeypatch.setattr(C, "check_conformer_eligibility", lambda smiles: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        monkeypatch.setattr(C, "pipeline_provenance", lambda: ("2.0.0", "test-clean-commit"))
        monkeypatch.setattr(
            C,
            "generate_conformers",
            lambda smiles, **kw: (
                [[("O", 0.0, 0.0, 0.0)]], [0.0], "MMFF94", [True]
            ),
        )
        return C

    @staticmethod
    def _tables():
        import pandas as pd

        full = pd.DataFrame([
            {"name": "Water", "cid": 1, "IsomericSMILES": "O"},
            {"name": "Ammonia", "cid": 2, "IsomericSMILES": "N"},
        ])
        return full, full.iloc[[0]].reset_index(drop=True)

    def test_append_false_subset_fails_before_mutation(self, tmp_path, monkeypatch):
        from pathlib import Path

        C = self._patch(monkeypatch)
        full, subset = self._tables()
        manifest_path = ensure_manifest(
            tmp_path,
            full,
            pipeline_version="2.0.0",
            pipeline_commit="test-clean-commit",
            rdkit_version="test-rdkit",
        )
        manifest_before = Path(manifest_path).read_bytes()
        prior_log = tmp_path / "conformer_log.csv"
        prior_log.write_bytes(b"prior log\n")
        prior_failure = tmp_path / "conformer_search_failed.csv"
        prior_failure.write_bytes(b"prior failure\n")

        with pytest.raises(ValueError, match="match the run manifest exactly"):
            C.search_conformers(
                subset,
                xyz_dir=str(tmp_path / "xyz"),
                log_csv=str(prior_log),
                failed_csv=str(prior_failure),
                manifest_path=manifest_path,
            )

        assert Path(manifest_path).read_bytes() == manifest_before
        assert prior_log.read_bytes() == b"prior log\n"
        assert prior_failure.read_bytes() == b"prior failure\n"
        assert not (tmp_path / "xyz").exists()

    def test_append_true_requires_current_plus_retained_complete_manifest(
        self, tmp_path, monkeypatch
    ):
        from pathlib import Path

        C = self._patch(monkeypatch)
        full, subset = self._tables()
        manifest_path = ensure_manifest(
            tmp_path,
            full,
            pipeline_version="2.0.0",
            pipeline_commit="test-clean-commit",
            rdkit_version="test-rdkit",
        )
        manifest_before = Path(manifest_path).read_bytes()

        with pytest.raises(ValueError, match="account for the complete run manifest"):
            C.search_conformers(
                subset,
                xyz_dir=str(tmp_path / "xyz"),
                log_csv=str(tmp_path / "conformer_log.csv"),
                failed_csv=str(tmp_path / "conformer_search_failed.csv"),
                append=True,
                manifest_path=manifest_path,
            )

        assert Path(manifest_path).read_bytes() == manifest_before
        assert not (tmp_path / "xyz").exists()
        assert not (tmp_path / "conformer_log.csv").exists()
        assert not (tmp_path / "conformer_search_failed.csv").exists()


class TestPackageBoundaryPreflight:
    """M-30: every v2 output destination is validated inside the run package
    before the first mutation, so an outside xyz_dir or authoritative conformer
    log fails atomically — no lineage removal, directory creation, failure-log
    deletion, XYZ write, or log rewrite occurs, and RDKit is never reached."""

    def _build_valid_run(self, tmp_path, monkeypatch):
        import pandas as pd

        import pipeline.conformers as C

        monkeypatch.setattr(C, "check_conformer_eligibility", lambda s: None)
        monkeypatch.setattr(C, "_rdkit_version", lambda: "test-rdkit")
        monkeypatch.setattr(
            C, "pipeline_provenance", lambda: ("2.0.0", "test-clean-commit")
        )
        coords = [[("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.96)]]
        monkeypatch.setattr(
            C, "generate_conformers",
            lambda smiles, **kw: (coords, [0.0], "MMFF94", [True]),
        )
        table = pd.DataFrame([{"name": "Water", "cid": 1, "IsomericSMILES": "O"}])
        manifest_path = ensure_manifest(
            tmp_path,
            table,
            pipeline_version="2.0.0",
            pipeline_commit="test-clean-commit",
            rdkit_version="test-rdkit",
        )
        xyz_dir = tmp_path / "xyz"
        log_csv = tmp_path / "conformer_log.csv"
        failed_csv = tmp_path / "conformer_search_failed.csv"
        C.search_conformers(
            table,
            xyz_dir=str(xyz_dir),
            log_csv=str(log_csv),
            failed_csv=str(failed_csv),
            manifest_path=manifest_path,
        )
        return C, table, manifest_path, xyz_dir, log_csv, failed_csv

    @pytest.mark.parametrize("target", ["xyz_dir", "log_csv"])
    def test_outside_package_destination_fails_atomically(
        self, tmp_path, monkeypatch, target
    ):
        from pathlib import Path

        C, table, manifest_path, xyz_dir, log_csv, failed_csv = (
            self._build_valid_run(tmp_path, monkeypatch)
        )
        # An existing failure log from a prior run must survive an aborted rerun.
        failed_csv.write_bytes(b"prior failure\n")

        manifest_before = Path(manifest_path).read_bytes()
        log_before = log_csv.read_bytes()
        failed_before = failed_csv.read_bytes()
        xyz_before = {p: p.read_bytes() for p in xyz_dir.glob("*.xyz")}
        assert xyz_before  # the valid run wrote at least one XYZ

        # Any generation attempt after preflight would be a bug: prove the stage
        # aborts during package-boundary preflight, before RDKit is reached.
        def _fail(*args, **kwargs):
            raise AssertionError("generate_conformers must not run after preflight")

        monkeypatch.setattr(C, "generate_conformers", _fail)

        outside = tmp_path.parent / f"m30_conformer_outside_{tmp_path.name}"
        if target == "xyz_dir":
            call = dict(xyz_dir=str(outside), log_csv=str(log_csv))
        else:
            call = dict(
                xyz_dir=str(xyz_dir),
                log_csv=str(outside / "conformer_log.csv"),
            )

        with pytest.raises(ValueError, match="inside the run package"):
            C.search_conformers(
                table,
                failed_csv=str(failed_csv),
                manifest_path=manifest_path,
                **call,
            )

        assert Path(manifest_path).read_bytes() == manifest_before
        assert log_csv.read_bytes() == log_before
        assert failed_csv.read_bytes() == failed_before
        assert {p: p.read_bytes() for p in xyz_dir.glob("*.xyz")} == xyz_before
        assert not outside.exists()
