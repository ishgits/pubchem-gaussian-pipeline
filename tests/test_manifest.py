"""Offline tests for the authoritative v2 run manifest contract."""

from __future__ import annotations

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess

import pandas as pd
import pytest

import pipeline.manifest as manifest_module
from pipeline.utils import pipeline_provenance

from pipeline.manifest import (
    artifact_abspath,
    canonical_json,
    configuration_hash,
    create_run_manifest,
    finalize_manifest,
    load_manifest,
    molecule_identity_hash,
    record_conformer_group,
    record_child_artifact,
    slurm_template_identity,
    stable_record_id,
    validate_manifest,
    verify_artifact,
)
from pipeline.slurm import DEFAULT_TEMPLATE


def _table(order=("Water", "Ammonia")):
    identities = {
        "Water": {"name": "Water", "cid": 962, "IsomericSMILES": "O"},
        "Ammonia": {"name": "Ammonia", "cid": 222, "IsomericSMILES": "N"},
    }
    return pd.DataFrame([identities[name] for name in order])


def _configs():
    return (
        {
            "method_policy": "MMFF94 preferred; UFF recorded fallback",
            "seed": 42,
            "n_generate": 20,
            "top_n": 3,
            "rmsd_prune": 0.5,
        },
        {
            "route_opt": "# opt b3lyp/6-31g(d)",
            "route_freq": "# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read",
            "title_suffix": "",
            "charge": 0,
            "multiplicity": 1,
            "nproc": 16,
            "link1": True,
        },
        {
            "account": "chem",
            "cpus": 16,
            "mem": "32G",
            "time": "24:00:00",
            "template_sha256": slurm_template_identity(DEFAULT_TEMPLATE),
        },
    )


def _create(tmp_path, table=None):
    conformer, gaussian, slurm = _configs()
    path = tmp_path / "run_manifest.json"
    create_run_manifest(
        _table() if table is None else table,
        conformer,
        gaussian,
        slurm,
        path=str(path),
        run_id="00000000-0000-4000-8000-000000000001",
        pipeline_version="2.0.0",
        pipeline_commit="abc1234",
        rdkit_version="2025.09.3",
    )
    return path


def _create_with_xyz(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = _create(
        tmp_path,
        table=pd.DataFrame([
            {"name": "Water", "cid": 962, "IsomericSMILES": "O"}
        ]),
    )
    manifest = load_manifest(str(path))
    molecule = next(
        record for record in manifest["molecules"]
        if record["molecule_name"] == "Water"
    )
    conformer_record_id = stable_record_id(
        manifest["run_id"],
        "conformer",
        f"{molecule['molecule_identity_hash']}:0",
    )
    artifact_id = stable_record_id(
        manifest["run_id"], "xyz", conformer_record_id
    )
    xyz = tmp_path / "conformer_xyz" / "water_c00.xyz"
    xyz.parent.mkdir(parents=True, exist_ok=True)
    xyz.write_text(
        "1\nlinked starting geometry\nO 0 0 0\n",
        encoding="utf-8",
    )
    record_conformer_group(
        str(path),
        name="Water",
        cid=962,
        smiles="O",
        conformers=[{
            "conformer_id": 0,
            "method": "MMFF94",
            "n_generated": 1,
            "n_kept": 1,
            "relative_energy_kcalmol": 0.0,
            "converged": True,
            "xyz_path": str(xyz),
            "artifact_id": artifact_id,
        }],
    )
    return path, xyz, artifact_id


def _write_unvalidated(path, manifest):
    path.write_text(
        json.dumps(manifest, sort_keys=True, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )


class TestCanonicalConfigurationHash:
    def test_identical_configuration_hashes_identically(self):
        value = {"b": 2, "a": [3, 1, 2], "nested": {"z": True, "x": "O"}}
        reordered = {"nested": {"x": "O", "z": True}, "a": [2, 3, 1], "b": 2}
        assert canonical_json(value, for_config_hash=True) == canonical_json(
            reordered, for_config_hash=True
        )
        assert configuration_hash(value) == configuration_hash(reordered)

    @pytest.mark.parametrize(
        "field,value",
        [("seed", 7), ("route_opt", "# opt m062x/def2tzvp"), ("charge", -1)],
    )
    def test_scientific_change_changes_hash(self, field, value):
        base = {"seed": 42, "route_opt": "# opt b3lyp/6-31g(d)", "charge": 0}
        changed = dict(base, **{field: value})
        assert configuration_hash(base) != configuration_hash(changed)

    def test_machine_paths_timestamps_and_output_hashes_are_excluded(self):
        base = {"seed": 42, "molecules": ["O"]}
        decorated = {
            **base,
            "created_at": "2026-07-13T12:00:00Z",
            "host_path": "/Users/alice/run",
            "sha256": "f" * 64,
        }
        assert configuration_hash(base) == configuration_hash(decorated)


class TestManifestValidation:
    @pytest.mark.parametrize("count", [0, 2])
    def test_provisional_manifest_requires_exactly_one_conformer(self, tmp_path, count):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        molecule = manifest["molecules"][0]
        molecule.update({
            "provenance_status": "provisional_undefined_stereo",
            "undefined_centers": "C1",
            "pubchem_smiles": "C[C@H](O)C",
            "arbitrated_smiles": "C[C@@H](O)C",
        })
        molecule["conformers"] = molecule["conformers"] * count
        with pytest.raises(ValueError, match="must contain exactly one conformer"):
            validate_manifest(manifest)

    def test_manifest_is_immutable_and_zero_job_valid(self, tmp_path):
        path = _create(tmp_path)
        manifest = finalize_manifest(str(path))
        assert manifest["run_id"] == "00000000-0000-4000-8000-000000000001"
        assert manifest["artifacts"] == []
        with pytest.raises(FileExistsError):
            _create(tmp_path)

    @pytest.mark.parametrize(
        "route_freq,missing_token",
        [
            ("# freq b3lyp/6-31g(d) Guess=Read", "Geom=AllChk"),
            ("# freq b3lyp/6-31g(d) Geom=AllChk", "Guess=Read"),
            ("# freq b3lyp/6-31g(d)", "Geom=AllChk"),
        ],
    )
    def test_manifest_creation_requires_link1_checkpoint_reads(
        self, tmp_path, route_freq, missing_token
    ):
        conformer, gaussian, slurm = _configs()
        gaussian["route_freq"] = route_freq
        with pytest.raises(ValueError, match=missing_token):
            create_run_manifest(
                _table(),
                conformer,
                gaussian,
                slurm,
                path=str(tmp_path / "bad_manifest.json"),
                run_id="00000000-0000-4000-8000-000000000001",
                pipeline_version="2.0.0",
                pipeline_commit="abc1234",
                rdkit_version="2025.09.3",
            )
        assert not (tmp_path / "bad_manifest.json").exists()

    def test_manifest_validation_rejects_tampered_link1_route(self, tmp_path):
        manifest = load_manifest(str(_create(tmp_path)))
        manifest["configuration"]["gaussian"]["route_freq"] = (
            "# freq b3lyp/6-31g(d)"
        )
        manifest["config_hash"] = configuration_hash(manifest["configuration"])
        with pytest.raises(ValueError, match="Geom=AllChk"):
            validate_manifest(manifest)

    def test_link1_checkpoint_keywords_are_case_insensitive(self, tmp_path):
        conformer, gaussian, slurm = _configs()
        gaussian["route_freq"] = (
            "# freq b3lyp/6-31g(d) geom = allchk guess = read"
        )
        path = tmp_path / "case_manifest.json"
        create_run_manifest(
            _table(),
            conformer,
            gaussian,
            slurm,
            path=str(path),
            run_id="00000000-0000-4000-8000-000000000001",
            pipeline_version="2.0.0",
            pipeline_commit="abc1234",
            rdkit_version="2025.09.3",
        )
        assert path.exists()

    def test_duplicate_molecule_record_rejected(self, tmp_path):
        manifest = load_manifest(str(_create(tmp_path)))
        manifest["molecules"].append(deepcopy(manifest["molecules"][0]))
        with pytest.raises(ValueError, match="Duplicate molecule"):
            validate_manifest(manifest)

    def test_manifest_creation_rejects_sanitized_basename_collision_before_write(
        self, tmp_path
    ):
        conformer, gaussian, slurm = _configs()
        table = pd.DataFrame([
            {"name": "Water", "cid": 962, "IsomericSMILES": "O"},
            {"name": "water", "cid": 962, "IsomericSMILES": "O"},
        ])
        path = tmp_path / "run_manifest.json"
        with pytest.raises(ValueError, match="both map to output basename"):
            create_run_manifest(
                table,
                conformer,
                gaussian,
                slurm,
                path=str(path),
                run_id="00000000-0000-4000-8000-000000000001",
                pipeline_version="2.0.0",
                pipeline_commit="abc1234",
                rdkit_version="2025.09.3",
            )
        assert not path.exists()

    def test_manifest_validation_rejects_tampered_basename_collision(self, tmp_path):
        manifest = load_manifest(str(_create(tmp_path)))
        first_name = manifest["molecules"][0]["molecule_name"]
        manifest["molecules"][1]["molecule_name"] = first_name.swapcase()
        with pytest.raises(ValueError, match="both map to output basename"):
            validate_manifest(manifest)

    def test_documented_manifest_does_not_dirty_clean_checkout(
        self, tmp_path, monkeypatch
    ):
        repo = tmp_path / "repo"
        repo.mkdir()
        project_root = Path(__file__).resolve().parents[1]
        (repo / ".gitignore").write_text(
            (project_root / ".gitignore").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
        subprocess.run(
            ["git", "config", "user.email", "tests@example.invalid"],
            cwd=repo,
            check=True,
        )
        subprocess.run(
            ["git", "config", "user.name", "Offline Tests"],
            cwd=repo,
            check=True,
        )
        subprocess.run(["git", "add", ".gitignore"], cwd=repo, check=True)
        subprocess.run(
            ["git", "commit", "-q", "-m", "test baseline"],
            cwd=repo,
            check=True,
        )
        before_version, before_commit = pipeline_provenance(cwd=str(repo))
        assert before_commit and not before_commit.endswith(".dirty")

        monkeypatch.setattr(
            manifest_module,
            "pipeline_provenance",
            lambda: pipeline_provenance(cwd=str(repo)),
        )
        conformer, gaussian, slurm = _configs()
        create_run_manifest(
            _table(),
            conformer,
            gaussian,
            slurm,
            path=str(repo / "run_manifest.json"),
            run_id="00000000-0000-4000-8000-000000000001",
            rdkit_version="2025.09.3",
        )

        after_version, after_commit = pipeline_provenance(cwd=str(repo))
        assert after_version == before_version
        assert after_commit == before_commit
        assert not after_commit.endswith(".dirty")

    def test_duplicate_artifact_id_and_path_rejected(self, tmp_path):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        manifest["artifacts"].append(deepcopy(manifest["artifacts"][0]))
        with pytest.raises(ValueError, match="Duplicate artifact"):
            validate_manifest(manifest)

    def test_artifact_path_cannot_escape_run_package(self, tmp_path):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        manifest["artifacts"][0]["relative_path"] = "../water.xyz"
        with pytest.raises(ValueError, match="inside the run package"):
            validate_manifest(manifest)

    def test_config_hash_drift_rejected(self, tmp_path):
        manifest = load_manifest(str(_create(tmp_path)))
        manifest["configuration"]["gaussian"]["charge"] = -1
        with pytest.raises(ValueError, match="config_hash"):
            validate_manifest(manifest)


class TestArtifactLineageAndHashes:
    def test_xyz_record_has_complete_lineage_and_verified_hash(self, tmp_path):
        path = _create(tmp_path)
        manifest = load_manifest(str(path))
        name, cid, smiles = "Water", 962, "O"
        conformer_record_id = stable_record_id(
            manifest["run_id"],
            "conformer",
            f"{molecule_identity_hash(name, cid, smiles)}:0",
        )
        artifact_id = stable_record_id(
            manifest["run_id"], "xyz", conformer_record_id
        )
        xyz = tmp_path / "conformer_xyz" / "water_c00.xyz"
        xyz.parent.mkdir()
        xyz.write_text("1\nlinked starting geometry\nO 0 0 0\n")
        recorded_id, digest = record_conformer_group(
            str(path),
            name=name,
            cid=cid,
            smiles=smiles,
            conformers=[{
                "conformer_id": 0,
                "method": "MMFF94",
                "n_generated": 1,
                "n_kept": 1,
                "relative_energy_kcalmol": 0.0,
                "converged": True,
                "xyz_path": str(xyz),
                "artifact_id": artifact_id,
            }],
        )[0]
        assert recorded_id == conformer_record_id
        artifact = verify_artifact(str(path), artifact_id)
        assert artifact["sha256"] == digest
        assert artifact["relative_path"] == "conformer_xyz/water_c00.xyz"

        xyz.write_text("1\ntampered\nO 1 0 0\n")
        with pytest.raises(ValueError, match="hash mismatch"):
            verify_artifact(str(path), artifact_id)


class TestCompleteConformerSchema:
    @pytest.mark.parametrize(
        "field,value", [("seed", 7), ("n_generate", 19), ("top_n", 2), ("rmsd_prune", 0.4)]
    )
    def test_conformer_search_knobs_must_match_authoritative_config(
        self, tmp_path, field, value
    ):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"][0][field] = value
        with pytest.raises(ValueError, match="configuration.conformer"):
            validate_manifest(manifest)

    def test_recording_rejects_nan_convergence_before_manifest_mutation(
        self, tmp_path
    ):
        path = _create(
            tmp_path,
            table=pd.DataFrame([
                {"name": "Water", "cid": 962, "IsomericSMILES": "O"}
            ]),
        )
        manifest = load_manifest(str(path))
        molecule = manifest["molecules"][0]
        conformer_record_id = stable_record_id(
            manifest["run_id"],
            "conformer",
            f"{molecule['molecule_identity_hash']}:0",
        )
        artifact_id = stable_record_id(
            manifest["run_id"], "xyz", conformer_record_id
        )
        xyz = tmp_path / "conformer_xyz" / "water_c00.xyz"
        xyz.parent.mkdir()
        xyz.write_text("1\nseed\nO 0 0 0\n", encoding="utf-8")
        before = path.read_bytes()

        with pytest.raises(ValueError, match="converged"):
            record_conformer_group(
                str(path),
                name="Water",
                cid=962,
                smiles="O",
                conformers=[{
                    "conformer_id": 0,
                    "method": "MMFF94",
                    "n_generated": 1,
                    "n_kept": 1,
                    "relative_energy_kcalmol": 0.0,
                    "converged": float("nan"),
                    "xyz_path": str(xyz),
                    "artifact_id": artifact_id,
                }],
            )
        assert path.read_bytes() == before

    @pytest.mark.parametrize(
        "field",
        [
            "conformer_record_id",
            "conformer_id",
            "method",
            "seed",
            "n_generate",
            "n_generated",
            "top_n",
            "n_kept",
            "rmsd_prune",
            "relative_energy_kcalmol",
            "converged",
            "xyz_artifact_id",
        ],
    )
    @pytest.mark.parametrize("operation", ["load", "finalize"])
    def test_missing_required_field_is_rejected_without_rewrite(
        self, tmp_path, field, operation
    ):
        package = tmp_path / f"{field}_{operation}"
        path, _, _ = _create_with_xyz(package)
        manifest = load_manifest(str(path))
        del manifest["molecules"][0]["conformers"][0][field]
        _write_unvalidated(path, manifest)
        before = path.read_bytes()

        with pytest.raises(ValueError, match="missing field"):
            if operation == "load":
                load_manifest(str(path))
            else:
                finalize_manifest(str(path))
        assert path.read_bytes() == before

    @pytest.mark.parametrize(
        "field,value,match",
        [
            ("conformer_record_id", "   ", "nonblank"),
            ("method", "", "nonblank"),
            ("xyz_artifact_id", " ", "nonblank"),
            ("conformer_id", True, "integer"),
            ("seed", 1.5, "integer"),
            ("n_generate", -1, "nonnegative"),
            ("n_generated", 21, "n_generated"),
            ("n_kept", 3, "n_kept"),
            ("top_n", -1, "nonnegative"),
            ("rmsd_prune", float("nan"), "finite"),
            ("rmsd_prune", float("inf"), "finite"),
            ("rmsd_prune", -0.1, "nonnegative"),
            ("relative_energy_kcalmol", float("nan"), "finite"),
            ("relative_energy_kcalmol", float("-inf"), "finite"),
            ("relative_energy_kcalmol", "0.0", "finite"),
            ("converged", "true", "JSON boolean"),
            ("converged", None, "converged"),
        ],
    )
    def test_invalid_conformer_values_are_rejected(
        self, tmp_path, field, value, match
    ):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"][0][field] = value
        with pytest.raises(ValueError, match=match):
            validate_manifest(manifest)


class TestExactConformerXyzLineage:
    def test_nonexistent_xyz_reference_is_rejected(self, tmp_path):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"][0]["xyz_artifact_id"] = (
            "xyz-does-not-exist"
        )
        with pytest.raises(ValueError, match="does not exist"):
            validate_manifest(manifest)

    def test_non_xyz_reference_is_rejected(self, tmp_path):
        path, _, xyz_artifact_id = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        conformer_record_id = manifest["molecules"][0]["conformers"][0][
            "conformer_record_id"
        ]
        com_id = stable_record_id(manifest["run_id"], "com", xyz_artifact_id)
        com_path = tmp_path / "gaussian_inputs" / "water.com"
        com_path.parent.mkdir()
        com_path.write_text("valid com placeholder for lineage test\n", encoding="utf-8")
        record_child_artifact(
            str(path),
            kind="com",
            artifact_id=com_id,
            parent_artifact_id=xyz_artifact_id,
            conformer_record_id=conformer_record_id,
            path=str(com_path),
        )
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"][0]["xyz_artifact_id"] = com_id
        with pytest.raises(ValueError, match="XYZ artifact"):
            validate_manifest(manifest)

    def test_two_conformers_cannot_share_one_xyz(self, tmp_path):
        path, _, first_xyz_id = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        molecule = manifest["molecules"][0]
        second_record_id = stable_record_id(
            manifest["run_id"],
            "conformer",
            f"{molecule['molecule_identity_hash']}:1",
        )
        second_xyz_id = stable_record_id(manifest["run_id"], "xyz", second_record_id)
        second_xyz = tmp_path / "conformer_xyz" / "water_c01.xyz"
        second_xyz.write_text("1\nsecond\nO 1 0 0\n", encoding="utf-8")
        first_xyz = tmp_path / "conformer_xyz" / "water_c00.xyz"
        record_conformer_group(
            str(path),
            name="Water",
            cid=962,
            smiles="O",
            conformers=[
                {
                    "conformer_id": 0,
                    "method": "MMFF94",
                    "n_generated": 2,
                    "n_kept": 2,
                    "relative_energy_kcalmol": 0.0,
                    "converged": True,
                    "xyz_path": str(first_xyz),
                    "artifact_id": first_xyz_id,
                },
                {
                    "conformer_id": 1,
                    "method": "MMFF94",
                    "n_generated": 2,
                    "n_kept": 2,
                    "relative_energy_kcalmol": 1.0,
                    "converged": True,
                    "xyz_path": str(second_xyz),
                    "artifact_id": second_xyz_id,
                },
            ],
        )
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"][1]["xyz_artifact_id"] = first_xyz_id
        with pytest.raises(ValueError, match="lineage|same XYZ"):
            validate_manifest(manifest)

    def test_orphan_xyz_is_rejected(self, tmp_path):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"] = []
        with pytest.raises(ValueError, match="conformer record|orphan"):
            validate_manifest(manifest)

    def test_second_xyz_for_one_conformer_is_rejected(self, tmp_path):
        path, _, _ = _create_with_xyz(tmp_path)
        manifest = load_manifest(str(path))
        duplicate = deepcopy(manifest["artifacts"][0])
        duplicate["artifact_id"] = "xyz-" + "f" * 24
        duplicate["relative_path"] = "conformer_xyz/duplicate.xyz"
        manifest["artifacts"].append(duplicate)
        with pytest.raises(ValueError, match="stable|multiple XYZ"):
            validate_manifest(manifest)


class TestResolvedArtifactContainment:
    def test_parent_directory_symlink_escape_fails_atomically(self, tmp_path):
        package = tmp_path / "run"
        path, xyz, artifact_id = _create_with_xyz(package)
        external_dir = tmp_path / "external_xyz"
        shutil.move(str(xyz.parent), external_dir)
        xyz.parent.symlink_to(external_dir, target_is_directory=True)
        external_file = external_dir / xyz.name
        external_before = external_file.read_bytes()
        manifest_before = path.read_bytes()

        for operation in (
            lambda: artifact_abspath(str(path), f"conformer_xyz/{xyz.name}"),
            lambda: verify_artifact(str(path), artifact_id),
            lambda: finalize_manifest(str(path)),
        ):
            with pytest.raises(ValueError, match="inside the run package"):
                operation()
        assert path.read_bytes() == manifest_before
        assert external_file.read_bytes() == external_before

    def test_direct_file_symlink_escape_is_rejected(self, tmp_path):
        package = tmp_path / "run"
        path, xyz, artifact_id = _create_with_xyz(package)
        external = tmp_path / "external.xyz"
        external.write_bytes(xyz.read_bytes())
        xyz.unlink()
        xyz.symlink_to(external)
        with pytest.raises(ValueError, match="inside the run package"):
            verify_artifact(str(path), artifact_id)

    def test_internal_symlink_remains_supported(self, tmp_path):
        package = tmp_path / "run"
        path, xyz, artifact_id = _create_with_xyz(package)
        internal = package / "stored" / xyz.name
        internal.parent.mkdir()
        shutil.move(str(xyz), internal)
        xyz.symlink_to(internal)
        assert artifact_abspath(str(path), f"conformer_xyz/{xyz.name}") == str(internal)
        assert verify_artifact(str(path), artifact_id)["artifact_id"] == artifact_id


class TestAtomicConformerGroupPublication:
    @staticmethod
    def _manifest(tmp_path):
        return _create(
            tmp_path,
            table=pd.DataFrame([
                {"name": "Water", "cid": 962, "IsomericSMILES": "O"}
            ]),
        )

    @staticmethod
    def _payload(path, tmp_path, *, count=3, prefix="new", staged=True):
        manifest = load_manifest(str(path))
        molecule = manifest["molecules"][0]
        final_dir = tmp_path / "conformer_xyz"
        final_dir.mkdir(parents=True, exist_ok=True)
        staging_dir = tmp_path / f".staging-{prefix}"
        if staged:
            staging_dir.mkdir(parents=True, exist_ok=True)
        payload = []
        for conformer_id in range(count):
            final_path = final_dir / f"water_c{conformer_id:02d}.xyz"
            source_path = (
                staging_dir / final_path.name if staged else final_path
            )
            source_path.write_text(
                f"1\n{prefix}-{conformer_id}\nO {conformer_id} 0 0\n",
                encoding="utf-8",
            )
            conformer_record_id = stable_record_id(
                manifest["run_id"],
                "conformer",
                f"{molecule['molecule_identity_hash']}:{conformer_id}",
            )
            artifact_id = stable_record_id(
                manifest["run_id"], "xyz", conformer_record_id
            )
            item = {
                "conformer_id": conformer_id,
                "method": "MMFF94",
                "n_generated": count,
                "n_kept": count,
                "relative_energy_kcalmol": float(conformer_id),
                "converged": True,
                "xyz_path": str(final_path),
                "artifact_id": artifact_id,
            }
            if staged:
                item["staged_xyz_path"] = str(source_path)
            payload.append(item)
        return payload

    def _publish(self, path, payload):
        return record_conformer_group(
            str(path),
            name="Water",
            cid=962,
            smiles="O",
            conformers=payload,
        )

    def test_three_conformers_publish_with_one_manifest_write(
        self, tmp_path, monkeypatch
    ):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=3)
        calls = []
        real_write = manifest_module.write_manifest

        def counted_write(manifest_path, manifest):
            calls.append(deepcopy(manifest))
            return real_write(manifest_path, manifest)

        monkeypatch.setattr(manifest_module, "write_manifest", counted_write)
        recorded = self._publish(path, payload)

        assert len(calls) == 1
        assert len(recorded) == 3
        final = load_manifest(str(path))
        assert len(final["molecules"][0]["conformers"]) == 3
        assert len([a for a in final["artifacts"] if a["kind"] == "xyz"]) == 3
        assert not list(tmp_path.glob(".staging-*"))

    def test_foreign_final_xyz_fails_before_publication(self, tmp_path):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=1)
        final_path = Path(payload[0]["xyz_path"])
        final_path.write_bytes(b"foreign sentinel\n")
        manifest_before = path.read_bytes()

        with pytest.raises(FileExistsError, match="not a tracked prior artifact"):
            self._publish(path, payload)

        assert final_path.read_bytes() == b"foreign sentinel\n"
        assert path.read_bytes() == manifest_before
        assert not list(tmp_path.rglob("*.backup-*"))

    @pytest.mark.parametrize(
        "ids",
        [[0, 2], [1, 2], [0, 1, 3], [0, 0]],
    )
    def test_invalid_complete_id_sets_fail_before_mutation(self, tmp_path, ids):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=len(ids))
        for item, conformer_id in zip(payload, ids):
            item["conformer_id"] = conformer_id
            manifest = load_manifest(str(path))
            molecule = manifest["molecules"][0]
            record_id = stable_record_id(
                manifest["run_id"],
                "conformer",
                f"{molecule['molecule_identity_hash']}:{conformer_id}",
            )
            item["artifact_id"] = stable_record_id(
                manifest["run_id"], "xyz", record_id
            )
        before = path.read_bytes()
        with pytest.raises(ValueError, match="contiguous|Duplicate"):
            self._publish(path, payload)
        assert path.read_bytes() == before
        assert not list(tmp_path.glob(".staging-*"))

    def test_incomplete_payload_fails_before_manifest_mutation(self, tmp_path):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=1)
        payload[0]["n_generated"] = 3
        payload[0]["n_kept"] = 3
        before = path.read_bytes()
        with pytest.raises(ValueError, match="incomplete"):
            self._publish(path, payload)
        assert path.read_bytes() == before

    def test_duplicate_staged_sources_fail_before_manifest_mutation(self, tmp_path):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=2)
        Path(payload[1]["staged_xyz_path"]).unlink()
        payload[1]["staged_xyz_path"] = payload[0]["staged_xyz_path"]
        before = path.read_bytes()
        with pytest.raises(ValueError, match="Duplicate staged"):
            self._publish(path, payload)
        assert path.read_bytes() == before
        assert not list(tmp_path.glob(".staging-*"))

    def test_staged_source_cannot_collide_with_another_final_path(self, tmp_path):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=2)
        colliding_path = Path(payload[1]["xyz_path"])
        colliding_path.write_text("1\nprotected\nO 0 0 0\n", encoding="utf-8")
        payload[0]["staged_xyz_path"] = str(colliding_path)
        before = path.read_bytes()
        protected = colliding_path.read_bytes()
        with pytest.raises(ValueError, match="collides"):
            self._publish(path, payload)
        assert path.read_bytes() == before
        assert colliding_path.read_bytes() == protected

    @pytest.mark.parametrize(
        "field,value",
        [
            ("method", "UFF"),
            ("seed", 7),
            ("n_generate", 19),
            ("n_generated", 2),
            ("top_n", 2),
            ("n_kept", 2),
            ("rmsd_prune", 0.4),
        ],
    )
    def test_inconsistent_group_metadata_is_rejected(
        self, tmp_path, field, value
    ):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=3, staged=False)
        self._publish(path, payload)
        manifest = load_manifest(str(path))
        manifest["molecules"][0]["conformers"][1][field] = value
        with pytest.raises(
            ValueError, match="inconsistent|configuration|n_kept|n_generated"
        ):
            validate_manifest(manifest)

    @pytest.mark.parametrize(
        "ids,flags,n_kept,accepted",
        [
            ([0, 1, 2], [True, True, True], 3, True),
            ([0], [False], 1, True),
            ([0, 1], [True, False], 2, False),
            ([0, 1], [False, False], 2, False),
            ([0], [False], 2, False),
            ([1], [False], 1, False),
        ],
    )
    def test_complete_group_convergence_policy(
        self, tmp_path, ids, flags, n_kept, accepted
    ):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=len(ids), staged=False)
        manifest = load_manifest(str(path))
        molecule = manifest["molecules"][0]
        for item, conformer_id, converged in zip(payload, ids, flags):
            item["conformer_id"] = conformer_id
            item["n_kept"] = n_kept
            item["n_generated"] = max(n_kept, len(ids))
            item["converged"] = converged
            record_id = stable_record_id(
                manifest["run_id"],
                "conformer",
                f"{molecule['molecule_identity_hash']}:{conformer_id}",
            )
            item["artifact_id"] = stable_record_id(
                manifest["run_id"], "xyz", record_id
            )
        if accepted:
            self._publish(path, payload)
            assert len(load_manifest(str(path))["molecules"][0]["conformers"]) == n_kept
        else:
            with pytest.raises(ValueError, match="convergence|incomplete|contiguous"):
                self._publish(path, payload)

    def test_move_failure_restores_old_group_and_files(self, tmp_path, monkeypatch):
        path = self._manifest(tmp_path)
        old_payload = self._payload(
            path, tmp_path, count=2, prefix="old", staged=False
        )
        self._publish(path, old_payload)
        manifest_before = path.read_bytes()
        old_bytes = {
            item["xyz_path"]: Path(item["xyz_path"]).read_bytes()
            for item in old_payload
        }
        new_payload = self._payload(path, tmp_path, count=2, prefix="new")
        staged_sources = {item["staged_xyz_path"] for item in new_payload}
        real_replace = os.replace
        moves = 0

        def fail_second_staged_move(source, destination):
            nonlocal moves
            if str(source) in staged_sources:
                moves += 1
                if moves == 2:
                    raise OSError("simulated placement failure")
            return real_replace(source, destination)

        monkeypatch.setattr(manifest_module.os, "replace", fail_second_staged_move)
        with pytest.raises(OSError, match="placement failure"):
            self._publish(path, new_payload)
        assert path.read_bytes() == manifest_before
        for final_path, expected in old_bytes.items():
            assert Path(final_path).read_bytes() == expected
        assert not list(tmp_path.glob(".staging-*"))

    def test_manifest_write_failure_restores_old_files(self, tmp_path, monkeypatch):
        path = self._manifest(tmp_path)
        old_payload = self._payload(
            path, tmp_path, count=2, prefix="old", staged=False
        )
        self._publish(path, old_payload)
        manifest_before = path.read_bytes()
        old_bytes = {
            item["xyz_path"]: Path(item["xyz_path"]).read_bytes()
            for item in old_payload
        }
        new_payload = self._payload(path, tmp_path, count=2, prefix="new")
        monkeypatch.setattr(
            manifest_module,
            "write_manifest",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                OSError("simulated manifest failure")
            ),
        )
        with pytest.raises(OSError, match="manifest failure"):
            self._publish(path, new_payload)
        assert path.read_bytes() == manifest_before
        for final_path, expected in old_bytes.items():
            assert Path(final_path).read_bytes() == expected
        assert not list(tmp_path.glob(".staging-*"))

    def test_successful_replacement_removes_old_descendant_lineage(self, tmp_path):
        path = self._manifest(tmp_path)
        old_payload = self._payload(
            path, tmp_path, count=2, prefix="old", staged=False
        )
        recorded = self._publish(path, old_payload)
        manifest = load_manifest(str(path))
        for item, (conformer_record_id, _digest) in zip(old_payload, recorded):
            com_path = tmp_path / "gaussian_inputs" / f"{item['conformer_id']}.com"
            com_path.parent.mkdir(exist_ok=True)
            com_path.write_text("old com\n", encoding="utf-8")
            com_id = stable_record_id(
                manifest["run_id"], "com", item["artifact_id"]
            )
            record_child_artifact(
                str(path),
                kind="com",
                artifact_id=com_id,
                parent_artifact_id=item["artifact_id"],
                conformer_record_id=conformer_record_id,
                path=str(com_path),
            )
            sh_path = tmp_path / "slurm_scripts" / f"{item['conformer_id']}.sh"
            sh_path.parent.mkdir(exist_ok=True)
            sh_path.write_text("#!/bin/bash\n", encoding="utf-8")
            sh_id = stable_record_id(manifest["run_id"], "sh", com_id)
            record_child_artifact(
                str(path),
                kind="sh",
                artifact_id=sh_id,
                parent_artifact_id=com_id,
                conformer_record_id=conformer_record_id,
                path=str(sh_path),
            )

        new_payload = self._payload(path, tmp_path, count=3, prefix="new")
        self._publish(path, new_payload)
        final = load_manifest(str(path))
        assert len(final["molecules"][0]["conformers"]) == 3
        assert {artifact["kind"] for artifact in final["artifacts"]} == {"xyz"}
        assert len(final["artifacts"]) == 3
        for item in new_payload:
            assert "new" in Path(item["xyz_path"]).read_text(encoding="utf-8")

    def test_successful_smaller_group_removes_obsolete_xyz_files(self, tmp_path):
        path = self._manifest(tmp_path)
        old_payload = self._payload(
            path, tmp_path, count=3, prefix="old", staged=False
        )
        self._publish(path, old_payload)
        new_payload = self._payload(path, tmp_path, count=1, prefix="new")
        self._publish(path, new_payload)

        assert Path(old_payload[0]["xyz_path"]).exists()
        assert not Path(old_payload[1]["xyz_path"]).exists()
        assert not Path(old_payload[2]["xyz_path"]).exists()
        final = load_manifest(str(path))
        assert len(final["molecules"][0]["conformers"]) == 1
        assert len(final["artifacts"]) == 1

    def test_internal_destination_symlink_updates_resolved_artifact(
        self, tmp_path
    ):
        path = self._manifest(tmp_path)
        payload = self._payload(path, tmp_path, count=1, prefix="old", staged=False)
        final_path = Path(payload[0]["xyz_path"])
        target = tmp_path / "stored" / final_path.name
        target.parent.mkdir()
        final_path.replace(target)
        final_path.symlink_to(target)
        self._publish(path, payload)

        replacement = self._payload(path, tmp_path, count=1, prefix="new")
        self._publish(path, replacement)

        assert final_path.is_symlink()
        assert "new-0" in target.read_text(encoding="utf-8")
        manifest = load_manifest(str(path))
        artifact_id = manifest["molecules"][0]["conformers"][0][
            "xyz_artifact_id"
        ]
        assert verify_artifact(str(path), artifact_id)["relative_path"] == (
            "stored/water_c00.xyz"
        )
