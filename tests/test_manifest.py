"""Offline tests for the authoritative v2 run manifest contract."""

from __future__ import annotations

from copy import deepcopy

import pandas as pd
import pytest

from pipeline.manifest import (
    canonical_json,
    configuration_hash,
    create_run_manifest,
    finalize_manifest,
    load_manifest,
    molecule_identity_hash,
    record_conformer_xyz,
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

    def test_duplicate_artifact_id_and_path_rejected(self, tmp_path):
        path = _create(tmp_path)
        manifest = load_manifest(str(path))
        molecule = manifest["molecules"][0]
        conformer_record_id = stable_record_id(
            manifest["run_id"],
            "conformer",
            f"{molecule['molecule_identity_hash']}:0",
        )
        molecule["conformers"].append({
            "conformer_record_id": conformer_record_id,
            "conformer_id": 0,
        })
        artifact = {
            "artifact_id": "xyz-duplicate",
            "kind": "xyz",
            "conformer_record_id": conformer_record_id,
            "relative_path": "xyz/water.xyz",
            "sha256": "a" * 64,
        }
        manifest["artifacts"] = [artifact, deepcopy(artifact)]
        with pytest.raises(ValueError, match="Duplicate artifact"):
            validate_manifest(manifest)

    def test_artifact_path_cannot_escape_run_package(self, tmp_path):
        manifest = load_manifest(str(_create(tmp_path)))
        molecule = manifest["molecules"][0]
        conformer_record_id = stable_record_id(
            manifest["run_id"],
            "conformer",
            f"{molecule['molecule_identity_hash']}:0",
        )
        molecule["conformers"].append({
            "conformer_record_id": conformer_record_id,
            "conformer_id": 0,
        })
        manifest["artifacts"].append({
            "artifact_id": stable_record_id(
                manifest["run_id"], "xyz", conformer_record_id
            ),
            "kind": "xyz",
            "conformer_record_id": conformer_record_id,
            "relative_path": "../water.xyz",
            "sha256": "a" * 64,
        })
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
        recorded_id, digest = record_conformer_xyz(
            str(path),
            name=name,
            cid=cid,
            smiles=smiles,
            conformer_id=0,
            method="MMFF94",
            n_generated=1,
            n_kept=1,
            relative_energy_kcalmol=0.0,
            converged=True,
            xyz_path=str(xyz),
            artifact_id=artifact_id,
        )
        assert recorded_id == conformer_record_id
        artifact = verify_artifact(str(path), artifact_id)
        assert artifact["sha256"] == digest
        assert artifact["relative_path"] == "conformer_xyz/water_c00.xyz"

        xyz.write_text("1\ntampered\nO 1 0 0\n")
        with pytest.raises(ValueError, match="hash mismatch"):
            verify_artifact(str(path), artifact_id)
