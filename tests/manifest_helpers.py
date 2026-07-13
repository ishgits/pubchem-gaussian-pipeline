"""Shared builders for strict v2 manifest-linked offline tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from pipeline.conformers import METHOD_POLICY
from pipeline.manifest import create_run_manifest, slurm_template_identity
from pipeline.manifest import (
    load_manifest,
    molecule_identity_hash,
    record_conformer_xyz,
    record_child_artifact,
    stable_record_id,
)
from pipeline.slurm import DEFAULT_TEMPLATE
from pipeline.utils import pipeline_provenance


ROUTE_OPT = "# opt b3lyp/6-31g(d)"
ROUTE_FREQ = "# freq b3lyp/6-31g(d) Geom=AllChk Guess=Read"


def ensure_manifest(
    tmp_path,
    molecule_table,
    *,
    seed=42,
    n_generate=20,
    top_n=3,
    rmsd_prune=0.5,
    route_opt=ROUTE_OPT,
    route_freq=ROUTE_FREQ,
    title_suffix="",
    charge=0,
    multiplicity=1,
    nproc=16,
    link1=True,
    account="myaccount",
    cpus=16,
    mem="32G",
    time="24:00:00",
    pipeline_version=None,
    pipeline_commit=None,
    rdkit_version=None,
):
    if pipeline_version is None or pipeline_commit is None:
        detected_version, detected_commit = pipeline_provenance()
        pipeline_version = pipeline_version or detected_version
        pipeline_commit = detected_commit if pipeline_commit is None else pipeline_commit
    if rdkit_version is None:
        import rdkit

        rdkit_version = rdkit.__version__

    conformer = {
        "method_policy": METHOD_POLICY,
        "seed": seed,
        "n_generate": n_generate,
        "top_n": top_n,
        "rmsd_prune": rmsd_prune,
    }
    gaussian = {
        "route_opt": route_opt,
        "route_freq": route_freq,
        "title_suffix": title_suffix,
        "charge": charge,
        "multiplicity": multiplicity,
        "nproc": nproc,
        "link1": link1,
    }
    slurm = {
        "account": account,
        "cpus": cpus,
        "mem": mem,
        "time": time,
        "template_sha256": slurm_template_identity(DEFAULT_TEMPLATE),
    }
    identity = {
        "rows": molecule_table[["name", "cid", "IsomericSMILES"]].to_dict("records"),
        "conformer": conformer,
        "gaussian": gaussian,
        "slurm": slurm,
        "pipeline_version": pipeline_version,
        "pipeline_commit": pipeline_commit,
        "rdkit_version": rdkit_version,
    }
    digest = hashlib.sha256(
        json.dumps(identity, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()[:16]
    path = tmp_path / f"run_manifest_{digest}.json"
    if not path.exists():
        create_run_manifest(
            molecule_table,
            conformer,
            gaussian,
            slurm,
            path=str(path),
            pipeline_version=pipeline_version,
            pipeline_commit=pipeline_commit,
            rdkit_version=rdkit_version,
        )
    return str(path)


def write_linked_conformer_log(
    tmp_path,
    rows,
    sample_xyz,
    *,
    route_opt=ROUTE_OPT,
    route_freq=ROUTE_FREQ,
    title_suffix="",
    charge=0,
    multiplicity=1,
    pipeline_version="2.0.0",
    pipeline_commit="abc1234",
    rdkit_version="2025.09.3",
):
    """Create an internally consistent manifest/XYZ/conformer-log test package."""
    import pandas as pd

    molecule_rows = []
    seen = set()
    for row in rows:
        name = row.get("name", "Ribose")
        if name not in seen:
            seen.add(name)
            molecule_rows.append({
                "name": name,
                "cid": row.get("cid", 5779),
                "IsomericSMILES": row.get("smiles", "C[C@H](O)C"),
            })
    table = pd.DataFrame(molecule_rows, columns=["name", "cid", "IsomericSMILES"])
    manifest_path = ensure_manifest(
        tmp_path,
        table,
        route_opt=route_opt,
        route_freq=route_freq,
        charge=charge,
        multiplicity=multiplicity,
        title_suffix=title_suffix,
        pipeline_version=pipeline_version,
        pipeline_commit=pipeline_commit,
        rdkit_version=rdkit_version,
    )
    manifest = load_manifest(manifest_path)
    source_text = open(sample_xyz, encoding="utf-8").read()
    output_rows = []
    grouped_counts = {}
    for row in rows:
        grouped_counts[row.get("name", "Ribose")] = grouped_counts.get(
            row.get("name", "Ribose"), 0
        ) + 1
    for row in rows:
        name = row.get("name", "Ribose")
        cid = row.get("cid", 5779)
        smiles = row.get("smiles", "C[C@H](O)C")
        conformer_id = int(row.get("conformer_id", 0))
        xyz_path = tmp_path / "xyz" / f"{name.lower()}_c{conformer_id:02d}.xyz"
        xyz_path.parent.mkdir(parents=True, exist_ok=True)
        xyz_path.write_text(source_text, encoding="utf-8")
        conformer_record_id = stable_record_id(
            manifest["run_id"],
            "conformer",
            f"{molecule_identity_hash(name, cid, smiles)}:{conformer_id}",
        )
        artifact_id = stable_record_id(
            manifest["run_id"], "xyz", conformer_record_id
        )
        recorded_id, digest = record_conformer_xyz(
            manifest_path,
            name=name,
            cid=cid,
            smiles=smiles,
            conformer_id=conformer_id,
            method=row.get("method", "MMFF94"),
            n_generated=row.get("n_generated", grouped_counts[name]),
            n_kept=row.get("n_kept", grouped_counts[name]),
            relative_energy_kcalmol=row.get("rel_energy_kcalmol", 0.0),
            converged=row.get("converged", True),
            xyz_path=str(xyz_path),
            artifact_id=artifact_id,
        )
        output = dict(row)
        output.update({
            "run_id": manifest["run_id"],
            "artifact_id": artifact_id,
            "config_hash": manifest["config_hash"],
            "name": name,
            "cid": cid,
            "smiles": smiles,
            "conformer_id": conformer_id,
            "rel_energy_kcalmol": row.get("rel_energy_kcalmol", 0.0),
            "method": row.get("method", "MMFF94"),
            "n_generated": row.get("n_generated", grouped_counts[name]),
            "n_kept": row.get("n_kept", grouped_counts[name]),
            "converged": row.get("converged", True),
            "xyz_path": str(xyz_path),
            "xyz_sha256": digest,
            "pipeline_version": row.get("pipeline_version", pipeline_version),
            "pipeline_commit": row.get("pipeline_commit", pipeline_commit),
            "rdkit_version": row.get("rdkit_version", rdkit_version),
            "conformer_record_id": recorded_id,
        })
        output_rows.append(output)
    log_path = tmp_path / "conformer_log.csv"
    pd.DataFrame(output_rows).to_csv(log_path, index=False)
    return str(log_path), manifest_path


def write_linked_com_log(tmp_path, specifications, sample_xyz):
    """Create arbitrary COM paths with valid XYZ→COM manifest lineage."""
    import pandas as pd

    conformer_rows = [
        {
            "name": spec.get("name", f"Mol{index}"),
            "cid": spec.get("cid", index + 1),
            "smiles": spec.get("smiles", "O"),
            "conformer_id": spec.get("conformer_id", 0),
        }
        for index, spec in enumerate(specifications)
    ]
    conformer_log, manifest_path = write_linked_conformer_log(
        tmp_path, conformer_rows, sample_xyz
    )
    source = pd.read_csv(conformer_log, dtype=str, keep_default_na=False)
    manifest = load_manifest(manifest_path)
    output = []
    for index, spec in enumerate(specifications):
        row = source.iloc[index]
        com_path = spec.get(
            "com_path", tmp_path / "gaussian_inputs" / f"job_{index}_F.com"
        )
        com_path = Path(com_path)
        com_path.parent.mkdir(parents=True, exist_ok=True)
        com_path.write_text(spec.get("content", "%chk=test.chk\n"), encoding="utf-8")
        artifact_id = stable_record_id(
            manifest["run_id"], "com", row["artifact_id"]
        )
        digest = record_child_artifact(
            manifest_path,
            kind="com",
            artifact_id=artifact_id,
            parent_artifact_id=row["artifact_id"],
            conformer_record_id=row["conformer_record_id"],
            path=str(com_path),
        )
        output.append({
            "run_id": manifest["run_id"],
            "artifact_id": artifact_id,
            "config_hash": manifest["config_hash"],
            "name": row["name"],
            "conformer_id": row["conformer_id"],
            "conformer_record_id": row["conformer_record_id"],
            "xyz_artifact_id": row["artifact_id"],
            "xyz_path": row["xyz_path"],
            "com_path": str(com_path),
            "com_sha256": digest,
            "pipeline_version": row["pipeline_version"],
            "pipeline_commit": row["pipeline_commit"],
            "rdkit_version": row["rdkit_version"],
        })
    log_path = tmp_path / "com_write_log.csv"
    pd.DataFrame(output).to_csv(log_path, index=False)
    return str(log_path), manifest_path


def direct_com_context(
    tmp_path,
    sample_xyz,
    *,
    route_opt=ROUTE_OPT,
    route_freq=ROUTE_FREQ,
    charge=0,
    multiplicity=1,
    conformer_id=0,
    rel_energy_kcalmol=0.0,
    converged=True,
    rdkit_version="2025.09.3",
    pipeline_commit="abc1234",
):
    """Return a linked XYZ path and required direct-v2 COM identity kwargs."""
    import pandas as pd

    log_path, manifest_path = write_linked_conformer_log(
        tmp_path,
        [{
            "name": "Ribose",
            "conformer_id": conformer_id,
            "rel_energy_kcalmol": rel_energy_kcalmol,
            "converged": converged,
        }],
        sample_xyz,
        route_opt=route_opt,
        route_freq=route_freq,
        charge=charge,
        multiplicity=multiplicity,
        pipeline_commit=pipeline_commit,
        rdkit_version=rdkit_version,
    )
    row = pd.read_csv(log_path, dtype=str, keep_default_na=False).iloc[0]
    manifest = load_manifest(manifest_path)
    artifact_id = stable_record_id(
        manifest["run_id"], "com", row["artifact_id"]
    )
    return row["xyz_path"], {
        "run_id": manifest["run_id"],
        "artifact_id": artifact_id,
        "config_hash": manifest["config_hash"],
        "manifest_path": manifest_path,
        "parent_artifact_id": row["artifact_id"],
        "conformer_record_id": row["conformer_record_id"],
    }
