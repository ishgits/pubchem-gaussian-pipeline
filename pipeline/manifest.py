"""Authoritative v2 run manifest and artifact-lineage utilities.

The manifest is the supported provenance authority for the RDKit conformer
path.  It deliberately stores only relative artifact paths and computes the
configuration hash from scientific/operational inputs, never from timestamps,
machine-specific absolute paths, run/artifact identifiers, or output hashes.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
import uuid
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import pipeline_provenance, sanitize_basename


MANIFEST_SCHEMA = "2.0"
_CONFIG_EXCLUDED_KEYS = {
    "artifact_id",
    "config_hash",
    "created_at",
    "file_hash",
    "run_id",
    "sha256",
    "timestamp",
    "updated_at",
}
_CONFORMER_CONFIG_KEYS = {
    "method_policy", "seed", "n_generate", "top_n", "rmsd_prune",
}
_GAUSSIAN_CONFIG_KEYS = {
    "route_opt", "route_freq", "title_suffix", "charge", "multiplicity",
    "nproc", "link1",
}
_SLURM_CONFIG_KEYS = {
    "account", "cpus", "mem", "time", "template_sha256",
}


def _require_link1_checkpoint_reads(route_freq: str) -> None:
    """Require the v2 Link1 frequency route to read checkpoint state.

    The generated Link1 frequency section deliberately omits title,
    charge/multiplicity, and coordinates.  It is therefore valid only when the
    route reads both geometry and wavefunction from the optimization
    checkpoint.  Match the required assignments case-insensitively while
    allowing ordinary Gaussian whitespace.
    """
    route = str(route_freq)
    required = {
        "Geom=AllChk": r"(?<![A-Za-z0-9_])geom\s*=\s*allchk(?![A-Za-z0-9_])",
        "Guess=Read": r"(?<![A-Za-z0-9_])guess\s*=\s*read(?![A-Za-z0-9_])",
    }
    missing = [
        label for label, pattern in required.items()
        if re.search(pattern, route, flags=re.IGNORECASE) is None
    ]
    if missing:
        raise ValueError(
            "The v2 Link1 frequency route must contain checkpoint-read "
            "keyword(s): " + ", ".join(missing) + "."
        )


def require_exact_artifact_id_set(
    manifest: dict,
    kind: str,
    observed_ids,
    *,
    source_label: str,
) -> None:
    """Require a stage index to cover exactly the manifest artifacts of *kind*.

    Stage CSVs are subordinate indexes under the frozen v2 contract.  Accepting
    a truncated CSV would silently drop downstream scientific or operational
    artifacts, so both missing and extra identifiers are fatal before mutation.
    Duplicate rows remain the responsibility of the caller's row-level checks.
    """
    if kind not in {"xyz", "com", "sh"}:
        raise ValueError(f"Unsupported manifest artifact kind: {kind!r}.")
    expected = {
        str(artifact["artifact_id"])
        for artifact in manifest.get("artifacts", [])
        if artifact.get("kind") == kind
    }
    observed = {
        str(value).strip()
        for value in observed_ids
        if value is not None and str(value).strip()
    }
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing {kind} artifact_id(s) {missing}")
        if extra:
            details.append(f"extra {kind} artifact_id(s) {extra}")
        raise ValueError(
            f"{source_label} does not exactly match manifest {kind} artifacts: "
            + "; ".join(details)
        )


def sha256_file(path: str | os.PathLike[str]) -> str:
    """Return the SHA-256 digest of a regular file's bytes."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    """Return the SHA-256 digest of UTF-8 *text*."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _json_scalar(value: Any) -> Any:
    """Normalize pandas/numpy-like scalars into strict, portable JSON values."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Manifest values must not contain NaN or infinity.")
        return float(value)
    if hasattr(value, "item"):
        return _json_scalar(value.item())
    return value


def _canonical_value(value: Any, *, strip_excluded: bool) -> Any:
    if isinstance(value, dict):
        normalized = {}
        for raw_key, raw_value in value.items():
            key = str(raw_key)
            if strip_excluded:
                lowered = key.lower()
                if (
                    key in _CONFIG_EXCLUDED_KEYS
                    or "timestamp" in lowered
                    or lowered.endswith("_at")
                    or (
                        isinstance(raw_value, str)
                        and os.path.isabs(os.path.expanduser(raw_value))
                    )
                ):
                    continue
            normalized[key] = _canonical_value(
                raw_value, strip_excluded=strip_excluded
            )
        return {key: normalized[key] for key in sorted(normalized)}
    if isinstance(value, (list, tuple, set)):
        items = [
            _canonical_value(item, strip_excluded=strip_excluded)
            for item in value
            if not (
                strip_excluded
                and isinstance(item, str)
                and os.path.isabs(os.path.expanduser(item))
            )
        ]
        # Configuration list order is not scientific identity in v2: molecule
        # request order, record order, and dictionary insertion order must hash
        # identically.  Sorting by each item's own canonical JSON is total and
        # deterministic for heterogeneous JSON values.
        return sorted(
            items,
            key=lambda item: json.dumps(
                item, sort_keys=True, separators=(",", ":"), ensure_ascii=False
            ),
        )
    return _json_scalar(value)


def canonical_json(value: Any, *, for_config_hash: bool = False) -> str:
    """Serialize *value* as canonical JSON.

    When ``for_config_hash`` is true, identifiers, timestamps, and output-file
    hashes are recursively excluded and list ordering is normalized.
    """
    normalized = _canonical_value(value, strip_excluded=for_config_hash)
    return json.dumps(
        normalized,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def configuration_hash(configuration: dict) -> str:
    """Hash one complete canonical scientific/operational configuration."""
    return sha256_text(canonical_json(configuration, for_config_hash=True))


def molecule_identity_hash(molecule_name: str, cid, smiles: str) -> str:
    """Return the stable identity digest used by molecule/conformer records."""
    identity = {
        "molecule_name": str(molecule_name),
        "CID": None if cid is None or pd.isna(cid) else int(float(cid)),
        "IsomericSMILES": "" if smiles is None or pd.isna(smiles) else str(smiles),
    }
    return sha256_text(canonical_json(identity))


def stable_record_id(run_id: str, kind: str, logical_key: str) -> str:
    """Return a deterministic, collision-resistant ID within one immutable run."""
    digest = sha256_text(f"{run_id}\0{kind}\0{logical_key}")
    return f"{kind}-{digest[:24]}"


def slurm_template_identity(template: str) -> str:
    """Return the exact template-text identity stored in SLURM configuration."""
    return sha256_text(template)


def _require_exact_config(section: str, value: dict, keys: set[str]) -> dict:
    if not isinstance(value, dict):
        raise ValueError(f"{section} configuration must be a dictionary.")
    missing = sorted(keys - set(value))
    if missing:
        raise ValueError(
            f"{section} configuration is missing required field(s): "
            + ", ".join(missing)
        )
    return deepcopy(value)


def _molecule_records(molecule_table) -> list[dict]:
    if isinstance(molecule_table, str):
        molecule_table = pd.read_csv(molecule_table)
    required = {"name", "cid", "IsomericSMILES"}
    missing = sorted(required - set(molecule_table.columns))
    if missing:
        raise ValueError(
            "Molecule table is missing manifest identity column(s): "
            + ", ".join(missing)
        )

    records = []
    names = set()
    identities = set()
    output_basenames: dict[str, str] = {}
    for _, row in molecule_table.iterrows():
        name = str(row["name"])
        if name in names:
            raise ValueError(f"Duplicate molecule record/name: {name!r}.")
        names.add(name)
        basename = sanitize_basename(name)
        if basename == "":
            raise ValueError(
                f"Molecule label {name!r} sanitizes to an empty filename; give "
                "it a name with at least one alphanumeric character."
            )
        if basename in output_basenames:
            previous = output_basenames[basename]
            raise ValueError(
                f"Molecule labels {previous!r} and {name!r} both map to output "
                f"basename {basename!r}. Use unique labels that remain distinct "
                "after filename sanitization."
            )
        output_basenames[basename] = name
        cid = None if pd.isna(row["cid"]) else int(float(row["cid"]))
        smiles = "" if pd.isna(row["IsomericSMILES"]) else str(row["IsomericSMILES"])
        identity_hash = molecule_identity_hash(name, cid, smiles)
        if identity_hash in identities:
            raise ValueError(f"Duplicate molecule identity record: {name!r}.")
        identities.add(identity_hash)
        records.append({
            "molecule_name": name,
            "CID": cid,
            "IsomericSMILES": smiles,
            "molecule_identity_hash": identity_hash,
            "conformers": [],
        })
    return records


def create_run_manifest(
    molecule_table,
    conformer: dict,
    gaussian: dict,
    slurm: dict,
    path: str = "run_manifest.json",
    *,
    run_id: str | None = None,
    pipeline_version: str | None = None,
    pipeline_commit: str | None = None,
    rdkit_version: str | None = None,
) -> dict:
    """Create the immutable identity/configuration shell for one v2 execution."""
    if os.path.exists(path):
        raise FileExistsError(
            f"Refusing to replace existing immutable run manifest: {path!r}."
        )
    if pipeline_version is None or pipeline_commit is None:
        detected_version, detected_commit = pipeline_provenance()
        pipeline_version = pipeline_version or detected_version
        pipeline_commit = detected_commit if pipeline_commit is None else pipeline_commit
    if rdkit_version is None:
        import rdkit

        rdkit_version = rdkit.__version__

    try:
        run_id = str(uuid.UUID(str(run_id))) if run_id is not None else str(uuid.uuid4())
    except (TypeError, ValueError, AttributeError):
        raise ValueError("run_id must be a valid collision-resistant UUID.")
    molecules = _molecule_records(molecule_table)
    conformer_config = _require_exact_config(
        "Conformer", conformer, _CONFORMER_CONFIG_KEYS
    )
    gaussian_config = _require_exact_config(
        "Gaussian", gaussian, _GAUSSIAN_CONFIG_KEYS
    )
    for route_field in ("route_opt", "route_freq"):
        if not isinstance(gaussian_config[route_field], str) or not gaussian_config[route_field].strip():
            raise ValueError(f"Gaussian {route_field} must be a nonblank route line.")
    if gaussian_config["link1"] is not True:
        raise ValueError("The v2 manifest requires the Link1 opt→freq contract.")
    _require_link1_checkpoint_reads(gaussian_config["route_freq"])
    slurm_config = _require_exact_config("SLURM", slurm, _SLURM_CONFIG_KEYS)

    configured_molecules = [
        {key: record[key] for key in (
            "molecule_name", "CID", "IsomericSMILES", "molecule_identity_hash"
        )}
        for record in molecules
    ]
    configuration = {
        "molecules": configured_molecules,
        "pipeline_version": str(pipeline_version),
        "pipeline_commit": str(pipeline_commit or ""),
        "rdkit_version": str(rdkit_version),
        "conformer": conformer_config,
        "gaussian": gaussian_config,
        "slurm": slurm_config,
    }
    manifest = {
        "manifest_schema": MANIFEST_SCHEMA,
        "run_id": run_id,
        "config_hash": configuration_hash(configuration),
        "pipeline_version": str(pipeline_version),
        "pipeline_commit": str(pipeline_commit or ""),
        "rdkit_version": str(rdkit_version),
        "configuration": configuration,
        "molecules": molecules,
        "artifacts": [],
    }
    write_manifest(path, manifest)
    return manifest


def _duplicate_values(records: list[dict], key: str) -> set[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for record in records:
        value = str(record.get(key, ""))
        if not value:
            continue
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return duplicates


def validate_manifest(manifest: dict) -> None:
    """Reject malformed manifests, duplicate records/IDs, and config drift."""
    required = {
        "manifest_schema", "run_id", "config_hash", "pipeline_version",
        "pipeline_commit", "rdkit_version", "configuration", "molecules",
        "artifacts",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError("Manifest is missing field(s): " + ", ".join(missing))
    if manifest["manifest_schema"] != MANIFEST_SCHEMA:
        raise ValueError(
            f"Unsupported manifest_schema {manifest['manifest_schema']!r}; "
            f"expected {MANIFEST_SCHEMA!r}."
        )
    try:
        if str(uuid.UUID(str(manifest["run_id"]))) != str(manifest["run_id"]):
            raise ValueError
    except (TypeError, ValueError, AttributeError):
        raise ValueError("Manifest run_id must be a canonical UUID.")
    if not re.fullmatch(r"[0-9a-f]{64}", str(manifest["config_hash"])):
        raise ValueError("Manifest config_hash must be 64 lowercase hex characters.")
    expected_hash = configuration_hash(manifest["configuration"])
    if manifest["config_hash"] != expected_hash:
        raise ValueError("Manifest config_hash does not match its configuration.")
    for field in ("pipeline_version", "pipeline_commit", "rdkit_version"):
        if manifest[field] != manifest["configuration"].get(field):
            raise ValueError(f"Manifest {field} disagrees with configuration.{field}.")
    for field in ("pipeline_version", "rdkit_version"):
        if not str(manifest[field]).strip():
            raise ValueError(f"Manifest {field} must be nonblank.")
    configuration = manifest["configuration"]
    for section, keys in (
        ("conformer", _CONFORMER_CONFIG_KEYS),
        ("gaussian", _GAUSSIAN_CONFIG_KEYS),
        ("slurm", _SLURM_CONFIG_KEYS),
    ):
        value = configuration.get(section)
        if not isinstance(value, dict) or not keys.issubset(value):
            raise ValueError(f"Manifest configuration.{section} is incomplete.")
    if configuration["gaussian"].get("link1") is not True:
        raise ValueError("Manifest violates the Link1 opt→freq contract.")
    for route_field in ("route_opt", "route_freq"):
        if not str(configuration["gaussian"].get(route_field, "")).strip():
            raise ValueError(f"Manifest Gaussian {route_field} must be nonblank.")
    _require_link1_checkpoint_reads(configuration["gaussian"]["route_freq"])

    molecules = manifest["molecules"]
    if not isinstance(molecules, list):
        raise ValueError("Manifest molecules must be a list.")
    for key in ("molecule_name", "molecule_identity_hash"):
        duplicates = _duplicate_values(molecules, key)
        if duplicates:
            raise ValueError(f"Duplicate molecule {key} record(s): {sorted(duplicates)}")
    output_basenames: dict[str, str] = {}
    for molecule in molecules:
        name = str(molecule.get("molecule_name", ""))
        basename = sanitize_basename(name)
        if basename == "":
            raise ValueError(
                f"Molecule label {name!r} sanitizes to an empty filename; give "
                "it a name with at least one alphanumeric character."
            )
        if basename in output_basenames:
            previous = output_basenames[basename]
            raise ValueError(
                f"Molecule labels {previous!r} and {name!r} both map to output "
                f"basename {basename!r}. Use unique labels that remain distinct "
                "after filename sanitization."
            )
        output_basenames[basename] = name

    configured_molecules = configuration.get("molecules")
    expected_configured_molecules = [
        {key: molecule[key] for key in (
            "molecule_name", "CID", "IsomericSMILES", "molecule_identity_hash"
        )}
        for molecule in molecules
    ]
    if canonical_json(configured_molecules) != canonical_json(expected_configured_molecules):
        raise ValueError("Manifest molecule records disagree with configuration.molecules.")

    conformers = []
    for molecule in molecules:
        for key in (
            "molecule_name", "CID", "IsomericSMILES", "molecule_identity_hash",
            "conformers",
        ):
            if key not in molecule:
                raise ValueError(f"Manifest molecule record is missing {key!r}.")
        expected_identity = molecule_identity_hash(
            molecule["molecule_name"], molecule["CID"], molecule["IsomericSMILES"]
        )
        if molecule["molecule_identity_hash"] != expected_identity:
            raise ValueError("Manifest molecule identity hash mismatch.")
        molecule_conformers = molecule["conformers"]
        conformer_ids = [record.get("conformer_id") for record in molecule_conformers]
        if len(conformer_ids) != len(set(conformer_ids)):
            raise ValueError(
                f"Duplicate conformer record for molecule {molecule['molecule_name']!r}."
            )
        for record in molecule_conformers:
            expected_record_id = stable_record_id(
                manifest["run_id"],
                "conformer",
                f"{molecule['molecule_identity_hash']}:{int(record['conformer_id'])}",
            )
            if record.get("conformer_record_id") != expected_record_id:
                raise ValueError("Manifest conformer record ID is not stable for its lineage.")
        conformers.extend(molecule_conformers)
    duplicates = _duplicate_values(conformers, "conformer_record_id")
    if duplicates:
        raise ValueError(f"Duplicate conformer record ID(s): {sorted(duplicates)}")

    artifacts = manifest["artifacts"]
    if not isinstance(artifacts, list):
        raise ValueError("Manifest artifacts must be a list.")
    for key in ("artifact_id", "relative_path"):
        duplicates = _duplicate_values(artifacts, key)
        if duplicates:
            raise ValueError(f"Duplicate artifact {key}(s): {sorted(duplicates)}")
    artifact_ids = {artifact.get("artifact_id") for artifact in artifacts}
    artifacts_by_id = {
        artifact.get("artifact_id"): artifact for artifact in artifacts
    }
    conformer_ids = {record.get("conformer_record_id") for record in conformers}
    for artifact in artifacts:
        for key in ("artifact_id", "kind", "relative_path", "sha256"):
            if not str(artifact.get(key, "")).strip():
                raise ValueError(f"Manifest artifact is missing nonblank {key!r}.")
        path = str(artifact["relative_path"])
        if os.path.isabs(path) or ".." in Path(path).parts:
            raise ValueError("Manifest artifact paths must stay inside the run package.")
        if not re.fullmatch(r"[0-9a-f]{64}", str(artifact["sha256"])):
            raise ValueError("Manifest artifact SHA-256 must be 64 lowercase hex characters.")
        parent = artifact.get("parent_artifact_id")
        if parent and parent not in artifact_ids:
            raise ValueError(f"Artifact parent does not exist: {parent!r}.")
        conformer_record_id = artifact.get("conformer_record_id")
        if conformer_record_id and conformer_record_id not in conformer_ids:
            raise ValueError(
                f"Artifact conformer record does not exist: {conformer_record_id!r}."
            )
        kind = artifact["kind"]
        if kind not in {"xyz", "com", "sh"}:
            raise ValueError(f"Unsupported manifest artifact kind: {kind!r}.")
        if kind == "xyz" and parent:
            raise ValueError("XYZ artifacts must not have a parent artifact.")
        expected_artifact_id = stable_record_id(
            manifest["run_id"],
            kind,
            conformer_record_id if kind == "xyz" else parent,
        )
        if artifact["artifact_id"] != expected_artifact_id:
            raise ValueError("Manifest artifact ID is not stable for its lineage.")
        if kind in {"com", "sh"}:
            expected_parent_kind = "xyz" if kind == "com" else "com"
            parent_record = artifacts_by_id.get(parent)
            if parent_record is None or parent_record.get("kind") != expected_parent_kind:
                raise ValueError(
                    f"{kind} artifact requires a {expected_parent_kind} parent."
                )
            if parent_record.get("conformer_record_id") != conformer_record_id:
                raise ValueError("Artifact parent and child disagree on conformer lineage.")


def write_manifest(path: str, manifest: dict) -> None:
    """Validate and atomically write canonical, human-readable manifest JSON."""
    validate_manifest(manifest)
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".run_manifest.", dir=parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(
                _canonical_value(manifest, strip_excluded=False),
                handle,
                sort_keys=True,
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            )
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def load_manifest(path: str = "run_manifest.json") -> dict:
    """Load and validate a run manifest."""
    with open(path, encoding="utf-8") as handle:
        manifest = json.load(handle)
    validate_manifest(manifest)
    return manifest


def relative_artifact_path(path: str, manifest_path: str) -> str:
    """Return a portable POSIX path relative to the manifest package root."""
    # M-01: resolve symlinks before the containment check. os.path.abspath only
    # normalizes text, so a symlink inside the package that targets an external
    # directory would pass commonpath while the bytes land outside the package.
    # realpath here must stay symmetric with artifact_abspath() below.
    package_root = os.path.realpath(os.path.dirname(os.path.abspath(manifest_path)))
    absolute_path = os.path.realpath(path)
    if os.path.commonpath((package_root, absolute_path)) != package_root:
        raise ValueError("Artifact path must stay inside the run package.")
    relative = os.path.relpath(absolute_path, package_root)
    return Path(relative).as_posix()


def artifact_abspath(manifest_path: str, relative_path: str) -> str:
    """Resolve a stored relative artifact path against the package root."""
    if os.path.isabs(relative_path):
        raise ValueError("Manifest artifact paths must be relative.")
    # M-01: resolve symlinks so this stays symmetric with relative_artifact_path();
    # otherwise a stored path (relative to the resolved root) would round-trip
    # against an unresolved root under a symlinked package.
    package_root = os.path.realpath(os.path.dirname(os.path.abspath(manifest_path)))
    return os.path.realpath(os.path.join(package_root, relative_path))


def find_artifact(manifest: dict, artifact_id: str) -> dict:
    matches = [
        artifact for artifact in manifest["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one manifest artifact {artifact_id!r}; found {len(matches)}."
        )
    return matches[0]


def find_conformer_record(manifest: dict, conformer_record_id: str) -> tuple[dict, dict]:
    """Return the unique ``(molecule, conformer)`` lineage record for an ID."""
    matches = []
    for molecule in manifest["molecules"]:
        for record in molecule["conformers"]:
            if record["conformer_record_id"] == conformer_record_id:
                matches.append((molecule, record))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one conformer record {conformer_record_id!r}; "
            f"found {len(matches)}."
        )
    return matches[0]


def verify_artifact(manifest_path: str, artifact_id: str) -> dict:
    """Verify one artifact's regular-file status, nonzero size, and SHA-256."""
    manifest = load_manifest(manifest_path)
    artifact = find_artifact(manifest, artifact_id)
    path = artifact_abspath(manifest_path, artifact["relative_path"])
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError(f"Manifest artifact is missing, irregular, or empty: {path!r}.")
    actual = sha256_file(path)
    if actual != artifact["sha256"]:
        raise ValueError(
            f"Manifest artifact hash mismatch for {artifact_id!r}: "
            f"expected {artifact['sha256']}, got {actual}."
        )
    return artifact


def planned_artifact_id(
    manifest_path: str, kind: str, logical_key: str
) -> str:
    manifest = load_manifest(manifest_path)
    return stable_record_id(manifest["run_id"], kind, logical_key)


def _find_molecule(manifest: dict, name: str, cid, smiles: str) -> dict:
    identity = molecule_identity_hash(name, cid, smiles)
    matches = [
        molecule for molecule in manifest["molecules"]
        if molecule["molecule_identity_hash"] == identity
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one manifest molecule for {name!r}; found {len(matches)}."
        )
    return matches[0]


def assert_stage_configuration(manifest_path: str, section: str, config: dict) -> dict:
    """Require the runtime stage configuration to match the frozen manifest."""
    manifest = load_manifest(manifest_path)
    expected = manifest["configuration"].get(section)
    if canonical_json(expected) != canonical_json(config):
        raise ValueError(
            f"Runtime {section} configuration disagrees with run manifest."
        )
    return manifest


def record_conformer_xyz(
    manifest_path: str,
    *,
    name: str,
    cid,
    smiles: str,
    conformer_id: int,
    method: str,
    n_generated: int,
    n_kept: int,
    relative_energy_kcalmol: float,
    converged: bool,
    xyz_path: str,
    artifact_id: str,
) -> tuple[str, str]:
    """Record one complete conformer row and its already-written XYZ artifact."""
    manifest = load_manifest(manifest_path)
    molecule = _find_molecule(manifest, name, cid, smiles)
    config = manifest["configuration"]["conformer"]
    conformer_record_id = stable_record_id(
        manifest["run_id"],
        "conformer",
        f"{molecule['molecule_identity_hash']}:{int(conformer_id)}",
    )
    if any(
        record.get("conformer_record_id") == conformer_record_id
        for record in molecule["conformers"]
    ):
        raise ValueError(f"Duplicate conformer record: {conformer_record_id!r}.")
    expected_artifact_id = stable_record_id(
        manifest["run_id"], "xyz", conformer_record_id
    )
    if artifact_id != expected_artifact_id:
        raise ValueError("XYZ artifact ID does not match its stable lineage key.")
    relative_path = relative_artifact_path(xyz_path, manifest_path)
    digest = sha256_file(xyz_path)
    record = {
        "conformer_record_id": conformer_record_id,
        "conformer_id": int(conformer_id),
        "method": str(method),
        "seed": int(config["seed"]),
        "n_generate": int(config["n_generate"]),
        "n_generated": int(n_generated),
        "top_n": int(config["top_n"]),
        "n_kept": int(n_kept),
        "rmsd_prune": float(config["rmsd_prune"]),
        "relative_energy_kcalmol": float(relative_energy_kcalmol),
        "converged": bool(converged),
        "xyz_artifact_id": artifact_id,
    }
    molecule["conformers"].append(record)
    manifest["artifacts"].append({
        "artifact_id": artifact_id,
        "kind": "xyz",
        "conformer_record_id": conformer_record_id,
        "relative_path": relative_path,
        "sha256": digest,
    })
    write_manifest(manifest_path, manifest)
    return conformer_record_id, digest


def record_child_artifact(
    manifest_path: str,
    *,
    kind: str,
    artifact_id: str,
    parent_artifact_id: str,
    conformer_record_id: str,
    path: str,
) -> str:
    """Record a COM or SH artifact after validating its exact parent lineage."""
    if kind not in {"com", "sh"}:
        raise ValueError(f"Unsupported child artifact kind: {kind!r}.")
    manifest = load_manifest(manifest_path)
    parent = find_artifact(manifest, parent_artifact_id)
    expected_parent_kind = "xyz" if kind == "com" else "com"
    if parent["kind"] != expected_parent_kind:
        raise ValueError(
            f"{kind} artifact requires a {expected_parent_kind} parent."
        )
    if parent.get("conformer_record_id") != conformer_record_id:
        raise ValueError("Child artifact and parent disagree on conformer lineage.")
    expected_id = stable_record_id(
        manifest["run_id"], kind, parent_artifact_id
    )
    if artifact_id != expected_id:
        raise ValueError(f"{kind} artifact ID does not match its lineage key.")
    relative_path = relative_artifact_path(path, manifest_path)
    digest = sha256_file(path)
    manifest["artifacts"].append({
        "artifact_id": artifact_id,
        "kind": kind,
        "parent_artifact_id": parent_artifact_id,
        "conformer_record_id": conformer_record_id,
        "relative_path": relative_path,
        "sha256": digest,
    })
    write_manifest(manifest_path, manifest)
    return digest


def remove_artifacts_by_kind(manifest_path: str, kind: str) -> None:
    """Remove one downstream artifact layer before deterministic regeneration."""
    manifest = load_manifest(manifest_path)
    removed = {
        artifact["artifact_id"]
        for artifact in manifest["artifacts"]
        if artifact["kind"] == kind
    }
    if not removed:
        return
    if any(
        artifact.get("parent_artifact_id") in removed
        for artifact in manifest["artifacts"]
        if artifact["kind"] != kind
    ):
        raise ValueError(
            f"Cannot remove {kind} manifest records while child artifacts remain."
        )
    manifest["artifacts"] = [
        artifact for artifact in manifest["artifacts"]
        if artifact["kind"] != kind
    ]
    write_manifest(manifest_path, manifest)


def remove_conformer_lineage(manifest_path: str, molecule_names: set[str]) -> None:
    """Remove complete XYZ→COM→SH lineage for molecules being regenerated."""
    if not molecule_names:
        return
    manifest = load_manifest(manifest_path)
    conformer_record_ids = set()
    for molecule in manifest["molecules"]:
        if molecule["molecule_name"] in molecule_names:
            conformer_record_ids.update(
                record["conformer_record_id"] for record in molecule["conformers"]
            )
            molecule["conformers"] = []
    manifest["artifacts"] = [
        artifact for artifact in manifest["artifacts"]
        if artifact.get("conformer_record_id") not in conformer_record_ids
    ]
    write_manifest(manifest_path, manifest)


def finalize_manifest(manifest_path: str = "run_manifest.json") -> dict:
    """Verify every recorded artifact and rewrite the validated final manifest."""
    manifest = load_manifest(manifest_path)
    for artifact in manifest["artifacts"]:
        verify_artifact(manifest_path, artifact["artifact_id"])
    write_manifest(manifest_path, manifest)
    return manifest
