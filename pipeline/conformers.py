"""
Conformer search stage (v2).

For each molecule the pipeline generates an RDKit conformer ensemble
(ETKDGv3 embed, RMSD-pruned), MMFF94-optimizes and ranks it, and carries the
top-N lowest-energy, distinct conformers forward as DFT starting geometries,
each with recorded provenance (seed, method, relative energy).

Scientific notes (see AGENTS.md §2 and docs/architecture.md v2):

- These are force-field *starting* geometries, NOT optimized minima. The DFT
  step makes the final call among the carried candidates.
- Conformer energies here are MMFF94 (or UFF fallback) in **kcal/mol**. They are
  labeled kcal/mol everywhere and must never be mixed with the DFT Hartree
  energies produced downstream.
- Distinctness is handled by RDKit's ``pruneRmsThresh`` at embed time, so rigid
  molecules collapse to a single conformer emergently (no special-case code).

Only :func:`select_top_n` is import-safe without RDKit; the embedding functions
import RDKit lazily so the pure ranking helper can be unit-tested with no RDKit
installed (matches the repo's no-network / no-heavy-dependency test rule).
"""

from __future__ import annotations

import math
import os
import shutil
import tempfile

import pandas as pd

from .manifest import (
    assert_stage_configuration,
    find_artifact,
    find_conformer_record,
    load_manifest,
    molecule_identity_hash,
    record_conformer_group,
    relative_artifact_path,
    remove_conformer_lineage,
    sha256_file,
    stable_record_id,
)
from .utils import (
    ensure_dir,
    normalize_cid,
    parse_strict_bool,
    pipeline_provenance,
    sanitize_basename,
)

# Locked defaults (docs/implementation-plan.md v2 — confirm at approval).
N_GENERATE = 20
TOP_N = 3
RMSD_PRUNE = 0.5  # Ångström
SEED = 42
METHOD_POLICY = "MMFF94 preferred; UFF recorded fallback"

# Force-field optimization iteration budgets (M-04). The first pass uses
# FF_MAXITERS; any conformer still flagged not-converged is retried once with
# FF_MAXITERS_RETRY before being judged failed.
FF_MAXITERS = 2000
FF_MAXITERS_RETRY = 10000

# Marker written into a Gaussian .com title when a molecule's only carried
# geometry is an unconverged best-effort FF seed (M-04 decision 2b).
UNCONVERGED_FF_SEED = "UNCONVERGED_FF_SEED"


# ---------------------------------------------------------------------------
# Pure ranking helper (no RDKit) — Task 1
# ---------------------------------------------------------------------------

def select_top_n(energies_kcal, n: int) -> list[int]:
    """
    Return the indices of the *n* lowest-energy conformers, lowest first.

    Parameters
    ----------
    energies_kcal : sequence of float
        Conformer energies in kcal/mol (MMFF94 or UFF). Order corresponds to the
        conformer list they were computed from.
    n : int
        Maximum number of conformers to keep. If fewer conformers exist than
        *n*, all are returned (ranked). ``n <= 0`` returns an empty list.

    Returns
    -------
    list[int]
        Indices into *energies_kcal*, sorted by ascending energy. Ties are
        broken by original index so the result is deterministic.
    """
    if n <= 0:
        return []
    indexed = list(enumerate(energies_kcal))
    # Sort by (energy, original index) → deterministic, stable tie-break.
    indexed.sort(key=lambda pair: (pair[1], pair[0]))
    return [i for i, _ in indexed[:n]]


# ---------------------------------------------------------------------------
# Convergence handling (M-04) — pure helpers, no RDKit
# ---------------------------------------------------------------------------

def _finalize_convergence(first, retry=None):
    """
    Merge first-pass and (optional) retry FF results into per-conformer
    ``(converged: bool, energy_kcal: float)`` tuples.

    Each input is a sequence of RDKit ``(not_converged, energy)`` tuples aligned
    by conformer index, where ``not_converged == 0`` means the optimization
    converged. Where the first pass converged, its result is kept; where it did
    not and a *retry* pass is supplied, the retry result is used (M-04 step 2).
    """
    out = []
    for i, (nc1, e1) in enumerate(first):
        if nc1 == 0 or retry is None:
            out.append((nc1 == 0, float(e1)))
        else:
            nc2, e2 = retry[i]
            out.append((nc2 == 0, float(e2)))
    return out


def select_converged_top_n(energies_kcal, converged, top_n: int):
    """
    Rank and select conformers, honoring convergence (M-04 decisions 1a / 2b).

    Parameters
    ----------
    energies_kcal : sequence of float
        FF energies (kcal/mol) for every conformer, aligned with *converged*.
    converged : sequence of bool
        Per-conformer convergence flag (from :func:`_finalize_convergence`).
    top_n : int
        Maximum number of converged conformers to keep.

    Returns
    -------
    (kept_indices, all_failed)
        ``kept_indices`` — indices into *energies_kcal*, lowest energy first.
        When at least one conformer converged, only converged conformers are
        eligible (1a). When **none** converged, exactly one best-effort geometry
        — the lowest FF energy overall — is returned and ``all_failed`` is True
        (2b); its energy is unreliable and the caller must flag it.
    """
    converged_idx = [i for i, c in enumerate(converged) if c]
    if converged_idx:
        conv_energies = [energies_kcal[i] for i in converged_idx]
        keep_pos = select_top_n(conv_energies, top_n)
        return [converged_idx[p] for p in keep_pos], False

    # No conformer converged (2b): carry exactly one lowest-energy best effort.
    if not list(energies_kcal):
        return [], True
    return select_top_n(energies_kcal, 1), True


# ---------------------------------------------------------------------------
# RDKit embed + rank core — Task 2
# ---------------------------------------------------------------------------

def _optimize_confs(mol, method: str, max_iters: int):
    """
    Force-field-optimize every conformer of *mol* in place, returning RDKit's
    per-conformer ``(not_converged, energy_kcal)`` list.

    Isolated at module level (with a lazy RDKit import) so the retry logic in
    :func:`generate_conformers` can be exercised deterministically in tests.
    """
    from rdkit.Chem import AllChem

    if method == "MMFF94":
        return AllChem.MMFFOptimizeMoleculeConfs(
            mol, mmffVariant="MMFF94", maxIters=max_iters
        )
    return AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_iters)


def _optimize_single_conf(mol, method: str, conf_id: int, max_iters: int):
    """
    Re-optimize a **single** conformer of *mol* in place and return its
    ``(not_converged, energy_kcal)`` — the same shape as one row of
    :func:`_optimize_confs`.

    Used by the M-04 retry so only the conformers that failed the first pass are
    touched (M-05): already-converged conformers keep their first-pass geometry
    *and* energy, so the ranked energy always describes the coordinates that get
    written. ``CalcEnergy`` returns the MMFF/UFF energy in kcal/mol, matching the
    batch optimizer. Isolated at module level (lazy RDKit) so the retry path is
    monkeypatchable in tests.
    """
    from rdkit.Chem import AllChem

    if method == "MMFF94":
        not_converged = AllChem.MMFFOptimizeMolecule(
            mol, mmffVariant="MMFF94", maxIters=max_iters, confId=conf_id
        )
        props = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    else:
        not_converged = AllChem.UFFOptimizeMolecule(
            mol, confId=conf_id, maxIters=max_iters
        )
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    return not_converged, ff.CalcEnergy()


def _rdkit_version() -> str:
    import rdkit

    return rdkit.__version__


def check_conformer_eligibility(smiles) -> str | None:
    """
    Decide whether *smiles* may enter the conformer search, or why it is skipped.

    Returns ``None`` if the molecule is eligible, otherwise a short skip-reason
    string (logged to ``conformer_search_failed.csv``). Skipping — never silent
    auto-correction — is the safe failure mode for ambiguous chemistry.

    Skip reasons:

    - ``"no IsomericSMILES"`` — empty/missing SMILES; nothing to embed.
    - ``"unparseable SMILES"`` — RDKit cannot parse the string.
    - ``"undefined stereochemistry"`` — the molecule has ≥1 *unspecified* stereo
      element (chiral center or double-bond geometry). Embedding would let RDKit
      pick a stereoisomer arbitrarily, changing the chemistry silently, so we
      skip instead. A molecule with **no** stereo elements (adenine, water) is
      eligible and returns ``None``.
    """
    if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)) or str(smiles).strip() == "":
        return "no IsomericSMILES"

    from rdkit import Chem

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return "unparseable SMILES"

    # FindPotentialStereo flags every stereogenic element and whether it is
    # specified; any Unspecified element means the SMILES leaves stereo open.
    for element in Chem.FindPotentialStereo(mol):
        if element.specified == Chem.StereoSpecified.Unspecified:
            return "undefined stereochemistry"
    return None


def _conf_coords(mol, conf_id: int) -> list[tuple[str, float, float, float]]:
    """Extract (element symbol, x, y, z) tuples for one conformer, in Ångström."""
    conf = mol.GetConformer(conf_id)
    coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append((atom.GetSymbol(), pos.x, pos.y, pos.z))
    return coords


def generate_conformers(
    smiles: str,
    n_generate: int = N_GENERATE,
    rmsd_prune: float = RMSD_PRUNE,
    seed: int = SEED,
):
    """
    Embed and force-field-rank a conformer ensemble for *smiles*.

    Uses RDKit ETKDGv3 ``EmbedMultipleConfs`` with ``pruneRmsThresh=rmsd_prune``
    (RMSD-based duplicate removal is a built-in embed parameter, not separate
    code) and a fixed ``randomSeed=seed`` for reproducibility. Conformers are
    optimized and scored with MMFF94; if MMFF parameters are unavailable for the
    molecule, UFF is used as a fallback with a logged warning.

    Parameters
    ----------
    smiles : str
        SMILES string (``IsomericSMILES`` from the molecule table; stereochemistry
        preserved).
    n_generate, rmsd_prune, seed
        See module-level locked defaults.

    Returns
    -------
    (coords_list, energies_kcal, method, converged)
        ``coords_list`` — list (one per surviving conformer) of
        ``[(symbol, x, y, z), ...]`` atom tuples in Ångström.
        ``energies_kcal`` — force-field energies in kcal/mol, same order.
        ``method`` — ``"MMFF94"`` or ``"UFF"`` (the FF actually used).
        ``converged`` — per-conformer ``bool`` FF-convergence flag after the
        retry pass, same order (M-04).

    Notes
    -----
    Returning ``method`` and ``converged`` extends the 2-tuple signature in the
    plan: provenance requires recording which FF ran (MMFF94 vs UFF fallback) and
    whether each conformer's optimization actually converged — a log that implies
    success when the optimizer did not converge would violate the provenance /
    "ran ≠ validated" invariants. Recorded as a deviation in
    docs/implementation-status-v2.md.

    Any conformer flagged not-converged on the first pass is re-optimized once
    with ``FF_MAXITERS_RETRY`` iterations before its flag is finalized (M-04).

    Raises
    ------
    ValueError
        If *smiles* cannot be parsed.
    RuntimeError
        If embedding produces no conformers.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = rmsd_prune
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_generate, params=params))
    if not conf_ids:
        raise RuntimeError(
            f"RDKit embedding produced no conformers for SMILES {smiles!r} "
            f"(n_generate={n_generate}, seed={seed})"
        )

    # MMFF94 preferred; UFF only if MMFF params are missing for this molecule.
    if AllChem.MMFFHasAllMoleculeParams(mol):
        method = "MMFF94"
    else:
        method = "UFF"
        print(
            f"WARNING: MMFF94 parameters unavailable for SMILES {smiles!r}; "
            f"falling back to UFF (logged, not silent)."
        )

    # First optimization pass over the whole ensemble. Then retry ONLY the
    # conformers that failed, one at a time, with more iterations (M-04). Retrying
    # just the failed conformers — rather than re-optimizing the whole molecule in
    # place — is required so already-converged conformers keep their first-pass
    # geometry AND energy: otherwise a whole-ensemble retry would move converged
    # geometries while _finalize_convergence kept their first-pass energies,
    # leaving ranked energies describing a different geometry than the XYZ/.com we
    # write (Codex round-02 M-05).
    first = _optimize_confs(mol, method, FF_MAXITERS)
    if any(nc != 0 for nc, _ in first):
        retry = list(first)
        for i, (nc, _e) in enumerate(first):
            if nc != 0:
                retry[i] = _optimize_single_conf(
                    mol, method, conf_ids[i], FF_MAXITERS_RETRY
                )
        results = _finalize_convergence(first, retry)
    else:
        results = _finalize_convergence(first)

    coords_list = []
    energies_kcal = []
    converged = []
    for conf_id, (conv, energy) in zip(conf_ids, results):
        coords_list.append(_conf_coords(mol, conf_id))
        energies_kcal.append(float(energy))
        converged.append(bool(conv))

    return coords_list, energies_kcal, method, converged


# ---------------------------------------------------------------------------
# XYZ writer
# ---------------------------------------------------------------------------

def _write_xyz(
    xyz_path: str,
    coords: list[tuple[str, float, float, float]],
    comment: str = "",
) -> None:
    """Write one conformer to an XYZ file (coordinates in Ångström)."""
    ensure_dir(os.path.dirname(xyz_path) or ".")
    lines = [str(len(coords)), comment]
    for sym, x, y, z in coords:
        lines.append(f"{sym:<2} {x:>16.8f} {y:>12.8f} {z:>12.8f}")
    with open(xyz_path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Batch driver — Task 3
# ---------------------------------------------------------------------------

# Column order for conformer_log.csv (one row per kept conformer).
_LOG_COLUMNS = [
    "run_id",
    "artifact_id",
    "config_hash",
    "name",
    "cid",
    "smiles",
    "conformer_id",
    "rel_energy_kcalmol",
    "xyz_path",
    "xyz_sha256",
    "rdkit_version",
    "pipeline_version",
    "pipeline_commit",
    "seed",
    "n_generate",
    "top_n",
    "method",
    "n_generated",
    "n_kept",
    "rmsd_prune",
    "converged",
]

# Recorded fields that define a run's conformer-search configuration. A resumed
# molecule may be skipped ONLY if its existing rows were produced with the same
# values for all of these; otherwise its rows are stale and it is regenerated
# (M-09). `seed`, `n_generate`, `top_n`, `rmsd_prune` are the requested search
# knobs; `pipeline_version`, `pipeline_commit`, and `rdkit_version` are run-level
# guards — a code change to the generation logic or a different RDKit
# (ETKDGv3+MMFF geometry is RDKit-version-dependent, B-02 call 1b) invalidates a
# cached geometry. Clean nonblank commits must match; dirty commits are never
# reusable because one ``sha.dirty`` marker cannot identify working-tree content
# (M-20). `cid`/`smiles` are per-molecule identity and are checked separately in
# `_resume_partition`, not here.
_RESUME_CONFIG_FIELDS = (
    "seed", "n_generate", "top_n", "rmsd_prune", "pipeline_version",
    "pipeline_commit", "rdkit_version",
)


def _commit_key(value) -> str:
    """Normalize best-effort commit provenance; missing/NaN becomes blank."""
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _row_config_matches(row: dict, config: dict) -> bool:
    """
    True iff a recorded conformer_log *row* was produced with *config*.

    Missing columns (a pre-provenance / older-schema log) or unparseable values
    count as a mismatch, so such rows are conservatively regenerated rather than
    trusted. `rmsd_prune` is compared with a small float tolerance; the ints
    (seed/n_generate/top_n) and strings (`pipeline_version`, `rdkit_version`) are
    compared exactly. When both source commits are available they must also match;
    a dirty commit on either side always invalidates reuse because the marker does
    not identify the uncommitted content (M-20). Either unavailable commit also
    disables reuse; there is no version-only fallback. Per-molecule identity
    (`cid`/`smiles`) is validated by the caller, not here.
    """
    try:
        for field in ("pipeline_version", "rdkit_version"):
            if str(row[field]) != str(config[field]):
                return False
        for field in ("seed", "n_generate", "top_n"):
            if int(row[field]) != int(config[field]):
                return False
        if not math.isclose(
            float(row["rmsd_prune"]), float(config["rmsd_prune"]), rel_tol=0.0, abs_tol=1e-9
        ):
            return False
    except (KeyError, TypeError, ValueError):
        return False

    row_commit = _commit_key(row.get("pipeline_commit"))
    config_commit = _commit_key(config.get("pipeline_commit"))
    if not row_commit or not config_commit:
        return False
    if row_commit.endswith(".dirty") or config_commit.endswith(".dirty"):
        return False
    if row_commit != config_commit:
        return False
    return True


def _row_manifest_matches(row: dict, manifest_path: str, manifest: dict) -> bool:
    """Validate one retained XYZ row against exact manifest identity and bytes."""
    try:
        if str(row["run_id"]) != manifest["run_id"]:
            return False
        if str(row["config_hash"]) != manifest["config_hash"]:
            return False
        artifact_id = str(row["artifact_id"])
        artifact = find_artifact(manifest, artifact_id)
        if artifact["kind"] != "xyz":
            return False
        molecule, conformer = find_conformer_record(
            manifest, artifact["conformer_record_id"]
        )
        if str(row.get("name")) != molecule["molecule_name"]:
            return False
        if not _row_identity_matches(
            row, molecule["CID"], molecule["IsomericSMILES"]
        ):
            return False
        if _integer_key(row.get("conformer_id")) != conformer["conformer_id"]:
            return False
        if str(row.get("method")) != conformer["method"]:
            return False
        if _integer_key(row.get("n_generated")) != conformer["n_generated"]:
            return False
        if _integer_key(row.get("n_kept")) != conformer["n_kept"]:
            return False
        if not math.isclose(
            float(row.get("rel_energy_kcalmol")),
            float(conformer["relative_energy_kcalmol"]),
            rel_tol=0.0,
            abs_tol=1e-9,
        ):
            return False
        converged_value = parse_strict_bool(
            row.get("converged"),
            field_name=f"Conformer row {row.get('conformer_id')!r} converged",
        )
        if converged_value is not conformer["converged"]:
            return False
        if str(row["xyz_sha256"]) != artifact["sha256"]:
            return False
        if relative_artifact_path(str(row["xyz_path"]), manifest_path) != artifact["relative_path"]:
            return False
        return sha256_file(str(row["xyz_path"])) == artifact["sha256"]
    except (KeyError, OSError, TypeError, ValueError):
        return False


def _cid_key(x):
    """Normalize a CID for identity comparison (int or None); never raises."""
    try:
        return normalize_cid(x)
    except (ValueError, TypeError):
        return None


def _smiles_key(x) -> str:
    """Normalize a SMILES for identity comparison; NaN/None → empty string."""
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return ""
    return str(x).strip()


def _row_identity_matches(row: dict, cid, smiles) -> bool:
    """
    True iff a recorded *row* was produced for the same molecule identity
    (``cid`` and ``smiles``) requested this run (B-02). A corrected structure —
    same molecule `name` but a different CID/SMILES — is therefore treated as
    stale and regenerated, never silently reused.
    """
    return (
        _cid_key(row.get("cid")) == _cid_key(cid)
        and _smiles_key(row.get("smiles")) == _smiles_key(smiles)
    )


def _row_xyz_present(row: dict) -> bool:
    """True iff a recorded row's ``xyz_path`` exists on disk and is non-empty (B-02)."""
    path = row.get("xyz_path")
    if path is None or (isinstance(path, float) and pd.isna(path)) or str(path).strip() == "":
        return False
    return os.path.exists(str(path)) and os.path.getsize(str(path)) > 0


def _integer_key(value):
    """Normalize an integer-like CSV value; return ``None`` when invalid."""
    if isinstance(value, bool):
        return None
    try:
        if value is None or pd.isna(value):
            return None
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


def _resume_group_is_complete(
    rows: list[dict], manifest_path: str | None = None, manifest: dict | None = None
) -> bool:
    """
    Validate that one molecule's resume rows form a complete conformer set.

    A resumable group must contain the number of rows declared by ``n_kept``;
    every row must agree on that positive count; conformer IDs must be the unique
    contiguous set ``0..n_kept-1``; and XYZ paths must be unique, present, and
    non-empty (M-12). This rejects truncated, duplicated, or manually damaged
    logs even when each surviving row is individually valid.
    """
    if not rows:
        return False

    n_kept_values = [_integer_key(row.get("n_kept")) for row in rows]
    if any(value is None or value < 1 for value in n_kept_values):
        return False
    if len(set(n_kept_values)) != 1:
        return False
    n_kept = n_kept_values[0]
    if len(rows) != n_kept:
        return False

    conformer_ids = [_integer_key(row.get("conformer_id")) for row in rows]
    if any(value is None or value < 0 for value in conformer_ids):
        return False
    if sorted(conformer_ids) != list(range(n_kept)):
        return False

    xyz_paths = []
    for row in rows:
        if not _row_xyz_present(row):
            return False
        xyz_paths.append(os.path.normcase(os.path.abspath(str(row["xyz_path"]))))
    if len(set(xyz_paths)) != len(xyz_paths):
        return False
    if manifest_path is not None:
        manifest = manifest or load_manifest(manifest_path)
        if not all(_row_manifest_matches(row, manifest_path, manifest) for row in rows):
            return False
    return True


def _group_identity_is_consistent(rows: list[dict]) -> bool:
    """True iff retained rows agree on one non-empty molecule identity (M-15)."""
    if not rows:
        return False

    names = {str(row.get("name")) for row in rows}
    cid_values = {_cid_key(row.get("cid")) for row in rows}
    smiles_values = {_smiles_key(row.get("smiles")) for row in rows}
    return (
        len(names) == 1
        and len(cid_values) == 1
        and None not in cid_values
        and len(smiles_values) == 1
        and "" not in smiles_values
    )


def _carry_forward_group_is_valid(
    rows: list[dict],
    config: dict,
    manifest_path: str | None = None,
    manifest: dict | None = None,
) -> bool:
    """Validate an unrequested group before append-mode carry-forward (M-15)."""
    return (
        _resume_group_is_complete(rows, manifest_path, manifest)
        and _group_identity_is_consistent(rows)
        and all(_row_config_matches(row, config) for row in rows)
        # The commit is best-effort and may be blank, but the field itself must
        # exist so an old schema cannot be mistaken for current provenance.
        and all("pipeline_commit" in row for row in rows)
    )


def validate_unique_output_basenames(labels: list[str]) -> None:
    """Reject distinct molecule labels that map to the same output basename.

    Repeated occurrences of the same label are allowed because a conformer log
    has one row per retained conformer. Distinct labels must remain distinct
    after :func:`sanitize_basename`, and every label must produce a non-empty
    basename (B-05/B-07).
    """
    seen_labels: set[str] = set()
    seen_basenames: dict[str, str] = {}
    for value in labels:
        label = str(value)
        if label in seen_labels:
            continue
        seen_labels.add(label)
        basename = sanitize_basename(label)
        if basename == "":
            raise ValueError(
                f"Molecule label {label!r} sanitizes to an empty filename; give "
                f"it a name with at least one alphanumeric character."
            )
        if basename in seen_basenames:
            previous = seen_basenames[basename]
            raise ValueError(
                f"Molecule labels {previous!r} and {label!r} both map to output "
                f"basename {basename!r}. Use unique labels that remain distinct "
                f"after filename sanitization."
            )
        seen_basenames[basename] = label


def _resume_partition(
    existing,
    config: dict,
    requested: dict,
    preserve_unrequested: bool = False,
    manifest_path: str | None = None,
    manifest: dict | None = None,
):
    """
    Split an existing conformer_log into resumable, carry-forward, and stale sets.

    *requested* maps each molecule name requested this run to its identity dict
    ``{"cid": ..., "smiles": ...}`` (B-02): identity is per-molecule, so it is
    matched against the requested molecule, not a run-level config.

    Groups rows by molecule name and, for each molecule **requested this run**,
    keeps it only if (a) its rows form the complete, self-consistent conformer
    group declared by ``n_kept`` (M-12), (b) every row matches the run-level
    *config* (search knobs, pipeline + RDKit version — M-09 / B-02), and (c)
    every row matches the requested molecule's ``cid``/``smiles`` identity.

    Molecules **not requested this run** are dropped from the new log by default
    (M-02 call 2a): the conformer log represents the molecules requested this run,
    so a changed molecule list yields a clean current run. When
    *preserve_unrequested* is True (the ``append=True`` opt-in), each retained
    group must pass complete-group, XYZ, current-config/version, internal-identity,
    and provenance checks (M-15). Invalid unrequested groups cannot be regenerated
    from the current molecule table, so they are reported to the caller instead
    of being silently kept or dropped.

    Returns ``(done_names, kept_rows, stale_names, invalid_retained)`` where
    ``done_names`` are skipped as complete, ``kept_rows`` are preserved in the
    new log, ``stale_names`` are requested groups that must regenerate, and
    ``invalid_retained`` maps invalid unrequested group names to failure reasons.
    """
    groups: dict[str, list] = {}
    for rec in existing.to_dict("records"):
        groups.setdefault(str(rec.get("name")), []).append(rec)

    done_names: set = set()
    kept_rows: list = []
    stale_names: set = set()
    invalid_retained: dict[str, str] = {}
    for name, rows in groups.items():
        if name not in requested:
            if preserve_unrequested:
                if _carry_forward_group_is_valid(
                    rows, config, manifest_path, manifest
                ):
                    kept_rows.extend(rows)
                else:
                    invalid_retained[name] = (
                        "incomplete group, missing XYZ, config/version mismatch, "
                        "missing provenance, or inconsistent identity"
                    )
            continue
        cid = requested[name].get("cid")
        smiles = requested[name].get("smiles")
        if _resume_group_is_complete(rows, manifest_path, manifest) and all(
            _row_config_matches(r, config)
            and _row_identity_matches(r, cid, smiles)
            for r in rows
        ):
            done_names.add(name)
            kept_rows.extend(rows)
        else:
            stale_names.add(name)  # config/version/identity drift — regenerate
    return done_names, kept_rows, stale_names, invalid_retained


def search_conformers(
    molecule_table,
    xyz_dir: str = "conformer_xyz",
    log_csv: str = "conformer_log.csv",
    failed_csv: str = "conformer_search_failed.csv",
    n_generate: int = N_GENERATE,
    top_n: int = TOP_N,
    rmsd_prune: float = RMSD_PRUNE,
    seed: int = SEED,
    append: bool = False,
    manifest_path: str = "run_manifest.json",
) -> pd.DataFrame:
    """
    Generate, rank, and record the top-*top_n* distinct conformers per molecule.

    *molecule_table* is the ``build_molecule_table`` results DataFrame (or a path
    to a CSV of it). Each row must carry a ``name`` and an ``IsomericSMILES``
    column (the SMILES already retrieved in step 1); rows without a usable SMILES
    are recorded as failures, never silently skipped.

    For each molecule: embed an ensemble, FF-optimize/rank, keep the top-N
    lowest-energy distinct conformers (distinctness via ``pruneRmsThresh``),
    write ``{base}_c{ii}.xyz`` per conformer, and append one provenance row per
    kept conformer to *log_csv*. ``rel_energy_kcalmol`` is ΔE from that
    molecule's lowest-energy carried conformer (so the minimum is 0.0), in
    kcal/mol.

    Convergence (M-04): only FF-converged conformers are eligible for the top-N
    (decision 1a); each kept row records a ``converged`` bool. If **no** conformer
    converges for a molecule (even after the retry pass in
    :func:`generate_conformers`), exactly one lowest-energy best-effort geometry
    is carried with ``converged=False`` and a logged warning (decision 2b); its
    XYZ comment carries the ``UNCONVERGED_FF_SEED`` marker and its FF energy is
    unreliable.

    Resume-safe (M-09 / B-02): a molecule already in *log_csv* is skipped **only**
    when its recorded run-level config (``seed``, ``n_generate``, ``top_n``,
    ``rmsd_prune``, ``pipeline_version``, ``rdkit_version``) *and* per-molecule
    identity (``cid``, ``smiles``) match this call; its rows form the complete set
    declared by ``n_kept`` with unique contiguous IDs and unique XYZ paths; and
    every recorded ``xyz_path`` still exists and is non-empty. If any condition
    fails — or the log predates those provenance columns — its rows are treated as
    stale, dropped, and the molecule is regenerated (with a warning), so
    downstream Gaussian inputs are never built from conformers produced under a
    different configuration, for a different structure, or from a damaged log.

    Manifest coverage (B-10): with ``append=False``, *molecule_table* must match
    the immutable manifest molecule set exactly. With ``append=True``, a subset is
    allowed only when valid retained rows account for every other manifest
    molecule; no configured molecule may disappear silently. Before any output
    mutation, append mode validates
    sanitized output basenames across the current labels and all retained labels
    (B-07), then validates every retained group for completeness, XYZ existence,
    current config/version, internal identity, and provenance (M-15). A collision
    or invalid retained group raises rather than overwriting, dropping, or
    silently preserving questionable chemistry.

    Returns the full conformer log as a DataFrame.
    """
    if isinstance(molecule_table, str):
        molecule_table = pd.read_csv(molecule_table)

    # Early parameter validation (MIN-03): fail loudly at entry rather than emit a
    # nonsensical or empty search silently.
    if n_generate < 1:
        raise ValueError(f"n_generate must be >= 1, got {n_generate}.")
    if top_n < 1:
        raise ValueError(f"top_n must be >= 1, got {top_n}.")
    if rmsd_prune < 0:
        raise ValueError(f"rmsd_prune must be >= 0, got {rmsd_prune}.")
    current_labels: list[str] = []
    seen_labels: set[str] = set()
    for _, row in molecule_table.iterrows():
        name = row.get("name")
        if name is None:
            continue
        label = str(name)
        if label in seen_labels:
            raise ValueError(
                f"Duplicate molecule label {label!r} in molecule_table; labels "
                f"must be unique so each molecule maps to a distinct output set."
            )
        seen_labels.add(label)
        current_labels.append(label)

    # B-07: append mode retains unrequested rows from the existing log, so output
    # identity must be validated over the UNION of current and retained labels.
    # Load that log before any directory creation, failure-log deletion, or file
    # write, and reuse the DataFrame later for resume partitioning. This also
    # rejects an already-corrupt append log containing colliding labels.
    existing = None
    labels_to_validate = list(current_labels)
    if append and os.path.exists(log_csv):
        existing = pd.read_csv(log_csv)
        current_label_set = set(current_labels)
        labels_to_validate.extend(
            str(name)
            for name in existing.get("name", pd.Series(dtype=object)).tolist()
            if str(name) not in current_label_set
        )
    validate_unique_output_basenames(labels_to_validate)

    # Provenance captured once per run (M-06): pipeline version + best-effort git
    # commit, so conformer_log rows and XYZ files tie back to the code revision.
    pipeline_version, pipeline_commit = pipeline_provenance()
    # RDKit version is captured up front (B-02): ETKDGv3+MMFF geometry is
    # RDKit-version-dependent, so the resume check must compare each recorded
    # rdkit_version against this run's before reusing a cached geometry.
    rdkit_ver = _rdkit_version()

    # This run's conformer-search configuration; resumed rows must match it.
    run_config = {
        "method_policy": METHOD_POLICY,
        "seed": seed,
        "n_generate": n_generate,
        "top_n": top_n,
        "rmsd_prune": rmsd_prune,
        "pipeline_version": pipeline_version,
        "pipeline_commit": pipeline_commit,
        "rdkit_version": rdkit_ver,
    }
    manifest_config = {
        key: run_config[key]
        for key in ("method_policy", "seed", "n_generate", "top_n", "rmsd_prune")
    }
    manifest = assert_stage_configuration(
        manifest_path, "conformer", manifest_config
    )
    for field, actual in (
        ("pipeline_version", pipeline_version),
        ("pipeline_commit", pipeline_commit),
        ("rdkit_version", rdkit_ver),
    ):
        if str(manifest[field]) != str(actual):
            raise ValueError(
                f"Runtime {field} disagrees with run manifest; create a new manifest."
            )
    # Per-molecule identity requested this run (B-02): name → {cid, smiles}. Resume
    # matches recorded rows against the requested molecule's identity, so a
    # corrected CID/SMILES under an unchanged name regenerates instead of reusing
    # a stale geometry.
    requested = {
        str(row.get("name")): {
            "cid": row.get("cid"),
            "smiles": row.get("IsomericSMILES"),
        }
        for _, row in molecule_table.iterrows()
        if row.get("name") is not None
    }
    requested_manifest_records = [
        {
            "molecule_name": name,
            "CID": _cid_key(identity["cid"]),
            "IsomericSMILES": _smiles_key(identity["smiles"]),
            "molecule_identity_hash": molecule_identity_hash(
                name, identity["cid"], identity["smiles"]
            ),
        }
        for name, identity in requested.items()
    ]
    configured_records = manifest["configuration"]["molecules"]
    configured_identities = {
        record["molecule_identity_hash"] for record in configured_records
    }
    requested_identities = {
        record["molecule_identity_hash"] for record in requested_manifest_records
    }
    if not requested_identities.issubset(configured_identities):
        raise ValueError(
            "Runtime molecule identities disagree with run manifest; create a new manifest."
        )
    if not append and requested_identities != configured_identities:
        missing = sorted(
            record["molecule_name"]
            for record in configured_records
            if record["molecule_identity_hash"] not in requested_identities
        )
        raise ValueError(
            "append=False requires the runtime molecule table to match the run "
            "manifest exactly; missing manifest molecule(s): " + ", ".join(missing)
        )

    # Resume-safe (M-09 / B-02 / M-15): only skip a requested molecule whose
    # existing rows were
    # produced with THIS run's config AND identity, and whose XYZ files still
    # exist. Rows from a different seed/n_generate/top_n/rmsd_prune, pipeline or
    # RDKit version, a different cid/smiles, or with a missing XYZ (or a
    # pre-provenance log) are stale — drop them and regenerate so downstream never
    # builds on outdated conformers.
    if os.path.exists(log_csv):
        if existing is None:
            existing = pd.read_csv(log_csv)
        done_names, log_rows, stale_names, invalid_retained = _resume_partition(
            existing,
            run_config,
            requested,
            preserve_unrequested=append,
            manifest_path=manifest_path,
            manifest=manifest,
        )
        if invalid_retained:
            details = "; ".join(
                f"{name!r}: {reason}"
                for name, reason in sorted(invalid_retained.items())
            )
            raise ValueError(
                "append=True cannot carry forward invalid existing conformer "
                "group(s): "
                + details
                + ". Rerun with those molecules included so they can be "
                "regenerated, repair the existing log/XYZ files, or use "
                "append=False."
            )
        for name in sorted(stale_names):
            print(
                f"WARNING: {name!r} in {log_csv} has stale config/identity, "
                f"missing geometry, or an incomplete conformer group; "
                f"regenerating (stale rows dropped) so downstream inputs are "
                f"not built on invalid conformers."
            )
    else:
        done_names = set()
        log_rows = []

    if append:
        accounted_names = set(requested)
        accounted_names.update(str(row.get("name")) for row in log_rows)
        configured_names = {
            str(record["molecule_name"]) for record in configured_records
        }
        if accounted_names != configured_names:
            missing = sorted(configured_names - accounted_names)
            extra = sorted(accounted_names - configured_names)
            details = []
            if missing:
                details.append("missing manifest molecule(s): " + ", ".join(missing))
            if extra:
                details.append("unexpected molecule(s): " + ", ".join(extra))
            raise ValueError(
                "append=True requires current rows plus valid retained rows to "
                "account for the complete run manifest; " + "; ".join(details)
            )

    # M-30: every v2 output destination must stay inside the manifest package,
    # validated BEFORE the first mutation. An xyz_dir or authoritative conformer
    # log outside the package must raise before any conformer-lineage removal,
    # directory creation, failure-log deletion, XYZ write, or log rewrite. Every
    # per-conformer XYZ path is built only as os.path.join(xyz_dir, base + ...),
    # and the basename is already validated, so directory-level checks suffice.
    relative_artifact_path(xyz_dir, manifest_path)
    relative_artifact_path(log_csv, manifest_path)

    # All append validation above is intentionally complete before the first
    # output mutation (M-15). An invalid retained group must leave the prior log,
    # XYZ files, and failure log byte-for-byte unchanged.
    ensure_dir(xyz_dir)

    # Clear any stale failure log from a prior run (MIN-02) so a later clean run
    # never leaves the notebook surfacing failures that no longer apply. Rewritten
    # at the end only if this run actually has failures.
    if os.path.exists(failed_csv):
        os.remove(failed_csv)

    failed = []

    for _, row in molecule_table.iterrows():
        name = row.get("name")
        if name is None or str(name) in done_names:
            continue

        # NB: this reads the molecule TABLE's own "IsomericSMILES" column
        # (populated by pubchem._isomeric_smiles from the live "SMILES" key), not
        # a PubChem property dict — so it carries real stereo SMILES, not the dead
        # PubChem key. See MOLECULE_TABLE_COLUMNS.
        smiles = row.get("IsomericSMILES")
        cid = row.get("cid")

        # Stereo/validity gate (remediation B-01, decision 2a): skip + log rather
        # than let RDKit auto-assign stereo or embed an empty/unparseable SMILES.
        skip_reason = check_conformer_eligibility(smiles)
        if skip_reason is not None:
            remove_conformer_lineage(manifest_path, {str(name)})
            failed.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "error": skip_reason,
            })
            continue

        try:
            coords_list, energies_kcal, method, converged = generate_conformers(
                str(smiles),
                n_generate=n_generate,
                rmsd_prune=rmsd_prune,
                seed=seed,
            )
        except Exception as e:  # noqa: BLE001 — log and continue, never crash the batch
            remove_conformer_lineage(manifest_path, {str(name)})
            failed.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "error": repr(e),
            })
            continue

        n_generated = len(energies_kcal)
        # Rank/select only converged conformers (M-04 1a); if none converged,
        # carry exactly one flagged best-effort seed (2b).
        keep, all_failed = select_converged_top_n(energies_kcal, converged, top_n)
        if all_failed:
            print(
                f"WARNING: no conformer converged for {name!r} after retry; "
                f"carrying 1 best-effort {UNCONVERGED_FF_SEED} geometry — its FF "
                f"energy is unreliable and it is an unminimized DFT start."
            )
        # ΔE reference is the lowest energy among the CARRIED conformers, so the
        # kept minimum is 0.0 kcal/mol. (In the all-failed branch the single
        # carried seed is trivially its own reference.)
        e_min = min(energies_kcal[i] for i in keep) if keep else 0.0
        n_kept = len(keep)
        base = sanitize_basename(str(name))

        staging_dir = tempfile.mkdtemp(prefix=f".staging-{base}-", dir=xyz_dir)
        group_payload = []
        group_log_rows = []
        try:
            for ii, conf_idx in enumerate(keep):
                filename = f"{base}_c{ii:02d}.xyz"
                xyz_path = os.path.join(xyz_dir, filename)
                staged_xyz_path = os.path.join(staging_dir, filename)
                rel_e = round(energies_kcal[conf_idx] - e_min, 6)
                is_converged = bool(converged[conf_idx])
                unconv_tag = "" if is_converged else f" {UNCONVERGED_FF_SEED}"
                molecule_hash = molecule_identity_hash(name, cid, smiles)
                conformer_record_id = stable_record_id(
                    manifest["run_id"], "conformer", f"{molecule_hash}:{ii}"
                )
                artifact_id = stable_record_id(
                    manifest["run_id"], "xyz", conformer_record_id
                )
                _write_xyz(
                    staged_xyz_path,
                    coords_list[conf_idx],
                    comment=(
                        f"run_id={manifest['run_id']} artifact_id={artifact_id} "
                        f"config_hash={manifest['config_hash']} conformer_id={ii} "
                        f"relative_energy_kcalmol={rel_e:.6f} method={method} "
                        f"pipeline_version={pipeline_version} rdkit_version={rdkit_ver}"
                        f"{unconv_tag}"
                    ),
                )
                group_payload.append({
                    "conformer_id": ii,
                    "method": method,
                    "n_generated": n_generated,
                    "n_kept": n_kept,
                    "relative_energy_kcalmol": rel_e,
                    "converged": is_converged,
                    "xyz_path": xyz_path,
                    "staged_xyz_path": staged_xyz_path,
                    "artifact_id": artifact_id,
                })
                group_log_rows.append({
                    "run_id": manifest["run_id"],
                    "artifact_id": artifact_id,
                    "config_hash": manifest["config_hash"],
                    "name": name,
                    "cid": cid,
                    "smiles": smiles,
                    "conformer_id": ii,
                    "rel_energy_kcalmol": rel_e,
                    "xyz_path": xyz_path,
                    "rdkit_version": rdkit_ver,
                    "pipeline_version": pipeline_version,
                    "pipeline_commit": pipeline_commit,
                    "seed": seed,
                    "n_generate": n_generate,
                    "top_n": top_n,
                    "method": method,
                    "n_generated": n_generated,
                    "n_kept": n_kept,
                    "rmsd_prune": rmsd_prune,
                    "converged": is_converged,
                })

            recorded = record_conformer_group(
                manifest_path,
                name=str(name),
                cid=cid,
                smiles=str(smiles),
                conformers=group_payload,
            )
            for expected, result, log_row in zip(
                group_payload, recorded, group_log_rows
            ):
                recorded_conformer_id, xyz_digest = result
                expected_record_id = stable_record_id(
                    manifest["run_id"],
                    "conformer",
                    f"{molecule_identity_hash(name, cid, smiles)}:"
                    f"{expected['conformer_id']}",
                )
                if recorded_conformer_id != expected_record_id:
                    raise ValueError(
                        "Manifest conformer ID changed during group recording."
                    )
                log_row["xyz_sha256"] = xyz_digest
            log_rows.extend(group_log_rows)
        except Exception as e:  # noqa: BLE001 — record publication failure and continue
            failed.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "error": repr(e),
            })
        finally:
            shutil.rmtree(staging_dir, ignore_errors=True)

    out_df = pd.DataFrame(log_rows, columns=_LOG_COLUMNS)
    out_df.to_csv(log_csv, index=False)

    if failed:
        pd.DataFrame(failed).to_csv(failed_csv, index=False)
        print(f"WARNING: {len(failed)} conformer search(es) failed — see {failed_csv}")
    else:
        print("All conformer searches succeeded.")

    print(f"Wrote: {log_csv}")
    return out_df
