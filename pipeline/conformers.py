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
    artifact_abspath,
    assert_stage_configuration,
    find_artifact,
    find_conformer_record,
    load_manifest,
    molecule_identity_hash,
    record_conformer_group,
    relative_artifact_path,
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
# Provisional undefined-stereo structure (v2.1 judgment #5)
# ---------------------------------------------------------------------------

# Status recorded in conformer_log.csv / the manifest molecule record and used
# downstream to mark provisional artifacts (dE=NA, PROVISIONAL marker).
PROVENANCE_NORMAL = "normal"
PROVENANCE_PROVISIONAL = "provisional_undefined_stereo"


def undefined_stereo_labels(mol) -> list[str]:
    """Return human labels for each *unspecified* stereo element of *mol*.

    Atom stereocentres are labeled ``atom <idx>`` and double-bond stereo elements
    ``bond <a>-<b>``. The count of returned labels is ``k`` — the number of
    undefined centres — so ``2**k`` is the number of possible stereoisomers.
    """
    from rdkit import Chem

    labels: list[str] = []
    for element in Chem.FindPotentialStereo(mol):
        if element.specified != Chem.StereoSpecified.Unspecified:
            continue
        if "Bond" in str(element.type):
            bond = mol.GetBondWithIdx(int(element.centeredOn))
            labels.append(
                f"bond {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}"
            )
        else:
            labels.append(f"atom {int(element.centeredOn)}")
    return labels


def generate_provisional_conformer(smiles: str, seed: int = SEED):
    """Embed ONE provisional structure for an undefined-stereo molecule (§9).

    This is explicitly **not** an ensemble conformer search. RDKit fixes the
    undefined centre(s) to an arbitrary configuration at embed time; the
    post-embed isomeric SMILES is read back as the *arbitrated* structure.

    Returns ``(coords, method, converged, arbitrated_smiles, undefined_labels)``:

    - ``coords`` — ``[(symbol, x, y, z), ...]`` in Ångström for the one structure;
    - ``method`` — ``"MMFF94"`` or ``"UFF"`` (the FF used for the light cleanup);
    - ``converged`` — whether the light FF cleanup converged;
    - ``arbitrated_smiles`` — isomeric SMILES read back from the embedded 3D
      geometry (the arbitrary configuration RDKit chose), distinct from the
      underspecified PubChem SMILES;
    - ``undefined_labels`` — labels of the undefined stereo centre(s).

    Raises ``ValueError`` if *smiles* cannot be parsed and ``RuntimeError`` if the
    single embed fails.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    undefined_labels = undefined_stereo_labels(mol)

    molh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    # A single ETKDG embed (fixed seed) — one structure, not an ensemble.
    conf_id = AllChem.EmbedMolecule(molh, params)
    if conf_id < 0:
        raise RuntimeError(
            f"RDKit embedding produced no provisional structure for SMILES "
            f"{smiles!r} (seed={seed})"
        )

    # Light FF cleanup. MMFF94 preferred; UFF recorded fallback (never silent).
    if AllChem.MMFFHasAllMoleculeParams(molh):
        method = "MMFF94"
        not_converged = AllChem.MMFFOptimizeMolecule(
            molh, mmffVariant="MMFF94", maxIters=FF_MAXITERS, confId=conf_id
        )
    else:
        method = "UFF"
        print(
            f"WARNING: MMFF94 parameters unavailable for provisional SMILES "
            f"{smiles!r}; using UFF for the light cleanup (logged, not silent)."
        )
        not_converged = AllChem.UFFOptimizeMolecule(
            molh, confId=conf_id, maxIters=FF_MAXITERS
        )
    converged = not_converged == 0

    coords = _conf_coords(molh, conf_id)

    # Read back the arbitrated configuration from the embedded 3D geometry.
    heavy = Chem.RemoveHs(molh)
    Chem.AssignStereochemistryFrom3D(heavy)
    arbitrated_smiles = Chem.MolToSmiles(heavy)

    return coords, method, converged, arbitrated_smiles, undefined_labels


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
    "provenance_status",
    "undefined_centers",
    "pubchem_smiles",
    "arbitrated_smiles",
]


def _manifest_conformer_log_rows(
    manifest_path: str,
    manifest: dict | None = None,
    *,
    verify_xyz: bool = False,
) -> list[dict]:
    """Reconstruct the canonical conformer-log rows from manifest authority.

    The CSV is a subordinate index, so every field that overlaps manifest
    authority is derived from the manifest rather than copied from a possibly
    stale prior CSV.  ``verify_xyz=True`` additionally requires every recorded
    XYZ to exist as a nonempty regular file whose bytes match the manifest.
    """
    manifest = manifest or load_manifest(manifest_path)
    artifacts = {
        artifact["artifact_id"]: artifact
        for artifact in manifest["artifacts"]
        if artifact["kind"] == "xyz"
    }
    rows = []
    for molecule in manifest["molecules"]:
        provenance_status = molecule.get("provenance_status", "normal")
        undefined_centers = molecule.get("undefined_centers", "") or ""
        pubchem_smiles = molecule.get("pubchem_smiles", "") or ""
        arbitrated_smiles = molecule.get("arbitrated_smiles", "") or ""
        for conformer in sorted(
            molecule["conformers"], key=lambda record: record["conformer_id"]
        ):
            artifact_id = conformer["xyz_artifact_id"]
            artifact = artifacts.get(artifact_id)
            if artifact is None:
                raise ValueError(
                    f"Manifest conformer is missing XYZ artifact {artifact_id!r}."
                )
            xyz_path = artifact_abspath(manifest_path, artifact["relative_path"])
            if verify_xyz:
                if not os.path.isfile(xyz_path) or os.path.getsize(xyz_path) == 0:
                    raise ValueError(
                        f"Manifest XYZ is missing, irregular, or empty: {xyz_path!r}."
                    )
                if sha256_file(xyz_path) != artifact["sha256"]:
                    raise ValueError(
                        f"Manifest XYZ hash mismatch: {artifact_id!r}."
                    )
            rows.append({
                "run_id": manifest["run_id"],
                "artifact_id": artifact_id,
                "config_hash": manifest["config_hash"],
                "name": molecule["molecule_name"],
                "cid": molecule["CID"],
                "smiles": molecule["IsomericSMILES"],
                "conformer_id": conformer["conformer_id"],
                "rel_energy_kcalmol": conformer["relative_energy_kcalmol"],
                "xyz_path": xyz_path,
                "xyz_sha256": artifact["sha256"],
                "rdkit_version": manifest["rdkit_version"],
                "pipeline_version": manifest["pipeline_version"],
                "pipeline_commit": manifest["pipeline_commit"],
                "seed": conformer["seed"],
                "n_generate": conformer["n_generate"],
                "top_n": conformer["top_n"],
                "method": conformer["method"],
                "n_generated": conformer["n_generated"],
                "n_kept": conformer["n_kept"],
                "rmsd_prune": conformer["rmsd_prune"],
                "converged": conformer["converged"],
                "provenance_status": provenance_status,
                "undefined_centers": undefined_centers,
                "pubchem_smiles": pubchem_smiles,
                "arbitrated_smiles": arbitrated_smiles,
            })
    return rows


def _stage_conformer_log(rows: list[dict], log_csv: str) -> str:
    """Write a complete candidate conformer log to a same-package temp file."""
    parent = os.path.dirname(os.path.abspath(log_csv))
    os.makedirs(parent, exist_ok=True)
    fd, staged_path = tempfile.mkstemp(
        prefix=".conformer_log.", suffix=".csv", dir=parent
    )
    os.close(fd)
    try:
        pd.DataFrame(rows, columns=_LOG_COLUMNS).to_csv(staged_path, index=False)
        if not os.path.isfile(staged_path) or os.path.getsize(staged_path) == 0:
            raise OSError("Staged conformer log was not written completely.")
    except Exception:
        try:
            os.unlink(staged_path)
        except OSError:
            pass
        raise
    return staged_path


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


def search_conformers(
    molecule_table,
    xyz_dir: str = "conformer_xyz",
    log_csv: str = "conformer_log.csv",
    failed_csv: str = "conformer_search_failed.csv",
    n_generate: int = N_GENERATE,
    top_n: int = TOP_N,
    rmsd_prune: float = RMSD_PRUNE,
    seed: int = SEED,
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

    Undefined stereo (v2.1 judgment #5): a molecule whose ``IsomericSMILES``
    leaves stereochemistry unspecified is no longer skipped. It takes the
    provisional path — one arbitrated structure embedded from the PubChem SMILES
    (:func:`generate_provisional_conformer`), ``dE=NA``,
    ``provenance_status=provisional_undefined_stereo``, ``undefined_centers`` /
    ``arbitrated_smiles`` recorded, a PROVISIONAL marker in the XYZ, and a loud
    console warning. It is never presented as the defined stereoisomer.

    No resume/append (v2.1, contract §7): every run is fresh and immutable.
    Pointing this at a run folder whose manifest is already populated with
    conformers raises rather than reusing or appending — start a new run instead.

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

    # §6 collision safety, now within the single run folder: distinct labels must
    # remain distinct after sanitization before any output mutation.
    validate_unique_output_basenames(current_labels)

    # Provenance captured once per run (M-06): pipeline version + best-effort git
    # commit, so conformer_log rows and XYZ files tie back to the code revision.
    pipeline_version, pipeline_commit = pipeline_provenance()
    rdkit_ver = _rdkit_version()

    manifest_config = {
        "method_policy": METHOD_POLICY,
        "seed": seed,
        "n_generate": n_generate,
        "top_n": top_n,
        "rmsd_prune": rmsd_prune,
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

    # v2.1 (contract §7): resume/append is removed. An already-populated run
    # folder must not be reused, appended to, or repaired — it raises so the user
    # starts a fresh run. A fresh manifest has no conformer records and no
    # artifacts, so this never fires in normal single-run use.
    if any(molecule["conformers"] for molecule in manifest["molecules"]) or manifest["artifacts"]:
        raise ValueError(
            "This run folder is already populated with conformer records; v2.1 "
            "has no resume/append path. Start a new run with a fresh run_id "
            "instead of reusing an existing run package."
        )

    # Per-molecule identity requested this run: name → {cid, smiles}. Must match
    # the immutable manifest molecule set exactly (no subset without append).
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
    if requested_identities != configured_identities:
        missing = sorted(
            record["molecule_name"]
            for record in configured_records
            if record["molecule_identity_hash"] not in requested_identities
        )
        raise ValueError(
            "The runtime molecule table must match the run manifest exactly; "
            "missing manifest molecule(s): " + ", ".join(missing)
        )

    log_rows: list[dict] = []

    # M-30: every v2 output destination must stay inside the manifest package,
    # validated BEFORE the first mutation. An xyz_dir or authoritative conformer
    # log outside the package must raise before any directory creation,
    # failure-log deletion, XYZ write, or log rewrite.
    relative_artifact_path(xyz_dir, manifest_path)
    relative_artifact_path(log_csv, manifest_path)

    ensure_dir(xyz_dir)

    # Clear any stale failure log from a prior run (MIN-02) so a later clean run
    # never leaves the notebook surfacing failures that no longer apply. Rewritten
    # at the end only if this run actually has failures.
    if os.path.exists(failed_csv):
        os.remove(failed_csv)

    failed = []

    for _, row in molecule_table.iterrows():
        name = row.get("name")
        if name is None:
            continue

        # NB: this reads the molecule TABLE's own "IsomericSMILES" column
        # (populated by pubchem._isomeric_smiles from the live "SMILES" key), not
        # a PubChem property dict — so it carries real stereo SMILES, not the dead
        # PubChem key. See MOLECULE_TABLE_COLUMNS.
        smiles = row.get("IsomericSMILES")
        cid = row.get("cid")

        skip_reason = check_conformer_eligibility(smiles)

        # Per-molecule normalized publication inputs. The provisional path (§9)
        # and the normal ensemble path both feed the SAME publication block below
        # — the only difference rides on `provisional` and the marker fields.
        provisional = False
        undefined_centers_label = ""
        arbitrated_smiles = ""
        pubchem_smiles = ""

        if skip_reason is None:
            try:
                coords_list, energies_kcal, method, converged = generate_conformers(
                    str(smiles),
                    n_generate=n_generate,
                    rmsd_prune=rmsd_prune,
                    seed=seed,
                )
            except Exception as e:  # noqa: BLE001 — log and continue, never crash the batch
                failed.append({
                    "name": name, "cid": cid, "smiles": smiles, "error": repr(e),
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
            # ΔE reference is the lowest energy among the CARRIED conformers, so
            # the kept minimum is 0.0 kcal/mol.
            e_min = min(energies_kcal[i] for i in keep) if keep else 0.0
        elif skip_reason == "undefined stereochemistry":
            # v2.1 judgment #5: embed ONE provisional, loudly-flagged structure
            # rather than skipping. Never a silent stereo guess.
            try:
                (
                    prov_coords, method, prov_converged, arbitrated_smiles,
                    undefined_labels,
                ) = generate_provisional_conformer(str(smiles), seed=seed)
            except Exception as e:  # noqa: BLE001 — log and continue
                failed.append({
                    "name": name, "cid": cid, "smiles": smiles, "error": repr(e),
                })
                continue
            provisional = True
            pubchem_smiles = str(smiles)
            coords_list = [prov_coords]
            energies_kcal = [0.0]
            converged = [prov_converged]
            keep = [0]
            n_generated = 1
            e_min = 0.0
            undefined_centers_label = (
                ", ".join(undefined_labels)
                if undefined_labels
                else "unspecified center(s)"
            )
            k = len(undefined_labels)
            isomer_note = (
                f" ({k} undefined centres → {2 ** k} isomers possible, one arbitrated)"
                if k > 1
                else ""
            )
            print(
                f"WARNING: {name!r} has undefined stereochemistry at "
                f"{undefined_centers_label}{isomer_note}; conformer search skipped. "
                f"Embedded ONE provisional structure from the PubChem SMILES with an "
                f"ARBITRARY choice at the undefined centre(s) — an unvalidated DFT "
                f"start, NOT the compound's real configuration. dE=NA."
            )
        else:
            # no IsomericSMILES / unparseable SMILES — skip + log (never guess).
            failed.append({
                "name": name, "cid": cid, "smiles": smiles, "error": skip_reason,
            })
            continue

        n_kept = len(keep)
        base = sanitize_basename(str(name))
        provenance_status = (
            PROVENANCE_PROVISIONAL if provisional else PROVENANCE_NORMAL
        )

        staging_dir = tempfile.mkdtemp(prefix=f".staging-{base}-", dir=xyz_dir)
        staged_log_path = None
        group_payload = []
        group_log_rows = []
        try:
            for ii, conf_idx in enumerate(keep):
                filename = f"{base}_c{ii:02d}.xyz"
                xyz_path = os.path.join(xyz_dir, filename)
                staged_xyz_path = os.path.join(staging_dir, filename)
                # Provisional: dE=NA in artifacts (no ensemble reference), 0.0 in
                # the manifest/log (its own reference within a group of one).
                rel_e = 0.0 if provisional else round(energies_kcal[conf_idx] - e_min, 6)
                is_converged = bool(converged[conf_idx])
                unconv_tag = "" if is_converged else f" {UNCONVERGED_FF_SEED}"
                molecule_hash = molecule_identity_hash(name, cid, smiles)
                conformer_record_id = stable_record_id(
                    manifest["run_id"], "conformer", f"{molecule_hash}:{ii}"
                )
                artifact_id = stable_record_id(
                    manifest["run_id"], "xyz", conformer_record_id
                )
                # v2.1 per-artifact metadata (contract §5): XYZ line 2 carries
                # only inline science (dE, method) + the one artifact_id.
                if provisional:
                    de_field = "dE=NA"
                    marker = f" PROVISIONAL: stereo arbitrated at {undefined_centers_label}"
                else:
                    de_field = f"dE={rel_e:.6f} kcal/mol"
                    marker = ""
                _write_xyz(
                    staged_xyz_path,
                    coords_list[conf_idx],
                    comment=(
                        f"{de_field} method={method} artifact_id={artifact_id}"
                        f"{unconv_tag}{marker}"
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
                    "provenance_status": provenance_status,
                    "undefined_centers": undefined_centers_label,
                    "pubchem_smiles": pubchem_smiles,
                    "arbitrated_smiles": arbitrated_smiles,
                })

            for expected, log_row in zip(group_payload, group_log_rows):
                log_row["xyz_sha256"] = sha256_file(
                    expected["staged_xyz_path"]
                )

            # Build the complete candidate subordinate index from current
            # manifest authority, replacing only this molecule's rows.  This
            # prevents a caught group failure from publishing a partial CSV and
            # guarantees the staged log indexes the complete candidate manifest.
            current_manifest = load_manifest(manifest_path)
            candidate_log_rows = [
                row
                for row in _manifest_conformer_log_rows(
                    manifest_path, current_manifest
                )
                if row["name"] != str(name)
            ]
            candidate_log_rows.extend(group_log_rows)
            molecule_order = {
                molecule["molecule_name"]: index
                for index, molecule in enumerate(current_manifest["molecules"])
            }
            candidate_log_rows.sort(
                key=lambda item: (
                    molecule_order.get(str(item["name"]), len(molecule_order)),
                    int(item["conformer_id"]),
                )
            )
            staged_log_path = _stage_conformer_log(
                candidate_log_rows, log_csv
            )

            recorded = record_conformer_group(
                manifest_path,
                name=str(name),
                cid=cid,
                smiles=str(smiles),
                conformers=group_payload,
                conformer_log_path=log_csv,
                staged_conformer_log_path=staged_log_path,
                provenance_status=provenance_status,
                undefined_centers=undefined_centers_label if provisional else None,
                pubchem_smiles=pubchem_smiles if provisional else None,
                arbitrated_smiles=arbitrated_smiles if provisional else None,
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
                if log_row["xyz_sha256"] != xyz_digest:
                    raise ValueError(
                        "Manifest XYZ digest changed during group recording."
                    )
            log_rows.extend(group_log_rows)
        except Exception as e:  # noqa: BLE001 — record publication failure and continue
            failed.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "error": repr(e),
            })
        finally:
            if staged_log_path is not None:
                try:
                    os.unlink(staged_log_path)
                except OSError:
                    pass
            shutil.rmtree(staging_dir, ignore_errors=True)

    # Return the committed package state, never the in-memory rows from attempted
    # publications.  Each successful group transaction has already replaced the
    # complete log; each failed transaction left the prior complete log intact.
    if os.path.isfile(log_csv):
        out_df = pd.read_csv(log_csv)
        if list(out_df.columns) != _LOG_COLUMNS:
            raise ValueError("Committed conformer log has an unexpected schema.")
        out_df["pipeline_commit"] = (
            out_df["pipeline_commit"].fillna("").astype(str)
        )
    else:
        # Zero-job run (M-11): no molecule produced a conformer group, so no
        # group transaction wrote the log. Write a header-only log so the
        # downstream Gaussian stage can consume a valid, empty index.
        out_df = pd.DataFrame(columns=_LOG_COLUMNS)
        out_df.to_csv(log_csv, index=False)

    if failed:
        pd.DataFrame(failed).to_csv(failed_csv, index=False)
        print(f"WARNING: {len(failed)} conformer search(es) failed — see {failed_csv}")
    else:
        print("All conformer searches succeeded.")

    print(f"Wrote: {log_csv}")
    return out_df
