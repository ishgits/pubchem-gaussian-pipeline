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

import os

import pandas as pd

from .utils import ensure_dir, sanitize_basename

# Locked defaults (docs/implementation-plan.md v2 — confirm at approval).
N_GENERATE = 20
TOP_N = 3
RMSD_PRUNE = 0.5  # Ångström
SEED = 42


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
# RDKit embed + rank core — Task 2
# ---------------------------------------------------------------------------

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
    (coords_list, energies_kcal, method)
        ``coords_list`` — list (one per surviving conformer) of
        ``[(symbol, x, y, z), ...]`` atom tuples in Ångström.
        ``energies_kcal`` — force-field energies in kcal/mol, same order.
        ``method`` — ``"MMFF94"`` or ``"UFF"`` (the FF actually used).

    Notes
    -----
    Returning ``method`` is a small, deliberate extension of the 2-tuple
    signature in the plan: provenance requires recording which FF actually ran
    (MMFF94 vs UFF fallback), and that decision is made here. Recorded as a
    deviation in docs/implementation-status-v2.md.

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
        results = AllChem.MMFFOptimizeMoleculeConfs(
            mol, mmffVariant="MMFF94", maxIters=2000
        )
    else:
        method = "UFF"
        print(
            f"WARNING: MMFF94 parameters unavailable for SMILES {smiles!r}; "
            f"falling back to UFF (logged, not silent)."
        )
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=2000)

    # results: list of (not_converged, energy_kcal_per_mol), aligned with conf_ids.
    coords_list = []
    energies_kcal = []
    for conf_id, (_not_converged, energy) in zip(conf_ids, results):
        coords_list.append(_conf_coords(mol, conf_id))
        energies_kcal.append(float(energy))

    return coords_list, energies_kcal, method


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
    "name",
    "cid",
    "smiles",
    "conformer_id",
    "rel_energy_kcalmol",
    "xyz_path",
    "rdkit_version",
    "seed",
    "method",
    "n_generated",
    "n_kept",
    "rmsd_prune",
]


def search_conformers(
    molecule_table,
    xyz_dir: str = "conformer_xyz",
    log_csv: str = "conformer_log.csv",
    failed_csv: str = "conformer_search_failed.csv",
    n_generate: int = N_GENERATE,
    top_n: int = TOP_N,
    rmsd_prune: float = RMSD_PRUNE,
    seed: int = SEED,
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
    molecule's lowest-energy conformer (so the minimum is 0.0), in kcal/mol.

    Resume-safe: molecules already present in an existing *log_csv* are skipped.

    Returns the full conformer log as a DataFrame.
    """
    if isinstance(molecule_table, str):
        molecule_table = pd.read_csv(molecule_table)

    ensure_dir(xyz_dir)

    # Resume-safe: preserve any existing log rows and skip those molecules.
    if os.path.exists(log_csv):
        existing = pd.read_csv(log_csv)
        done_names = set(existing["name"].astype(str)) if "name" in existing else set()
        log_rows = existing.to_dict("records")
    else:
        done_names = set()
        log_rows = []

    failed = []
    rdkit_ver = None

    for _, row in molecule_table.iterrows():
        name = row.get("name")
        if name is None or str(name) in done_names:
            continue

        smiles = row.get("IsomericSMILES")
        cid = row.get("cid")

        # Stereo/validity gate (remediation B-01, decision 2a): skip + log rather
        # than let RDKit auto-assign stereo or embed an empty/unparseable SMILES.
        skip_reason = check_conformer_eligibility(smiles)
        if skip_reason is not None:
            failed.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "error": skip_reason,
            })
            continue

        try:
            if rdkit_ver is None:
                rdkit_ver = _rdkit_version()
            coords_list, energies_kcal, method = generate_conformers(
                str(smiles),
                n_generate=n_generate,
                rmsd_prune=rmsd_prune,
                seed=seed,
            )
        except Exception as e:  # noqa: BLE001 — log and continue, never crash the batch
            failed.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "error": repr(e),
            })
            continue

        n_generated = len(energies_kcal)
        keep = select_top_n(energies_kcal, top_n)
        e_min = min(energies_kcal)  # ΔE reference for this molecule (kcal/mol)
        n_kept = len(keep)
        base = sanitize_basename(str(name))

        for ii, conf_idx in enumerate(keep):
            xyz_path = os.path.join(xyz_dir, f"{base}_c{ii:02d}.xyz")
            rel_e = energies_kcal[conf_idx] - e_min
            _write_xyz(
                xyz_path,
                coords_list[conf_idx],
                comment=(
                    f"{base} c{ii:02d} dE={rel_e:.4f} kcal/mol "
                    f"method={method} seed={seed}"
                ),
            )
            log_rows.append({
                "name": name,
                "cid": cid,
                "smiles": smiles,
                "conformer_id": ii,
                "rel_energy_kcalmol": round(rel_e, 6),
                "xyz_path": xyz_path,
                "rdkit_version": rdkit_ver,
                "seed": seed,
                "method": method,
                "n_generated": n_generated,
                "n_kept": n_kept,
                "rmsd_prune": rmsd_prune,
            })

    out_df = pd.DataFrame(log_rows, columns=_LOG_COLUMNS)
    out_df.to_csv(log_csv, index=False)

    if failed:
        pd.DataFrame(failed).to_csv(failed_csv, index=False)
        print(f"WARNING: {len(failed)} conformer search(es) failed — see {failed_csv}")
    else:
        print("All conformer searches succeeded.")

    print(f"Wrote: {log_csv}")
    return out_df
