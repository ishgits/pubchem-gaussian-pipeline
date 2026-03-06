"""
PubChem PUG-REST helpers: name resolution, candidate scoring, fallback
queries, and SDF download.

All network calls use retry + exponential backoff + optional on-disk caching
to handle flaky connections and respect PubChem rate limits.
"""

from __future__ import annotations

import json
import os
import re
import time
from urllib.parse import quote

import pandas as pd
import requests

from .utils import ensure_dir, normalize_cid, sanitize_basename

# ---------------------------------------------------------------------------
# Module-level defaults (overridden at runtime by notebook config)
# ---------------------------------------------------------------------------
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_DEFAULT_HEADERS = {
    "User-Agent": "gaussian-input-pipeline/1.0 (research use)"
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _cache_path(cache_dir: str, key: str) -> str:
    safe = re.sub(r"[^a-zA-Z0-9_.-]+", "_", key)[:180]
    return os.path.join(cache_dir, safe + ".json")


def _get_json(
    url: str,
    *,
    cache_dir: str = ".pubchem_cache",
    cache_key: str | None = None,
    headers: dict | None = None,
    timeout_s: int = 60,
    max_retries: int = 4,
    backoff_s: float = 0.8,
    min_delay_s: float = 0.2,
) -> dict:
    """
    GET JSON from *url* with retries, exponential backoff, optional on-disk
    caching, and rate-limit politeness.
    """
    headers = headers or _DEFAULT_HEADERS

    if cache_key:
        ensure_dir(cache_dir)
        cp = _cache_path(cache_dir, cache_key)
        if os.path.exists(cp):
            with open(cp, "r") as f:
                return json.load(f)

    last_err = None
    for attempt in range(max_retries):
        try:
            time.sleep(min_delay_s)
            r = requests.get(url, headers=headers, timeout=timeout_s)
            if r.status_code == 200:
                data = r.json()
                if cache_key:
                    with open(_cache_path(cache_dir, cache_key), "w") as f:
                        json.dump(data, f)
                return data
            last_err = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:
            last_err = e
        time.sleep(backoff_s * (2 ** attempt))

    raise last_err  # type: ignore[misc]


# ---------------------------------------------------------------------------
# PubChem query functions
# ---------------------------------------------------------------------------

def get_cids_by_name(name_query: str, **kwargs) -> list[int]:
    """Resolve a molecule name to a list of PubChem CIDs."""
    q = quote(name_query)
    url = f"{PUBCHEM_BASE}/compound/name/{q}/cids/JSON"
    data = _get_json(url, cache_key=f"cids__{name_query}", **kwargs)
    return data.get("IdentifierList", {}).get("CID", [])


def get_props_by_cids(cids: list[int], **kwargs) -> list[dict]:
    """Fetch a standard set of properties for one or more CIDs."""
    props = ",".join([
        "IsomericSMILES",
        "CanonicalSMILES",
        "InChI",
        "InChIKey",
        "MolecularFormula",
        "MolecularWeight",
        "IUPACName",
        "Title",
    ])
    cid_str = ",".join(str(c) for c in cids)
    url = f"{PUBCHEM_BASE}/compound/cid/{cid_str}/property/{props}/JSON"
    data = _get_json(url, cache_key=f"props__{cid_str}", **kwargs)
    return data.get("PropertyTable", {}).get("Properties", [])


# ---------------------------------------------------------------------------
# Candidate scoring heuristic
# ---------------------------------------------------------------------------

def score_candidate(
    prop: dict,
    expected_formula: str | None = None,
    prefer_stereo: bool = True,
    keyword_boost: list[str] | None = None,
) -> int:
    """
    Score a PubChem property record for relevance.

    Higher score = better match. Factors considered:

    - Formula match with expected formula (if provided)
    - Presence of stereochemistry in the SMILES
    - CID magnitude (lower CIDs tend to be more canonical)
    - Keyword matches in IUPAC name or title
    """
    score = 0
    formula = prop.get("MolecularFormula", "")
    iso = prop.get("IsomericSMILES", "")
    iupac = (prop.get("IUPACName") or "").lower()
    title = (prop.get("Title") or "").lower()
    cid = prop.get("CID", 999_999_999)

    # Formula match
    if expected_formula and formula == expected_formula:
        score += 100

    # Stereo bonus (indicates well-defined 3D structure)
    if prefer_stereo and ("@" in iso or "/" in iso):
        score += 20

    # Lower CID = more canonical record
    if cid < 10_000:
        score += 10
    elif cid < 100_000:
        score += 5

    # Keyword boosting
    if keyword_boost:
        for kw in keyword_boost:
            kw_lower = kw.lower()
            if kw_lower in iupac or kw_lower in title:
                score += 15

    return score


# ---------------------------------------------------------------------------
# Full resolution pipeline
# ---------------------------------------------------------------------------

def resolve_pubchem_record(
    label: str,
    query: str,
    expected_formula: str | None = None,
    keyword_boost: list[str] | None = None,
    **kwargs,
) -> tuple[dict | None, dict]:
    """
    Resolve *query* via PubChem, score candidates, and return the best match.

    Returns
    -------
    (best_property_dict | None, diagnostics_dict)
    """
    info: dict = {
        "label": label,
        "query": query,
        "status": "UNKNOWN",
        "n_cids": 0,
        "selected_cid": None,
        "selected_reason": "",
        "warnings": [],
    }

    try:
        cids = get_cids_by_name(query, **kwargs)
    except Exception as e:
        info["status"] = "QUERY_FAILED"
        info["warnings"].append(str(e))
        return None, info

    if not cids:
        info["status"] = "NO_CIDS"
        return None, info

    info["n_cids"] = len(cids)

    # Fetch properties for up to the top 5 CIDs
    try:
        props = get_props_by_cids(cids[:5], **kwargs)
    except Exception as e:
        info["status"] = "PROPS_FAILED"
        info["warnings"].append(str(e))
        return None, info

    if not props:
        info["status"] = "NO_PROPS"
        return None, info

    # Score and rank
    scored = [
        (score_candidate(p, expected_formula=expected_formula, keyword_boost=keyword_boost), p)
        for p in props
    ]
    scored.sort(key=lambda x: x[0], reverse=True)
    best_score, best_prop = scored[0]

    info["status"] = "OK"
    info["selected_cid"] = best_prop.get("CID")
    info["selected_reason"] = f"score={best_score}"

    if expected_formula:
        actual = best_prop.get("MolecularFormula", "")
        if actual != expected_formula:
            info["warnings"].append(f"formula_mismatch: expected={expected_formula}, got={actual}")

    return best_prop, info


def resolve_with_fallback(
    label: str,
    primary_query: str,
    fallback_queries: dict[str, list[str]] | None = None,
    expected_formula: str | None = None,
    **kwargs,
) -> tuple[dict | None, dict]:
    """
    Try *primary_query* first, then iterate through any fallback queries
    defined for *label*.
    """
    prop, info = resolve_pubchem_record(label, primary_query, expected_formula=expected_formula, **kwargs)
    if prop is not None:
        return prop, info

    if fallback_queries:
        for q in fallback_queries.get(label, []):
            prop2, info2 = resolve_pubchem_record(label, q, expected_formula=expected_formula, **kwargs)
            if prop2 is not None:
                info2["selected_reason"] = f"fallback_query={q}"
                return prop2, info2

    return None, info


# ---------------------------------------------------------------------------
# Table builder (Step 1)
# ---------------------------------------------------------------------------

def build_molecule_table(
    molecules: list[str],
    alias: dict[str, str] | None = None,
    fallback_queries: dict[str, list[str]] | None = None,
    expected_formulas: dict[str, str] | None = None,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Resolve every molecule in the list and return (results_df, diagnostics_df).

    Parameters
    ----------
    molecules : list[str]
        Human-readable molecule labels.
    alias : dict
        Maps label → PubChem query string (for names that differ from common usage).
    fallback_queries : dict
        Maps label → list of alternative query strings to try on failure.
    expected_formulas : dict
        Maps label → expected molecular formula (used for scoring).
    """
    alias = alias or {}
    fallback_queries = fallback_queries or {}
    expected_formulas = expected_formulas or {}

    rows = []
    diagnostics = []

    for label in molecules:
        query = alias.get(label, label)
        exp_formula = expected_formulas.get(label, None)

        prop, info = resolve_with_fallback(
            label, query,
            fallback_queries=fallback_queries,
            expected_formula=exp_formula,
            **kwargs,
        )
        diagnostics.append(info)

        if prop is None:
            rows.append({
                "name": label,
                "pubchem_query": query,
                "cid": None,
                "formula": None,
                "iupac_name": None,
                "title": None,
                "status": info["status"],
                "warnings": "; ".join(info.get("warnings", [])),
            })
        else:
            rows.append({
                "name": label,
                "pubchem_query": query,
                "cid": prop.get("CID"),
                "formula": prop.get("MolecularFormula"),
                "iupac_name": prop.get("IUPACName"),
                "title": prop.get("Title"),
                "status": info["status"],
                "warnings": "; ".join(info.get("warnings", [])),
            })

    return pd.DataFrame(rows), pd.DataFrame(diagnostics)


# ---------------------------------------------------------------------------
# SDF download (Step 2)
# ---------------------------------------------------------------------------

def download_pubchem_sdf(
    cid: int,
    outpath: str,
    prefer_3d: bool = True,
    headers: dict | None = None,
    timeout_s: int = 60,
    min_delay_s: float = 0.2,
) -> str:
    """
    Download SDF from PubChem by CID.

    Tries the 3D record first (if *prefer_3d*); falls back to 2D.
    Returns the URL that succeeded.
    """
    headers = headers or _DEFAULT_HEADERS
    ensure_dir(os.path.dirname(outpath) or ".")

    urls = []
    if prefer_3d:
        urls.append(f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d")
    urls.append(f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF")

    last_err = None
    for url in urls:
        try:
            time.sleep(min_delay_s)
            r = requests.get(url, headers=headers, timeout=timeout_s)
            if r.status_code == 200 and len(r.text) > 0:
                with open(outpath, "w") as f:
                    f.write(r.text)
                return url
            last_err = f"HTTP {r.status_code}: {r.text[:200]}"
        except Exception as e:
            last_err = repr(e)

    raise RuntimeError(f"Failed to download SDF for CID {cid}. Last error: {last_err}")


def download_sdfs(
    input_csv: str,
    sdf_dir: str = "pubchem_sdf",
    log_csv: str = "sdf_download_log.csv",
    **kwargs,
) -> pd.DataFrame:
    """
    Batch-download SDFs for every resolved CID in *input_csv*.

    Resume-safe: skips files that already exist.
    """
    df = pd.read_csv(input_csv)
    df = df[df["cid"].notna()].copy()
    df["cid_int"] = df["cid"].apply(normalize_cid)

    ensure_dir(sdf_dir)
    download_log = []
    failed = []

    for _, row in df.iterrows():
        name = row["name"]
        cid = row["cid_int"]
        base = sanitize_basename(name)
        sdf_path = os.path.join(sdf_dir, f"{base}.sdf")

        # Resume-safe
        if os.path.exists(sdf_path) and os.path.getsize(sdf_path) > 0:
            download_log.append({
                "name": name, "cid": cid,
                "sdf_path": sdf_path, "source_url": "SKIPPED_EXISTS",
            })
            continue

        try:
            src_url = download_pubchem_sdf(cid, sdf_path, **kwargs)
            download_log.append({
                "name": name, "cid": cid,
                "sdf_path": sdf_path, "source_url": src_url,
            })
        except Exception as e:
            failed.append({"name": name, "cid": cid, "error": repr(e)})

    log_df = pd.DataFrame(download_log)
    log_df.to_csv(log_csv, index=False)

    if failed:
        fail_df = pd.DataFrame(failed)
        fail_df.to_csv("sdf_download_failed.csv", index=False)
        print(f"WARNING: {len(failed)} downloads failed — see sdf_download_failed.csv")
    else:
        print("All SDF downloads succeeded.")

    print(f"Wrote: {log_csv}")
    return log_df
