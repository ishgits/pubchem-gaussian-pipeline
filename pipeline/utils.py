"""
Shared utility functions used across the pipeline.
"""

import os
import re

import pandas as pd


def sanitize_basename(name: str) -> str:
    """
    Convert a molecule name into a filesystem-safe, lowercase basename.

    Examples
    --------
    >>> sanitize_basename("Adenine + Ribose")
    'adenine_ribose'
    >>> sanitize_basename("2,6-Diaminopurine")
    '26_diaminopurine'
    """
    s = name.strip().lower()
    s = s.replace(" + ", "_").replace("+", "_").replace(",", "")
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^a-z0-9_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def ensure_dir(path: str) -> str:
    """Create directory (and parents) if it doesn't exist. Returns the path."""
    os.makedirs(path, exist_ok=True)
    return path


def normalize_cid(x) -> int | None:
    """Coerce a CID value (possibly float-string from CSV) to int, or None."""
    if pd.isna(x):
        return None
    return int(float(str(x).strip()))
