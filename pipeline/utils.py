"""
Shared utility functions used across the pipeline.
"""

import os
import math
import numbers
import re
import subprocess

import pandas as pd


def parse_strict_bool(value, *, field_name: str) -> bool:
    """Parse an explicit boolean without falling back to Python truthiness.

    This accepts the representations emitted by JSON and ordinary CSV readers,
    while rejecting missing values, non-finite values, and arbitrary nonzero
    numbers.  ``field_name`` is included in every error so callers can identify
    the damaged provenance field or row.
    """
    try:
        if value is None or pd.isna(value):
            raise ValueError(f"{field_name} is missing; expected true or false.")
    except (TypeError, ValueError):
        # Array-like values are never valid scalar boolean provenance.
        raise ValueError(
            f"{field_name} has invalid value {value!r}; expected true or false."
        )

    if isinstance(value, bool):
        return value
    if hasattr(value, "item"):
        scalar = value.item()
        if scalar is not value:
            return parse_strict_bool(scalar, field_name=field_name)
    if isinstance(value, numbers.Integral) and value in (0, 1):
        return bool(value)
    if isinstance(value, numbers.Real):
        numeric = float(value)
        if math.isfinite(numeric) and numeric in (0.0, 1.0):
            return bool(int(numeric))
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no"}:
            return False
    raise ValueError(
        f"{field_name} has invalid value {value!r}; expected true or false."
    )


def git_short_sha(cwd: str | None = None) -> str:
    """
    Best-effort full git SHA of the pipeline code's HEAD, for provenance (M-06).

    Appends ``.dirty`` when the working tree has uncommitted changes (so an output
    produced from a modified tree is visibly *not* reproducible from the commit
    alone). Returns ``""`` on **any** failure — git absent, not a repository,
    subprocess error, or timeout. It never raises and never touches the network,
    so it is safe in the offline test suite (AGENTS.md §4) and on installed / HPC
    copies with no ``.git`` directory.

    Parameters
    ----------
    cwd : str, optional
        Directory to run git in. Defaults to the pipeline package directory so the
        recorded commit identifies the *code*, regardless of the process's working
        directory. A caller (or test) may point this elsewhere.
    """
    root = cwd or os.path.dirname(os.path.abspath(__file__))
    try:
        rev = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root, capture_output=True, text=True, timeout=5,
        )
        if rev.returncode != 0:
            return ""
        commit = rev.stdout.strip()
        if not commit:
            return ""
        status = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=root, capture_output=True, text=True, timeout=5,
        )
        if status.returncode == 0 and status.stdout.strip():
            commit += ".dirty"
        return commit
    except Exception:
        # Any failure (FileNotFoundError, timeout, etc.) → no commit recorded.
        return ""


def pipeline_provenance(cwd: str | None = None) -> tuple[str, str]:
    """
    Return ``(pipeline_version, git_commit)`` for provenance logging (M-06).

    ``pipeline_version`` is ``pipeline.__version__`` (a manually bumped string);
    ``git_commit`` is :func:`git_short_sha` (the legacy helper name now returns
    the full SHA; best-effort, possibly ``""``). Read
    the version lazily to avoid an import cycle with the package ``__init__``.
    """
    from . import __version__

    return __version__, git_short_sha(cwd=cwd)


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
