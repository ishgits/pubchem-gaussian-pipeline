"""Tests for pipeline.utils"""

import subprocess
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pipeline
from pipeline.utils import (
    sanitize_basename,
    normalize_cid,
    git_short_sha,
    pipeline_provenance,
)
import pipeline.utils as _utils


class TestSanitizeBasename:
    def test_simple_name(self):
        assert sanitize_basename("Adenine") == "adenine"

    def test_plus_separator(self):
        assert sanitize_basename("Adenine + Ribose") == "adenine_ribose"

    def test_comma_and_numbers(self):
        # Comma is stripped, dash becomes underscore
        assert sanitize_basename("2,6-Diaminopurine") == "26_diaminopurine"

    def test_whitespace_collapse(self):
        assert sanitize_basename("  barbituric   acid  ") == "barbituric_acid"

    def test_special_chars(self):
        result = sanitize_basename("β-D-ribofuranose (test)")
        # Should strip non-alphanumeric chars and not have double underscores
        assert "__" not in result
        assert result == result.lower()

    def test_empty_string(self):
        assert sanitize_basename("") == ""


class TestNormalizeCid:
    def test_int(self):
        assert normalize_cid(5793) == 5793

    def test_float_string(self):
        assert normalize_cid("5793.0") == 5793

    def test_none(self):
        assert normalize_cid(None) is None

    def test_nan(self):
        import math
        assert normalize_cid(float("nan")) is None


def _fake_run(rev_stdout="", rev_rc=0, status_stdout="", status_rc=0, raise_exc=None):
    """Build a fake subprocess.run that answers rev-parse and status calls."""
    def run(cmd, **kwargs):
        if raise_exc is not None:
            raise raise_exc
        if "rev-parse" in cmd:
            return subprocess.CompletedProcess(cmd, rev_rc, stdout=rev_stdout, stderr="")
        return subprocess.CompletedProcess(cmd, status_rc, stdout=status_stdout, stderr="")
    return run


class TestGitShortSha:
    """M-06: best-effort, offline-safe git provenance — never raises, no network."""

    def test_returns_empty_when_git_absent(self, monkeypatch):
        # git binary missing → FileNotFoundError, swallowed → "".
        monkeypatch.setattr(_utils.subprocess, "run",
                            _fake_run(raise_exc=FileNotFoundError("no git")))
        assert git_short_sha() == ""

    def test_returns_empty_on_nonzero_exit(self, monkeypatch):
        # Not a git repository → non-zero rev-parse → "".
        monkeypatch.setattr(_utils.subprocess, "run", _fake_run(rev_rc=128))
        assert git_short_sha() == ""

    def test_returns_empty_on_timeout(self, monkeypatch):
        monkeypatch.setattr(_utils.subprocess, "run",
                            _fake_run(raise_exc=subprocess.TimeoutExpired("git", 5)))
        assert git_short_sha() == ""

    def test_clean_tree_returns_sha(self, monkeypatch):
        monkeypatch.setattr(_utils.subprocess, "run",
                            _fake_run(rev_stdout="abc1234\n", status_stdout=""))
        assert git_short_sha() == "abc1234"

    def test_dirty_tree_appends_marker(self, monkeypatch):
        monkeypatch.setattr(_utils.subprocess, "run",
                            _fake_run(rev_stdout="abc1234\n", status_stdout=" M pipeline/x.py\n"))
        assert git_short_sha() == "abc1234.dirty"


class TestPipelineProvenance:
    """M-06: returns (version, commit); version is always the package version."""

    def test_version_and_empty_commit_when_git_fails(self, monkeypatch):
        monkeypatch.setattr(_utils, "git_short_sha", lambda cwd=None: "")
        version, commit = pipeline_provenance()
        assert version == pipeline.__version__
        assert commit == ""

    def test_version_and_commit_when_git_ok(self, monkeypatch):
        monkeypatch.setattr(_utils, "git_short_sha", lambda cwd=None: "deadbee")
        version, commit = pipeline_provenance()
        assert version == pipeline.__version__
        assert commit == "deadbee"
