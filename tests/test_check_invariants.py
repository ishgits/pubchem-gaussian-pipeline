"""Tests for scripts/check_invariants.py (the mechanical floor).

check_invariants.py is a standalone script, not a package module, so it is loaded
here by file path. The status-doc drift guard (M-07) is tested through its pure
helper `_status_doc_problems`, which takes text — no file IO, fully offline.
"""

import importlib.util
import pathlib

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "check_invariants.py"
_spec = importlib.util.spec_from_file_location("check_invariants", SCRIPT)
check_invariants = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_invariants)


TEMPLATE = """# implementation-status.md

**PR:** #<n>   **Branch:** <name>   **Round:** <n>

## 1. What was implemented
-

## 3. Deviations from architecture.md
- <none | describe each deviation and the reason>

## Provenance
- pipeline version: <>
"""

POPULATED = """# implementation-status.md

**PR:** #3   Branch: feat/conformer-search-v2   Round: 3

## 1. What was implemented
- Added pipeline_version + pipeline_commit provenance columns.

## 3. Deviations from architecture.md
- generate_conformers returns a 4-tuple (adds method + converged).

## Provenance
- pipeline version: 0.2.0
"""


class TestStatusDocDriftGuard:
    def test_fires_on_template(self):
        problems = check_invariants._status_doc_problems(TEMPLATE)
        assert problems  # non-empty: placeholders + empty bullets detected
        joined = " ".join(problems)
        assert "template placeholder" in joined
        assert "empty bullet" in joined

    def test_passes_on_populated(self):
        assert check_invariants._status_doc_problems(POPULATED) == []

    def test_detects_pr_placeholder(self):
        # The "#<n>" PR marker is caught as an angle-bracket placeholder.
        problems = check_invariants._status_doc_problems("**PR:** #<n>")
        assert any("<n>" in p for p in problems)

    def test_detects_empty_bullet_only(self):
        problems = check_invariants._status_doc_problems("## X\n- \n- real item\n")
        assert len(problems) == 1
        assert "empty bullet" in problems[0]
