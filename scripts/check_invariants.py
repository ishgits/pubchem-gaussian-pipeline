#!/usr/bin/env python3
"""
Mechanical scientific-invariant checks for pubchem-gaussian-pipeline.

The OBJECTIVE floor that runs in CI before any human/agent review, so reviewers
spend attention on judgment, not on "did you leave a TODO." Grow this file as
invariants become code-checkable — it is meant to accrete.

Design rule for this file: NEVER false-fire on clean code. A flaky red floor is
worse than no floor, because the whole workflow keys off "green = review-ready."

Exit non-zero on any violation.
"""
from __future__ import annotations

import ast
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "pipeline"
STATUS_DOC = ROOT / "docs" / "implementation-status.md"
GAUSSIAN = PIPELINE / "gaussian.py"

# 1) Placeholder science: real API results / coordinates / energies must never
#    be replaced by hardcoded, illustrative, random, or mocked values in a
#    non-test code path. (tests/ is exempt; fixtures are fine there.)
PLACEHOLDER_PATTERNS = [
    r"\bPLACEHOLDER\b",
    r"\bDUMMY_DATA\b",
    r"\bFIXME_SCIENCE\b",
    r"illustrative[_ ]only",
    r"#\s*fake\b",
    r"np\.random\.",
    r"numpy\.random\.",
]

# 2) Route-line integrity — enforced mechanically ONLY where route lines live in
#    code (constants named ROUTE_OPT / ROUTE_FREQ in the package). Today they are
#    passed in from user config, so this is a no-op that activates automatically
#    if/when route construction moves into pipeline/. Until then, route-line
#    integrity is a REVIEW guideline (AGENTS.md §2/§6), not a CI gate.
REQUIRED_OPT_TOKENS = ["b3lyp", "scrf", "6-311"]
REQUIRED_FREQ_TOKENS = ["b3lyp", "scrf", "6-311", "temperature"]


def check_placeholders() -> list[str]:
    problems = []
    for py in PIPELINE.rglob("*.py"):
        text = py.read_text(encoding="utf-8", errors="replace")
        # ignore docstring example lines (kept honest by review, not CI)
        for pat in PLACEHOLDER_PATTERNS:
            for m in re.finditer(pat, text):
                line = text[: m.start()].count("\n") + 1
                problems.append(
                    f"[placeholder] {py.relative_to(ROOT)}:{line} matches /{pat}/"
                )
    return problems


def check_route_constants() -> list[str]:
    problems = []
    for py in PIPELINE.rglob("*.py"):
        for i, raw in enumerate(py.read_text(errors="replace").splitlines(), 1):
            m = re.match(r"\s*(ROUTE_OPT|ROUTE_FREQ)\s*=\s*(.+)", raw)
            if not m:
                continue
            name, value = m.group(1), m.group(2).lower()
            required = REQUIRED_FREQ_TOKENS if name == "ROUTE_FREQ" else REQUIRED_OPT_TOKENS
            missing = [t for t in required if t not in value]
            if missing:
                problems.append(
                    f"[route-const] {py.relative_to(ROOT)}:{i} {name} missing "
                    f"physics token(s): {', '.join(missing)}"
                )
    return problems


# 3) Canonical status-doc drift guard (M-07). AGENTS.md §5 names
#    docs/implementation-status.md as the merge gate. It must never sit as the
#    empty template while the real status lives in implementation-status-v2.md, so
#    reviewers/maintainers following the mandated gate can't miss v2 deviations.
#    This turns "canonical went stale" from a reviewer catch into a floor failure.
#    (The v2 working-history drafts are archived under docs/review-history/v2/.)
def _status_doc_problems(text: str) -> list[str]:
    """Return template-marker problems in a status-doc string (pure, testable)."""
    problems = []
    # Angle-bracket template placeholders: <n>, <name>, <none | ...>, <>, <e.g...>.
    # The populated status doc must contain no bare angle-bracket tokens.
    for m in re.finditer(r"<[^>\n]*>", text):
        problems.append(f"template placeholder {m.group(0)!r}")
    # Empty bullets: a list line that is just "-" / "- " with nothing after it.
    for i, line in enumerate(text.splitlines(), 1):
        if re.match(r"^\s*-\s*$", line):
            problems.append(f"empty bullet at line {i}")
    return problems


def check_status_doc() -> list[str]:
    if not STATUS_DOC.exists():
        return [f"[status-doc] {STATUS_DOC.name} is missing"]
    text = STATUS_DOC.read_text(encoding="utf-8", errors="replace")
    return [
        f"[status-doc] docs/implementation-status.md still has {p} "
        f"(populate the canonical merge-gate doc — AGENTS.md §5 gate)"
        for p in _status_doc_problems(text)
    ]


# 4) V2 Gaussian provenance threading (M-14). Every conformer-derived COM must
#    carry pipeline/commit/RDKit identity from its conformer-log row. Use the AST
#    so this guard checks semantic structure without depending on formatting.
_GAUSSIAN_PROVENANCE_FIELDS = (
    "pipeline_version",
    "pipeline_commit",
    "rdkit_version",
)


def _function_runtime_literals(function: ast.FunctionDef) -> str:
    """Return string literals from a function body, excluding its docstring."""
    body = function.body
    if body and isinstance(body[0], ast.Expr):
        value = body[0].value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            body = body[1:]
    literals = []
    for statement in body:
        for node in ast.walk(statement):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                literals.append(node.value)
    return "".join(literals)


def _gaussian_provenance_problems(text: str) -> list[str]:
    """Return missing M-14 writer/threading semantics in Gaussian source text."""
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        return [f"cannot parse pipeline/gaussian.py: {exc}"]

    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    writer = functions.get("write_gaussian_com")
    batch = functions.get("write_gaussian_coms_from_conformers")
    problems = []
    if writer is None:
        problems.append("write_gaussian_com is missing")
    if batch is None:
        problems.append("write_gaussian_coms_from_conformers is missing")
    if writer is None or batch is None:
        return problems

    writer_args = {
        arg.arg
        for arg in (
            writer.args.posonlyargs + writer.args.args + writer.args.kwonlyargs
        )
    }
    for field in _GAUSSIAN_PROVENANCE_FIELDS:
        if field not in writer_args:
            problems.append(f"write_gaussian_com missing parameter {field}")

    literals = _function_runtime_literals(writer)
    for token in ("provenance ", "pipeline=", "commit=", "rdkit="):
        if token not in literals:
            problems.append(f"write_gaussian_com title missing token {token!r}")

    row_reads = set()
    forwarded = set()
    for node in ast.walk(batch):
        if not isinstance(node, ast.Call):
            continue
        if (
            isinstance(node.func, ast.Attribute)
            and node.func.attr == "get"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            row_reads.add(node.args[0].value)
        if isinstance(node.func, ast.Name) and node.func.id == "write_gaussian_com":
            forwarded.update(keyword.arg for keyword in node.keywords if keyword.arg)

    for field in _GAUSSIAN_PROVENANCE_FIELDS:
        if field not in row_reads:
            problems.append(
                f"write_gaussian_coms_from_conformers does not read {field}"
            )
        if field not in forwarded:
            problems.append(
                f"write_gaussian_coms_from_conformers does not forward {field}"
            )
    return problems


def check_gaussian_provenance() -> list[str]:
    if not GAUSSIAN.exists():
        return ["[gaussian-provenance] pipeline/gaussian.py is missing"]
    text = GAUSSIAN.read_text(encoding="utf-8", errors="replace")
    return [
        f"[gaussian-provenance] {problem}"
        for problem in _gaussian_provenance_problems(text)
    ]


def main() -> int:
    problems = (
        check_placeholders()
        + check_route_constants()
        + check_status_doc()
        + check_gaussian_provenance()
    )
    if problems:
        print("SCIENTIFIC INVARIANT CHECK FAILED:\n")
        for p in problems:
            print("  " + p)
        print(f"\n{len(problems)} violation(s).")
        return 1
    print("Scientific invariant checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
