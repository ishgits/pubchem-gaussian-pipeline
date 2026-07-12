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

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "pipeline"

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


def main() -> int:
    problems = check_placeholders() + check_route_constants()
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
