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
CONFORMERS = PIPELINE / "conformers.py"
MANIFEST = PIPELINE / "manifest.py"
SLURM = PIPELINE / "slurm.py"

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


# 4) V2 Gaussian provenance enforcement (M-14/M-16/M-17). Every
#    conformer-derived COM must carry pipeline/commit/RDKit identity, whether it
#    comes through the batch API or a direct conformer-specific writer call. Use
#    the AST so this guard checks semantic structure without depending on
#    formatting.
_GAUSSIAN_PROVENANCE_FIELDS = (
    "pipeline_version",
    "pipeline_commit",
    "rdkit_version",
)
_GAUSSIAN_LINKAGE_FIELDS = ("run_id", "artifact_id", "config_hash")


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
    """Return missing M-14/M-16/M-17 provenance semantics in Gaussian source."""
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        return [f"cannot parse pipeline/gaussian.py: {exc}"]

    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    writer = functions.get("write_gaussian_com")
    batch = functions.get("write_gaussian_coms_from_conformers")
    validator = functions.get("_validate_required_conformer_provenance")
    direct_validator = functions.get("_validate_direct_conformer_provenance")
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
    for field in _GAUSSIAN_PROVENANCE_FIELDS + _GAUSSIAN_LINKAGE_FIELDS:
        if field not in writer_args:
            problems.append(f"write_gaussian_com missing parameter {field}")

    literals = _function_runtime_literals(writer)
    for token in ("provenance ", "pipeline=", "commit=", "rdkit="):
        if token not in literals:
            problems.append(f"write_gaussian_com title missing token {token!r}")

    # M-17: a direct call with `conformer_id` is still a v2 scientific output and
    # must not bypass the protected batch boundary. The helper must condition on
    # that identifier, reject blank pipeline/RDKit values, and run before the
    # writer creates its output directory or opens a file.
    if direct_validator is None:
        problems.append("_validate_direct_conformer_provenance is missing")
    else:
        direct_args = {
            arg.arg
            for arg in (
                direct_validator.args.posonlyargs
                + direct_validator.args.args
                + direct_validator.args.kwonlyargs
            )
        }
        for field in (
            "conformer_id", "pipeline_version", "rdkit_version",
            "run_id", "artifact_id", "config_hash", "manifest_path",
            "parent_artifact_id", "conformer_record_id",
        ):
            if field not in direct_args:
                problems.append(
                    f"direct conformer provenance validator missing parameter {field}"
                )

        conditional_on_conformer = any(
            isinstance(node, ast.If)
            and any(
                isinstance(child, ast.Name) and child.id == "conformer_id"
                for child in ast.walk(node.test)
            )
            for node in ast.walk(direct_validator)
        )
        if not conditional_on_conformer:
            problems.append(
                "direct conformer provenance validation is not conditional on "
                "conformer_id"
            )

        normalized_names = set()
        for node in ast.walk(direct_validator):
            if not (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "_optional_text"
                and node.args
                and isinstance(node.args[0], ast.Name)
            ):
                continue
            normalized_names.add(node.args[0].id)
        for field in (
            "pipeline_version", "rdkit_version", "run_id", "artifact_id",
            "config_hash", "manifest_path", "parent_artifact_id",
            "conformer_record_id",
        ):
            if field not in normalized_names:
                problems.append(
                    f"direct conformer provenance does not require nonblank {field}"
                )
        if not any(isinstance(node, ast.Raise) for node in ast.walk(direct_validator)):
            problems.append("direct conformer provenance validator never raises")

    direct_validation_lines = []
    direct_mutation_lines = []
    for node in ast.walk(writer):
        if not isinstance(node, ast.Call):
            continue
        if (
            isinstance(node.func, ast.Name)
            and node.func.id == "_validate_direct_conformer_provenance"
        ):
            direct_validation_lines.append(node.lineno)
            forwarded_names = {
                arg.id for arg in node.args if isinstance(arg, ast.Name)
            }
            for field in (
                "conformer_id", "pipeline_version", "rdkit_version",
                "run_id", "artifact_id", "config_hash", "manifest_path",
                "parent_artifact_id", "conformer_record_id",
            ):
                if field not in forwarded_names:
                    problems.append(
                        f"write_gaussian_com does not pass {field} to direct "
                        "provenance validation"
                    )
        if isinstance(node.func, ast.Name) and node.func.id in {"ensure_dir", "open"}:
            direct_mutation_lines.append(node.lineno)
    if not direct_validation_lines:
        problems.append("write_gaussian_com does not validate direct conformer provenance")
    elif direct_mutation_lines and min(direct_validation_lines) >= min(direct_mutation_lines):
        problems.append("direct conformer provenance validation occurs after mutation")

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

    # M-16: the downstream writer must reject nonempty external logs that lack
    # the source versions. Checking this here prevents a future refactor from
    # preserving title-token threading while accidentally removing the boundary
    # validation that makes those tokens trustworthy.
    required_columns = set()
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name)
            and target.id == "_REQUIRED_CONFORMER_PROVENANCE_COLUMNS"
            for target in node.targets
        ):
            continue
        if isinstance(node.value, (ast.Tuple, ast.List)):
            required_columns = {
                elt.value
                for elt in node.value.elts
                if isinstance(elt, ast.Constant) and isinstance(elt.value, str)
            }
    for field in (
        "run_id", "artifact_id", "config_hash", "pipeline_version",
        "rdkit_version", "xyz_sha256",
    ):
        if field not in required_columns:
            problems.append(f"required conformer provenance omits {field}")

    if validator is None:
        problems.append("_validate_required_conformer_provenance is missing")
    else:
        validator_calls = {
            node.func.id
            for node in ast.walk(validator)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        if "_optional_text" not in validator_calls:
            problems.append(
                "required conformer provenance does not reject blank/NaN values"
            )
        if not any(isinstance(node, ast.Raise) for node in ast.walk(validator)):
            problems.append("required conformer provenance validator never raises")

    read_lines = []
    validation_lines = []
    mutation_lines = []
    for node in ast.walk(batch):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute) and node.func.attr == "read_csv":
            read_lines.append(node.lineno)
        if (
            isinstance(node.func, ast.Name)
            and node.func.id == "_validate_required_conformer_provenance"
        ):
            validation_lines.append(node.lineno)
        if (
            isinstance(node.func, ast.Name)
            and node.func.id == "write_gaussian_com"
        ) or (
            isinstance(node.func, ast.Attribute)
            and node.func.attr in {"remove", "to_csv"}
        ):
            mutation_lines.append(node.lineno)
    if not validation_lines:
        problems.append(
            "write_gaussian_coms_from_conformers does not validate required provenance"
        )
    elif not read_lines or min(validation_lines) <= min(read_lines):
        problems.append("required conformer provenance validation is not after CSV read")
    elif mutation_lines and min(validation_lines) >= min(mutation_lines):
        problems.append("required conformer provenance validation occurs after mutation")
    return problems


def check_gaussian_provenance() -> list[str]:
    if not GAUSSIAN.exists():
        return ["[gaussian-provenance] pipeline/gaussian.py is missing"]
    text = GAUSSIAN.read_text(encoding="utf-8", errors="replace")
    return [
        f"[gaussian-provenance] {problem}"
        for problem in _gaussian_provenance_problems(text)
    ]


# 5) Conformer reuse/output provenance (M-15/M-19/M-20). Unrequested groups
# cannot be regenerated from the current molecule table, so append mode must
# validate and reject them before any output mutation rather than copying them
# untouched. Generated XYZ comments must identify RDKit, and cached clean commit
# provenance must match while dirty worktrees are never reusable.
def _append_integrity_problems(text: str) -> list[str]:
    """Return missing M-15 append-boundary semantics in conformer source text."""
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        return [f"cannot parse pipeline/conformers.py: {exc}"]

    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    carry = functions.get("_carry_forward_group_is_valid")
    config_match = functions.get("_row_config_matches")
    partition = functions.get("_resume_partition")
    search = functions.get("search_conformers")
    problems = []
    for name, function in (
        ("_carry_forward_group_is_valid", carry),
        ("_row_config_matches", config_match),
        ("_resume_partition", partition),
        ("search_conformers", search),
    ):
        if function is None:
            problems.append(f"{name} is missing")
    if carry is None or config_match is None or partition is None or search is None:
        return problems

    config_literals = _function_runtime_literals(config_match)
    commit_getters = {
        node.func.value.id
        for node in ast.walk(config_match)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "get"
        and isinstance(node.func.value, ast.Name)
        and node.args
        and isinstance(node.args[0], ast.Constant)
        and node.args[0].value == "pipeline_commit"
    }
    if not {"row", "config"}.issubset(commit_getters):
        problems.append("resume config matching omits pipeline_commit")
    blank_commit_guards = {
        node.operand.id
        for node in ast.walk(config_match)
        if isinstance(node, ast.UnaryOp)
        and isinstance(node.op, ast.Not)
        and isinstance(node.operand, ast.Name)
        and node.operand.id in {"row_commit", "config_commit"}
    }
    if blank_commit_guards != {"row_commit", "config_commit"}:
        problems.append("resume config matching permits a missing commit")
    if ".dirty" not in config_literals:
        problems.append("resume config matching omits dirty-commit rejection")

    carry_calls = {
        node.func.id
        for node in ast.walk(carry)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    for required in (
        "_resume_group_is_complete",
        "_group_identity_is_consistent",
        "_row_config_matches",
    ):
        if required not in carry_calls:
            problems.append(f"carry-forward validation omits {required}")
    if "pipeline_commit" not in _function_runtime_literals(carry):
        problems.append("carry-forward validation omits pipeline_commit field presence")

    search_literals = _function_runtime_literals(search)
    if "rdkit_version=" not in search_literals:
        problems.append("generated XYZ provenance omits rdkit_version= token")

    run_config_keys = set()
    for node in ast.walk(search):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Dict):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "run_config"
            for target in node.targets
        ):
            continue
        run_config_keys.update(
            key.value
            for key in node.value.keys
            if isinstance(key, ast.Constant) and isinstance(key.value, str)
        )
    if "pipeline_commit" not in run_config_keys:
        problems.append("search_conformers run_config omits pipeline_commit")

    partition_calls = {
        node.func.id
        for node in ast.walk(partition)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    partition_names = {
        node.id for node in ast.walk(partition) if isinstance(node, ast.Name)
    }
    if "_carry_forward_group_is_valid" not in partition_calls:
        problems.append("_resume_partition does not validate carry-forward groups")
    if "preserve_unrequested" not in partition_names:
        problems.append("_resume_partition lacks append carry-forward gating")
    if "invalid_retained" not in partition_names:
        problems.append("_resume_partition does not report invalid retained groups")

    partition_lines = []
    mutation_lines = []
    for node in ast.walk(search):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Name) and node.func.id == "_resume_partition":
            partition_lines.append(node.lineno)
        if isinstance(node.func, ast.Name) and node.func.id == "ensure_dir":
            mutation_lines.append(node.lineno)
        if isinstance(node.func, ast.Attribute) and node.func.attr in {"remove", "to_csv"}:
            mutation_lines.append(node.lineno)

    guard_lines = []
    for node in ast.walk(search):
        if not isinstance(node, ast.If):
            continue
        test_names = {
            child.id for child in ast.walk(node.test) if isinstance(child, ast.Name)
        }
        if "invalid_retained" in test_names and any(
            isinstance(child, ast.Raise) for child in ast.walk(node)
        ):
            guard_lines.append(node.lineno)
    if not partition_lines:
        problems.append("search_conformers does not partition the existing log")
    if not guard_lines:
        problems.append("search_conformers does not raise for invalid retained groups")
    elif mutation_lines and min(guard_lines) >= min(mutation_lines):
        problems.append("invalid retained groups are checked after output mutation")
    elif partition_lines and min(partition_lines) >= min(guard_lines):
        problems.append("invalid retained groups are checked before partitioning")
    return problems


def check_append_integrity() -> list[str]:
    if not CONFORMERS.exists():
        return ["[append-integrity] pipeline/conformers.py is missing"]
    text = CONFORMERS.read_text(encoding="utf-8", errors="replace")
    return [
        f"[append-integrity] {problem}"
        for problem in _append_integrity_problems(text)
    ]


# 6) Frozen manifest-centric v2.0 matrix. This intentionally checks only the
# approved per-artifact linkage fields and authoritative manifest mechanisms;
# it must not grow into a requirement to duplicate every search knob into XYZ
# comments or COM title blocks.
def _frozen_matrix_problems(
    manifest_text: str,
    conformer_text: str,
    gaussian_text: str,
    slurm_text: str,
) -> list[str]:
    problems = []
    try:
        conformer_tree = ast.parse(conformer_text)
        gaussian_tree = ast.parse(gaussian_text)
        slurm_tree = ast.parse(slurm_text)
        conformer_functions = {
            node.name: node for node in conformer_tree.body
            if isinstance(node, ast.FunctionDef)
        }
        gaussian_functions = {
            node.name: node for node in gaussian_tree.body
            if isinstance(node, ast.FunctionDef)
        }
        slurm_functions = {
            node.name: node for node in slurm_tree.body
            if isinstance(node, ast.FunctionDef)
        }
        xyz_runtime = _function_runtime_literals(conformer_functions["search_conformers"])
        com_runtime = _function_runtime_literals(gaussian_functions["write_gaussian_com"])
        sh_runtime = _function_runtime_literals(slurm_functions["write_slurm_script"])
    except (KeyError, SyntaxError) as exc:
        return [f"cannot parse frozen stage writers: {exc}"]
    manifest_tokens = (
        'MANIFEST_SCHEMA = "2.0"',
        "def canonical_json(",
        "def configuration_hash(",
        "def stable_record_id(",
        "def sha256_file(",
        "def require_exact_artifact_id_set(",
        "def _require_link1_checkpoint_reads(",
        "def validate_manifest(",
        "def record_conformer_xyz(",
        "def record_child_artifact(",
        '"Geom=AllChk"',
        '"Guess=Read"',
        '"molecule_identity_hash"',
        '"conformer_record_id"',
        '"parent_artifact_id"',
        '"relative_path"',
        '"sha256"',
    )
    for token in manifest_tokens:
        if token not in manifest_text:
            problems.append(f"manifest contract omits {token!r}")

    xyz_fields = (
        "run_id=", "artifact_id=", "config_hash=", "conformer_id=",
        "relative_energy_kcalmol=", "method=", "pipeline_version=",
        "rdkit_version=",
    )
    for token in xyz_fields:
        if token not in xyz_runtime:
            problems.append(f"XYZ linkage omits {token!r}")

    com_fields = (
        "run_id=", "artifact_id=", "config_hash=", "conformer_id=",
        "relative_energy_kcalmol=", "pipeline_version=", "rdkit_version=",
    )
    for token in com_fields:
        if token not in com_runtime:
            problems.append(f"COM linkage omits {token!r}")
    for token in ("manifest_path", "record_child_artifact", "sha256_file"):
        if token not in gaussian_text:
            problems.append(f"Gaussian manifest boundary omits {token!r}")
    if "require_exact_artifact_id_set(" not in gaussian_text:
        problems.append("Gaussian boundary omits exact manifest XYZ-set validation")

    sh_fields = (
        "# run_id=", "# artifact_id=", "# source_com_relative_path=",
        "# source_com_sha256=",
    )
    for token in sh_fields:
        if token not in sh_runtime:
            problems.append(f"SLURM linkage omits {token!r}")
    for token in (
        "zero-byte com_path", "duplicate normalized com_path",
        "collapse to script basename", "record_child_artifact",
        "require_exact_artifact_id_set(",
    ):
        if token not in slurm_text:
            problems.append(f"SLURM one-to-one validation omits {token!r}")
    return problems


def check_frozen_matrix() -> list[str]:
    paths = (MANIFEST, CONFORMERS, GAUSSIAN, SLURM)
    missing = [path for path in paths if not path.exists()]
    if missing:
        return [f"[frozen-v2-matrix] missing {path.relative_to(ROOT)}" for path in missing]
    problems = _frozen_matrix_problems(
        MANIFEST.read_text(encoding="utf-8", errors="replace"),
        CONFORMERS.read_text(encoding="utf-8", errors="replace"),
        GAUSSIAN.read_text(encoding="utf-8", errors="replace"),
        SLURM.read_text(encoding="utf-8", errors="replace"),
    )
    return [f"[frozen-v2-matrix] {problem}" for problem in problems]


def main() -> int:
    problems = (
        check_placeholders()
        + check_route_constants()
        + check_status_doc()
        + check_gaussian_provenance()
        + check_append_integrity()
        + check_frozen_matrix()
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
