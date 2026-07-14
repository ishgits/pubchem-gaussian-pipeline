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


class TestGeneratedArtifactIgnoreGuard:
    def test_current_gitignore_passes(self):
        text = (SCRIPT.parents[1] / ".gitignore").read_text()
        assert check_invariants._generated_artifact_ignore_problems(text) == []

    def test_detects_unignored_manifest(self):
        problems = check_invariants._generated_artifact_ignore_problems(
            "conformer_xyz/\ngaussian_inputs/\n"
        )
        assert any("run_manifest.json" in problem for problem in problems)


class TestGaussianProvenanceGuard:
    """M-14/M-16/M-17: enforce provenance across batch and direct v2 paths."""

    @staticmethod
    def _source():
        return (SCRIPT.parents[1] / "pipeline" / "gaussian.py").read_text()

    def test_current_gaussian_source_passes(self):
        assert check_invariants._gaussian_provenance_problems(self._source()) == []

    def test_detects_missing_writer_parameter(self):
        broken = self._source().replace(
            "    pipeline_version: str | None = None,\n", "", 1
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("missing parameter pipeline_version" in p for p in problems)

    def test_detects_missing_row_read(self):
        broken = self._source().replace(
            'row.get("rdkit_version")', 'row.get("other_version")', 1
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("does not read rdkit_version" in p for p in problems)

    def test_detects_missing_forwarding(self):
        broken = self._source().replace(
            "                pipeline_commit=pipeline_commit,\n", "", 1
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("does not forward pipeline_commit" in p for p in problems)

    def test_detects_missing_required_source_version(self):
        broken = self._source().replace(
            '    "rdkit_version",\n    "xyz_sha256",',
            '    "xyz_sha256",',
            1,
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("required conformer provenance omits rdkit_version" in p for p in problems)

    def test_detects_missing_batch_boundary_validation(self):
        broken = self._source().replace(
            "    _validate_required_conformer_provenance(conf_log)\n", "", 1
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("does not validate required provenance" in p for p in problems)

    def test_detects_missing_direct_boundary_validation(self):
        broken = self._source().replace(
            "    _validate_direct_conformer_provenance(\n",
            "    _ignored_direct_conformer_provenance(\n",
            1,
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("does not validate direct conformer provenance" in p for p in problems)

    def test_detects_direct_validation_without_conformer_condition(self):
        broken = self._source().replace(
            "    if conformer_id is None:\n", "    if False:\n", 1
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("not conditional on conformer_id" in p for p in problems)

    def test_detects_missing_direct_rdkit_requirement(self):
        broken = self._source().replace(
            "    if _optional_text(rdkit_version) is None:\n",
            "    if False:\n",
            1,
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("does not require nonblank rdkit_version" in p for p in problems)

    def test_detects_direct_validation_after_mutation(self):
        marker = "    # M-17: `conformer_id` selects the v2 scientific-output path"
        broken = self._source().replace(
            marker,
            "    ensure_dir(outdir)\n\n" + marker,
            1,
        )
        problems = check_invariants._gaussian_provenance_problems(broken)
        assert any("direct conformer provenance validation occurs after mutation" in p for p in problems)


class TestNoResumeGuard:
    """v2.1 (contract §7): the resume/append subsystem is removed; the guard
    proves the path is GONE, not merely refactored, and holds the reduced XYZ
    provenance tokens."""

    @staticmethod
    def _source():
        return (SCRIPT.parents[1] / "pipeline" / "conformers.py").read_text()

    def test_current_conformer_source_passes(self):
        assert check_invariants._no_resume_problems(self._source()) == []

    def test_detects_reintroduced_resume_helper(self):
        broken = self._source().replace(
            "def search_conformers(",
            "def _resume_partition(existing):\n    return None\n\n\ndef search_conformers(",
            1,
        )
        problems = check_invariants._no_resume_problems(broken)
        assert any("_resume_partition was not removed" in p for p in problems)

    def test_detects_reintroduced_append_parameter(self):
        broken = self._source().replace(
            "    seed: int = SEED,\n    manifest_path: str = \"run_manifest.json\",",
            "    seed: int = SEED,\n    append: bool = False,\n    manifest_path: str = \"run_manifest.json\",",
            1,
        )
        problems = check_invariants._no_resume_problems(broken)
        assert any("still exposes an append parameter" in p for p in problems)

    def test_detects_missing_populated_run_guard(self):
        broken = self._source().replace("already populated", "some-other-text")
        problems = check_invariants._no_resume_problems(broken)
        assert any("already-populated run folder" in p for p in problems)

    def test_detects_missing_xyz_metadata_token(self):
        broken = self._source().replace("method={method}", "ff={method}", 1)
        problems = check_invariants._no_resume_problems(broken)
        assert any("omits 'method='" in p for p in problems)


class TestFrozenManifestMatrixGuard:
    @staticmethod
    def _sources():
        root = SCRIPT.parents[1] / "pipeline"
        return tuple(
            (root / name).read_text()
            for name in ("manifest.py", "conformers.py", "gaussian.py", "slurm.py")
        )

    def test_current_sources_pass(self):
        assert check_invariants._frozen_matrix_problems(*self._sources()) == []

    def test_detects_missing_atomic_conformer_group_writer(self):
        manifest, conformers, gaussian, slurm = self._sources()
        manifest = manifest.replace(
            "def record_conformer_group(",
            "def removed_record_conformer_group(",
            1,
        )
        problems = check_invariants._frozen_matrix_problems(
            manifest, conformers, gaussian, slurm
        )
        assert any("record_conformer_group" in problem for problem in problems)

    def test_detects_missing_xyz_linkage_field(self):
        manifest, conformers, gaussian, slurm = self._sources()
        conformers = conformers.replace("method={method}", "ff={method}", 1)
        problems = check_invariants._frozen_matrix_problems(
            manifest, conformers, gaussian, slurm
        )
        assert any("XYZ linkage omits" in problem for problem in problems)

    def test_detects_missing_zero_byte_guard(self):
        manifest, conformers, gaussian, slurm = self._sources()
        slurm = slurm.replace("zero-byte com_path", "empty input", 1)
        problems = check_invariants._frozen_matrix_problems(
            manifest, conformers, gaussian, slurm
        )
        assert any("zero-byte com_path" in problem for problem in problems)

    def test_detects_missing_exact_xyz_set_guard(self):
        manifest, conformers, gaussian, slurm = self._sources()
        gaussian = gaussian.replace(
            "require_exact_artifact_id_set(", "removed_exact_set_guard(", 1
        )
        problems = check_invariants._frozen_matrix_problems(
            manifest, conformers, gaussian, slurm
        )
        assert any("exact manifest XYZ-set" in problem for problem in problems)

    def test_detects_missing_link1_checkpoint_guard(self):
        manifest, conformers, gaussian, slurm = self._sources()
        manifest = manifest.replace(
            "def _require_link1_checkpoint_reads(",
            "def removed_link1_checkpoint_guard(",
            1,
        )
        problems = check_invariants._frozen_matrix_problems(
            manifest, conformers, gaussian, slurm
        )
        assert any("_require_link1_checkpoint_reads" in problem for problem in problems)

    def test_detects_missing_manifest_basename_collision_guard(self):
        manifest, conformers, gaussian, slurm = self._sources()
        manifest = manifest.replace("sanitize_basename(name)", "str(name)")
        problems = check_invariants._frozen_matrix_problems(
            manifest, conformers, gaussian, slurm
        )
        assert any("sanitize_basename" in problem for problem in problems)
