"""Tests for pipeline.pubchem scoring heuristic (offline — no network calls)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.pubchem import (
    MOLECULE_TABLE_COLUMNS,
    _isomeric_smiles,
    _resolved_row,
    score_candidate,
)


def _mock_prop(cid=100, formula="C5H5N5", iso_smiles="c1ncc2[nH]cnc2n1",
               iupac="7H-purin-6-amine", title="adenine"):
    return {
        "CID": cid,
        "MolecularFormula": formula,
        "IsomericSMILES": iso_smiles,
        "IUPACName": iupac,
        "Title": title,
    }


class TestScoreCandidate:
    def test_formula_match_boosts_score(self):
        prop = _mock_prop(formula="C5H5N5")
        score_match = score_candidate(prop, expected_formula="C5H5N5")
        score_no_match = score_candidate(prop, expected_formula="C6H6N6")
        assert score_match > score_no_match

    def test_stereo_bonus(self):
        prop_stereo = _mock_prop(iso_smiles="C[C@@H](O)N")
        prop_flat = _mock_prop(iso_smiles="CC(O)N")
        assert score_candidate(prop_stereo) > score_candidate(prop_flat)

    def test_low_cid_bonus(self):
        prop_low = _mock_prop(cid=500)
        prop_high = _mock_prop(cid=500_000)
        assert score_candidate(prop_low) > score_candidate(prop_high)

    def test_keyword_boost(self):
        prop = _mock_prop(iupac="adenosine", title="adenosine")
        score_with = score_candidate(prop, keyword_boost=["adenosine"])
        score_without = score_candidate(prop, keyword_boost=[])
        assert score_with > score_without

    def test_no_formula_no_crash(self):
        prop = _mock_prop()
        # Should work fine with no expected formula
        score = score_candidate(prop)
        assert isinstance(score, int)


class TestIsomericSmiles:
    """B-01: read the stereo-bearing SMILES across PubChem's 2025 key rename."""

    def test_current_smiles_key_preferred(self):
        # Post-2025 PubChem: stereo SMILES arrives under "SMILES".
        prop = {"SMILES": "C[C@@H](N)C(=O)O", "ConnectivitySMILES": "CC(N)C(=O)O"}
        assert _isomeric_smiles(prop) == "C[C@@H](N)C(=O)O"

    def test_legacy_isomericsmiles_fallback(self):
        # Old caches used "IsomericSMILES".
        prop = {"IsomericSMILES": "C[C@@H](N)C(=O)O"}
        assert _isomeric_smiles(prop) == "C[C@@H](N)C(=O)O"

    def test_never_uses_connectivity_smiles(self):
        # ConnectivitySMILES drops stereo — must not be returned as a fallback.
        prop = {"ConnectivitySMILES": "CC(N)C(=O)O"}
        assert _isomeric_smiles(prop) == ""

    def test_missing_returns_empty(self):
        assert _isomeric_smiles({}) == ""


class TestResolvedRow:
    """B-01: IsomericSMILES must be carried into the resolved molecule table."""

    _info = {"status": "OK", "warnings": []}

    def test_resolved_row_contains_isomeric_smiles(self):
        prop = {"CID": 5950, "SMILES": "C[C@@H](N)C(=O)O", "MolecularFormula": "C3H7NO2"}
        row = _resolved_row("Alanine", "alanine", prop, self._info)
        assert "IsomericSMILES" in row
        assert row["IsomericSMILES"] == "C[C@@H](N)C(=O)O"

    def test_resolved_row_uses_stereo_not_connectivity(self):
        # ConnectivitySMILES (stereo-free) must NOT be substituted.
        prop = {"CID": 5950, "SMILES": "C[C@@H](N)C(=O)O", "ConnectivitySMILES": "CC(N)C(=O)O"}
        row = _resolved_row("Alanine", "alanine", prop, self._info)
        assert row["IsomericSMILES"] == "C[C@@H](N)C(=O)O"

    def test_unresolved_row_has_empty_smiles_column(self):
        # Column present (empty) even when resolution failed → stable schema.
        info = {"status": "NO_CIDS", "warnings": ["no hits"]}
        row = _resolved_row("Nope", "nope", None, info)
        assert row["IsomericSMILES"] == ""
        assert row["cid"] is None

    def test_schema_matches_declared_columns(self):
        prop = _mock_prop()
        row = _resolved_row("Adenine", "adenine", prop, self._info)
        assert set(row.keys()) == set(MOLECULE_TABLE_COLUMNS)

    def test_downstream_consumer_columns_present(self):
        # Regression: download_sdfs needs name+cid; conformers needs IsomericSMILES.
        prop = _mock_prop()
        row = _resolved_row("Adenine", "adenine", prop, self._info)
        for col in ("name", "cid", "IsomericSMILES"):
            assert col in row
