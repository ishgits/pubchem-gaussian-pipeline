"""Tests for pipeline.pubchem scoring heuristic (offline — no network calls)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.pubchem import score_candidate


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
