"""Tests for pipeline.pubchem scoring heuristic (offline — no network calls)."""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.pubchem import (
    MOLECULE_TABLE_COLUMNS,
    _cache_path,
    _cache_request,
    _get_json,
    _isomeric_smiles,
    _resolved_row,
    build_molecule_table,
    score_candidate,
)
import pipeline.pubchem as pubchem_module
import pytest


class _Response:
    status_code = 200
    text = ""

    def __init__(self, url):
        self.url = url

    def json(self):
        return {"requested_url": self.url}


class TestRequestBoundCache:
    @pytest.mark.parametrize(
        "key_a,url_a,key_b,url_b",
        [
            ("cids__A B", "https://example.test/A%20B", "cids__A_B", "https://example.test/A_B"),
            ("cids__A/B", "https://example.test/A%2FB", "cids__A_B", "https://example.test/A_B"),
            (
                "cids__" + "x" * 200 + "A",
                "https://example.test/" + "x" * 200 + "A",
                "cids__" + "x" * 200 + "B",
                "https://example.test/" + "x" * 200 + "B",
            ),
        ],
    )
    def test_lossy_readable_keys_cannot_share_entries(
        self, tmp_path, monkeypatch, key_a, url_a, key_b, url_b
    ):
        calls = []

        def get(url, **_kwargs):
            calls.append(url)
            return _Response(url)

        monkeypatch.setattr(pubchem_module.requests, "get", get)
        kwargs = dict(cache_dir=str(tmp_path), max_retries=1, min_delay_s=0)

        first_a = _get_json(url_a, cache_key=key_a, **kwargs)
        first_b = _get_json(url_b, cache_key=key_b, **kwargs)
        assert _cache_path(str(tmp_path), key_a, url_a) != _cache_path(
            str(tmp_path), key_b, url_b
        )
        assert first_a != first_b
        assert calls == [url_a, url_b]

        assert _get_json(url_a, cache_key=key_a, **kwargs) == first_a
        assert _get_json(url_b, cache_key=key_b, **kwargs) == first_b
        assert calls == [url_a, url_b]

    def test_tampered_request_envelope_is_a_cache_miss(self, tmp_path, monkeypatch):
        url = "https://example.test/requested"
        key = "cids__requested"
        path = _cache_path(str(tmp_path), key, url)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "cache_schema": 2,
                    "request": _cache_request("https://example.test/different"),
                    "response": {"wrong": True},
                },
                handle,
            )

        calls = []

        def get(requested_url, **_kwargs):
            calls.append(requested_url)
            return _Response(requested_url)

        monkeypatch.setattr(pubchem_module.requests, "get", get)
        result = _get_json(
            url,
            cache_key=key,
            cache_dir=str(tmp_path),
            max_retries=1,
            min_delay_s=0,
        )
        assert result == {"requested_url": url}
        assert calls == [url]
        with open(path, encoding="utf-8") as handle:
            repaired = json.load(handle)
        assert repaired["request"] == _cache_request(url)

    def test_unverifiable_legacy_payload_is_a_cache_miss(self, tmp_path, monkeypatch):
        url = "https://example.test/legacy"
        key = "cids__legacy"
        path = _cache_path(str(tmp_path), key, url)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump({"IdentifierList": {"CID": [1]}}, handle)
        monkeypatch.setattr(
            pubchem_module.requests, "get", lambda requested_url, **_kwargs: _Response(requested_url)
        )
        assert _get_json(
            url,
            cache_key=key,
            cache_dir=str(tmp_path),
            max_retries=1,
            min_delay_s=0,
        ) == {"requested_url": url}


class TestFallbackQueryProvenance:
    def test_table_records_the_query_that_succeeded(self, monkeypatch):
        selected = _mock_prop(cid=321)
        calls = []

        def resolve(label, query, **_kwargs):
            calls.append(query)
            if query == "primary":
                return None, {
                    "label": label,
                    "query": query,
                    "status": "NO_CIDS",
                    "warnings": [],
                }
            return selected, {
                "label": label,
                "query": query,
                "status": "OK",
                "warnings": [],
            }

        monkeypatch.setattr(pubchem_module, "resolve_pubchem_record", resolve)
        table, diagnostics = build_molecule_table(
            ["Molecule"],
            alias={"Molecule": "primary"},
            fallback_queries={"Molecule": ["successful fallback"]},
        )
        assert calls == ["primary", "successful fallback"]
        assert table.iloc[0]["pubchem_query"] == "successful fallback"
        assert int(table.iloc[0]["cid"]) == 321
        assert diagnostics.iloc[0]["query"] == "successful fallback"


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


class TestScoreCandidateCurrentSchema:
    """M-03: scoring must read stereo from the current PubChem ``SMILES`` key,
    not the dead ``IsomericSMILES`` key. These records use the post-2025 schema
    only (no legacy key), so they fail under the old ``prop.get("IsomericSMILES")``
    read and pass once scoring routes through ``_isomeric_smiles``."""

    def _prop(self, cid, smiles, formula="C5H10O5"):
        return {
            "CID": cid,
            "MolecularFormula": formula,
            "SMILES": smiles,                          # stereo-bearing (current)
            "ConnectivitySMILES": smiles.replace("@", ""),  # flat (current)
            "IUPACName": "",
            "Title": "",
        }

    def test_stereo_bonus_from_current_smiles_key(self):
        # D-ribose-like: stereo under "SMILES". With the dead-key bug both read
        # "" and tie; via the helper the stereo record earns the +20 bonus.
        stereo = self._prop(5793, "C1[C@H]([C@H](C(C(O1)O)O)O)O")
        flat = self._prop(5793, "C1C(C(C(C(O1)O)O)O)O")
        assert score_candidate(stereo) > score_candidate(flat)

    def test_stereo_candidate_wins_tie_over_lower_cid_achiral(self):
        # Correct chiral record (higher CID) vs a lower-CID achiral candidate,
        # same formula. The stereo bonus (+20) must outweigh the low-CID bonus
        # (+10) so the resolver picks the stereodefined CID — this is exactly the
        # wrong-CID selection M-03 fixes.
        chiral = self._prop(5793, "C1[C@H]([C@H](C(C(O1)O)O)O)O")
        achiral_low_cid = self._prop(500, "C1C(C(C(C(O1)O)O)O)O")
        assert score_candidate(chiral) > score_candidate(achiral_low_cid)

    def test_directional_bond_stereo_markers_receive_bonus(self):
        flat = self._prop(50_000, "FC=CF", formula="C2H2F2")
        for smiles in ("F/C=C/F", r"F\C=C\F"):
            stereo = self._prop(50_000, smiles, formula="C2H2F2")
            assert score_candidate(stereo) > score_candidate(flat)

    def test_backslash_stereo_beats_lower_cid_stereo_free_candidate(self):
        stereo = self._prop(50_000, r"F\C=C\F", formula="C2H2F2")
        stereo_free = self._prop(1_000, "FC=CF", formula="C2H2F2")
        assert score_candidate(stereo) > score_candidate(stereo_free)

    def test_connectivity_smiles_never_supplies_stereo_bonus(self):
        prop = self._prop(50_000, "FC=CF", formula="C2H2F2")
        prop["ConnectivitySMILES"] = r"F\C=C\F"
        assert score_candidate(prop, prefer_stereo=True) == score_candidate(
            prop, prefer_stereo=False
        )


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
