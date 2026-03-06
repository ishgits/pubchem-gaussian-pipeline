"""Tests for pipeline.utils"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.utils import sanitize_basename, normalize_cid


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
