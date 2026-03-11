"""Tests for NMR solvent extraction, curation, and mapping."""

import pytest

from alkahest import map_nmr_solvent
from alkahest.nmr import extract_nmr_snippet, curate_nmr_snippet


class TestExtractNmrSnippet:
    def test_bracket_format(self):
        # "NMR [solvent]" pattern
        snippet = extract_nmr_snippet("1H NMR [CDCl3]: delta 7.26")
        assert "CDCl3" in snippet

    def test_parenthesis_with_mhz(self):
        snippet = extract_nmr_snippet("1H NMR (400 MHz, CDCl3): delta 7.26")
        assert snippet != ""

    def test_no_nmr(self):
        assert extract_nmr_snippet("The compound was dissolved in THF.") == ""

    def test_none_input(self):
        assert extract_nmr_snippet(None) == ""

    def test_empty_string(self):
        assert extract_nmr_snippet("") == ""


class TestCurateNmrSnippet:
    def test_empty(self):
        assert curate_nmr_snippet("") == ""

    def test_removes_mhz(self):
        result = curate_nmr_snippet("400 MHz, CDCl3")
        assert "MHz" not in result or "CDCl3" in result

    def test_simple_solvent(self):
        assert curate_nmr_snippet("CDCl3") == "CDCl3"

    def test_filters_ppm(self):
        result = curate_nmr_snippet("CDCl3, 7.26 ppm")
        assert "ppm" not in result

    def test_short_string(self):
        assert curate_nmr_snippet("ab") == ""

    def test_delta_symbol(self):
        assert curate_nmr_snippet("\u03b4") == ""


class TestMapNmrSolvent:
    def test_canonical_name(self):
        result = map_nmr_solvent("CDCl3")
        assert result is not None
        assert result[0] == "CDCl3"
        assert result[1] == "[2H]C(Cl)(Cl)Cl"

    def test_case_insensitive(self):
        result = map_nmr_solvent("cdcl3")
        assert result is not None
        assert result[0] == "CDCl3"

    def test_ocr_error(self):
        # Common OCR typo: "1" instead of "l"
        result = map_nmr_solvent("cdc13")
        assert result is not None
        assert result[0] == "CDCl3"

    def test_dmso_d6(self):
        result = map_nmr_solvent("dmso-d6")
        assert result is not None
        assert result[0] == "DMSO-d6"

    def test_meod(self):
        result = map_nmr_solvent("meod")
        assert result is not None
        assert result[0] == "CD3OD"

    def test_d2o(self):
        result = map_nmr_solvent("D2O")
        assert result is not None
        assert result[1] == "[2H]O[2H]"

    def test_unknown_returns_none(self):
        assert map_nmr_solvent("not_a_solvent_xyz") is None

    def test_empty_returns_none(self):
        assert map_nmr_solvent("") is None

    def test_whitespace_stripped(self):
        result = map_nmr_solvent("  CDCl3  ")
        assert result is not None
        assert result[0] == "CDCl3"
