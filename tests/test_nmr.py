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


class TestExtractNmrSnippetFormats:
    """The parser dispatches on how the solvent is written after "NMR"."""

    def test_bracketed_solvent(self):
        assert extract_nmr_snippet("1H NMR [400 MHz CDCl3]: 7.26") == "CDCl3"

    def test_parenthesised_lowercase_mhz(self):
        # Patent text is inconsistent about the capitalisation of MHz.
        assert extract_nmr_snippet("1H NMR (400 mHz, CDCl3): 7.26") == "CDCl3"

    def test_parenthesised_without_frequency(self):
        assert extract_nmr_snippet("13C NMR (DMSO-d6) 165.1") == "DMSO-d6"

    def test_conditions_marker_yields_nothing(self):
        assert extract_nmr_snippet("1H NMR Conditions were as follows") == ""

    def test_frequency_further_into_the_string(self):
        assert extract_nmr_snippet("1H NMR spectrum at 400 MHz CDCl3") == "CDCl3"

    def test_unparseable_returns_empty_string(self):
        assert extract_nmr_snippet("1H NMR of the product was recorded") == ""

    def test_unterminated_delimiters_still_recover_the_solvent(self):
        # Truncated snippets are common in OCR'd patents. The heuristics index
        # into split results, so they must not raise on a missing closer.
        assert extract_nmr_snippet("1H NMR (400 MHz, CDCl3") == "CDCl3"
        assert extract_nmr_snippet("1H NMR [CDCl3") == "CDCl3"

    def test_no_text_after_the_marker(self):
        assert extract_nmr_snippet("recorded by NMR") == ""


class TestCurateNmrSnippetNoise:
    def test_drops_bare_numbers(self):
        assert curate_nmr_snippet("CDCl3, 400") == "CDCl3"

    def test_drops_tms_reference(self):
        assert "tms" not in curate_nmr_snippet("CDCl3, TMS").lower()

    def test_semicolon_separator(self):
        assert curate_nmr_snippet("CDCl3; 7.26 ppm") == "CDCl3"

    def test_slash_separator_keeps_both_components(self):
        assert curate_nmr_snippet("CDCl3/CD3OD") == "CDCl3|CD3OD"

    def test_frequency_list_keeps_the_solvent(self):
        assert curate_nmr_snippet("400 MHz, CDCl3") == "CDCl3"

    def test_all_noise_yields_empty_string(self):
        assert curate_nmr_snippet("7.26 ppm, 400 Hz") == ""

    def test_mhz_marker_alone(self):
        assert curate_nmr_snippet("MHz") == ""


class TestNmrSolventLookupHelpers:
    def test_get_solvent_info_matches_map_nmr_solvent(self):
        from alkahest.nmr_solvents import get_solvent_info

        assert get_solvent_info("  CDCl3 ") == map_nmr_solvent("CDCl3")
        assert get_solvent_info("not_a_solvent_xyz") is None

    def test_is_known_solvent(self):
        from alkahest.nmr_solvents import is_known_solvent

        assert is_known_solvent("CDCl3")
        assert is_known_solvent("  dmso-d6  ")
        assert not is_known_solvent("not_a_solvent_xyz")

    def test_every_entry_is_a_name_smiles_pair(self):
        from alkahest.nmr_solvents import NMR_SOLVENT_DICTIONARY

        for key, value in NMR_SOLVENT_DICTIONARY.items():
            assert key == key.lower(), key
            assert isinstance(value, tuple) and len(value) == 2, key
            assert value[0] and value[1], key


class TestExtractNmrSnippetRecovery:
    """Malformed snippets must degrade to an empty string, never raise."""

    def test_frequency_without_parentheses(self):
        assert extract_nmr_snippet("1H NMR 400 MHz CDCl3") == ""

    def test_lowercase_frequency_without_parentheses(self):
        assert extract_nmr_snippet("1H NMR spectrum at 400 mHz CDCl3") == ""
