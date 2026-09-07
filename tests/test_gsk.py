"""Tests for GSK Solvent Selection Guide scoring."""

import pandas as pd
import pytest

from alkahest import score_solvent, get_gsk_guide


class TestGetGskGuide:
    def test_returns_dataframe(self):
        guide = get_gsk_guide()
        assert isinstance(guide, pd.DataFrame)

    def test_expected_columns(self):
        guide = get_gsk_guide()
        assert "Solvent" in guide.columns
        assert "SMILES" in guide.columns
        assert "RAG" in guide.columns

    def test_has_entries(self):
        guide = get_gsk_guide()
        assert len(guide) > 100

    def test_rag_values(self):
        guide = get_gsk_guide()
        assert set(guide["RAG"].unique()).issubset({"G", "A", "R"})

    def test_commas_in_solvent_names(self):
        # Regression: CSV parser must handle commas in names like "1,4-dioxane"
        guide = get_gsk_guide()
        names = guide["Solvent"].str.lower().tolist()
        assert any("dioxane" in n for n in names)

    def test_canonical_smiles_column_present(self):
        # Lookups are made on canonical SMILES; the column must ship with the data
        # so that scoring works without RDKit installed.
        guide = get_gsk_guide()
        assert "Canonical SMILES" in guide.columns
        assert guide["Canonical SMILES"].notna().all()
        assert (guide["Canonical SMILES"].str.len() > 0).all()


class TestScoreSolvent:
    def test_green_solvent(self):
        assert score_solvent("CCOC(C)=O") == "G"  # EtOAc

    def test_amber_solvent(self):
        assert score_solvent("CS(C)=O") == "A"  # DMSO

    def test_red_solvent(self):
        assert score_solvent("C1CCOC1") == "R"  # THF
        assert score_solvent("ClCCl") == "R"  # DCM

    def test_unknown_solvent(self):
        assert score_solvent("C#CC#C") == "Unknown"

    def test_empty_input(self):
        assert score_solvent("") == "Unknown"
        assert score_solvent("solvent-free") == "Unknown"
        assert score_solvent("not-reported") == "Unknown"

    def test_mixture_worst_score(self):
        # EtOAc (G) + THF (R) -> R
        assert score_solvent("CCOC(C)=O.C1CCOC1") == "R"

    def test_deuterated_mapped_to_parent(self):
        # CDCl3 -> CHCl3 -> R
        assert score_solvent("[2H]C(Cl)(Cl)Cl") == "R"

    def test_deuterated_dmso(self):
        # DMSO-d6 -> DMSO -> A
        assert score_solvent("[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]") == "A"

    def test_water(self):
        assert score_solvent("O") == "G"

    def test_non_canonical_guide_entries_are_scored(self):
        # Regression: the guide stores many SMILES in a non-canonical form, so
        # exact string matching used to return "Unknown" for ~40% of the
        # solvent dictionary. Lookups now go through canonical SMILES.
        assert score_solvent("COc1ccccc1") == "G"  # anisole
        assert score_solvent("OCc1ccccc1") == "G"  # benzyl alcohol
        assert score_solvent("O=C1CCCCC1") == "G"  # cyclohexanone
        assert score_solvent("Clc1ccccc1") == "R"  # chlorobenzene
        assert score_solvent("O=S1(=O)CCCC1") == "R"  # sulfolane
        assert score_solvent("NC=O") == "R"  # formamide

    def test_whole_solvent_dictionary_is_covered(self):
        # Only solvents genuinely absent from the 2016 guide may be Unknown.
        from alkahest.solvents import SOLVENT_DICT

        unknown = {s for s in set(SOLVENT_DICT.values()) if score_solvent(s) == "Unknown"}
        assert unknown == {
            "CCN(C(C)C)C(C)C",  # DIPEA
            "COCCOCCOCCOC",  # triglyme
            "Cc1cccc(C)c1",  # m-xylene
            "Cc1ccccc1C",  # o-xylene
            "FC(F)(F)C(=O)O",  # trifluoroacetic acid
        }

    def test_every_deuterated_solvent_is_scored(self):
        from alkahest.solvents import DEUTERATED_SOLVENTS

        for smiles in set(DEUTERATED_SOLVENTS.values()):
            assert score_solvent(smiles) != "Unknown", smiles

    def test_deuterated_dioxane_maps_to_parent(self):
        # Regression: dioxane-d8 was missing from the deuterated parent map.
        assert score_solvent("[2H]C1([2H])OC([2H])([2H])C([2H])([2H])OC1([2H])[2H]") == "R"

    def test_mixture_ignores_unknown_component(self):
        # A known-red component still dominates an unrecognised one.
        assert score_solvent("C#CC#C.C1CCOC1") == "R"

    def test_guide_is_cached(self):
        assert get_gsk_guide() is get_gsk_guide()

    def test_sentinels_inside_a_mixture_are_skipped(self):
        assert score_solvent("solvent-free.C1CCOC1") == "R"
        assert score_solvent("C1CCOC1..O") == "R"
        assert score_solvent("solvent-free.not-reported") == "Unknown"


class TestScoringWithoutRdkit:
    """Scoring must work on the shipped canonical column alone."""

    @pytest.fixture
    def no_rdkit(self, monkeypatch):
        import sys

        monkeypatch.setitem(sys.modules, "rdkit", None)
        monkeypatch.setitem(sys.modules, "rdkit.Chem", None)

    def test_canonical_column_carries_the_lookup(self, no_rdkit):
        from alkahest.gsk import _canonicalize, _strip_isotopes

        assert _canonicalize("COc1ccccc1") is None
        assert _strip_isotopes("[2H]C(Cl)(Cl)Cl") is None
        # Still scored, because the guide ships canonical SMILES.
        assert score_solvent("COc1ccccc1") == "G"
        assert score_solvent("Clc1ccccc1") == "R"

    def test_deuterated_falls_back_to_the_parent_map(self, no_rdkit):
        assert score_solvent("[2H]C(Cl)(Cl)Cl") == "R"
        assert score_solvent("[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]") == "A"

    def test_unlisted_notation_is_unknown_without_rdkit(self, no_rdkit):
        # THF written ring-first is valid but is not one of the spellings the
        # guide lists, so without RDKit there is nothing to match it against.
        assert score_solvent("O1CCCC1") == "Unknown"


class TestScoringWithRdkit:
    def test_arbitrary_notation_resolves(self):
        pytest.importorskip("rdkit")
        assert score_solvent("O1CCCC1") == "R"  # THF, written ring-first
        assert score_solvent("OC(=O)C") == "A"  # acetic acid, reversed

    def test_isotopes_are_stripped_generically(self):
        pytest.importorskip("rdkit")
        # Chlorobenzene-d5 is not in the deuterated parent map, so it can only
        # be scored by stripping the isotope labels down to chlorobenzene.
        from alkahest.gsk import _DEUTERATED_TO_PARENT

        chlorobenzene_d5 = "[2H]c1c([2H])c([2H])c(Cl)c([2H])c1[2H]"
        assert chlorobenzene_d5 not in _DEUTERATED_TO_PARENT
        assert score_solvent(chlorobenzene_d5) == "R"

    def test_deuterated_solvent_absent_from_the_guide_is_unknown(self):
        pytest.importorskip("rdkit")
        assert score_solvent("[2H]C([2H])([2H])C([2H])([2H])[2H]") == "Unknown"

    def test_invalid_smiles_is_unknown_not_an_error(self):
        pytest.importorskip("rdkit")
        assert score_solvent("not_a_smiles") == "Unknown"
        assert score_solvent("C1CCOC1.not_a_smiles") == "R"
