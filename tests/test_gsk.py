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
