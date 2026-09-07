"""Tests for solvent extraction and classification."""

import pandas as pd
import pytest

from alkahest import CategorizedSolvents, extract_solvents_categorized, extract_solvents_batch


class TestCategorizedSolvents:
    def test_empty(self):
        cs = CategorizedSolvents()
        assert cs.reaction == {}
        assert cs.workup == {}
        assert cs.purification == {}
        assert cs.analytical == {}

    def test_to_smiles_dict_empty(self):
        cs = CategorizedSolvents()
        d = cs.to_smiles_dict()
        assert d["reaction"] == "solvent-free"
        assert d["workup"] == "solvent-free"
        assert d["purification"] == "solvent-free"
        assert d["analytical"] == "solvent-free"

    def test_to_smiles_dict_populated(self):
        cs = CategorizedSolvents(reaction={"thf": "C1CCOC1"})
        d = cs.to_smiles_dict()
        assert d["reaction"] == "C1CCOC1"

    def test_get_smiles_accessors(self):
        cs = CategorizedSolvents(
            reaction={"thf": "C1CCOC1", "tetrahydrofuran": "C1CCOC1"},
            workup={"water": "O"},
        )
        # Duplicates collapsed
        assert cs.get_reaction_smiles() == ["C1CCOC1"]
        assert cs.get_workup_smiles() == ["O"]
        assert cs.get_purification_smiles() == []

    def test_repr(self):
        cs = CategorizedSolvents(reaction={"thf": "C1CCOC1"})
        r = repr(cs)
        assert "reaction=" in r
        assert "thf" in r


class TestExtractSolventsCategorized:
    def test_empty_input(self):
        assert extract_solvents_categorized("") == CategorizedSolvents()
        assert extract_solvents_categorized(None) == CategorizedSolvents()

    def test_reaction_solvent(self):
        result = extract_solvents_categorized(
            "The compound was dissolved in THF (20 mL) and stirred."
        )
        assert "C1CCOC1" in result.get_reaction_smiles()

    def test_workup_solvent(self):
        result = extract_solvents_categorized(
            "After completion, the mixture was quenched with water "
            "and extracted with ethyl acetate."
        )
        assert "O" in result.get_workup_smiles()
        assert "CCOC(C)=O" in result.get_workup_smiles()

    def test_purification_solvent(self):
        result = extract_solvents_categorized(
            "Purification by column chromatography on silica gel "
            "(hexanes/ethyl acetate 4:1) gave the product."
        )
        purif = result.get_purification_smiles()
        assert len(purif) > 0

    def test_analytical_deuterated(self):
        result = extract_solvents_categorized(
            "1H NMR (400 MHz, CDCl3): delta 7.26."
        )
        assert len(result.get_analytical_smiles()) > 0

    def test_deuterated_always_analytical(self):
        # Even if a deuterated solvent appears early in the text,
        # it should be classified as analytical
        result = extract_solvents_categorized(
            "CDCl3 was used. The compound was dissolved in THF."
        )
        analytical_smiles = result.get_analytical_smiles()
        assert "[2H]C(Cl)(Cl)Cl" in analytical_smiles

    def test_full_procedure(self):
        procedure = """
        The starting material (1.0 g) was dissolved in THF (20 mL) and cooled to
        0 degrees C. The reaction was stirred at room temperature for 2 hours.
        The mixture was quenched with water and extracted with ethyl acetate
        (3 x 30 mL). The combined organic layers were washed with brine.
        Purification by column chromatography (silica gel, hexanes/ethyl acetate
        4:1) gave the product. 1H NMR (400 MHz, CDCl3): delta 7.45 (d, 2H).
        """
        result = extract_solvents_categorized(procedure)
        assert len(result.get_reaction_smiles()) > 0
        assert len(result.get_workup_smiles()) > 0
        assert len(result.get_analytical_smiles()) > 0

    def test_no_false_positive_in_compound_name(self):
        # "ether" inside "diethylether" context shouldn't create spurious matches
        result = extract_solvents_categorized(
            "trimethylamine (2 eq) was added to the solution"
        )
        # "trimethylamine" should not produce a match for "amine"
        assert result.reaction == {} or all(
            "amine" not in k for k in result.reaction
        )


class TestExtractSolventsBatch:
    def test_batch_processing(self):
        df = pd.DataFrame({
            "procedure": [
                "Dissolved in THF and stirred.",
                "The mixture was extracted with ethyl acetate.",
            ]
        })
        result = extract_solvents_batch(df, procedure_column="procedure")
        assert "SOLV_RXN" in result.columns
        assert "SOLV_WORKUP" in result.columns
        assert "SOLV_PURIF" in result.columns
        assert "SOLV_ANAL" in result.columns
        assert len(result) == 2

    def test_custom_prefix(self):
        df = pd.DataFrame({"text": ["Dissolved in THF."]})
        result = extract_solvents_batch(df, procedure_column="text", prefix="MY_")
        assert "MY_RXN" in result.columns

    def test_empty_dataframe(self):
        # Regression: progress_apply on an empty frame returns an empty
        # DataFrame with no columns, which used to raise KeyError.
        df = pd.DataFrame({"procedure": []})
        result = extract_solvents_batch(df, procedure_column="procedure")
        assert len(result) == 0
        for column in ("SOLV_RXN", "SOLV_WORKUP", "SOLV_PURIF", "SOLV_ANAL"):
            assert column in result.columns

    def test_single_row(self):
        df = pd.DataFrame({"procedure": ["Dissolved in DMSO and stirred."]})
        result = extract_solvents_batch(df, procedure_column="procedure")
        assert len(result) == 1
        assert "SOLV_RXN" in result.columns
        assert result["SOLV_RXN"].iloc[0] != "solvent-free"
