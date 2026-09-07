"""
Alkahest: Solvent extraction and classification from chemical procedure texts.

Provides tools for:
- Extracting solvents from experimental procedure texts using rule-based NLP
- Classifying solvents by role: reaction, workup, purification, analytical
- Mapping NMR solvent mentions (including typos/OCR errors) to canonical SMILES
- Scoring solvents using the GSK Solvent Selection Guide (red/amber/green)
"""

from alkahest.extraction import (
    CategorizedSolvents,
    extract_solvents_categorized,
    extract_solvents_batch,
    find_phase_boundaries,
)
from alkahest.gsk import score_solvent, get_gsk_guide
from alkahest.nmr import map_nmr_solvent, extract_nmr_snippet, curate_nmr_snippet
from alkahest.nmr_solvents import NMR_SOLVENT_DICTIONARY
from alkahest.solvents import SOLVENT_DICT, DEUTERATED_SOLVENTS

__version__ = "0.1.0"

__all__ = [
    # Extraction
    "CategorizedSolvents",
    "extract_solvents_categorized",
    "extract_solvents_batch",
    "find_phase_boundaries",
    # GSK scoring
    "score_solvent",
    "get_gsk_guide",
    # NMR
    "map_nmr_solvent",
    "extract_nmr_snippet",
    "curate_nmr_snippet",
    # Dictionaries
    "SOLVENT_DICT",
    "DEUTERATED_SOLVENTS",
    "NMR_SOLVENT_DICTIONARY",
]
