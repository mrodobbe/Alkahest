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
)
from alkahest.gsk import score_solvent, get_gsk_guide
from alkahest.nmr import map_nmr_solvent

__version__ = "0.1.0"

__all__ = [
    "CategorizedSolvents",
    "extract_solvents_categorized",
    "extract_solvents_batch",
    "score_solvent",
    "get_gsk_guide",
    "map_nmr_solvent",
]
