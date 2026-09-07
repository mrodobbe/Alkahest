"""
Example: Extract solvents from ORD data.

This example shows the full pipeline from Open Reaction Database (ORD) files
through Rxn-INSIGHT analysis to solvent extraction and GSK scoring.
"""

import pandas as pd

# =====================================================================
# Step 1: Load ORD data using Rxn-INSIGHT
# =====================================================================
# Rxn-INSIGHT provides ORDDatabase for reading ORD protocol buffer files.
#
#   import rxn_insight as ri
#
#   db = ri.ORDDatabase("dataset.pb.gz")
#   df = db.analyze()
#   df.to_parquet("analyzed_reactions.parquet")
#
# The resulting DataFrame contains columns like REACTION (SMILES),
# REACTANTS, PRODUCTS, REAGENT, CATALYST, SOLVENT, procedure, etc.

# =====================================================================
# Step 2: Extract solvents from procedure texts
# =====================================================================
from alkahest import extract_solvents_categorized, extract_solvents_batch

# Single procedure text
procedure = """
The starting material (1.0 g) was dissolved in THF (20 mL) and cooled to
0 degrees C. NaH (60% in mineral oil, 0.5 g) was added portionwise. The
reaction was stirred at room temperature for 2 hours. The mixture was
quenched with water and extracted with ethyl acetate (3 x 30 mL). The
combined organic layers were washed with brine, dried over Na2SO4, and
concentrated. Purification by column chromatography (silica gel, hexanes/
ethyl acetate 4:1) gave the product as a white solid (0.85 g, 78%),
mp 112-114 degrees C. 1H NMR (400 MHz, CDCl3): delta 7.45 (d, J=8.0 Hz, 2H).
"""
result = extract_solvents_categorized(procedure)
print(result)
# CategorizedSolvents(
#   reaction=['thf'],
#   workup=['ethyl acetate', 'water'],
#   purification=['ethyl acetate', 'hexanes'],
#   analytical=['cdcl3']
# )

# SMILES output
print(result.to_smiles_dict())

# =====================================================================
# Step 3: Batch processing a DataFrame
# =====================================================================
# df = pd.read_parquet("analyzed_reactions.parquet")
# df = extract_solvents_batch(df, procedure_column="procedure", prefix="SOLV_")
# This adds columns: SOLV_RXN, SOLV_WORKUP, SOLV_PURIF, SOLV_ANAL

# =====================================================================
# Step 4: Score solvents with the GSK Solvent Selection Guide
# =====================================================================
from alkahest import score_solvent

print(score_solvent("C1CCOC1"))   # THF -> R (Red)
print(score_solvent("CCOC(C)=O")) # EtOAc -> G (Green)
print(score_solvent("ClCCl"))     # DCM -> R (Red)
print(score_solvent("CS(C)=O"))   # DMSO -> A (Amber)

# Deuterated solvents are mapped to parents for scoring
print(score_solvent("[2H]C(Cl)(Cl)Cl"))  # CDCl3 -> scored as CHCl3 -> R
print(score_solvent("COc1ccccc1"))       # Anisole -> G
print(score_solvent("CCN(C(C)C)C(C)C"))  # DIPEA: not in the guide -> Unknown

# =====================================================================
# Step 5: Map NMR solvent names to canonical SMILES
# =====================================================================
from alkahest import map_nmr_solvent

print(map_nmr_solvent("CDCl3"))     # ('CDCl3', '[2H]C(Cl)(Cl)Cl')
print(map_nmr_solvent("cdc13"))     # OCR typo -> still maps to CDCl3
print(map_nmr_solvent("dmso-d6"))   # ('DMSO-d6', '[2H]C([2H])...')
print(map_nmr_solvent("meod"))      # ('CD3OD', '[2H]C([2H])([2H])O[2H]')
