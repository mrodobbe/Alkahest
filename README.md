# Alkahest

[![Tests](https://github.com/mrodobbe/alkahest/actions/workflows/tests.yml/badge.svg)](https://github.com/mrodobbe/alkahest/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/mrodobbe/alkahest/branch/main/graph/badge.svg)](https://codecov.io/gh/mrodobbe/alkahest)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)

Solvent extraction and classification from chemical procedure texts using rule-based NLP.

This package accompanies the paper: *Solvent abuse: The persistence of toxic solvents in chemical synthesis* by Maarten R. Dobbelaere and Helen F. Sneddon.

## Features

- **Solvent extraction**: Identify solvents from experimental procedure texts using a dictionary of 350+ name variants (IUPAC, common names, abbreviations)
- **Phase classification**: Classify solvents by experimental role (reaction, workup, purification, analytical) using phase boundary detection and context-aware rules
- **NMR solvent mapping**: Map 600+ NMR solvent string variants (including typos, OCR errors) to canonical SMILES
- **GSK scoring**: Score solvents as Green/Amber/Red using the 2016 GSK Solvent Selection Guide ([Alder et al., *Green Chem.*, 2016, 18, 3879-3890](https://doi.org/10.1039/C6GC00611F))

## Installation

```bash
pip install .
```

For reading Open Reaction Database (ORD) files:
```bash
pip install ".[ord]"
```

For Rxn-INSIGHT integration (reaction classification):
```bash
pip install ".[rxn-insight]"
```

## Quick Start

```python
from alkahest import extract_solvents_categorized, score_solvent, map_nmr_solvent

# Extract solvents from a procedure text
result = extract_solvents_categorized("""
    The compound was dissolved in THF (20 mL) and stirred for 2 h.
    The mixture was quenched with water and extracted with ethyl acetate.
    Purification by column chromatography (silica, hexanes/EtOAc).
    1H NMR (CDCl3): delta 7.45 (d, 2H).
""")
print(result)
# CategorizedSolvents(
#   reaction=['thf'],
#   workup=['water', 'ethyl acetate'],
#   purification=['hexanes', 'ethyl acetate'],
#   analytical=['cdcl3']
# )

# Score solvents with the GSK guide
score_solvent("C1CCOC1")   # THF -> "R" (Red)
score_solvent("CCOC(C)=O") # EtOAc -> "G" (Green)

# Map NMR solvent names (handles typos and OCR errors)
map_nmr_solvent("cdc13")   # ('CDCl3', '[2H]C(Cl)(Cl)Cl')
```

### Batch Processing

```python
import pandas as pd
from alkahest import extract_solvents_batch

df = pd.read_parquet("reactions.parquet")
df = extract_solvents_batch(df, procedure_column="procedure")
# Adds: SOLV_RXN, SOLV_WORKUP, SOLV_PURIF, SOLV_ANAL
```

## HPC Scripts

For processing large datasets in parallel:

```bash
# Rxn-INSIGHT reaction analysis
python scripts/run_rxn_insight.py --input data.parquet --output_dir output/ -i 0 -n 1000 -c 6

# Solvent extraction
python scripts/extract_solvents.py --input data.parquet --output_dir output/ --job_index 0 --n_chunks 130
```

## License

This software is available for academic and non-commercial use only. See [LICENSE](LICENSE) for details. For commercial licensing, contact mrodobbe.Dobbelaere@UGent.be.

## Citation

If you use this software, please cite:

```
M.R. Dobbelaere and H.F. Sneddon, "Solvent abuse: The persistence of toxic solvents in chemical synthesis", 2026.
```

If you use the GSK solvent scoring functionality (`score_solvent`, `get_gsk_guide`), please also cite the original guide:

```
C.M. Alder, J.D. Hayler, R.K. Henderson, A.M. Redman, L. Shukla, L.E. Shuster, H.F. Sneddon,
"Updating and further expanding GSK's solvent sustainability guide",
Green Chem., 2016, 18, 3879-3890. DOI: 10.1039/C6GC00611F
```
