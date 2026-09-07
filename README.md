# Alkahest

[![Tests](https://github.com/mrodobbe/Alkahest/actions/workflows/tests.yml/badge.svg)](https://github.com/mrodobbe/Alkahest/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/mrodobbe/Alkahest/branch/master/graph/badge.svg)](https://codecov.io/gh/mrodobbe/Alkahest)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)

Solvent extraction and classification from chemical procedure texts using rule-based NLP.

This package accompanies the paper: *The Stubborn Persistence of Toxic Solvents in Chemical Synthesis* by Maarten R. Dobbelaere and Helen F. Sneddon, accepted in *Angewandte Chemie International Edition*.

An interactive companion to the study, which lets you look up the solvents actually used for a given reaction type across 1.34M USPTO reactions (1976-2016), is available at **[solventexplorer.com](https://solventexplorer.com)**.

## Features

- **Solvent extraction**: Identify solvents from experimental procedure texts using a dictionary of 350+ name variants (IUPAC, common names, abbreviations)
- **Phase classification**: Classify solvents by experimental role (reaction, workup, purification, analytical) using phase boundary detection and context-aware rules
- **NMR solvent mapping**: Map 600+ NMR solvent string variants (including typos, OCR errors) to canonical SMILES
- **GSK scoring**: Score solvents as Green/Amber/Red using the 2016 GSK Solvent Selection Guide ([Alder et al., *Green Chem.*, 2016, 18, 3879-3890](https://doi.org/10.1039/C6GC00611F))

## Installation

```bash
pip install alkahest-chem
```

> **Note on the package name.** The distribution is named `alkahest-chem`; the
> import name is `alkahest`. Do **not** run `pip install alkahest` — that name
> belongs to an unrelated computer algebra system on PyPI.

For reading Open Reaction Database (ORD) files:
```bash
pip install "alkahest-chem[ord]"
```

For Rxn-INSIGHT integration (reaction classification):
```bash
pip install "alkahest-chem[rxn-insight]"
```

For everything:
```bash
pip install "alkahest-chem[all]"
```

To install the development version from source:

```bash
git clone https://github.com/mrodobbe/Alkahest.git
cd Alkahest
pip install .
```

## Quick Start

```python
from alkahest import extract_solvents_categorized, score_solvent, map_nmr_solvent

# Extract solvents from a procedure text
result = extract_solvents_categorized("""
    The compound (1.0 g) was dissolved in THF (20 mL) and stirred for 2 h.
    The mixture was quenched with water and extracted with ethyl acetate.
    The combined organic layers were washed with brine and concentrated.
    Purification by column chromatography (silica gel, hexanes/EtOAc 4:1)
    gave the product as a white solid (0.85 g, 78%), mp 112-114 C.
    1H NMR (400 MHz, CDCl3): delta 7.45 (d, 2H).
""")
print(result)
# CategorizedSolvents(
#   reaction=['thf'],
#   workup=['ethyl acetate', 'water'],
#   purification=['hexanes', 'etoac'],
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

## Scope and limitations

The extractor is a dictionary and regex pipeline tuned for the high-volume
patterns of USPTO patent prose, not a general-purpose chemical named-entity
recogniser. It is deliberately recall-oriented, and the published analysis
relies on aggregate trends across 1.34M reactions rather than on any single
procedure being parsed perfectly. Two consequences are worth knowing before
reusing it on other corpora:

- Short dictionary keys collide with common non-solvent tokens. `EDC` (the
  coupling reagent) matches 1,2-dichloroethane, `DEC` (as in a decomposition
  melting point) matches diethyl carbonate, and `DMA` matches
  N,N-dimethylaniline rather than dimethylacetamide.
- Names are matched independently rather than consumed longest-first, so a
  span may be counted twice. "Petroleum ether" also registers diethyl ether.

Solvent identities are reported as SMILES; role assignment near a phase
boundary is heuristic and is the least reliable part of the output.

## License

This software is available for academic and non-commercial use only. See [LICENSE](LICENSE) for details. For commercial licensing, contact mrodobbe.Dobbelaere@UGent.be.

## Citation

If you use this software, please cite:

```
M. R. Dobbelaere and H. F. Sneddon, "The Stubborn Persistence of Toxic Solvents
in Chemical Synthesis", Angew. Chem. Int. Ed., 2026, accepted.
```

If you use the GSK solvent scoring functionality (`score_solvent`, `get_gsk_guide`), please also cite the original guide:

```
C.M. Alder, J.D. Hayler, R.K. Henderson, A.M. Redman, L. Shukla, L.E. Shuster, H.F. Sneddon,
"Updating and further expanding GSK's solvent sustainability guide",
Green Chem., 2016, 18, 3879-3890. DOI: 10.1039/C6GC00611F
```
