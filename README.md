# Alkahest

[![Tests](https://github.com/mrodobbe/Alkahest/actions/workflows/tests.yml/badge.svg)](https://github.com/mrodobbe/Alkahest/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/mrodobbe/Alkahest/branch/master/graph/badge.svg)](https://codecov.io/gh/mrodobbe/Alkahest)
[![PyPI](https://img.shields.io/pypi/v/alkahest-chem.svg)](https://pypi.org/project/alkahest-chem/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)

Find out which solvents a written experimental procedure actually used, and what
role each one played.

Alkahest reads free-text procedures of the kind found in patents and papers, and
returns the solvents as SMILES, sorted into the stage of the experiment they
belong to. It is a rule-based pipeline: a dictionary of name variants, phase
boundary detection, and local context rules. There is no model to train and no
network call, so the same text always gives the same answer.

This package accompanies the paper *The Stubborn Persistence of Toxic Solvents in
Chemical Synthesis* by Maarten R. Dobbelaere and Helen F. Sneddon, accepted in
*Angewandte Chemie International Edition*.

An interactive companion to the study, which lets you look up the solvents
actually used for a given reaction type across 1.34M USPTO reactions
(1976-2016), is available at **[solventexplorer.com](https://solventexplorer.com)**.

## Installation

```bash
pip install alkahest-chem
```

> **Note on the name.** The distribution is `alkahest-chem`; the import name is
> `alkahest`. Do not run `pip install alkahest` — that name belongs to an
> unrelated computer algebra system on PyPI.

Optional extras:

| Extra | Adds | Install |
|---|---|---|
| `rdkit` | SMILES canonicalisation, so GSK scoring accepts any valid notation | `pip install "alkahest-chem[rdkit]"` |
| `ord` | Reading Open Reaction Database protocol buffer files | `pip install "alkahest-chem[ord]"` |
| `rxn-insight` | Reaction classification via Rxn-INSIGHT | `pip install "alkahest-chem[rxn-insight]"` |
| `all` | All of the above | `pip install "alkahest-chem[all]"` |

## Quick Start

```python
from alkahest import extract_solvents_categorized

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

print(result.to_smiles_dict())
# {'reaction': 'C1CCOC1',
#  'workup': 'CCOC(C)=O.O',
#  'purification': 'CCCCCC.CCOC(C)=O',
#  'analytical': '[2H]C(Cl)(Cl)Cl'}
```

Every solvent is placed in exactly one of four roles:

| Role | Meaning | Typical cue in the text |
|---|---|---|
| `reaction` | the reaction medium itself | "dissolved in", "a solution of", anything before the first phase boundary |
| `workup` | post-reaction processing | "quenched with", "extracted with", "washed with" |
| `purification` | chromatography, recrystallisation | "column", "silica", "eluted with" |
| `analytical` | NMR, MS and other measurements | after an "NMR" marker; deuterated solvents always land here |

A category with no solvents reads `"solvent-free"` in the SMILES dict, and
several solvents in one category are joined with `.` into a single SMILES string.

### Scoring solvents

`score_solvent` grades a solvent with the 2016 GSK Solvent Selection Guide.
Mixtures take the worst component's score, and deuterated solvents are scored as
their non-deuterated parent.

```python
from alkahest import score_solvent

score_solvent("C1CCOC1")            # THF   -> 'R'  (red: major issues)
score_solvent("CS(C)=O")            # DMSO  -> 'A'  (amber: some issues)
score_solvent("CCOC(C)=O")          # EtOAc -> 'G'  (green: few issues)
score_solvent("CCOC(C)=O.C1CCOC1")  # worst of the two -> 'R'
score_solvent("[2H]C(Cl)(Cl)Cl")    # CDCl3, scored as chloroform -> 'R'
score_solvent("CCN(C(C)C)C(C)C")    # DIPEA, not in the guide -> 'Unknown'
```

### NMR solvents

Patent text spells NMR solvents in hundreds of ways, including OCR damage.
`map_nmr_solvent` resolves them to a canonical name and SMILES, and returns
`None` for anything it does not recognise.

```python
from alkahest import map_nmr_solvent

map_nmr_solvent("CDCl3")     # ('CDCl3', '[2H]C(Cl)(Cl)Cl')
map_nmr_solvent("cdc13")     # OCR damage, digit 1 for letter l -> same result
map_nmr_solvent("dmso-d6")   # ('DMSO-d6', '[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]')
map_nmr_solvent("meod")      # ('CD3OD', '[2H]C([2H])([2H])O[2H]')
map_nmr_solvent("xyzzy")     # None
```

### Whole DataFrames

`extract_solvents_batch` adds four columns to a DataFrame of procedures. It
modifies the frame in place and returns it, so pass a copy if you need the
original untouched.

```python
import pandas as pd
from alkahest import extract_solvents_batch

df = pd.read_parquet("reactions.parquet")
df = extract_solvents_batch(df, procedure_column="procedure")
# Adds SOLV_RXN, SOLV_WORKUP, SOLV_PURIF, SOLV_ANAL
```

## Command line

Installing the package provides two commands. Both split their input into
chunks, so a large corpus can run as an array job on a cluster.

```bash
# Solvent extraction: chunk 0 of 130
alkahest-extract --input reactions.parquet --output_dir output/ \
    --job_index 0 --n_chunks 130

# Rxn-INSIGHT reaction analysis, 6 cores
alkahest-rxn-insight --input reactions.parquet --output_dir output/ \
    -i 0 -n 1000 -c 6
```

Chunk `i` of `n` covers rows `i * len(df) // n` to `(i + 1) * len(df) // n`, so
the chunks tile the input exactly whatever the row count. Under SLURM:

```bash
#SBATCH --array=0-129
alkahest-extract --input reactions.parquet --output_dir output/ \
    --job_index $SLURM_ARRAY_TASK_ID --n_chunks 130
```

## API

| Function | Returns |
|---|---|
| `extract_solvents_categorized(text)` | `CategorizedSolvents` with a dict per role |
| `extract_solvents_batch(df, ...)` | the DataFrame, with four solvent columns added |
| `find_phase_boundaries(text)` | character offset where each phase begins, or `None` |
| `score_solvent(smiles)` | `'G'`, `'A'`, `'R'` or `'Unknown'` |
| `get_gsk_guide()` | the guide as a DataFrame |
| `map_nmr_solvent(text)` | `(name, smiles)` or `None` |
| `extract_nmr_snippet(procedure)` | the raw solvent snippet following "NMR" |
| `curate_nmr_snippet(snippet)` | that snippet with frequencies and shifts stripped |

The dictionaries are importable too: `SOLVENT_DICT`, `DEUTERATED_SOLVENTS` and
`NMR_SOLVENT_DICTIONARY`.

## Development

```bash
git clone https://github.com/mrodobbe/Alkahest.git
cd Alkahest
pip install -e ".[test]"
pytest --cov=alkahest --cov-report=term-missing
```

The HTML documentation is built with Sphinx:

```bash
pip install ".[docs]"
sphinx-build -b html docs docs/_build/html
```

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
