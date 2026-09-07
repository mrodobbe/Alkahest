User Guide
==========

This guide explains how to use Alkahest to extract and classify solvents
from experimental procedure texts commonly found in patents and journal articles.


Solvent Extraction
------------------

The core function :func:`~alkahest.extract_solvents_categorized` takes a
procedure text and returns solvents classified into four categories:

- **reaction**: solvents used during the chemical reaction itself
- **workup**: solvents used in post-reaction processing (extraction, washing, quenching)
- **purification**: solvents used in column chromatography, recrystallization, etc.
- **analytical**: solvents used for NMR, MS, and other analytical measurements

.. code-block:: python

   from alkahest import extract_solvents_categorized

   result = extract_solvents_categorized("""
       A solution of the amine (1.5 g) in dichloromethane (30 mL) was treated
       with K2CO3 (2 eq) and acyl chloride (1.1 eq) at 0 C. After 2 h,
       the reaction was quenched with sat. NaHCO3 and extracted with DCM
       (3 x 20 mL). The combined organics were dried over MgSO4 and concentrated.
       The residue was purified by flash chromatography (hexanes -> 50% EtOAc
       in hexanes) to afford the product as a colourless oil (1.2 g, 76%).
       Data for the product: 1H NMR (500 MHz, DMSO-d6): delta 8.21 (s, 1H).
   """)

   print(result)
   # CategorizedSolvents(
   #   reaction=['dichloromethane'],
   #   workup=['dcm'],
   #   purification=['hexanes', 'etoac'],
   #   analytical=['dmso-d6']
   # )


How classification works
^^^^^^^^^^^^^^^^^^^^^^^^

Alkahest uses a rule-based pipeline:

1. **Dictionary lookup**: A dictionary of 350+ solvent name variants (IUPAC names,
   common names, abbreviations) is searched using a longest-match-first strategy
   to avoid partial matches.

2. **Phase boundary detection**: Regex patterns identify where the text transitions
   between experimental phases. For example, words like "quenched" or "extracted"
   mark the start of workup; "column" or "chromatography" mark purification;
   "NMR" or "ppm" mark analytical sections.

3. **Context-aware classification**: Each solvent mention is classified based on its
   position relative to detected phase boundaries, combined with a local context
   window (60-80 characters) that catches phrases like "dissolved in" (reaction)
   or "washed with" (workup).

4. **Deuterated solvent handling**: Deuterated solvents (CDCl3, DMSO-d6, etc.) are
   always classified as analytical, regardless of their position in the text.

5. **False positive filtering**: A heuristic filter removes matches that appear
   inside IUPAC compound names (e.g., "methyl" inside "trimethylamine").


Working with SMILES output
^^^^^^^^^^^^^^^^^^^^^^^^^^

Each category in the result maps solvent names to their SMILES. Use
:meth:`~alkahest.CategorizedSolvents.to_smiles_dict` to get dot-separated SMILES
strings per category:

.. code-block:: python

   smiles = result.to_smiles_dict()
   print(smiles["reaction"])
   # "ClCCl"  (DCM)

   # Individual category accessors
   result.get_reaction_smiles()       # list of unique SMILES
   result.get_workup_smiles()
   result.get_purification_smiles()
   result.get_analytical_smiles()

Categories without solvents return ``"solvent-free"`` in the SMILES dict.


Batch Processing
----------------

For DataFrames (e.g., from ORD or USPTO data), use
:func:`~alkahest.extract_solvents_batch`:

.. code-block:: python

   import pandas as pd
   from alkahest import extract_solvents_batch

   df = pd.read_parquet("reactions.parquet")
   df = extract_solvents_batch(df, procedure_column="procedure")

This adds four new columns to the DataFrame:

- ``SOLV_RXN`` — reaction solvent SMILES (dot-separated)
- ``SOLV_WORKUP`` — workup solvent SMILES
- ``SOLV_PURIF`` — purification solvent SMILES
- ``SOLV_ANAL`` — analytical solvent SMILES

You can customize the column prefix:

.. code-block:: python

   df = extract_solvents_batch(df, procedure_column="procedure", prefix="MY_")
   # Creates: MY_RXN, MY_WORKUP, MY_PURIF, MY_ANAL


Using ORD Data with Rxn-INSIGHT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To extract solvents from `Open Reaction Database <https://open-reaction-database.org/>`_
(ORD) files, first process them with
`Rxn-INSIGHT <https://github.com/mrodobbe/Rxn-INSIGHT>`_:

.. code-block:: python

   import rxn_insight as ri

   db = ri.ORDDatabase("dataset.pb.gz")
   df = db.analyze()
   df.to_parquet("analyzed_reactions.parquet")

Then extract solvents from the resulting DataFrame:

.. code-block:: python

   from alkahest import extract_solvents_batch

   df = pd.read_parquet("analyzed_reactions.parquet")
   df = extract_solvents_batch(df, procedure_column="procedure")


GSK Solvent Scoring
-------------------

.. note::

   The GSK scoring data is based on the 2016 GSK Solvent Selection Guide.
   If you use ``score_solvent`` or ``get_gsk_guide``, please cite the original paper:

   Alder, C. M.; Hayler, J. D.; Henderson, R. K.; Redman, A. M.; Shukla, L.;
   Shuster, L. E.; Sneddon, H. F. *Updating and further expanding GSK's solvent
   sustainability guide.* Green Chem. **2016**, 18, 3879--3890.
   `DOI: 10.1039/C6GC00611F <https://doi.org/10.1039/C6GC00611F>`_

Score solvents using the 2016 GSK Solvent Selection Guide with
:func:`~alkahest.score_solvent`:

.. code-block:: python

   from alkahest import score_solvent

   score_solvent("C1CCOC1")    # THF -> "R" (Red: major issues)
   score_solvent("CCOC(C)=O")  # EtOAc -> "G" (Green: few issues)
   score_solvent("ClCCl")      # DCM -> "R" (Red)
   score_solvent("CS(C)=O")    # DMSO -> "A" (Amber: some issues)

Scores:

- **G** (Green): Few issues — recommended
- **A** (Amber): Some issues — usable but consider alternatives
- **R** (Red): Major issues — avoid if possible

Mixtures and deuterated solvents
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For solvent mixtures, pass dot-separated SMILES. The worst score is returned:

.. code-block:: python

   score_solvent("CCOC(C)=O.C1CCOC1")  # EtOAc + THF -> "R" (worst of G and R)

Deuterated solvents are automatically mapped to their non-deuterated parent
compounds before scoring:

.. code-block:: python

   score_solvent("[2H]C(Cl)(Cl)Cl")  # CDCl3 -> scored as CHCl3 -> "R"

Use :func:`~alkahest.get_gsk_guide` to access the full guide as a DataFrame:

.. code-block:: python

   from alkahest import get_gsk_guide

   guide = get_gsk_guide()
   print(guide.head())
   #            Solvent        SMILES Alternative SMILES Canonical SMILES RAG
   # 0            Water             O                NaN                O   G
   # 1      Lactic acid   CC(C(=O)O)O                NaN     CC(O)C(=O)O   G
   # ...

The ``Canonical SMILES`` column is pre-computed with RDKit so that lookups
succeed on canonical input without RDKit installed. When RDKit *is* available,
the query SMILES is canonicalized too, so any valid notation resolves.


NMR Solvent Mapping
-------------------

Map NMR solvent names (including typos and OCR errors) to canonical names and
SMILES with :func:`~alkahest.map_nmr_solvent`:

.. code-block:: python

   from alkahest import map_nmr_solvent

   map_nmr_solvent("CDCl3")       # ('CDCl3', '[2H]C(Cl)(Cl)Cl')
   map_nmr_solvent("dmso-d6")     # ('DMSO-d6', '[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]')
   map_nmr_solvent("cdc13")       # OCR typo: "1" instead of "l" -> still maps to CDCl3
   map_nmr_solvent("meod")        # ('CD3OD', '[2H]C([2H])([2H])O[2H]')
   map_nmr_solvent("unknown")     # None (not in dictionary)

The dictionary contains 600+ variants covering:

- Standard names (CDCl3, DMSO-d6, D2O, etc.)
- Common abbreviations (meod, aceton-d6, etc.)
- OCR errors from digitized patents (cdc13, cdcl 3, dmso d6, etc.)
- Typos and alternative spellings

The function returns ``None`` for unrecognized solvents.


Conventions
-----------

- SMILES strings are stored in RDKit canonical form
- Multiple solvents within a category are joined with ``.`` (dot-separated SMILES)
- Categories without solvents use the string ``"solvent-free"``
- NMR solvent mixtures use ``|`` as separator in curated strings
