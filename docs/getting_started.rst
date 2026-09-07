Getting Started
===============

Installation
------------

Install the core package:

.. code-block:: bash

   pip install alkahest-chem

Or install the development version from a local clone:

.. code-block:: bash

   git clone https://github.com/mrodobbe/Alkahest.git
   cd Alkahest
   pip install .

.. warning::

   The distribution is named ``alkahest-chem`` and the import name is
   ``alkahest``. Do **not** run ``pip install alkahest``: that name belongs to
   an unrelated computer algebra system on PyPI.

Optional dependencies
^^^^^^^^^^^^^^^^^^^^^

For RDKit support (SMILES canonicalization):

.. code-block:: bash

   pip install "alkahest-chem[rdkit]"

For reading `Open Reaction Database <https://open-reaction-database.org/>`_ (ORD) files:

.. code-block:: bash

   pip install "alkahest-chem[ord]"

For `Rxn-INSIGHT <https://github.com/mrodobbe/Rxn-INSIGHT>`_ integration
(reaction classification from ORD data):

.. code-block:: bash

   pip install "alkahest-chem[rxn-insight]"

For everything:

.. code-block:: bash

   pip install "alkahest-chem[all]"


Quick Example
-------------

Extract solvents from a procedure text, score them, and map NMR solvents:

.. code-block:: python

   from alkahest import extract_solvents_categorized, score_solvent, map_nmr_solvent

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

   # Score with GSK Solvent Selection Guide
   score_solvent("C1CCOC1")    # THF -> "R" (Red)
   score_solvent("CCOC(C)=O")  # EtOAc -> "G" (Green)

   # Map NMR solvent names (handles typos and OCR errors)
   map_nmr_solvent("cdc13")    # ('CDCl3', '[2H]C(Cl)(Cl)Cl')
