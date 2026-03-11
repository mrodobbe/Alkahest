Getting Started
===============

Installation
------------

Install the core package:

.. code-block:: bash

   pip install alkahest

Or install from source:

.. code-block:: bash

   git clone https://github.com/mrodobbe/alkahest.git
   cd alkahest
   pip install .

Optional dependencies
^^^^^^^^^^^^^^^^^^^^^

For RDKit support (SMILES canonicalization):

.. code-block:: bash

   pip install "alkahest[rdkit]"

For reading `Open Reaction Database <https://open-reaction-database.org/>`_ (ORD) files:

.. code-block:: bash

   pip install "alkahest[ord]"

For `Rxn-INSIGHT <https://github.com/mrodobbe/Rxn-INSIGHT>`_ integration
(reaction classification from ORD data):

.. code-block:: bash

   pip install "alkahest[rxn-insight]"

For everything:

.. code-block:: bash

   pip install "alkahest[all]"


Quick Example
-------------

Extract solvents from a procedure text, score them, and map NMR solvents:

.. code-block:: python

   from alkahest import extract_solvents_categorized, score_solvent, map_nmr_solvent

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

   # Score with GSK Solvent Selection Guide
   score_solvent("C1CCOC1")    # THF -> "R" (Red)
   score_solvent("CCOC(C)=O")  # EtOAc -> "G" (Green)

   # Map NMR solvent names (handles typos and OCR errors)
   map_nmr_solvent("cdc13")    # ('CDCl3', '[2H]C(Cl)(Cl)Cl')
