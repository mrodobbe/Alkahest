HPC Batch Processing
====================

For processing large datasets (e.g., the full USPTO corpus) on an HPC cluster,
two batch scripts are provided in the ``scripts/`` directory.


Solvent Extraction
------------------

.. code-block:: bash

   python scripts/extract_solvents.py \
       --input reactions.parquet \
       --output_dir output/ \
       --job_index 0 \
       --n_chunks 130 \
       --procedure_column procedure \
       --prefix SOLV_

Arguments:

- ``--input``: Path to input Parquet file
- ``--output_dir``: Directory for output Parquet chunks
- ``--job_index``: Zero-based index of the current chunk to process
- ``--n_chunks``: Total number of chunks to split into
- ``--procedure_column``: Name of the text column (default: ``procedure``)
- ``--prefix``: Column prefix for output (default: ``SOLV_``)

Each job processes one chunk and writes a separate Parquet file. Submit as an
array job on your cluster scheduler.


Rxn-INSIGHT Analysis
--------------------

.. code-block:: bash

   python scripts/run_rxn_insight.py \
       --input data.parquet \
       --output_dir output/ \
       -i 0 \
       -n 1000 \
       -c 6

Arguments:

- ``--input``: Path to input Parquet file
- ``--output_dir``: Directory for output files
- ``-i``: Zero-based index of the current chunk
- ``-n``: Total number of chunks
- ``-c``: Number of CPU cores

This command needs the Rxn-INSIGHT extra:
``pip install "alkahest-chem[rxn-insight]"``.


Example SLURM Submission
------------------------

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=solvent_extract
   #SBATCH --array=0-129
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=4G
   #SBATCH --time=02:00:00

   python scripts/extract_solvents.py \
       --input data/reactions.parquet \
       --output_dir output/ \
       --job_index $SLURM_ARRAY_TASK_ID \
       --n_chunks 130
