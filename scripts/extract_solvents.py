"""
HPC batch script for solvent extraction.

Kept for backwards compatibility with existing cluster job scripts. Installing
the package also provides this as the ``alkahest-extract`` command.

Usage:
    python scripts/extract_solvents.py --input data.parquet --output_dir output/ --job_index 0 --n_chunks 130
"""

import sys

from alkahest.cli import extract_solvents

if __name__ == "__main__":
    sys.exit(extract_solvents())
