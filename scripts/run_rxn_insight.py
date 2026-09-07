"""
HPC batch script for Rxn-INSIGHT reaction analysis.

Kept for backwards compatibility with existing cluster job scripts. Installing
the package also provides this as the ``alkahest-rxn-insight`` command.

Usage:
    python scripts/run_rxn_insight.py --input data.parquet --output_dir output/ -i 0 -n 1000 -c 6
"""

import sys

from alkahest.cli import run_rxn_insight

if __name__ == "__main__":
    sys.exit(run_rxn_insight())
