"""
HPC batch script for solvent extraction.

Splits a large Parquet DataFrame into chunks and extracts solvents
from procedure texts in parallel across multiple jobs.

Usage:
    python scripts/extract_solvents.py --input data.parquet --output_dir output/ --job_index 0 --n_chunks 130
"""

import argparse

import pandas as pd

from alkahest import extract_solvents_batch


def main():
    parser = argparse.ArgumentParser(
        description="Extract solvents from procedure texts (HPC batch mode)"
    )
    parser.add_argument(
        "--input", type=str, required=True, help="Path to input Parquet file"
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Directory for output Parquet files"
    )
    parser.add_argument(
        "--job_index", type=int, required=True, help="Index of the current job (0 to N-1)"
    )
    parser.add_argument(
        "--n_chunks", type=int, default=130, help="Total number of chunks"
    )
    parser.add_argument(
        "--procedure_column",
        type=str,
        default="procedure",
        help="Name of the procedure text column",
    )
    parser.add_argument(
        "--prefix", type=str, default="SOLV_", help="Prefix for output solvent columns"
    )
    args = parser.parse_args()

    if args.job_index >= args.n_chunks:
        raise ValueError(
            f"Job index {args.job_index} must be less than n_chunks {args.n_chunks}"
        )

    df = pd.read_parquet(args.input)
    total_rows = len(df)
    chunk_size = total_rows // args.n_chunks

    start_idx = args.job_index * chunk_size
    end_idx = start_idx + chunk_size
    if args.job_index == args.n_chunks - 1:
        end_idx = total_rows

    df_subset = df.iloc[start_idx:end_idx]
    df_subset = extract_solvents_batch(
        df_subset, args.procedure_column, prefix=args.prefix
    )
    df_subset.to_parquet(f"{args.output_dir}/solvents_{args.job_index}.parquet")
    print(f"Processed rows {start_idx}-{end_idx} ({end_idx - start_idx} rows)")


if __name__ == "__main__":
    main()
