"""
HPC batch script for Rxn-INSIGHT reaction analysis.

Splits a large Parquet DataFrame into chunks and runs Rxn-INSIGHT
classification in parallel across multiple jobs.

Usage:
    python scripts/run_rxn_insight.py --input data.parquet --output_dir output/ --index 0 --num_scripts 1000 --num_cores 6
"""

import argparse

import pandas as pd
import rxn_insight as ri


def main():
    parser = argparse.ArgumentParser(description="Run Rxn-INSIGHT analysis (HPC batch)")
    parser.add_argument(
        "--input", type=str, required=True, help="Path to input Parquet file"
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Directory for output Parquet files"
    )
    parser.add_argument(
        "-i", "--index", type=int, required=True, help="Job index (0 to N-1)"
    )
    parser.add_argument(
        "-n", "--num_scripts", type=int, required=True, help="Total number of jobs"
    )
    parser.add_argument(
        "-c", "--num_cores", type=int, required=True, help="Number of parallel cores"
    )
    args = parser.parse_args()

    if not 0 <= args.index <= args.num_scripts - 1:
        raise ValueError(
            f"Index must be between 0 and {args.num_scripts - 1}, got {args.index}"
        )

    df = pd.read_parquet(args.input)

    start_idx = int(args.index * len(df.index) / args.num_scripts)
    end_idx = int((args.index + 1) * len(df.index) / args.num_scripts)

    df = df.iloc[start_idx:end_idx]
    db = ri.Database()
    res = db.create_database_from_df(df=df, reaction_column="REACTION", n_jobs=args.num_cores)
    res.to_parquet(f"{args.output_dir}/analyzed_{args.index}.parquet")


if __name__ == "__main__":
    main()
