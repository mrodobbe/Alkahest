"""
Command-line entry points.

Installing the package provides two commands:

``alkahest-extract``
    Extract and classify solvents from the procedure texts in a Parquet file.

``alkahest-rxn-insight``
    Run Rxn-INSIGHT reaction analysis over a Parquet file.

Both split their input into chunks so that a large corpus can be processed as
an array job on an HPC cluster. Chunk ``i`` of ``n`` covers the rows from
``i * len(df) // n`` up to ``(i + 1) * len(df) // n``, so the chunks tile the
input exactly, whatever the row count.
"""

import argparse
from typing import Optional, Sequence

import pandas as pd

from alkahest.extraction import extract_solvents_batch


def _chunk_bounds(n_rows: int, index: int, n_chunks: int) -> tuple[int, int]:
    """Return the half-open row range covered by chunk ``index`` of ``n_chunks``."""
    if n_chunks < 1:
        raise ValueError(f"Number of chunks must be at least 1, got {n_chunks}")
    if not 0 <= index < n_chunks:
        raise ValueError(
            f"Chunk index must be between 0 and {n_chunks - 1}, got {index}"
        )
    start = index * n_rows // n_chunks
    end = (index + 1) * n_rows // n_chunks
    return start, end


def extract_solvents(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point for ``alkahest-extract``."""
    parser = argparse.ArgumentParser(
        prog="alkahest-extract",
        description="Extract solvents from procedure texts in a Parquet file.",
    )
    parser.add_argument("--input", required=True, help="Path to input Parquet file")
    parser.add_argument(
        "--output_dir", required=True, help="Directory for output Parquet files"
    )
    parser.add_argument(
        "--job_index", type=int, default=0, help="Index of this job (0 to n_chunks-1)"
    )
    parser.add_argument(
        "--n_chunks", type=int, default=1, help="Total number of chunks"
    )
    parser.add_argument(
        "--procedure_column",
        default="procedure",
        help="Name of the procedure text column",
    )
    parser.add_argument(
        "--prefix", default="SOLV_", help="Prefix for the output solvent columns"
    )
    args = parser.parse_args(argv)

    df = pd.read_parquet(args.input)
    start, end = _chunk_bounds(len(df), args.job_index, args.n_chunks)

    subset = df.iloc[start:end].copy()
    subset = extract_solvents_batch(subset, args.procedure_column, prefix=args.prefix)

    destination = f"{args.output_dir}/solvents_{args.job_index}.parquet"
    subset.to_parquet(destination)
    print(f"Processed rows {start}-{end} ({end - start} rows) -> {destination}")
    return 0


def run_rxn_insight(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point for ``alkahest-rxn-insight``."""
    parser = argparse.ArgumentParser(
        prog="alkahest-rxn-insight",
        description="Run Rxn-INSIGHT reaction analysis over a Parquet file.",
    )
    parser.add_argument("--input", required=True, help="Path to input Parquet file")
    parser.add_argument(
        "--output_dir", required=True, help="Directory for output Parquet files"
    )
    parser.add_argument(
        "-i", "--index", type=int, default=0, help="Index of this job (0 to n-1)"
    )
    parser.add_argument(
        "-n", "--num_scripts", type=int, default=1, help="Total number of jobs"
    )
    parser.add_argument(
        "-c", "--num_cores", type=int, default=1, help="Number of parallel cores"
    )
    parser.add_argument(
        "--reaction_column", default="REACTION", help="Name of the reaction column"
    )
    args = parser.parse_args(argv)

    try:
        import rxn_insight as ri
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise SystemExit(
            "Rxn-INSIGHT is required for this command but is not installed. "
            'Install it with:  pip install "alkahest-chem[rxn-insight]"'
        ) from exc

    df = pd.read_parquet(args.input)
    start, end = _chunk_bounds(len(df), args.index, args.num_scripts)

    subset = df.iloc[start:end]
    database = ri.Database()
    analysed = database.create_database_from_df(
        df=subset, reaction_column=args.reaction_column, n_jobs=args.num_cores
    )

    destination = f"{args.output_dir}/analyzed_{args.index}.parquet"
    analysed.to_parquet(destination)
    print(f"Analysed rows {start}-{end} ({end - start} rows) -> {destination}")
    return 0
