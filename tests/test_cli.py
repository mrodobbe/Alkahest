"""Tests for the command-line entry points."""

import pandas as pd
import pytest

from alkahest.cli import _chunk_bounds, extract_solvents, run_rxn_insight

def _has_parquet_engine() -> bool:
    try:
        pd.io.parquet.get_engine("auto")
    except Exception:
        return False
    return True


# The chunking tests are pure arithmetic; only the round trips need an engine.
needs_parquet = pytest.mark.skipif(
    not _has_parquet_engine(), reason="no parquet engine installed"
)

PROCEDURES = [
    "The material was dissolved in THF (20 mL) and stirred for 2 h.",
    "The mixture was quenched with water and extracted with ethyl acetate.",
    "Purification by column chromatography (silica gel, hexanes/EtOAc 4:1).",
    "1H NMR (400 MHz, CDCl3): delta 7.45 (d, 2H).",
    "Dissolved in DMSO and heated to 80 C.",
]


@pytest.fixture
def parquet_input(tmp_path):
    path = tmp_path / "reactions.parquet"
    pd.DataFrame({"procedure": PROCEDURES}).to_parquet(path)
    return path


class TestChunkBounds:
    @pytest.mark.parametrize("n_rows", [0, 1, 5, 7, 1338762])
    @pytest.mark.parametrize("n_chunks", [1, 3, 130])
    def test_chunks_tile_the_input_exactly(self, n_rows, n_chunks):
        bounds = [_chunk_bounds(n_rows, i, n_chunks) for i in range(n_chunks)]
        assert bounds[0][0] == 0
        assert bounds[-1][1] == n_rows
        assert all(bounds[i][1] == bounds[i + 1][0] for i in range(n_chunks - 1))
        assert sum(end - start for start, end in bounds) == n_rows

    def test_more_chunks_than_rows_is_allowed(self):
        # Regression: integer division used to give every chunk a size of zero
        # and hand the whole frame to the last job.
        bounds = [_chunk_bounds(3, i, 10) for i in range(10)]
        assert sum(end - start for start, end in bounds) == 3

    def test_rejects_index_out_of_range(self):
        with pytest.raises(ValueError, match="Chunk index"):
            _chunk_bounds(10, 5, 5)
        with pytest.raises(ValueError, match="Chunk index"):
            _chunk_bounds(10, -1, 5)

    def test_rejects_non_positive_chunk_count(self):
        with pytest.raises(ValueError, match="at least 1"):
            _chunk_bounds(10, 0, 0)


@needs_parquet
class TestExtractSolventsCommand:
    def test_writes_expected_columns(self, parquet_input, tmp_path):
        assert extract_solvents(
            ["--input", str(parquet_input), "--output_dir", str(tmp_path)]
        ) == 0
        result = pd.read_parquet(tmp_path / "solvents_0.parquet")
        assert len(result) == len(PROCEDURES)
        for column in ("SOLV_RXN", "SOLV_WORKUP", "SOLV_PURIF", "SOLV_ANAL"):
            assert column in result.columns
        assert result["SOLV_RXN"].iloc[0] == "C1CCOC1"

    def test_chunks_cover_every_row(self, parquet_input, tmp_path):
        for index in range(3):
            extract_solvents([
                "--input", str(parquet_input), "--output_dir", str(tmp_path),
                "--job_index", str(index), "--n_chunks", "3",
            ])
        chunks = [pd.read_parquet(tmp_path / f"solvents_{i}.parquet") for i in range(3)]
        assert sum(len(c) for c in chunks) == len(PROCEDURES)
        combined = pd.concat(chunks)
        assert combined["procedure"].tolist() == PROCEDURES

    def test_custom_prefix_and_column(self, tmp_path):
        source = tmp_path / "custom.parquet"
        pd.DataFrame({"text": ["Dissolved in DMSO."]}).to_parquet(source)
        extract_solvents([
            "--input", str(source), "--output_dir", str(tmp_path),
            "--procedure_column", "text", "--prefix", "MY_",
        ])
        result = pd.read_parquet(tmp_path / "solvents_0.parquet")
        assert "MY_RXN" in result.columns

    def test_does_not_mutate_the_input_file(self, parquet_input, tmp_path):
        extract_solvents(["--input", str(parquet_input), "--output_dir", str(tmp_path)])
        assert pd.read_parquet(parquet_input).columns.tolist() == ["procedure"]

    def test_rejects_bad_chunk_index(self, parquet_input, tmp_path):
        with pytest.raises(ValueError):
            extract_solvents([
                "--input", str(parquet_input), "--output_dir", str(tmp_path),
                "--job_index", "4", "--n_chunks", "2",
            ])


@needs_parquet
class TestRxnInsightCommand:
    def test_points_at_the_extra_when_rxn_insight_is_missing(
        self, parquet_input, tmp_path, monkeypatch
    ):
        import sys

        monkeypatch.setitem(sys.modules, "rxn_insight", None)
        with pytest.raises(SystemExit, match="rxn-insight"):
            run_rxn_insight(
                ["--input", str(parquet_input), "--output_dir", str(tmp_path)]
            )

    def test_analyses_the_requested_chunk(self, parquet_input, tmp_path, monkeypatch):
        import sys
        import types

        seen = {}

        class FakeDatabase:
            def create_database_from_df(self, df, reaction_column, n_jobs):
                seen["rows"] = len(df)
                seen["reaction_column"] = reaction_column
                seen["n_jobs"] = n_jobs
                return pd.DataFrame({"CLASS": ["acylation"] * len(df)})

        fake = types.ModuleType("rxn_insight")
        fake.Database = FakeDatabase
        monkeypatch.setitem(sys.modules, "rxn_insight", fake)

        assert run_rxn_insight([
            "--input", str(parquet_input), "--output_dir", str(tmp_path),
            "-i", "1", "-n", "2", "-c", "4", "--reaction_column", "procedure",
        ]) == 0

        assert seen == {"rows": 3, "reaction_column": "procedure", "n_jobs": 4}
        written = pd.read_parquet(tmp_path / "analyzed_1.parquet")
        assert len(written) == 3
