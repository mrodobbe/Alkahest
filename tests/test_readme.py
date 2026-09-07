"""The README is executable documentation, so it is tested like code.

Every documented output in this project was wrong at one point. These tests run
the README's own snippets and compare against the values it advertises, so the
examples cannot drift away from the behaviour again.
"""

import io
import pathlib
import re
import contextlib

import pytest

from alkahest import map_nmr_solvent, score_solvent

README = pathlib.Path(__file__).resolve().parent.parent / "README.md"

pytestmark = pytest.mark.skipif(
    not README.exists(), reason="README.md is not part of the installed package"
)


def _python_blocks() -> list[str]:
    text = README.read_text(encoding="utf-8")
    return re.findall(r"```python\n(.*?)```", text, re.S)


def test_readme_has_examples():
    assert len(_python_blocks()) >= 3


@pytest.mark.parametrize("index", range(4))
def test_every_snippet_runs(index):
    blocks = _python_blocks()
    if index >= len(blocks):
        pytest.skip("no such block")
    block = blocks[index]
    if "read_parquet" in block:
        pytest.skip("snippet needs a data file the repository does not ship")
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(block, f"README block {index}", "exec"), {})


def test_quick_start_prints_what_it_claims():
    block = next(b for b in _python_blocks() if "extract_solvents_categorized" in b)
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        exec(compile(block, "README quick start", "exec"), {})
    printed = buffer.getvalue()

    # The commented-out expectation in the README, line for line.
    assert "reaction=['thf']" in printed
    assert "workup=['ethyl acetate', 'water']" in printed
    assert "purification=['hexanes', 'etoac']" in printed
    assert "analytical=['cdcl3']" in printed
    assert "'reaction': 'C1CCOC1'" in printed
    assert "'workup': 'CCOC(C)=O.O'" in printed
    assert "'purification': 'CCCCCC.CCOC(C)=O'" in printed
    assert "'analytical': '[2H]C(Cl)(Cl)Cl'" in printed


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCOC1", "R"),
        ("CS(C)=O", "A"),
        ("CCOC(C)=O", "G"),
        ("CCOC(C)=O.C1CCOC1", "R"),
        ("[2H]C(Cl)(Cl)Cl", "R"),
        ("CCN(C(C)C)C(C)C", "Unknown"),
    ],
)
def test_documented_scores_are_correct(smiles, expected):
    assert score_solvent(smiles) == expected


@pytest.mark.parametrize(
    "text,expected",
    [("CDCl3", "CDCl3"), ("cdc13", "CDCl3"), ("dmso-d6", "DMSO-d6"), ("meod", "CD3OD")],
)
def test_documented_nmr_lookups_are_correct(text, expected):
    result = map_nmr_solvent(text)
    assert result is not None and result[0] == expected


def test_documented_unknown_nmr_lookup():
    assert map_nmr_solvent("xyzzy") is None


def test_install_instructions_never_name_the_wrong_project():
    # "pip install alkahest" would fetch an unrelated computer algebra system.
    text = README.read_text(encoding="utf-8")
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("pip install alkahest"):
            assert stripped.startswith("pip install alkahest-chem"), line
