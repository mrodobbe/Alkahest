"""
Solvent extraction and classification from experimental procedure texts.

Uses a rule-based NLP pipeline with a comprehensive solvent dictionary to
identify solvents and classify them by role (reaction, workup, purification,
analytical) based on phase boundaries and contextual cues.
"""

import re
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from alkahest.solvents import SOLVENT_DICT, DEUTERATED_SOLVENTS


@dataclass
class CategorizedSolvents:
    """Container for solvents classified by experimental role."""

    reaction: dict[str, str] = field(default_factory=dict)
    workup: dict[str, str] = field(default_factory=dict)
    purification: dict[str, str] = field(default_factory=dict)
    analytical: dict[str, str] = field(default_factory=dict)

    def get_reaction_smiles(self) -> list[str]:
        return list(set(self.reaction.values()))

    def get_workup_smiles(self) -> list[str]:
        return list(set(self.workup.values()))

    def get_purification_smiles(self) -> list[str]:
        return list(set(self.purification.values()))

    def get_analytical_smiles(self) -> list[str]:
        return list(set(self.analytical.values()))

    def to_smiles_dict(self) -> dict[str, str]:
        """Convert to dict of category -> dot-separated SMILES string."""
        return {
            "reaction": ".".join(sorted(self.get_reaction_smiles())) or "solvent-free",
            "workup": ".".join(sorted(self.get_workup_smiles())) or "solvent-free",
            "purification": ".".join(sorted(self.get_purification_smiles()))
            or "solvent-free",
            "analytical": ".".join(sorted(self.get_analytical_smiles()))
            or "solvent-free",
        }

    def __repr__(self):
        return (
            f"CategorizedSolvents(\n"
            f"  reaction={list(self.reaction.keys())},\n"
            f"  workup={list(self.workup.keys())},\n"
            f"  purification={list(self.purification.keys())},\n"
            f"  analytical={list(self.analytical.keys())}\n)"
        )


def find_phase_boundaries(text: str) -> dict[str, Optional[int]]:
    """Find where workup/purification/analytical sections start in procedure text."""
    text_lower = text.lower()

    workup_patterns = [
        r"quench",
        r"work.?up",
        r"poured\s+(onto|into)",
        r"layers?\s+were",
    ]
    purif_patterns = [
        r"column",
        r"chromatograph",
        r"silica",
        r"combiflash",
        r"recrystalli",
    ]
    anal_patterns = [r"\bnmr\b", r"δ\s*ppm", r"m/[ez]\s*\d+"]

    def first_match(patterns):
        pos = len(text)
        for p in patterns:
            m = re.search(p, text_lower)
            if m and m.start() < pos:
                pos = m.start()
        return pos if pos < len(text) else None

    return {
        "workup": first_match(workup_patterns),
        "purification": first_match(purif_patterns),
        "analytical": first_match(anal_patterns),
    }


def _is_inside_compound_name(text: str, pos: int) -> bool:
    """Check if a match position is likely inside an IUPAC compound name."""
    before = text[max(0, pos - 80) : pos]
    open_p, close_p = before.rfind("("), before.rfind(")")
    if open_p > close_p:
        if re.search(r"yl\)|phenyl|amino|hydroxy|\d-\w+yl", before[open_p:], re.I):
            return True
    if re.search(r"(\)|yl|ol|ine|ide|ate)-?\s*$", before, re.I):
        return True
    return False


def _classify_solvent(
    text: str, pos: int, name: str, boundaries: dict
) -> str:
    """Classify a solvent mention by position and local context."""
    ctx_before = text[max(0, pos - 60) : pos].lower()
    ctx_after = text[pos : min(len(text), pos + 80)].lower()
    ctx = ctx_before + ctx_after

    # 1. ANALYTICAL: near NMR/MS data
    if re.search(r"(nmr|δ\s*ppm|δ\s*=|\d+\s*hz|m/[ez])", ctx):
        return "analytical"

    # 2. Check for immediate workup indicators
    if re.search(r"(added\s+)?to\s+quench|quench.*" + re.escape(name), ctx_after):
        return "workup"
    if re.search(r"diluted?\s+with|extracted?\s+with|washed?\s+with", ctx_before):
        return "workup"

    # 3. Determine phase by position
    workup_start = boundaries.get("workup")
    purif_start = boundaries.get("purification")
    anal_start = boundaries.get("analytical")

    if anal_start and pos >= anal_start:
        return "analytical"

    if purif_start and pos >= purif_start:
        if re.search(r"column|silica|elut|gradient|\d+%", ctx):
            return "purification"

    if workup_start and pos >= workup_start:
        if re.search(
            r"in\s+anhydrous|anhydrous\s+\w+\s*\(|stirred\s+solution\s+in", ctx
        ):
            return "reaction"
        return "workup"

    # 4. Before any phase boundary = reaction zone
    if re.search(
        r"in\s+anhydrous|dissolved?\s+in|stirred?\s+in|solution\s+of|suspend", ctx
    ):
        return "reaction"

    return "reaction"


def extract_solvents_categorized(text: str) -> CategorizedSolvents:
    """
    Extract and classify solvents from a procedure text.

    Uses longest-match-first strategy over a comprehensive dictionary.
    Deuterated solvents are always classified as analytical.

    Parameters
    ----------
    text : str
        Experimental procedure text.

    Returns
    -------
    CategorizedSolvents
        Solvents categorized by role (reaction, workup, purification, analytical).
    """
    if not isinstance(text, str) or not text:
        return CategorizedSolvents()

    result = CategorizedSolvents()
    text_lower = text.lower()
    boundaries = find_phase_boundaries(text)
    seen = {
        "reaction": set(),
        "workup": set(),
        "purification": set(),
        "analytical": set(),
    }

    # Deuterated solvents (always analytical)
    for name in sorted(DEUTERATED_SOLVENTS.keys(), key=len, reverse=True):
        pattern = r"(?<![a-z-])" + re.escape(name) + r"(?![a-z-])"
        if re.search(pattern, text_lower):
            smiles = DEUTERATED_SOLVENTS[name]
            if smiles not in seen["analytical"]:
                result.analytical[name] = smiles
                seen["analytical"].add(smiles)

    # Regular solvents
    for name in sorted(SOLVENT_DICT.keys(), key=len, reverse=True):
        pattern = r"(?<![a-z-])" + re.escape(name) + r"(?![a-z-])"

        for match in re.finditer(pattern, text_lower):
            pos = match.start()

            if _is_inside_compound_name(text, pos):
                continue

            smiles = SOLVENT_DICT[name]
            category = _classify_solvent(text, pos, name, boundaries)

            if smiles not in seen[category]:
                getattr(result, category)[name] = smiles
                seen[category].add(smiles)

    return result


def _extract_all_solvents_for_row(
    row: pd.Series, procedure_column: str = "procedure"
) -> pd.Series:
    """Extract all 4 solvent categories from a single DataFrame row."""
    text = row.get(procedure_column, "")
    result = extract_solvents_categorized(text)
    smiles = result.to_smiles_dict()
    return pd.Series(
        [smiles["reaction"], smiles["workup"], smiles["purification"], smiles["analytical"]]
    )


def extract_solvents_batch(
    df: pd.DataFrame,
    procedure_column: str = "procedure",
    prefix: str = "SOLV_",
) -> pd.DataFrame:
    """
    Extract all solvent categories for an entire DataFrame.

    Adds 4 new columns: {prefix}RXN, {prefix}WORKUP, {prefix}PURIF, {prefix}ANAL.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a procedure text column.
    procedure_column : str
        Name of the column containing procedure text.
    prefix : str
        Prefix for new solvent columns.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with 4 new solvent columns added.
    """
    from tqdm import tqdm

    tqdm.pandas(desc="Extracting solvents")

    new_cols = df.progress_apply(
        lambda row: _extract_all_solvents_for_row(row, procedure_column),
        axis=1,
    )

    df[f"{prefix}RXN"] = new_cols[0]
    df[f"{prefix}WORKUP"] = new_cols[1]
    df[f"{prefix}PURIF"] = new_cols[2]
    df[f"{prefix}ANAL"] = new_cols[3]

    return df
