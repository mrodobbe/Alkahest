"""
GSK Solvent Selection Guide scoring.

Scores solvents as Green (G), Amber (A), or Red (R) based on the 2016 GSK
Solvent Selection Guide. For solvent mixtures, the worst score is used
(e.g., G + A = A).

Deuterated solvents are mapped to their non-deuterated parent compounds
before scoring.

If you use the GSK scoring functionality, please cite the original guide:

    Alder, C. M.; Hayler, J. D.; Henderson, R. K.; Redman, A. M.; Shukla, L.;
    Shuster, L. E.; Sneddon, H. F. "Updating and further expanding GSK's
    solvent sustainability guide." Green Chem. 2016, 18, 3879-3890.
    DOI: 10.1039/C6GC00611F
"""

from importlib import resources
from typing import Optional

import pandas as pd

_GSK_GUIDE: Optional[pd.DataFrame] = None

# Deuterated -> non-deuterated SMILES mapping for GSK scoring
_DEUTERATED_TO_PARENT = {
    "[2H]C(Cl)(Cl)Cl": "ClC(Cl)Cl",  # CDCl3 -> CHCl3
    "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]": "CS(C)=O",  # DMSO-d6 -> DMSO
    "[2H]C([2H])([2H])O[2H]": "CO",  # CD3OD -> MeOH
    "[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]": "CC(C)=O",  # Acetone-d6 -> Acetone
    "[2H]O[2H]": "O",  # D2O -> H2O
    "[2H]C([2H])([2H])C#N": "CC#N",  # CD3CN -> MeCN
    "[2H]C([2H])(Cl)Cl": "ClCCl",  # CD2Cl2 -> DCM
    "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]": "c1ccccc1",  # C6D6 -> Benzene
    "[2H]c1c([2H])c([2H])nc([2H])c1[2H]": "c1ccncc1",  # Pyridine-d5 -> Pyridine
    "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]": "C1CCOC1",  # THF-d8
    "[2H]C([2H])([2H])N(C([2H])([2H])[2H])C(=O)[2H]": "CN(C)C=O",  # DMF-d7
    "[2H]C([2H])([2H])c1c([2H])c([2H])c([2H])c([2H])c1[2H]": "Cc1ccccc1",  # Toluene-d8
}


def get_gsk_guide() -> pd.DataFrame:
    """
    Load the GSK Solvent Selection Guide data.

    If you use this data, please cite:
    Alder et al., Green Chem., 2016, 18, 3879-3890. DOI: 10.1039/C6GC00611F

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: Solvent, SMILES, Alternative SMILES, RAG.
    """
    global _GSK_GUIDE
    if _GSK_GUIDE is None:
        data_path = resources.files("alkahest") / "data" / "rag_gsk.csv"
        # CSV has unquoted commas in solvent names (e.g. "1,3-propanediol").
        # Parse by splitting from the right: RAG is last, Alt SMILES second-to-last,
        # SMILES third-to-last, everything else is the solvent name.
        rows = []
        with open(data_path) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split(",")
                rag = parts[-1]
                alt_smiles = parts[-2] if len(parts) > 3 else ""
                smiles = parts[-3] if len(parts) > 3 else parts[-2]
                name = ",".join(parts[: -3 if len(parts) > 3 else -2])
                rows.append(
                    {"Solvent": name, "SMILES": smiles, "Alternative SMILES": alt_smiles, "RAG": rag}
                )
        _GSK_GUIDE = pd.DataFrame(rows)
    return _GSK_GUIDE


def _build_smiles_to_rag() -> dict[str, str]:
    """Build lookup from SMILES -> RAG score."""
    guide = get_gsk_guide()
    mapping = {}
    for _, row in guide.iterrows():
        smiles = row["SMILES"]
        rag = row["RAG"]
        mapping[smiles] = rag
        if pd.notna(row.get("Alternative SMILES")):
            mapping[row["Alternative SMILES"]] = rag
    return mapping


def score_solvent(smiles: str) -> str:
    """
    Score a solvent SMILES using the GSK Solvent Selection Guide.

    For mixtures (dot-separated SMILES), returns the worst score.
    Deuterated solvents are mapped to their parent compounds.

    If you use this function, please cite:
    Alder et al., Green Chem., 2016, 18, 3879-3890. DOI: 10.1039/C6GC00611F

    Parameters
    ----------
    smiles : str
        Solvent SMILES string. Use dot-separated SMILES for mixtures.

    Returns
    -------
    str
        "G" (green/few issues), "A" (amber/some issues), "R" (red/major issues),
        or "Unknown" if the solvent is not in the guide.
    """
    if not smiles or smiles in ("solvent-free", "not-reported"):
        return "Unknown"

    lookup = _build_smiles_to_rag()
    severity = {"G": 0, "A": 1, "R": 2}

    components = smiles.split(".")
    worst = -1
    all_known = True

    for component in components:
        # Map deuterated to parent
        parent = _DEUTERATED_TO_PARENT.get(component, component)
        rag = lookup.get(parent)
        if rag is None:
            all_known = False
        else:
            worst = max(worst, severity.get(rag, -1))

    if worst == -1:
        return "Unknown"

    inv_severity = {v: k for k, v in severity.items()}
    return inv_severity[worst]
