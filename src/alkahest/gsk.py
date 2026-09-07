"""
GSK Solvent Selection Guide scoring.

Scores solvents as Green (G), Amber (A), or Red (R) based on the 2016 GSK
Solvent Selection Guide. For solvent mixtures, the worst score is used
(e.g., G + A = A).

Deuterated solvents are mapped to their non-deuterated parent compounds
before scoring.

Lookups are made on canonical SMILES. The shipped guide carries a
pre-computed ``Canonical SMILES`` column so that lookups succeed without
RDKit installed; when RDKit is available the query SMILES is canonicalized
too, so arbitrary input notations also resolve.

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
_SMILES_TO_RAG: Optional[dict[str, str]] = None

# Severity ordering used to pick the worst component of a mixture.
_SEVERITY = {"G": 0, "A": 1, "R": 2}
_INV_SEVERITY = {v: k for k, v in _SEVERITY.items()}

# Sentinel values that carry no solvent identity.
_NON_SOLVENT = {"", "solvent-free", "not-reported", "NEAT"}

# Deuterated -> non-deuterated SMILES mapping for GSK scoring. Used as the
# fallback when RDKit is not installed; with RDKit, isotope labels are
# stripped generically so any deuterated solvent resolves to its parent.
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
    "[2H]C([2H])([2H])N(C([2H])([2H])[2H])C([2H])=O": "CN(C)C=O",  # DMF-d7 (alt)
    "[2H]C([2H])([2H])c1c([2H])c([2H])c([2H])c([2H])c1[2H]": "Cc1ccccc1",  # Toluene-d8
    "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])OC1([2H])[2H]": "C1COCCO1",  # Dioxane-d8
    "[2H]C([2H])([2H])C(=O)O[2H]": "CC(=O)O",  # Acetic-acid-d4
    "[2H]C([2H])([2H])C([2H])([2H])O[2H]": "CCO",  # EtOD-d6 -> EtOH
    "[2H]C([2H])([2H])[N+](=O)[O-]": "C[N+](=O)[O-]",  # Nitromethane-d3
    "FC(F)(F)C(=O)O[2H]": "O=C(O)C(F)(F)F",  # TFA-d
}


def get_gsk_guide() -> pd.DataFrame:
    """
    Load the GSK Solvent Selection Guide data.

    If you use this data, please cite:
    Alder et al., Green Chem., 2016, 18, 3879-3890. DOI: 10.1039/C6GC00611F

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: Solvent, SMILES, Alternative SMILES,
        Canonical SMILES, RAG.
    """
    global _GSK_GUIDE
    if _GSK_GUIDE is None:
        data_path = resources.files("alkahest") / "data" / "rag_gsk.csv"
        with data_path.open(encoding="utf-8") as f:
            _GSK_GUIDE = pd.read_csv(f)
    return _GSK_GUIDE


def _canonicalize(smiles: str) -> Optional[str]:
    """Canonical SMILES via RDKit, or None if RDKit is absent or parsing fails."""
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def _strip_isotopes(smiles: str) -> Optional[str]:
    """Canonical SMILES with isotope labels removed (CDCl3 -> CHCl3)."""
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def _build_smiles_to_rag() -> dict[str, str]:
    """Build (and cache) the lookup from every known SMILES spelling -> RAG score."""
    global _SMILES_TO_RAG
    if _SMILES_TO_RAG is None:
        guide = get_gsk_guide()
        mapping: dict[str, str] = {}
        columns = ["SMILES", "Alternative SMILES", "Canonical SMILES"]
        for _, row in guide.iterrows():
            rag = row["RAG"]
            for column in columns:
                smiles = row.get(column)
                if pd.notna(smiles) and smiles:
                    mapping[smiles] = rag
        _SMILES_TO_RAG = mapping
    return _SMILES_TO_RAG


def _score_component(component: str, lookup: dict[str, str]) -> Optional[str]:
    """Resolve one SMILES component to a RAG score, or None if unknown."""
    rag = lookup.get(component)
    if rag is not None:
        return rag

    parent = _DEUTERATED_TO_PARENT.get(component)
    if parent is not None and parent in lookup:
        return lookup[parent]

    canonical = _canonicalize(component)
    if canonical is not None and canonical in lookup:
        return lookup[canonical]

    stripped = _strip_isotopes(component)
    if stripped is not None and stripped in lookup:
        return lookup[stripped]

    return None


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
    if not smiles or smiles in _NON_SOLVENT:
        return "Unknown"

    lookup = _build_smiles_to_rag()
    worst = -1

    for component in smiles.split("."):
        if not component or component in _NON_SOLVENT:
            continue
        rag = _score_component(component, lookup)
        if rag is not None:
            worst = max(worst, _SEVERITY.get(rag, -1))

    if worst == -1:
        return "Unknown"

    return _INV_SEVERITY[worst]
