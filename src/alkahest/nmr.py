"""
NMR solvent extraction and curation from procedure texts.

Extracts NMR solvent mentions from experimental procedures, curates raw
snippets by removing noise (frequencies, chemical shifts, metadata),
and maps them to canonical SMILES using a comprehensive dictionary of
600+ variants including typos and OCR errors.
"""

from typing import Optional

from alkahest.nmr_solvents import NMR_SOLVENT_DICTIONARY


def extract_nmr_snippet(procedure: str) -> str:
    """
    Extract the NMR solvent snippet from a procedure text.

    Looks for text immediately following "NMR" and extracts the solvent
    mention using heuristic rules based on common formatting patterns
    in patent procedures (brackets, parentheses, MHz markers).

    Parameters
    ----------
    procedure : str
        Full experimental procedure text.

    Returns
    -------
    str
        Raw solvent snippet, or empty string if not found.
    """
    if not isinstance(procedure, str) or "NMR" not in procedure:
        return ""

    try:
        nmr_part = procedure.split("NMR")[1]
        if "[" in nmr_part[:5]:
            return nmr_part.split("]")[0].split("[")[1].split(" ")[-1]
        elif "MHz" in nmr_part[:10]:
            return nmr_part.split(")")[0].split("(")[1].split(" ")[-1]
        elif "mHz" in nmr_part[:10]:
            return nmr_part.split(")")[0].split("(")[1].split(" ")[-1]
        elif "Conditions" in nmr_part[:20]:
            return ""
        elif "(" in nmr_part[:5] and "]" not in nmr_part[:10]:
            return nmr_part.split(")")[0].split("(")[1]
        elif "MHz" in nmr_part[:20]:
            return nmr_part.split("MHz")[1].split(" ")[-1]
        elif "mHz" in nmr_part[:20]:
            return nmr_part.split("MHz")[1].split(" ")[-1]
        else:
            return ""
    except IndexError:
        return ""


def curate_nmr_snippet(snippet: str) -> str:
    """
    Clean a raw NMR solvent snippet by removing noise.

    Strips out frequency data (MHz/Hz), chemical shift values (ppm),
    temperature references, TMS references, and other metadata that
    commonly co-occur with solvent mentions in procedure texts.

    Parameters
    ----------
    snippet : str
        Raw NMR solvent snippet from extract_nmr_snippet().

    Returns
    -------
    str
        Cleaned solvent string(s), pipe-separated if multiple.
        Empty string if no solvent found.
    """
    if not snippet:
        return ""

    noise_terms = ("ppm", "hz", "mer", "base", "tms", "\u00b0 C")

    def _is_noise(part: str) -> bool:
        lower = part.lower()
        return any(term in lower for term in noise_terms)

    def _is_number(part: str) -> bool:
        try:
            float(part.strip())
            return True
        except ValueError:
            return False

    if "mhz" in snippet.lower() and "," in snippet:
        solvent = []
        for part in snippet.split(","):
            if not _is_noise(part):
                solvent.append(part.strip())
        return "|".join(solvent)

    for sep in (",", ";", "/"):
        if sep in snippet:
            solvent = []
            for part in snippet.split(sep):
                if _is_number(part) or _is_noise(part):
                    continue
                solvent.append(part.strip())
            return "|".join(solvent) if solvent else ""

    if snippet in ("\u03b4", "MHz") or len(snippet) < 3:
        return ""

    return snippet


def map_nmr_solvent(text: str) -> Optional[tuple[str, str]]:
    """
    Map an NMR solvent string to its canonical name and SMILES.

    Uses a dictionary of 600+ variants including common typos,
    OCR errors, and alternative notations found in patent texts.

    Parameters
    ----------
    text : str
        Solvent string (e.g., "CDCl3", "dmso-d6", "chloroform-d").

    Returns
    -------
    tuple[str, str] or None
        (canonical_name, SMILES) if found, None otherwise.

    Examples
    --------
    >>> map_nmr_solvent("CDCl3")
    ('CDCl3', '[2H]C(Cl)(Cl)Cl')
    >>> map_nmr_solvent("dmso-d6")
    ('DMSO-d6', '[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]')
    >>> map_nmr_solvent("cdc13")  # OCR typo: 1 instead of l
    ('CDCl3', '[2H]C(Cl)(Cl)Cl')
    """
    if not text:
        return None
    key = text.lower().strip()
    return NMR_SOLVENT_DICTIONARY.get(key)
