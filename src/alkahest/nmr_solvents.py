"""
Comprehensive NMR Solvent Dictionary
Based on extraction from 1M+ chemical procedures

This dictionary maps all known NMR solvent variants (including typos, OCR errors,
and alternative notations) to their canonical SMILES representations.

Structure:
- Key: lowercase solvent string as found in text
- Value: tuple of (canonical_name, SMILES)
"""

# SMILES for common NMR solvents
SMILES = {
    # Deuterated solvents
    "CDCl3": "[2H]C(Cl)(Cl)Cl",
    "DMSO-d6": "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]",
    "CD3OD": "[2H]C([2H])([2H])O[2H]",
    "D2O": "[2H]O[2H]",
    "Acetone-d6": "[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]",
    "Benzene-d6": "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]",
    "CD3CN": "[2H]C([2H])([2H])C#N",
    "THF-d8": "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]",
    "Toluene-d8": "[2H]C([2H])([2H])c1c([2H])c([2H])c([2H])c([2H])c1[2H]",
    "CD2Cl2": "[2H]C([2H])(Cl)Cl",
    "Pyridine-d5": "[2H]c1c([2H])c([2H])nc([2H])c1[2H]",
    "Acetic-acid-d4": "[2H]C([2H])([2H])C(=O)O[2H]",
    "DMF-d7": "[2H]C([2H])([2H])N(C([2H])([2H])[2H])C([2H])=O",
    "Dioxane-d8": "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])OC1([2H])[2H]",
    "TFA-d": "FC(F)(F)C(=O)O[2H]",
    "EtOD-d6": "[2H]C([2H])([2H])C([2H])([2H])O[2H]",
    "Nitromethane-d3": "[2H]C([2H])([2H])[N+](=O)[O-]",
    
    # Non-deuterated solvents
    "DMSO": "CS(=O)C",
    "MeOH": "CO",
    "EtOH": "CCO",
    "Acetone": "CC(=O)C",
    "CHCl3": "ClC(Cl)Cl",
    "H2O": "O",
    "Acetic-acid": "CC(=O)O",
    "DMF": "CN(C)C=O",
    "THF": "C1CCOC1",
    "Acetonitrile": "CC#N",
    "Pyridine": "c1ccncc1",
    "Benzene": "c1ccccc1",
    "Toluene": "Cc1ccccc1",
    "CCl4": "ClC(Cl)(Cl)Cl",
    "Dioxane": "C1COCCO1",
    "TFA": "FC(F)(F)C(=O)O",
    
    # References
    "TMS": "[Si](C)(C)(C)C",
    "neat": "NEAT",
}


# Complete solvent dictionary mapping variants to (canonical_name, SMILES)
NMR_SOLVENT_DICTIONARY = {
    # =========================================================================
    # CDCl3 (Chloroform-d) - Most common NMR solvent
    # =========================================================================
    "cdcl3": ("CDCl3", SMILES["CDCl3"]),
    "chloroform-d": ("CDCl3", SMILES["CDCl3"]),
    "chloroform-d1": ("CDCl3", SMILES["CDCl3"]),
    "deuterochloroform": ("CDCl3", SMILES["CDCl3"]),
    "deuteriochloroform": ("CDCl3", SMILES["CDCl3"]),
    "deuterated chloroform": ("CDCl3", SMILES["CDCl3"]),
    "d-chloroform": ("CDCl3", SMILES["CDCl3"]),
    "d-chroloform": ("CDCl3", SMILES["CDCl3"]),  # typo
    "dr-chloroform": ("CDCl3", SMILES["CDCl3"]),  # typo
    "d1-chloroform": ("CDCl3", SMILES["CDCl3"]),
    "choroform-d": ("CDCl3", SMILES["CDCl3"]),  # typo
    "chloroform-d3": ("CDCl3", SMILES["CDCl3"]),  # wrong but appears
    "chloroform-d6": ("CDCl3", SMILES["CDCl3"]),  # wrong but appears
    "chloroform-a": ("CDCl3", SMILES["CDCl3"]),  # typo
    "chloroform-": ("CDCl3", SMILES["CDCl3"]),  # truncated
    "chloroform": ("CHCl3", SMILES["CHCl3"]),  # non-deuterated
    "heavy chloroform": ("CDCl3", SMILES["CDCl3"]),
    # Common typos/OCR errors
    "cdc13": ("CDCl3", SMILES["CDCl3"]),  # 1 instead of l
    "cdci3": ("CDCl3", SMILES["CDCl3"]),  # i instead of l
    "cdcl 3": ("CDCl3", SMILES["CDCl3"]),  # space
    "cdcl-3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl_3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3-d": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3-d1": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3-d3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3-d6": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3-d8": ("CDCl3", SMILES["CDCl3"]),
    "d-cdcl3": ("CDCl3", SMILES["CDCl3"]),
    "d1-cdcl3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl-d": ("CDCl3", SMILES["CDCl3"]),
    "cdcl": ("CDCl3", SMILES["CDCl3"]),  # truncated
    "cdc": ("CDCl3", SMILES["CDCl3"]),  # truncated
    "cdcl1": ("CDCl3", SMILES["CDCl3"]),
    "cdcl2": ("CDCl3", SMILES["CDCl3"]),
    "cdcl4": ("CDCl3", SMILES["CDCl3"]),
    "cdcl6": ("CDCl3", SMILES["CDCl3"]),
    "cdcl9": ("CDCl3", SMILES["CDCl3"]),
    "cdcl13": ("CDCl3", SMILES["CDCl3"]),
    "cdc3": ("CDCl3", SMILES["CDCl3"]),
    "cl3cd": ("CDCl3", SMILES["CDCl3"]),
    "dccl3": ("CDCl3", SMILES["CDCl3"]),
    "dcdl3": ("CDCl3", SMILES["CDCl3"]),
    "dcl3": ("CDCl3", SMILES["CDCl3"]),
    "cdl3": ("CDCl3", SMILES["CDCl3"]),
    "cdci": ("CDCl3", SMILES["CDCl3"]),
    "cdcls": ("CDCl3", SMILES["CDCl3"]),
    "cdcla": ("CDCl3", SMILES["CDCl3"]),
    "cdclb": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3l": ("CDCl3", SMILES["CDCl3"]),
    "cdcl31": ("CDCl3", SMILES["CDCl3"]),
    "cdcl32": ("CDCl3", SMILES["CDCl3"]),
    "cdcb": ("CDCl3", SMILES["CDCl3"]),
    "cdch": ("CDCl3", SMILES["CDCl3"]),
    "cdco3": ("CDCl3", SMILES["CDCl3"]),
    "cdcq3": ("CDCl3", SMILES["CDCl3"]),
    "cdct3": ("CDCl3", SMILES["CDCl3"]),
    "cdcfe": ("CDCl3", SMILES["CDCl3"]),
    "cdcp3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl_": ("CDCl3", SMILES["CDCl3"]),
    "cdccl3": ("CDCl3", SMILES["CDCl3"]),
    "cddl3": ("CDCl3", SMILES["CDCl3"]),
    "cecl3": ("CDCl3", SMILES["CDCl3"]),
    "cldcl3": ("CDCl3", SMILES["CDCl3"]),
    "clcl3": ("CDCl3", SMILES["CDCl3"]),
    "ccdl3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl33": ("CDCl3", SMILES["CDCl3"]),
    "cd13": ("CDCl3", SMILES["CDCl3"]),
    "cdi3": ("CDCl3", SMILES["CDCl3"]),
    "cocl3": ("CDCl3", SMILES["CDCl3"]),
    "ccdcl3": ("CDCl3", SMILES["CDCl3"]),
    "cdcdl3": ("CDCl3", SMILES["CDCl3"]),
    "dcc13": ("CDCl3", SMILES["CDCl3"]),
    "3cdcl": ("CDCl3", SMILES["CDCl3"]),
    ".cdcl3": ("CDCl3", SMILES["CDCl3"]),
    "ccl3d": ("CDCl3", SMILES["CDCl3"]),
    "cd1c3": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3.": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3'": ("CDCl3", SMILES["CDCl3"]),
    "<cdcl3>": ("CDCl3", SMILES["CDCl3"]),
    "chcl3-d": ("CDCl3", SMILES["CDCl3"]),
    "chcl3-d1": ("CDCl3", SMILES["CDCl3"]),
    "chcl3-d6": ("CDCl3", SMILES["CDCl3"]),
    "chcl-d": ("CDCl3", SMILES["CDCl3"]),
    "d-chcl3": ("CDCl3", SMILES["CDCl3"]),
    "cdl3od": ("CDCl3", SMILES["CDCl3"]),
    "cdcl3□": ("CDCl3", SMILES["CDCl3"]),
    
    # =========================================================================
    # DMSO-d6 (Dimethyl sulfoxide-d6)
    # =========================================================================
    "dmso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmsod6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6 dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "(d6)dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "[d6]dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d-6-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-d6-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d-6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d--6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d—6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso−d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso—d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso.d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso_d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dimethylsulfoxide-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dimethylsulfoxide-de": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dimethyl sulfoxide-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dimethyl-d6 sulfoxide": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dimethysulfoxide-d6": ("DMSO-d6", SMILES["DMSO-d6"]),  # typo
    "dimethyosulfoxide-d6": ("DMSO-d6", SMILES["DMSO-d6"]),  # typo
    "diemthylsulfoxide-d6": ("DMSO-d6", SMILES["DMSO-d6"]),  # typo
    "deuterodimethylsulfoxide": ("DMSO-d6", SMILES["DMSO-d6"]),
    "deuteriomethylsulfoxide": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dimethylsulphoxide": ("DMSO", SMILES["DMSO"]),
    "dimethylsulfoxide": ("DMSO", SMILES["DMSO"]),
    "me2so-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "me2so--d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "me2sod6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-me2so": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6me2so": ("DMSO-d6", SMILES["DMSO-d6"]),
    "(cd3)2so": ("DMSO-d6", SMILES["DMSO-d6"]),
    "cd3socd3": ("DMSO-d6", SMILES["DMSO-d6"]),
    "cd3so2cd3": ("DMSO-d6", SMILES["DMSO-d6"]),
    "c2d6so": ("DMSO-d6", SMILES["DMSO-d6"]),
    "c2d6os": ("DMSO-d6", SMILES["DMSO-d6"]),
    "sulfoxide-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    # Typos
    "dmso-do": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d0": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d1": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d3": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d4": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d5": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d8": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-db": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dc": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-da": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dd": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-de": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dg": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dh": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dr": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-ds": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dv": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d$": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d□": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d.6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-<d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-dd6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d61": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d63": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d66": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-d51": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-δ": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-δ6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-β6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-17": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-26": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-66": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-4": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-c": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-h6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso3": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso46": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmsod": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmsod-6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmsn": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmxo-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmao-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmos-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmbo-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmco-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmiso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmmso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dnmso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmoso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "pmso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "omso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "cmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dazo-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmo-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dsmo-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dsmo": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dsmo-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "mso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "mso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dm50-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dms0-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dms0": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dms-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dms": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dhso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dhso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dhsod6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dm-so-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dcmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "cdmso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-dms": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-dsmo": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d3-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d4-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d8-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d0-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "db-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "de-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "do-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d-6-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-dmso3": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-dmso-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "d6-cd3od": ("CD3OD", SMILES["CD3OD"]),  # mixed notation
    "dmdo": ("DMSO-d6", SMILES["DMSO-d6"]),  # typo
    "dmso--d6": ("DMSO-d6", SMILES["DMSO-d6"]),  # double dash
    "dmso-d": ("DMSO-d6", SMILES["DMSO-d6"]),  # truncated
    "me2so": ("DMSO", SMILES["DMSO"]),  # non-deuterated
    "d6mdso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "ddmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "6-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "66-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "4-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "<dmso>": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso.": ("DMSO-d6", SMILES["DMSO-d6"]),
    "[d6-dmso": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-6d": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmso-4-d6": ("DMSO-d6", SMILES["DMSO-d6"]),
    "dmsod6+cd3cood": ("DMSO-d6", SMILES["DMSO-d6"]),  # mixture
    
    # =========================================================================
    # CD3OD (Methanol-d4)
    # =========================================================================
    "cd3od": ("CD3OD", SMILES["CD3OD"]),
    "methanol-d4": ("CD3OD", SMILES["CD3OD"]),
    "methanol d4": ("CD3OD", SMILES["CD3OD"]),
    "methanol-d3": ("CD3OD", SMILES["CD3OD"]),
    "methanol-d6": ("CD3OD", SMILES["CD3OD"]),
    "methanol-d": ("CD3OD", SMILES["CD3OD"]),
    "methanol-d4": ("CD3OD", SMILES["CD3OD"]),
    "metanol-d4": ("CD3OD", SMILES["CD3OD"]),  # typo
    "deuterated methanol": ("CD3OD", SMILES["CD3OD"]),
    "deuteromethanol": ("CD3OD", SMILES["CD3OD"]),
    "tetradeuteriomethanol": ("CD3OD", SMILES["CD3OD"]),
    "d4-methanol": ("CD3OD", SMILES["CD3OD"]),
    "d4-meoh": ("CD3OD", SMILES["CD3OD"]),
    "d4-meod": ("CD3OD", SMILES["CD3OD"]),
    "d4-cd3od": ("CD3OD", SMILES["CD3OD"]),
    "d-methanol": ("CD3OD", SMILES["CD3OD"]),
    "d1-methanol": ("CD3OD", SMILES["CD3OD"]),
    "meod": ("CD3OD", SMILES["CD3OD"]),
    "meod-d4": ("CD3OD", SMILES["CD3OD"]),
    "meod-d3": ("CD3OD", SMILES["CD3OD"]),
    "meod-d6": ("CD3OD", SMILES["CD3OD"]),
    "meod.": ("CD3OD", SMILES["CD3OD"]),
    "meod3": ("CD3OD", SMILES["CD3OD"]),
    "meod4": ("CD3OD", SMILES["CD3OD"]),
    "meodd4": ("CD3OD", SMILES["CD3OD"]),
    "meoh-d4": ("CD3OD", SMILES["CD3OD"]),
    "meoh-d3": ("CD3OD", SMILES["CD3OD"]),
    "meoh-d6": ("CD3OD", SMILES["CD3OD"]),
    "meoh-d": ("CD3OD", SMILES["CD3OD"]),
    "meoh--d4": ("CD3OD", SMILES["CD3OD"]),
    "meoh—d4": ("CD3OD", SMILES["CD3OD"]),
    "meoh-4": ("CD3OD", SMILES["CD3OD"]),
    "meoh4": ("CD3OD", SMILES["CD3OD"]),
    "meohd": ("CD3OD", SMILES["CD3OD"]),
    "meohd4": ("CD3OD", SMILES["CD3OD"]),
    "d4meoh": ("CD3OD", SMILES["CD3OD"]),
    "d4-meoh": ("CD3OD", SMILES["CD3OD"]),
    "d6-meoh": ("CD3OD", SMILES["CD3OD"]),
    "meoh[d4": ("CD3OD", SMILES["CD3OD"]),
    "meoh-d4-d4": ("CD3OD", SMILES["CD3OD"]),
    "cd3od-d4": ("CD3OD", SMILES["CD3OD"]),
    "cd3od-d6": ("CD3OD", SMILES["CD3OD"]),
    "cd3od-do": ("CD3OD", SMILES["CD3OD"]),
    "cd3od.": ("CD3OD", SMILES["CD3OD"]),
    "cd3od3": ("CD3OD", SMILES["CD3OD"]),
    "cd-3od": ("CD3OD", SMILES["CD3OD"]),
    "cd3-od": ("CD3OD", SMILES["CD3OD"]),
    "cd30d": ("CD3OD", SMILES["CD3OD"]),  # zero instead of O
    "cd300d": ("CD3OD", SMILES["CD3OD"]),
    "cd30od": ("CD3OD", SMILES["CD3OD"]),
    "cd3o": ("CD3OD", SMILES["CD3OD"]),
    "cd3oh": ("CD3OD", SMILES["CD3OD"]),
    "cd3oc": ("CD3OD", SMILES["CD3OD"]),
    "cd3oi": ("CD3OD", SMILES["CD3OD"]),
    "cd3oδ": ("CD3OD", SMILES["CD3OD"]),
    "cdod": ("CD3OD", SMILES["CD3OD"]),
    "cdod3": ("CD3OD", SMILES["CD3OD"]),
    "cd2od": ("CD3OD", SMILES["CD3OD"]),
    "cd4od": ("CD3OD", SMILES["CD3OD"]),
    "cd6od": ("CD3OD", SMILES["CD3OD"]),
    "c3od": ("CD3OD", SMILES["CD3OD"]),
    "co3od": ("CD3OD", SMILES["CD3OD"]),
    "c03od": ("CD3OD", SMILES["CD3OD"]),
    "d3cod": ("CD3OD", SMILES["CD3OD"]),
    "d3cocd3": ("Acetone-d6", SMILES["Acetone-d6"]),  # actually acetone
    "d3cood": ("CD3OD", SMILES["CD3OD"]),
    "ch3od": ("CD3OD", SMILES["CD3OD"]),
    "ch30d": ("CD3OD", SMILES["CD3OD"]),
    "ch3oh-d4": ("CD3OD", SMILES["CD3OD"]),
    "ch3oh-d6": ("CD3OD", SMILES["CD3OD"]),
    "ch3oh-d": ("CD3OD", SMILES["CD3OD"]),
    "me-d3-od": ("CD3OD", SMILES["CD3OD"]),
    "d3-meod": ("CD3OD", SMILES["CD3OD"]),
    "d-meod": ("CD3OD", SMILES["CD3OD"]),
    "d-meoh": ("CD3OD", SMILES["CD3OD"]),
    "<cd3od>": ("CD3OD", SMILES["CD3OD"]),
    "methanol-δ4": ("CD3OD", SMILES["CD3OD"]),
    "methonal-d4": ("CD3OD", SMILES["CD3OD"]),  # typo
    "od-d4": ("CD3OD", SMILES["CD3OD"]),
    "cd6so": ("DMSO-d6", SMILES["DMSO-d6"]),  # actually DMSO
    "cd6oo": ("CD3OD", SMILES["CD3OD"]),
    
    # =========================================================================
    # D2O (Deuterium oxide)
    # =========================================================================
    "d2o": ("D2O", SMILES["D2O"]),
    "d20": ("D2O", SMILES["D2O"]),  # zero instead of O
    "d2q": ("D2O", SMILES["D2O"]),
    "d3o": ("D2O", SMILES["D2O"]),
    "d6o": ("D2O", SMILES["D2O"]),
    "d2o.": ("D2O", SMILES["D2O"]),
    "deuterium oxide": ("D2O", SMILES["D2O"]),
    "d2o-d2": ("D2O", SMILES["D2O"]),
    
    # =========================================================================
    # Acetone-d6
    # =========================================================================
    "acetone-d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "acetone d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "acetoned-d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "acetoned6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "acetoned": ("Acetone-d6", SMILES["Acetone-d6"]),
    "aceton-d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "aceton": ("Acetone", SMILES["Acetone"]),  # non-deuterated
    "acetone-d": ("Acetone-d6", SMILES["Acetone-d6"]),
    "acetone-d5": ("Acetone-d6", SMILES["Acetone-d6"]),
    "acetone-": ("Acetone-d6", SMILES["Acetone-d6"]),
    "d6-acetone": ("Acetone-d6", SMILES["Acetone-d6"]),
    "do-acetone": ("Acetone-d6", SMILES["Acetone-d6"]),
    "d5-acetone": ("Acetone-d6", SMILES["Acetone-d6"]),
    "deuterioacetone": ("Acetone-d6", SMILES["Acetone-d6"]),
    "(cd3)2co": ("Acetone-d6", SMILES["Acetone-d6"]),
    "cd3cocd3": ("Acetone-d6", SMILES["Acetone-d6"]),
    "cd3cocd3.": ("Acetone-d6", SMILES["Acetone-d6"]),
    "c3d6o": ("Acetone-d6", SMILES["Acetone-d6"]),
    "c2d6o": ("Acetone-d6", SMILES["Acetone-d6"]),
    "me2co-d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "me2co--d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "me2co": ("Acetone", SMILES["Acetone"]),
    "cd3co": ("Acetone-d6", SMILES["Acetone-d6"]),
    "cocd3": ("Acetone-d6", SMILES["Acetone-d6"]),
    "co-d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "ace-d6": ("Acetone-d6", SMILES["Acetone-d6"]),
    "4-acetone": ("Acetone-d6", SMILES["Acetone-d6"]),
    "cd6co": ("Acetone-d6", SMILES["Acetone-d6"]),
    "cd3c(o": ("Acetone-d6", SMILES["Acetone-d6"]),
    
    # =========================================================================
    # Benzene-d6
    # =========================================================================
    "benzene-d6": ("Benzene-d6", SMILES["Benzene-d6"]),
    "benzene d6": ("Benzene-d6", SMILES["Benzene-d6"]),
    "c6d6": ("Benzene-d6", SMILES["Benzene-d6"]),
    "d6-benzene": ("Benzene-d6", SMILES["Benzene-d6"]),
    "6d6": ("Benzene-d6", SMILES["Benzene-d6"]),
    "c6d": ("Benzene-d6", SMILES["Benzene-d6"]),
    
    # =========================================================================
    # Acetonitrile/CD3CN
    # =========================================================================
    "cd3cn": ("CD3CN", SMILES["CD3CN"]),
    "cd3cn-d3": ("CD3CN", SMILES["CD3CN"]),
    "acetonitrile-d3": ("CD3CN", SMILES["CD3CN"]),
    "acetonitrile-d6": ("CD3CN", SMILES["CD3CN"]),
    "acetonitrile": ("Acetonitrile", SMILES["Acetonitrile"]),
    "mecn": ("Acetonitrile", SMILES["Acetonitrile"]),
    "mecn-d3": ("CD3CN", SMILES["CD3CN"]),
    "ch3cn": ("Acetonitrile", SMILES["Acetonitrile"]),
    "ch3cn-d3": ("CD3CN", SMILES["CD3CN"]),
    "d3-ch3cn": ("CD3CN", SMILES["CD3CN"]),
    "cn-d3": ("CD3CN", SMILES["CD3CN"]),
    "cd3c1": ("CD3CN", SMILES["CD3CN"]),
    "cd3c3": ("CD3CN", SMILES["CD3CN"]),
    "cd3cl3": ("CDCl3", SMILES["CDCl3"]),  # actually CDCl3
    
    # =========================================================================
    # THF-d8 (Tetrahydrofuran-d8)
    # =========================================================================
    "thf-d8": ("THF-d8", SMILES["THF-d8"]),
    "d8-thf": ("THF-d8", SMILES["THF-d8"]),
    "tetrahydrofuran-d8": ("THF-d8", SMILES["THF-d8"]),
    "d8-tetrahydrofuran": ("THF-d8", SMILES["THF-d8"]),
    "d-thf": ("THF-d8", SMILES["THF-d8"]),
    "thf": ("THF", SMILES["THF"]),
    
    # =========================================================================
    # Toluene-d8
    # =========================================================================
    "toluene-d8": ("Toluene-d8", SMILES["Toluene-d8"]),
    "d8-toluene": ("Toluene-d8", SMILES["Toluene-d8"]),
    "toluene": ("Toluene", SMILES["Toluene"]),
    "c6d5cd3": ("Toluene-d8", SMILES["Toluene-d8"]),
    "c7d8": ("Toluene-d8", SMILES["Toluene-d8"]),
    "tol": ("Toluene", SMILES["Toluene"]),
    
    # =========================================================================
    # Dichloromethane-d2 (CD2Cl2)
    # =========================================================================
    "cd2cl2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "cd2cl2-d2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "cd2cl2-d4": ("CD2Cl2", SMILES["CD2Cl2"]),
    "dichloromethane-d2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "methylene chloride-d2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "chloride-d2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "dcm-d2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "cd2-cl2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "cd2cl1": ("CD2Cl2", SMILES["CD2Cl2"]),
    "ccl2d2": ("CD2Cl2", SMILES["CD2Cl2"]),
    "ch2cl2": ("CH2Cl2", "ClCCl"),  # non-deuterated
    "dcm": ("CH2Cl2", "ClCCl"),
    
    # =========================================================================
    # Pyridine-d5
    # =========================================================================
    "pyridine-d5": ("Pyridine-d5", SMILES["Pyridine-d5"]),
    "pyridine d5": ("Pyridine-d5", SMILES["Pyridine-d5"]),
    "pyridine-d5.": ("Pyridine-d5", SMILES["Pyridine-d5"]),
    "c5d5n": ("Pyridine-d5", SMILES["Pyridine-d5"]),
    "c5d5": ("Pyridine-d5", SMILES["Pyridine-d5"]),
    "pyr-d5": ("Pyridine-d5", SMILES["Pyridine-d5"]),
    "pyr": ("Pyridine", SMILES["Pyridine"]),
    "pyridine": ("Pyridine", SMILES["Pyridine"]),
    
    # =========================================================================
    # DMF-d7
    # =========================================================================
    "dmf-d7": ("DMF-d7", SMILES["DMF-d7"]),
    "dmf-d6": ("DMF-d7", SMILES["DMF-d7"]),
    "dmf-da": ("DMF-d7", SMILES["DMF-d7"]),
    "d7-dmf": ("DMF-d7", SMILES["DMF-d7"]),
    "d-dmf": ("DMF-d7", SMILES["DMF-d7"]),
    "dimethylformamide-d7": ("DMF-d7", SMILES["DMF-d7"]),
    "n-dimethylformamide": ("DMF", SMILES["DMF"]),
    "dmf": ("DMF", SMILES["DMF"]),
    
    # =========================================================================
    # Acetic acid-d4
    # =========================================================================
    "acetic acid-d4": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "acetic-acid-d4": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "acid-d4": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "cd3cood": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "cd3co2d": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "cd3cod": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "acod": ("Acetic-acid-d4", SMILES["Acetic-acid-d4"]),
    "acetic acid": ("Acetic-acid", SMILES["Acetic-acid"]),
    "acetic": ("Acetic-acid", SMILES["Acetic-acid"]),
    "acoh": ("Acetic-acid", SMILES["Acetic-acid"]),
    "ch3cood": ("Acetic-acid", SMILES["Acetic-acid"]),
    
    # =========================================================================
    # Dioxane-d8
    # =========================================================================
    "dioxane-d8": ("Dioxane-d8", SMILES["Dioxane-d8"]),
    "d8-dioxane": ("Dioxane-d8", SMILES["Dioxane-d8"]),
    "dioxane": ("Dioxane", SMILES["Dioxane"]),
    "1,4-dioxane": ("Dioxane", SMILES["Dioxane"]),
    
    # =========================================================================
    # TFA / Trifluoroacetic acid-d
    # =========================================================================
    "tfa-d": ("TFA-d", SMILES["TFA-d"]),
    "tfa-d1": ("TFA-d", SMILES["TFA-d"]),
    "tfa-dl": ("TFA-d", SMILES["TFA-d"]),
    "tfad": ("TFA-d", SMILES["TFA-d"]),
    "d1-tfa": ("TFA-d", SMILES["TFA-d"]),
    "dtfa": ("TFA-d", SMILES["TFA-d"]),
    "cf3cood": ("TFA-d", SMILES["TFA-d"]),
    "cf3co2d": ("TFA-d", SMILES["TFA-d"]),
    "trifluoroacetic": ("TFA", SMILES["TFA"]),
    "tfa": ("TFA", SMILES["TFA"]),
    "cf3cooh": ("TFA", SMILES["TFA"]),
    "cf3co2h": ("TFA", SMILES["TFA"]),
    "tfd": ("TFA-d", SMILES["TFA-d"]),
    "cfcood": ("TFA-d", SMILES["TFA-d"]),
    "d-trifluoroacetic": ("TFA-d", SMILES["TFA-d"]),
    
    # =========================================================================
    # Other deuterated solvents
    # =========================================================================
    # Ethanol-d6
    "etod-d6": ("EtOD-d6", SMILES["EtOD-d6"]),
    "cd3cd2od": ("EtOD-d6", SMILES["EtOD-d6"]),
    "ethanol-d6": ("EtOD-d6", SMILES["EtOD-d6"]),
    
    # Nitromethane-d3
    "nitromethane-d3": ("Nitromethane-d3", SMILES["Nitromethane-d3"]),
    
    # Glyme-d10
    "glyme-d10": ("Glyme-d10", "[2H]C([2H])([2H])OC([2H])([2H])C([2H])([2H])OC([2H])([2H])[2H]"),
    
    # Freon
    "freon": ("Freon", "ClC(F)(Cl)F"),
    "freon®": ("Freon", "ClC(F)(Cl)F"),
    "cfcl3": ("CFCl3", "ClC(Cl)(Cl)F"),
    "ccl3f": ("CFCl3", "ClC(Cl)(Cl)F"),
    "cl3cf": ("CFCl3", "ClC(Cl)(Cl)F"),
    
    # =========================================================================
    # Non-deuterated solvents
    # =========================================================================
    "dmso": ("DMSO", SMILES["DMSO"]),
    "meoh": ("MeOH", SMILES["MeOH"]),
    "methanol": ("MeOH", SMILES["MeOH"]),
    "mcoh": ("MeOH", SMILES["MeOH"]),  # typo
    "etoh": ("EtOH", SMILES["EtOH"]),
    "ethanol": ("EtOH", SMILES["EtOH"]),
    "acetone": ("Acetone", SMILES["Acetone"]),
    "chcl3": ("CHCl3", SMILES["CHCl3"]),
    "h2o": ("H2O", SMILES["H2O"]),
    "water": ("H2O", SMILES["H2O"]),
    "benzene": ("Benzene", SMILES["Benzene"]),
    "c6h6": ("Benzene", SMILES["Benzene"]),
    "ccl4": ("CCl4", SMILES["CCl4"]),
    "cc14": ("CCl4", SMILES["CCl4"]),  # typo
    "carbon tetrachloride": ("CCl4", SMILES["CCl4"]),
    "formamide": ("Formamide", "NC=O"),
    "et2o": ("Diethyl ether", "CCOCC"),
    "diethyl ether": ("Diethyl ether", "CCOCC"),
    "ethyl ether": ("Diethyl ether", "CCOCC"),
    "triethylamine": ("Triethylamine", "CCN(CC)CC"),
    
    # =========================================================================
    # References
    # =========================================================================
    "tms": ("TMS", SMILES["TMS"]),
    "tetramethylsilane": ("TMS", SMILES["TMS"]),
    "me4si": ("TMS", SMILES["TMS"]),
    "tmsi": ("TMS", SMILES["TMS"]),
    "dtms": ("TMS", SMILES["TMS"]),
    "tsp": ("TSP", "C[Si](C)(C)CCC[S](=O)(=O)[O-].[Na+]"),
    "dss": ("DSS", "C[Si](C)(C)CCCS(=O)(=O)[O-].[Na+]"),
    "hmds": ("HMDS", "C[Si](C)(C)N[Si](C)(C)C"),
    "hmdso": ("HMDSO", "C[Si](C)(C)O[Si](C)(C)C"),
    
    # =========================================================================
    # Special markers
    # =========================================================================
    "neat": ("neat", "NEAT"),
}


def get_solvent_info(text: str):
    """
    Look up a solvent string and return (canonical_name, SMILES) or None.
    """
    key = text.lower().strip()
    return NMR_SOLVENT_DICTIONARY.get(key)


def is_known_solvent(text: str) -> bool:
    """Check if a string is a known NMR solvent."""
    return text.lower().strip() in NMR_SOLVENT_DICTIONARY


# Statistics
if __name__ == "__main__":
    # Count unique canonical solvents
    canonical_names = set(v[0] for v in NMR_SOLVENT_DICTIONARY.values())
    print(f"Dictionary contains {len(NMR_SOLVENT_DICTIONARY)} variant entries")
    print(f"Mapping to {len(canonical_names)} unique canonical solvents")
    print(f"\nCanonical solvents:")
    for name in sorted(canonical_names):
        variants = [k for k, v in NMR_SOLVENT_DICTIONARY.items() if v[0] == name]
        print(f"  {name}: {len(variants)} variants")
