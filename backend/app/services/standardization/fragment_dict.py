"""
Fragment Dictionary for Standardization Provenance.

Provides a curated COUNTERION_NAMES dictionary and classify_fragment() function
for identifying and categorizing removed fragments in the parent extraction stage.
"""

from typing import Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Curated dictionary of common counterions, salts, and solvents.
# Keys are canonical SMILES. Values contain name and role.
# Roles: "counterion" (charged species), "salt" (full salt fragment),
#        "solvent" (neutral solvent), "unknown" (not in dictionary)
COUNTERION_NAMES: dict[str, dict[str, str]] = {
    # Halide counterions (anions)
    "[Cl-]": {"name": "chloride", "role": "counterion"},
    "[Br-]": {"name": "bromide", "role": "counterion"},
    "[I-]": {"name": "iodide", "role": "counterion"},
    "[F-]": {"name": "fluoride", "role": "counterion"},
    # Metal cations
    "[Na+]": {"name": "sodium", "role": "counterion"},
    "[K+]": {"name": "potassium", "role": "counterion"},
    "[Li+]": {"name": "lithium", "role": "counterion"},
    "[Ca+2]": {"name": "calcium", "role": "counterion"},
    "[Mg+2]": {"name": "magnesium", "role": "counterion"},
    "[Zn+2]": {"name": "zinc", "role": "counterion"},
    "[NH4+]": {"name": "ammonium", "role": "counterion"},
    # Neutral metal atoms (sometimes appear as disconnected fragments)
    "[Na]": {"name": "sodium", "role": "counterion"},
    "[K]": {"name": "potassium", "role": "counterion"},
    "[Li]": {"name": "lithium", "role": "counterion"},
    # Hydroxide/oxide
    "[OH-]": {"name": "hydroxide", "role": "counterion"},
    "[O-2]": {"name": "oxide", "role": "counterion"},
    # Acids commonly used as salt formers (salt role)
    "Cl": {"name": "hydrochloric acid", "role": "salt"},
    "Br": {"name": "hydrobromic acid", "role": "salt"},
    "I": {"name": "hydroiodic acid", "role": "salt"},
    "CC(=O)O": {"name": "acetic acid", "role": "salt"},
    "OC(=O)C(F)(F)F": {"name": "trifluoroacetic acid", "role": "salt"},
    "OC(=O)c1ccccc1": {"name": "benzoic acid", "role": "salt"},
    "OS(=O)(=O)O": {"name": "sulfuric acid", "role": "salt"},
    "OS(=O)(=O)c1ccc(C)cc1": {"name": "p-toluenesulfonic acid", "role": "salt"},
    "OP(=O)(O)O": {"name": "phosphoric acid", "role": "salt"},
    "OC(=O)O": {"name": "carbonic acid", "role": "salt"},
    "OC(=O)CC(O)(CC(=O)O)C(=O)O": {"name": "citric acid", "role": "salt"},
    "OC(O)=O": {"name": "formic acid", "role": "salt"},
    "OC(=O)C(O)=O": {"name": "oxalic acid", "role": "salt"},
    "OC(=O)CCC(=O)O": {"name": "succinic acid", "role": "salt"},
    "OC(=O)C=CC(=O)O": {"name": "maleic acid", "role": "salt"},
    "OC(=O)/C=C/C(=O)O": {"name": "fumaric acid", "role": "salt"},
    "OC(=O)[C@@H](O)C(=O)O": {"name": "tartaric acid", "role": "salt"},
    "OC(=O)[C@H](O)C(=O)O": {"name": "tartaric acid", "role": "salt"},
    "OC(=O)C(O)C(=O)O": {"name": "tartaric acid", "role": "salt"},
    "OS(=O)(=O)c1ccccc1": {"name": "benzenesulfonic acid", "role": "salt"},
    "OS(=O)(=O)CCCC": {"name": "butanesulfonic acid", "role": "salt"},
    # Solvents
    "O": {"name": "water", "role": "solvent"},
    "CO": {"name": "methanol", "role": "solvent"},
    "CCO": {"name": "ethanol", "role": "solvent"},
    "CC(C)O": {"name": "isopropanol", "role": "solvent"},
    "CS(C)=O": {"name": "DMSO", "role": "solvent"},
    "CN(C)C=O": {"name": "DMF", "role": "solvent"},
    "CC#N": {"name": "acetonitrile", "role": "solvent"},
    "CC(C)=O": {"name": "acetone", "role": "solvent"},
    "ClCCl": {"name": "dichloromethane", "role": "solvent"},
    "ClC(Cl)Cl": {"name": "chloroform", "role": "solvent"},
    "CCOC(C)=O": {"name": "ethyl acetate", "role": "solvent"},
    "C1CCOC1": {"name": "THF", "role": "solvent"},
    "ClCCCl": {"name": "1,2-dichloroethane", "role": "solvent"},
    "CCOCC": {"name": "diethyl ether", "role": "solvent"},
    "c1ccncc1": {"name": "pyridine", "role": "solvent"},
    "C1CCCCC1": {"name": "cyclohexane", "role": "solvent"},
    "Cc1ccccc1": {"name": "toluene", "role": "solvent"},
    "c1ccccc1": {"name": "benzene", "role": "solvent"},
    "CCCCCC": {"name": "hexane", "role": "solvent"},
    "CCCCC": {"name": "pentane", "role": "solvent"},
}


def classify_fragment(smiles: str) -> dict:
    """
    Classify a SMILES fragment using the COUNTERION_NAMES dictionary.

    Canonicalizes the input SMILES and looks it up in the dictionary.
    Falls back to role="unknown" for unrecognized fragments.
    Computes molecular weight for all fragments.

    Args:
        smiles: SMILES string of the fragment to classify.

    Returns:
        Dictionary with keys: smiles (canonical), name (str or None),
        role ("salt"|"solvent"|"counterion"|"unknown"), mw (float).
    """
    # Attempt to canonicalize the SMILES
    canonical_smiles = smiles
    mw = 0.0
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
            canonical_smiles = Chem.MolToSmiles(mol)
            mw = float(rdMolDescriptors.CalcExactMolWt(mol))
        else:
            # Invalid SMILES — return graceful fallback
            return {"smiles": smiles, "name": None, "role": "unknown", "mw": 0.0}
    except Exception:
        return {"smiles": smiles, "name": None, "role": "unknown", "mw": 0.0}

    # Look up in dictionary
    entry = COUNTERION_NAMES.get(canonical_smiles)
    if entry:
        return {
            "smiles": canonical_smiles,
            "name": entry["name"],
            "role": entry["role"],
            "mw": mw,
        }

    # Unknown fragment
    return {"smiles": canonical_smiles, "name": None, "role": "unknown", "mw": mw}


def get_fragment_name(smiles: str) -> Optional[str]:
    """
    Get the human-readable name for a known fragment.

    Args:
        smiles: SMILES string of the fragment.

    Returns:
        Name string if found, None otherwise.
    """
    result = classify_fragment(smiles)
    return result.get("name")
