"""
Identifier type detection for Universal Identifier Resolution.

Detects whether an input string is a SMILES, InChI, InChIKey, PubChem CID,
ChEMBL ID, CAS number, DrugBank ID, ChEBI ID, UNII, Wikipedia link, or
compound name.
"""

import re

from rdkit import Chem

# Ordered from most specific to least specific
_PATTERNS: list[tuple[str, re.Pattern]] = [
    ("inchi", re.compile(r"^InChI=", re.IGNORECASE)),
    ("inchikey", re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")),
    ("chembl_id", re.compile(r"^CHEMBL\d+$", re.IGNORECASE)),
    ("drugbank_id", re.compile(r"^DB\d{5}$", re.IGNORECASE)),
    ("chebi_id", re.compile(r"^CHEBI:\d+$", re.IGNORECASE)),
    ("cas", re.compile(r"^\d{2,7}-\d{2}-\d$")),
    ("pubchem_cid", re.compile(r"^(?:CID:)?\d+$", re.IGNORECASE)),
    ("wikipedia", re.compile(r"wikipedia\.org/wiki/", re.IGNORECASE)),
    ("unii", re.compile(r"^[A-Z0-9]{10}$")),
]

# Characters that strongly indicate SMILES
_SMILES_CHARS = re.compile(r"[()[\]=#@/\\]")

_VALID_TYPES = {
    "auto", "name", "smiles", "inchi", "inchikey", "pubchem_cid",
    "chembl_id", "cas", "drugbank_id", "unii", "chebi_id", "wikipedia",
}


def detect_identifier_type(identifier: str, identifier_type: str = "auto") -> str:
    """Detect the type of a chemical identifier.

    Args:
        identifier: The raw input string.
        identifier_type: Explicit type hint. If not "auto", returns it directly
                         (after validation).

    Returns:
        One of: "smiles", "inchi", "inchikey", "pubchem_cid", "chembl_id",
        "cas", "drugbank_id", "chebi_id", "unii", "wikipedia", "name".
    """
    if identifier_type != "auto" and identifier_type in _VALID_TYPES:
        return identifier_type

    stripped = identifier.strip()
    if not stripped:
        return "name"

    # Check regex patterns (most specific first)
    for type_name, pattern in _PATTERNS:
        if pattern.search(stripped):
            return type_name

    # Check if it parses as SMILES via RDKit
    if _SMILES_CHARS.search(stripped):
        return "smiles"

    # Try parsing as SMILES with RDKit (handles simple SMILES like "CCO")
    try:
        mol = Chem.MolFromSmiles(stripped, sanitize=False)
        if mol and mol.GetNumAtoms() > 0:
            # Pure alphanumeric short strings that parse as SMILES: prefer SMILES
            # unless they look like names (4+ alpha chars only)
            if re.match(r"^[a-zA-Z]{4,}$", stripped):
                return "name"
            return "smiles"
    except Exception:
        pass

    return "name"
