"""
Kazius Mutagenicity Toxicophore SMARTS Patterns

29 SMARTS-based toxicophore rules derived from the Kazius/McGuire/Saeh (2005)
mutagenicity study. These patterns identify structural features associated with
Ames mutagenicity. All alerts are WARNING level — presence prompts investigation
but does not imply definitive mutagenicity.

Reference:
    Kazius J, McGuire R, Saeh R. Derivation and Validation of Toxicophores for
    Mutagenicity Prediction. J Med Chem 48 (2005) 312-320.
    DOI: 10.1021/jm040835a
"""

from __future__ import annotations

from typing import Optional

from rdkit import Chem

from .alert_manager import AlertResult, AlertSeverity

# ---------------------------------------------------------------------------
# Pattern catalogue — 29 Kazius toxicophores
# ---------------------------------------------------------------------------
KAZIUS_TOXICOPHORES = [
    {
        "name": "aromatic_amine",
        "smarts": "[NX3;H2,H1;!$(NC=O)]c1ccccc1",
        "description": "Aromatic amine — Kazius mutagenicity toxicophore",
    },
    {
        "name": "aromatic_nitro",
        "smarts": "[NX3](=[OX1])([OX1])c1ccccc1",
        "description": "Aromatic nitro compound — Kazius mutagenicity toxicophore",
    },
    {
        "name": "polycyclic_aromatic_hydrocarbon",
        "smarts": "c1ccc2c(c1)ccc1ccccc12",
        "description": "Polycyclic aromatic hydrocarbon — Kazius mutagenicity toxicophore",
    },
    {
        "name": "epoxide_kazius",
        "smarts": "[OX2r3]1[#6r3][#6r3]1",
        "description": "Epoxide — Kazius mutagenicity toxicophore",
    },
    {
        "name": "aziridine_kazius",
        "smarts": "[NX3r3]1[#6r3][#6r3]1",
        "description": "Aziridine — Kazius mutagenicity toxicophore",
    },
    {
        "name": "aliphatic_halide",
        "smarts": "[CX4][F,Cl,Br,I]",
        "description": "Aliphatic halide — Kazius mutagenicity toxicophore",
    },
    {
        "name": "alkyl_nitrite",
        "smarts": "[NX2](=[OX1])[OX2][CX4]",
        "description": "Alkyl nitrite — Kazius mutagenicity toxicophore",
    },
    {
        "name": "n_nitroso",
        "smarts": "[NX3][NX2]=[OX1]",
        "description": "N-nitroso compound — Kazius mutagenicity toxicophore",
    },
    {
        "name": "diazo",
        "smarts": "[NX2]=[NX2]",
        "description": "Diazo compound — Kazius mutagenicity toxicophore",
    },
    {
        "name": "triazene",
        "smarts": "[NX3][NX2]=[NX2]",
        "description": "Triazene — Kazius mutagenicity toxicophore",
    },
    {
        "name": "nitrogen_mustard",
        "smarts": "[NX3]([CH2][CH2][Cl,Br,I])[CH2][CH2][Cl,Br,I]",
        "description": "Nitrogen mustard — Kazius mutagenicity toxicophore",
    },
    {
        "name": "azide",
        "smarts": "[NX1]=[NX2]=[NX1]",
        "description": "Azide — Kazius mutagenicity toxicophore",
    },
    {
        "name": "hydrazine",
        "smarts": "[NX3][NX3]",
        "description": "Hydrazine — Kazius mutagenicity toxicophore",
    },
    {
        "name": "acyl_halide_kazius",
        "smarts": "[CX3](=[OX1])[Cl,Br,I]",
        "description": "Acyl halide — Kazius mutagenicity toxicophore",
    },
    {
        "name": "michael_acceptor",
        "smarts": "[CX3]=[CX3][CX3](=[OX1])",
        "description": "Michael acceptor — Kazius mutagenicity toxicophore",
    },
    {
        "name": "unsubstituted_heteroatom_bonded_heteroatom",
        "smarts": "[#7,#8,#16]~[#7,#8,#16]",
        "description": (
            "Unsubstituted heteroatom bonded to heteroatom — "
            "Kazius mutagenicity toxicophore"
        ),
    },
    {
        "name": "aldehyde",
        "smarts": "[CX3H1](=[OX1])",
        "description": "Aldehyde — Kazius mutagenicity toxicophore",
    },
    {
        "name": "stilbene",
        "smarts": "c1ccccc1/C=C/c1ccccc1",
        "description": "Stilbene — Kazius mutagenicity toxicophore",
    },
    {
        "name": "azo_compound",
        "smarts": "[#6]/[NX2]=[NX2]/[#6]",
        "description": "Azo compound — Kazius mutagenicity toxicophore",
    },
    {
        "name": "beta_propiolactone",
        "smarts": "[OX2r4]1[CX3](=[OX1])[CX4][CX4]1",
        "description": "Beta-propiolactone — Kazius mutagenicity toxicophore",
    },
    {
        "name": "three_membered_heterocycle",
        "smarts": "[#7X3,#8X2,#16X2]1~[#6]~[#6]1",
        "description": "Three-membered heterocycle — Kazius mutagenicity toxicophore",
    },
    {
        "name": "nitro_compound",
        "smarts": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "description": "Nitro compound — Kazius mutagenicity toxicophore",
    },
    {
        "name": "aromatic_nitroso",
        "smarts": "[NX2](=[OX1])c",
        "description": "Aromatic nitroso — Kazius mutagenicity toxicophore",
    },
    {
        "name": "polycyclic_system",
        "smarts": "c1ccc2c(c1)cc1ccc3ccccc3c1c2",
        "description": "Polycyclic aromatic system — Kazius mutagenicity toxicophore",
    },
    {
        "name": "quinone",
        "smarts": "[#6]1(=[OX1])[#6]=[#6][#6](=[OX1])[#6]=[#6]1",
        "description": "Quinone — Kazius mutagenicity toxicophore",
    },
    {
        "name": "aliphatic_nitro",
        "smarts": "[CX4][NX3](=[OX1])=[OX1]",
        "description": "Aliphatic nitro — Kazius mutagenicity toxicophore",
    },
    {
        "name": "hydroxylamine",
        "smarts": "[NX3;H1,H0][OX2H]",
        "description": "Hydroxylamine — Kazius mutagenicity toxicophore",
    },
    {
        "name": "amino_phenol",
        "smarts": "[NX3;H2,H1]c1ccc([OX2H])cc1",
        "description": "Amino phenol — Kazius mutagenicity toxicophore",
    },
    {
        "name": "thiocarbonyl",
        "smarts": "[#6]=[SX1]",
        "description": "Thiocarbonyl — Kazius mutagenicity toxicophore",
    },
]

# ---------------------------------------------------------------------------
# Module-level compiled singleton
# ---------------------------------------------------------------------------
_COMPILED_KAZIUS: Optional[list] = None


def get_kazius_patterns() -> list:
    """
    Return the compiled KAZIUS_TOXICOPHORES list, compiling on first call.

    Each entry is a dict with keys: name, smarts, description, and mol
    (compiled RDKit query molecule from SMARTS).

    Returns:
        List of pattern dicts with compiled RDKit query molecules.
    """
    global _COMPILED_KAZIUS
    if _COMPILED_KAZIUS is not None:
        return _COMPILED_KAZIUS

    compiled = []
    for pattern in KAZIUS_TOXICOPHORES:
        mol_query = Chem.MolFromSmarts(pattern["smarts"])
        if mol_query is None:
            raise ValueError(
                f"Failed to compile Kazius SMARTS pattern '{pattern['name']}': "
                f"{pattern['smarts']}"
            )
        compiled.append({**pattern, "mol": mol_query})

    _COMPILED_KAZIUS = compiled
    return _COMPILED_KAZIUS


def screen_kazius(mol: Chem.Mol) -> list[AlertResult]:
    """
    Screen a molecule against all 29 Kazius mutagenicity toxicophores.

    All matches are returned as WARNING-level alerts — toxicophore hits prompt
    investigation but do not imply definitive mutagenicity.

    Args:
        mol: RDKit molecule to screen.

    Returns:
        List of AlertResult objects for each matching toxicophore pattern.
    """
    if mol is None:
        return []

    results: list[AlertResult] = []
    for pattern in get_kazius_patterns():
        matches = mol.GetSubstructMatches(pattern["mol"])
        if matches:
            atom_indices = sorted({atom for match in matches for atom in match})
            results.append(
                AlertResult(
                    pattern_name=pattern["name"],
                    description=pattern["description"],
                    severity=AlertSeverity.WARNING,
                    matched_atoms=atom_indices,
                    catalog_source="KAZIUS",
                    smarts=pattern["smarts"],
                    concern_group="Mutagenicity Toxicophore",
                )
            )

    return results
