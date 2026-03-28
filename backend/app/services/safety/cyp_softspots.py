"""
CYP Metabolism Soft-Spot Prediction (SAFE-01)

SMARTS-based prediction of CYP450 metabolism soft-spots with atom-level indices
and reaction type annotation.
"""
from __future__ import annotations

from typing import Optional

from rdkit import Chem

# 11 CYP soft-spot SMARTS patterns with reaction type annotations
CYP_SOFTSPOT_SMARTS: list[tuple[str, str, str]] = [
    ("benzylic_CH", "[CX4H2,CX4H1]c", "C-hydroxylation"),
    ("allylic_CH", "[CX4H2,CX4H1][CX3]=[CX3]", "C-hydroxylation"),
    ("omega_CH3", "[CH3][CH2][CH2]", "omega-hydroxylation"),
    ("N_dealkylation", "[NX3;!$(NC=O)]([CH3,CH2])", "N-dealkylation"),
    ("O_dealkylation", "[OX2]([CH3,CH2])[c,C]", "O-dealkylation"),
    ("S_oxidation", "[SX2][#6]", "S-oxidation"),
    ("ester_hydrolysis", "[CX3](=[OX1])[OX2][#6]", "ester-hydrolysis"),
    ("amide_hydrolysis", "[CX3](=[OX1])[NX3][#6]", "amide-hydrolysis"),
    ("N_oxidation", "[nX2]", "N-oxidation"),
    ("aromatic_para_H", "[cH]1[cH][c]([*])[cH][cH][c]1", "aromatic-hydroxylation"),
    ("omega1_CH2", "[CH2]([CH2][CH3])", "omega-1-hydroxylation"),
]

# Module-level compiled singleton
_CYP_COMPILED: Optional[list[tuple[str, Chem.Mol, str]]] = None


def get_cyp_patterns() -> list[tuple[str, Chem.Mol, str]]:
    """Return lazily compiled CYP SMARTS patterns as (name, compiled_mol, reaction_type) tuples.

    Compilation is performed once and cached in module-level singleton. Patterns that
    fail to compile are silently filtered out.

    Returns:
        List of (site_name, compiled_pattern, reaction_type) tuples.
    """
    global _CYP_COMPILED
    if _CYP_COMPILED is not None:
        return _CYP_COMPILED

    compiled = []
    for name, smarts, reaction_type in CYP_SOFTSPOT_SMARTS:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            compiled.append((name, pattern, reaction_type))

    _CYP_COMPILED = compiled
    return _CYP_COMPILED


def screen_cyp_softspots(mol: Chem.Mol) -> list[dict]:
    """Screen a molecule for CYP metabolism soft-spots using SMARTS matching.

    Each matched instance of a pattern is returned as a separate entry.
    The same site_name can appear multiple times for multiple matches within
    the same molecule.

    Args:
        mol: RDKit molecule to screen.

    Returns:
        List of dicts, each with keys:
          - site_name (str): e.g. "benzylic_CH"
          - reaction_type (str): e.g. "C-hydroxylation"
          - matched_atoms (list[int]): atom indices from the substructure match
    """
    results: list[dict] = []
    for site_name, pattern, reaction_type in get_cyp_patterns():
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            results.append(
                {
                    "site_name": site_name,
                    "reaction_type": reaction_type,
                    "matched_atoms": list(match),
                }
            )
    return results
