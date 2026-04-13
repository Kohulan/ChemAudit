"""
Custom SMARTS Structural Alert Patterns

21 curated SMARTS patterns covering reactive warheads, metal chelators,
redox-active groups, fluorescent interferents, and phospholipidosis risk.
Patterns are compiled once as a module-level singleton for performance.

Sources: ChemAudit BACKEND-IMPLEMENTATION-GUIDE (authoritative).
"""

from __future__ import annotations

from typing import Optional

from rdkit import Chem

from .alert_manager import AlertResult, AlertSeverity

# ---------------------------------------------------------------------------
# Pattern catalogue — 21 entries
# ---------------------------------------------------------------------------
CUSTOM_SMARTS_PATTERNS = [
    # --- Reactive Warheads (CRITICAL) ---
    {
        "name": "acyl_halide",
        "smarts": "[CX3](=[OX1])[F,Cl,Br,I]",
        "description": "Acyl halide — highly reactive electrophile",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "sulfonyl_halide",
        "smarts": "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]",
        "description": "Sulfonyl halide — reactive electrophile",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "vinyl_sulfone",
        "smarts": "[SX4](=[OX1])(=[OX1])/[CX3]=[CX3]",
        "description": "Vinyl sulfone — Michael acceptor",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "epoxide",
        "smarts": "[OX2r3]1[#6r3][#6r3]1",
        "description": "Epoxide — strained ring electrophile",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "aziridine",
        "smarts": "[NX3r3]1[#6r3][#6r3]1",
        "description": "Aziridine — strained ring electrophile",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "isocyanate",
        "smarts": "[NX2]=[CX2]=[OX1]",
        "description": "Isocyanate — electrophilic carbon",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "isothiocyanate",
        "smarts": "[NX2]=[CX2]=[SX1]",
        "description": "Isothiocyanate — electrophilic carbon",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "propargyl_halide",
        "smarts": "[F,Cl,Br,I]CC#C",
        "description": "Propargyl halide — SN2 reactive",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "alpha_haloketone",
        "smarts": "[CX3](=[OX1])[CX4][F,Cl,Br,I]",
        "description": "Alpha-haloketone — alkylating agent",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "acyl_hydrazide",
        "smarts": "[CX3](=[OX1])[NX3][NX3]",
        "description": "Acyl hydrazide — reactive nucleophile",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    {
        "name": "beta_lactam",
        "smarts": "[CX3]1(=[OX1])[NX3][CX4][CX4]1",
        "description": "Beta-lactam — strained ring",
        "concern_group": "Reactive Warheads",
        "severity": AlertSeverity.CRITICAL,
    },
    # --- Metal Chelators (WARNING) ---
    {
        "name": "hydroxamic_acid",
        "smarts": "[CX3](=[OX1])[NX3][OX2H]",
        "description": "Hydroxamic acid — metal chelator",
        "concern_group": "Metal Chelators",
        "severity": AlertSeverity.WARNING,
    },
    {
        "name": "catechol",
        "smarts": "[OX2H]c1ccccc1[OX2H]",
        "description": "Catechol — metal chelator and redox active",
        "concern_group": "Metal Chelators",
        "severity": AlertSeverity.WARNING,
    },
    {
        "name": "dithiocarbamate",
        "smarts": "[SX2][CX3](=[SX1])[#7]",
        "description": "Dithiocarbamate — metal chelator",
        "concern_group": "Metal Chelators",
        "severity": AlertSeverity.WARNING,
    },
    {
        "name": "phosphonic_acid",
        "smarts": "[PX4](=[OX1])([OX2H])([OX2H])",
        "description": "Phosphonic acid — strong metal chelator",
        "concern_group": "Metal Chelators",
        "severity": AlertSeverity.WARNING,
    },
    # --- Redox-Active Groups (WARNING) ---
    {
        "name": "quinone_redox",
        "smarts": "[#6]1(=[OX1])[#6]=[#6][#6](=[OX1])[#6]=[#6]1",
        "description": "Quinone — redox cycling",
        "concern_group": "Redox-Active Groups",
        "severity": AlertSeverity.WARNING,
    },
    {
        "name": "nitroquinoline",
        "smarts": "[nX2]1c([NX3](=[OX1])=[OX1])cc2ccccc12",
        "description": "Nitroquinoline — redox-active mutagen",
        "concern_group": "Redox-Active Groups",
        "severity": AlertSeverity.WARNING,
    },
    # --- Fluorescent Interferents (INFO) ---
    {
        "name": "coumarin",
        "smarts": "[#8r6]1[#6](=[OX1])[#6][#6]c2ccccc12",
        "description": "Coumarin — fluorescent interferent",
        "concern_group": "Fluorescent Interferents",
        "severity": AlertSeverity.INFO,
    },
    {
        "name": "fluorescein",
        "smarts": "[OX2H]c1ccc2c(c1)Oc1cc(O)ccc1C2=O",
        "description": "Fluorescein scaffold — fluorescent interferent",
        "concern_group": "Fluorescent Interferents",
        "severity": AlertSeverity.INFO,
    },
    {
        "name": "rhodamine_scaffold",
        "smarts": "c1ccc2c(c1)c1cc(N)ccc1[C@@H]2c1ccc(N)cc1",
        "description": "Rhodamine scaffold — fluorescent interferent",
        "concern_group": "Fluorescent Interferents",
        "severity": AlertSeverity.INFO,
    },
    # --- Phospholipidosis Risk (WARNING) ---
    {
        "name": "cad_motif",
        "smarts": "[NX3;H2,H1,H0;!$(NC=O)]CCCCCCCC",
        "description": (
            "CAD motif — phospholipidosis risk "
            "(amphiphilic with basic N and long chain)"
        ),
        "concern_group": "Phospholipidosis Risk",
        "severity": AlertSeverity.WARNING,
    },
]

# ---------------------------------------------------------------------------
# Module-level compiled singleton
# ---------------------------------------------------------------------------
_COMPILED_CUSTOM: Optional[list] = None


def get_custom_patterns() -> list:
    """
    Return the compiled CUSTOM_SMARTS_PATTERNS list, compiling on first call.

    Each entry is a dict with keys: name, smarts, description, concern_group,
    severity, and mol (compiled RDKit molecule from SMARTS).

    Returns:
        List of pattern dicts with compiled RDKit query molecules.
    """
    global _COMPILED_CUSTOM
    if _COMPILED_CUSTOM is not None:
        return _COMPILED_CUSTOM

    compiled = []
    for pattern in CUSTOM_SMARTS_PATTERNS:
        mol_query = Chem.MolFromSmarts(pattern["smarts"])
        if mol_query is None:
            raise ValueError(
                f"Failed to compile custom SMARTS pattern '{pattern['name']}': "
                f"{pattern['smarts']}"
            )
        compiled.append({**pattern, "mol": mol_query})

    _COMPILED_CUSTOM = compiled
    return _COMPILED_CUSTOM


def screen_custom_smarts(mol: Chem.Mol) -> list[AlertResult]:
    """
    Screen a molecule against all 21 custom SMARTS patterns.

    Args:
        mol: RDKit molecule to screen.

    Returns:
        List of AlertResult objects for each matching pattern.
    """
    if mol is None:
        return []

    results: list[AlertResult] = []
    for pattern in get_custom_patterns():
        matches = mol.GetSubstructMatches(pattern["mol"])
        if matches:
            # Flatten atom indices from all match instances, deduplicate
            atom_indices = sorted({atom for match in matches for atom in match})
            results.append(
                AlertResult(
                    pattern_name=pattern["name"],
                    description=pattern["description"],
                    severity=pattern["severity"],
                    matched_atoms=atom_indices,
                    catalog_source="CUSTOM",
                    smarts=pattern["smarts"],
                    concern_group=pattern["concern_group"],
                )
            )

    return results
