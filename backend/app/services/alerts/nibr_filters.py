"""
NIBR Novartis Screening Deck Filters

Provides access to the Novartis Institute for BioMedical Research (NIBR)
compound screening filters via the medchem library. The RDKit FilterCatalog
from medchem NIBRFilters contains 444 entries covering various problematic
substructures and functional groups used in industrial drug discovery.

Integration: Hybrid singleton pattern (D-23) — load catalog from medchem
NIBRFilters at startup; screen using the catalog's own GetMatches().

Reference:
    Schuffenhauer A, et al. Evolution of Novartis' Small Molecule Screening
    Deck Design. J Med Chem 63 (2020) 14425-14447.
    DOI: 10.1021/acs.jmedchem.0c01332
"""

from __future__ import annotations

from typing import Optional

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog

from .alert_manager import AlertResult, AlertSeverity

# ---------------------------------------------------------------------------
# Module-level NIBR catalog singleton
# ---------------------------------------------------------------------------
_NIBR_SINGLETON: Optional[FilterCatalog] = None

# Severity code mapping from NIBR catalog description strings
_NIBR_SEVERITY_MAP: dict[int, AlertSeverity] = {
    0: AlertSeverity.INFO,
    1: AlertSeverity.WARNING,
    2: AlertSeverity.CRITICAL,
}

_NIBR_CONCERN_MAP: dict[int, str] = {
    0: "NIBR Annotations",
    1: "NIBR Screening Flags",
    2: "NIBR Excluded Patterns",
}


def get_nibr_catalog() -> FilterCatalog:
    """
    Return the NIBR RDKit FilterCatalog, initialising on first call.

    The catalog is loaded from medchem.structural.NIBRFilters and contains
    444 entries. Pre-warming at FastAPI startup is strongly recommended to
    avoid a 5-10s first-call spike (see Plan 03 lifespan integration).

    Returns:
        Initialised RDKit FilterCatalog with NIBR entries.
    """
    global _NIBR_SINGLETON
    if _NIBR_SINGLETON is not None:
        return _NIBR_SINGLETON

    from medchem.structural import NIBRFilters as _MedchemNIBR  # noqa: PLC0415

    _NIBR_SINGLETON = _MedchemNIBR().catalog
    return _NIBR_SINGLETON


def screen_nibr(mol: Chem.Mol) -> list[AlertResult]:
    """
    Screen a molecule against the NIBR Novartis structural alert catalog.

    Iterates catalog.GetMatches(mol), parses the description field for
    NIBR metadata, and maps severity codes to AlertSeverity levels:
      - 0 → INFO  (NIBR Annotations)
      - 1 → WARNING (NIBR Screening Flags)
      - 2 → CRITICAL (NIBR Excluded Patterns)

    Description format: ``NIBR||<name>||<severity_code>||<n_covalent>||<special_mol>``

    Args:
        mol: RDKit molecule to screen.

    Returns:
        List of AlertResult objects for each NIBR catalog match.
    """
    if mol is None:
        return []

    catalog = get_nibr_catalog()
    results: list[AlertResult] = []

    matches = catalog.GetMatches(mol)
    for entry in matches:
        raw_desc = entry.GetDescription() or ""

        # Parse NIBR pipe-separated description
        parts = raw_desc.split("||")

        if len(parts) >= 3:
            name = parts[1].strip() if len(parts) > 1 else raw_desc
            try:
                severity_code = int(parts[2].strip())
            except (ValueError, IndexError):
                severity_code = 1  # default to WARNING
        else:
            name = raw_desc
            severity_code = 1

        severity = _NIBR_SEVERITY_MAP.get(severity_code, AlertSeverity.WARNING)
        concern_group = _NIBR_CONCERN_MAP.get(severity_code, "NIBR Screening Flags")

        # Extract matched atom indices via FilterMatch atom pairs
        atom_indices: list[int] = []
        try:
            filter_matches = entry.GetFilterMatches(mol)
            for fm in filter_matches:
                for pair in fm.atomPairs:
                    atom_indices.append(int(pair[1]))
        except Exception:  # noqa: BLE001
            pass

        atom_indices = sorted(set(atom_indices))

        results.append(
            AlertResult(
                pattern_name=name,
                description=f"NIBR filter: {name}",
                severity=severity,
                matched_atoms=atom_indices,
                catalog_source="NIBR",
                catalog_description="NIBR Novartis Screening Deck Filters",
                concern_group=concern_group,
            )
        )

    return results
