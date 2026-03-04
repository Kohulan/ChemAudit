"""
Cross-database comparator service.

Fetches compound representations from PubChem, ChEMBL, and COCONUT,
then compares them field-by-field using RDKit canonicalization.
"""

import asyncio
import logging
from typing import Optional

from rdkit import Chem

from app.schemas.integrations import (
    ChEMBLRequest,
    ChEMBLResult,
    COCONUTRequest,
    COCONUTResult,
    ConsistencyResult,
    DatabaseEntry,
    PropertyComparison,
    PubChemRequest,
    PubChemResult,
)
from app.services.integrations.chembl import get_bioactivity
from app.services.integrations.coconut import lookup_natural_product
from app.services.integrations.pubchem import get_compound_info

logger = logging.getLogger(__name__)


def _canonical_smiles(smiles: Optional[str]) -> Optional[str]:
    """Canonicalize a SMILES string via RDKit for fair comparison."""
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return smiles


def _compare_field(
    prop_name: str,
    values: dict[str, Optional[str]],
    tolerance: Optional[float] = None,
) -> PropertyComparison:
    """Compare a single property across databases.

    Returns a PropertyComparison with status:
    - "match": all non-None values are identical
    - "minor_difference": values differ but within tolerance or tautomeric
    - "mismatch": values are structurally different
    - "missing": fewer than 2 databases have data
    """
    non_null = {db: v for db, v in values.items() if v is not None}

    if len(non_null) < 2:
        return PropertyComparison(
            property_name=prop_name,
            values=values,
            status="missing",
            detail=f"Only {len(non_null)} database(s) have this property",
        )

    unique_values = set(non_null.values())

    if len(unique_values) == 1:
        return PropertyComparison(
            property_name=prop_name,
            values=values,
            status="match",
        )

    # For numeric fields: check within tolerance
    if tolerance is not None:
        try:
            nums = [float(v) for v in non_null.values()]
            if max(nums) - min(nums) <= tolerance:
                return PropertyComparison(
                    property_name=prop_name,
                    values=values,
                    status="match",
                    detail=f"Within tolerance ({tolerance})",
                )
        except (ValueError, TypeError):
            pass

    # For SMILES: check if canonicalized forms match
    if prop_name == "canonical_smiles":
        canon_values = {db: _canonical_smiles(v) for db, v in non_null.items()}
        unique_canon = set(v for v in canon_values.values() if v is not None)
        if len(unique_canon) == 1:
            return PropertyComparison(
                property_name=prop_name,
                values=values,
                status="match",
                detail="Identical after RDKit canonicalization",
            )
        return PropertyComparison(
            property_name=prop_name,
            values=values,
            status="mismatch",
            detail="Canonical SMILES differ across databases",
        )

    return PropertyComparison(
        property_name=prop_name,
        values=values,
        status="mismatch",
        detail=f"Values differ: {', '.join(str(v) for v in unique_values)}",
    )


async def compare_across_databases(
    smiles: Optional[str] = None,
    inchikey: Optional[str] = None,
) -> ConsistencyResult:
    """Compare compound representation across PubChem, ChEMBL, and COCONUT.

    Args:
        smiles: SMILES string of the compound.
        inchikey: InChIKey of the compound.

    Returns:
        ConsistencyResult with per-database entries, field comparisons, and verdict.
    """
    if not smiles and not inchikey:
        return ConsistencyResult(overall_verdict="no_data", summary="No input provided")

    # Generate InChIKey from SMILES if not provided
    if smiles and not inchikey:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                inchikey = Chem.MolToInchiKey(mol)
        except Exception:
            pass

    # Query all databases in parallel
    pc_result, chembl_result, coconut_result = await asyncio.gather(
        get_compound_info(PubChemRequest(smiles=smiles, inchikey=inchikey)),
        get_bioactivity(ChEMBLRequest(smiles=smiles, inchikey=inchikey)),
        lookup_natural_product(COCONUTRequest(smiles=smiles, inchikey=inchikey)),
        return_exceptions=True,
    )

    # Handle exceptions gracefully
    if isinstance(pc_result, Exception):
        pc_result = PubChemResult(found=False)
    if isinstance(chembl_result, Exception):
        chembl_result = ChEMBLResult(found=False)
    if isinstance(coconut_result, Exception):
        coconut_result = COCONUTResult(found=False)

    # Build database entries
    entries: list[DatabaseEntry] = []

    if pc_result.found:
        entries.append(DatabaseEntry(
            database="PubChem",
            found=True,
            canonical_smiles=pc_result.canonical_smiles,
            inchi=pc_result.inchi,
            inchikey=pc_result.inchikey,
            molecular_formula=pc_result.molecular_formula,
            molecular_weight=pc_result.molecular_weight,
            name=pc_result.iupac_name,
            url=pc_result.url,
        ))
    else:
        entries.append(DatabaseEntry(database="PubChem", found=False))

    if chembl_result.found:
        entries.append(DatabaseEntry(
            database="ChEMBL",
            found=True,
            molecular_formula=chembl_result.molecular_formula,
            molecular_weight=chembl_result.molecular_weight,
            name=chembl_result.pref_name,
            url=chembl_result.url,
        ))
    else:
        entries.append(DatabaseEntry(database="ChEMBL", found=False))

    if coconut_result.found:
        entries.append(DatabaseEntry(
            database="COCONUT",
            found=True,
            canonical_smiles=coconut_result.smiles,
            inchikey=coconut_result.inchikey,
            molecular_formula=coconut_result.molecular_formula,
            molecular_weight=coconut_result.molecular_weight,
            name=coconut_result.name,
            url=coconut_result.url,
        ))
    else:
        entries.append(DatabaseEntry(database="COCONUT", found=False))

    # Check how many databases found it
    found_entries = [e for e in entries if e.found]
    if len(found_entries) == 0:
        return ConsistencyResult(
            entries=entries,
            overall_verdict="no_data",
            summary="Compound not found in any database",
        )

    if len(found_entries) == 1:
        db_name = found_entries[0].database
        return ConsistencyResult(
            entries=entries,
            overall_verdict="consistent",
            summary=f"Found only in {db_name} — no cross-database comparison possible",
        )

    # Build comparisons for common properties
    comparisons: list[PropertyComparison] = []

    # Molecular formula
    formula_values = {e.database: e.molecular_formula for e in found_entries}
    comparisons.append(_compare_field("molecular_formula", formula_values))

    # Molecular weight (tolerance 0.1 Da)
    mw_values = {
        e.database: str(e.molecular_weight) if e.molecular_weight is not None else None
        for e in found_entries
    }
    comparisons.append(_compare_field("molecular_weight", mw_values, tolerance=0.1))

    # SMILES (only if available)
    smiles_values = {e.database: e.canonical_smiles for e in found_entries}
    if any(v is not None for v in smiles_values.values()):
        comparisons.append(_compare_field("canonical_smiles", smiles_values))

    # InChIKey
    ik_values = {e.database: e.inchikey for e in found_entries}
    if any(v is not None for v in ik_values.values()):
        comparisons.append(_compare_field("inchikey", ik_values))

    # Determine overall verdict
    statuses = [c.status for c in comparisons]
    if "mismatch" in statuses:
        verdict = "major_discrepancies"
        summary = "Structural differences found across databases"
    elif "minor_difference" in statuses:
        verdict = "minor_differences"
        summary = "Minor differences detected (tautomers or stereochemistry)"
    else:
        verdict = "consistent"
        summary = "All compared properties match across databases"

    return ConsistencyResult(
        entries=entries,
        comparisons=comparisons,
        overall_verdict=verdict,
        summary=summary,
    )
