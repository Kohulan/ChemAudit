"""
Cross-database comparator service.

Fetches compound representations from PubChem, ChEMBL, COCONUT, and Wikidata,
then compares them field-by-field using RDKit canonicalization.
"""

import asyncio
import logging
from typing import Optional

from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

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
from app.services.integrations.wikidata import WikidataClient

logger = logging.getLogger(__name__)

RDKIT_VERSION = rdBase.rdkitVersion


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


def _tautomer_canonical_smiles(smiles: str) -> Optional[str]:
    """Canonicalize SMILES via tautomer enumeration for tautomer-invariant comparison."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        canon = rdMolStandardize.TautomerEnumerator()
        canon_mol = canon.Canonicalize(mol)
        return Chem.MolToSmiles(canon_mol)
    except Exception:
        return None


def _compare_smiles(values: dict[str, Optional[str]]) -> PropertyComparison:
    """Compare SMILES across databases with rich detail."""
    non_null = {db: v for db, v in values.items() if v is not None}
    toolkit = f"RDKit {RDKIT_VERSION}"

    if len(non_null) < 2:
        return PropertyComparison(
            property_name="canonical_smiles",
            values=values,
            status="missing",
            detail=(
                f"Only {len(non_null)} database(s) returned SMILES. "
                f"At least 2 are needed for comparison."
            ),
        )

    # Step 1: direct string comparison
    unique_raw = set(non_null.values())
    if len(unique_raw) == 1:
        return PropertyComparison(
            property_name="canonical_smiles",
            values=values,
            status="match",
            detail=(
                f"All databases return identical SMILES strings. "
                f"Canonicalized with {toolkit}."
            ),
        )

    # Step 2: RDKit canonical comparison
    canon_values = {db: _canonical_smiles(v) for db, v in non_null.items()}
    unique_canon = set(v for v in canon_values.values() if v is not None)
    if len(unique_canon) == 1:
        return PropertyComparison(
            property_name="canonical_smiles",
            values=values,
            status="match",
            detail=(
                f"Raw SMILES strings differ between databases but resolve to the "
                f"same canonical form after {toolkit} canonicalization. "
                f"Each database uses its own toolkit to generate SMILES, "
                f"so surface-level text differences are expected."
            ),
        )

    # Step 3: tautomer-invariant comparison
    tautomer_values = {
        db: _tautomer_canonical_smiles(v) for db, v in non_null.items()
    }
    unique_tautomer = set(v for v in tautomer_values.values() if v is not None)
    if len(unique_tautomer) == 1:
        return PropertyComparison(
            property_name="canonical_smiles",
            values=values,
            status="minor_difference",
            detail=(
                f"SMILES differ after {toolkit} canonicalization but resolve to the "
                f"same form after tautomer enumeration ({toolkit} "
                f"TautomerEnumerator). These are tautomeric forms of the same "
                f"compound — databases store different protonation or keto/enol states."
            ),
        )

    # Step 4: genuine mismatch — explain possible causes
    differing_dbs = [
        f"{db}: {canon_values.get(db, v)}" for db, v in non_null.items()
    ]
    return PropertyComparison(
        property_name="canonical_smiles",
        values=values,
        status="mismatch",
        detail=(
            f"Canonical SMILES differ even after {toolkit} canonicalization and "
            f"tautomer normalization. Possible causes: different stereoisomers "
            f"(E/Z or R/S), different charge or salt forms, or genuinely different "
            f"structures registered under the same identifier. "
            f"Canonical forms: {'; '.join(differing_dbs)}."
        ),
    )


def _compare_inchikey(values: dict[str, Optional[str]]) -> PropertyComparison:
    """Compare InChIKeys across databases with layer-level analysis."""
    non_null = {db: v for db, v in values.items() if v is not None}

    if len(non_null) < 2:
        return PropertyComparison(
            property_name="inchikey",
            values=values,
            status="missing",
            detail=(
                f"Only {len(non_null)} database(s) returned an InChIKey. "
                f"At least 2 are needed for comparison."
            ),
        )

    unique_values = set(non_null.values())
    if len(unique_values) == 1:
        return PropertyComparison(
            property_name="inchikey",
            values=values,
            status="match",
            detail=(
                "InChIKeys are identical across all databases. "
                "An InChIKey has 3 layers: connectivity (first 14 chars), "
                "stereochemistry + charge (next 8 chars), and version (last char). "
                "Full match confirms the same molecular structure."
            ),
        )

    # Check layer-by-layer: connectivity (first 14), stereo (chars 15-24)
    connectivity = {db: v[:14] for db, v in non_null.items()}
    unique_connectivity = set(connectivity.values())

    if len(unique_connectivity) == 1:
        # Same connectivity, different stereo/charge layers
        return PropertyComparison(
            property_name="inchikey",
            values=values,
            status="minor_difference",
            detail=(
                "The first 14 characters (connectivity layer) match across all "
                "databases — same atom connectivity. The remaining characters "
                "(stereochemistry and charge layers) differ, indicating different "
                "stereoisomers (R/S, E/Z), protonation states, or salt forms of "
                "the same base structure."
            ),
        )

    return PropertyComparison(
        property_name="inchikey",
        values=values,
        status="mismatch",
        detail=(
            "InChIKeys differ in the connectivity layer (first 14 characters), "
            "meaning the databases have different molecular graphs (atom "
            "connectivity) for this compound. This is a significant discrepancy — "
            "the entries may represent different chemical entities."
        ),
    )


def _compare_inchi(values: dict[str, Optional[str]]) -> PropertyComparison:
    """Compare InChI strings across databases with layer analysis."""
    non_null = {db: v for db, v in values.items() if v is not None}

    if len(non_null) < 2:
        return PropertyComparison(
            property_name="inchi",
            values=values,
            status="missing",
            detail=(
                f"Only {len(non_null)} database(s) returned an InChI string. "
                f"InChI is not always available from all sources — PubChem and "
                f"Wikidata typically provide it, while ChEMBL and COCONUT may not."
            ),
        )

    unique_values = set(non_null.values())
    if len(unique_values) == 1:
        return PropertyComparison(
            property_name="inchi",
            values=values,
            status="match",
            detail=(
                "InChI strings are identical across all databases. InChI encodes "
                "molecular structure in layers: formula, connectivity, hydrogen, "
                "charge, and stereochemistry. Full match confirms identical "
                "structure representation."
            ),
        )

    # Parse InChI layers to find where they diverge
    detail_parts = []
    layers = {}
    for db, inchi in non_null.items():
        parts = inchi.split("/") if inchi else []
        layers[db] = parts

    # Check if formula layer (index 1) matches
    formulas = {db: parts[1] if len(parts) > 1 else None for db, parts in layers.items()}
    unique_formulas = set(v for v in formulas.values() if v is not None)
    if len(unique_formulas) == 1:
        detail_parts.append(
            "Molecular formula layer matches across databases"
        )
        # Check connectivity layer (index 2, starts with 'c')
        conn = {
            db: parts[2] if len(parts) > 2 and parts[2].startswith("c") else None
            for db, parts in layers.items()
        }
        unique_conn = set(v for v in conn.values() if v is not None)
        if len(unique_conn) <= 1:
            detail_parts.append(
                "Connectivity layer also matches — differences are in "
                "stereochemistry, hydrogen, or charge layers"
            )
        else:
            detail_parts.append(
                "Connectivity layers differ — different atom bonding"
            )
    else:
        detail_parts.append(
            "Molecular formula layers differ across databases"
        )

    status = "minor_difference" if len(unique_formulas) == 1 else "mismatch"

    return PropertyComparison(
        property_name="inchi",
        values=values,
        status=status,
        detail=(
            f"InChI strings differ. {'. '.join(detail_parts)}. "
            f"InChI layers encode: /formula/connectivity/hydrogen/charge/"
            f"stereochemistry — differences in later layers typically indicate "
            f"stereoisomers or protonation variants of the same base structure."
        ),
    )


def _build_pubchem_entry(result: PubChemResult) -> DatabaseEntry:
    """Build a DatabaseEntry from a PubChem result."""
    if not result.found:
        return DatabaseEntry(database="PubChem", found=False)
    return DatabaseEntry(
        database="PubChem",
        found=True,
        canonical_smiles=result.canonical_smiles,
        inchi=result.inchi,
        inchikey=result.inchikey,
        molecular_formula=result.molecular_formula,
        molecular_weight=result.molecular_weight,
        name=result.iupac_name,
        url=result.url,
    )


def _build_chembl_entry(result: ChEMBLResult) -> DatabaseEntry:
    """Build a DatabaseEntry from a ChEMBL result."""
    if not result.found:
        return DatabaseEntry(database="ChEMBL", found=False)
    return DatabaseEntry(
        database="ChEMBL",
        found=True,
        canonical_smiles=result.canonical_smiles,
        inchi=result.inchi,
        inchikey=result.inchikey,
        molecular_formula=result.molecular_formula,
        molecular_weight=result.molecular_weight,
        name=result.pref_name,
        url=result.url,
    )


def _build_coconut_entry(result: COCONUTResult) -> DatabaseEntry:
    """Build a DatabaseEntry from a COCONUT result."""
    if not result.found:
        return DatabaseEntry(database="COCONUT", found=False)
    return DatabaseEntry(
        database="COCONUT",
        found=True,
        canonical_smiles=result.smiles,
        inchi=result.inchi,
        inchikey=result.inchikey,
        molecular_formula=result.molecular_formula,
        molecular_weight=result.molecular_weight,
        name=result.name,
        url=result.url,
    )


async def _lookup_wikidata(inchikey: Optional[str] = None) -> Optional[dict]:
    """Lookup compound in Wikidata by InChIKey.

    Returns:
        Dict with SMILES, formula, MW, InChIKey, name, or None on failure.
    """
    if not inchikey:
        return None
    try:
        return await WikidataClient().resolve_from_inchikey(inchikey)
    except Exception:
        return None


def _build_wikidata_entry(result: Optional[dict]) -> DatabaseEntry:
    """Build a DatabaseEntry from a Wikidata result dict."""
    if not result:
        return DatabaseEntry(database="Wikidata", found=False)
    return DatabaseEntry(
        database="Wikidata",
        found=True,
        canonical_smiles=result.get("smiles"),
        inchi=result.get("inchi"),
        inchikey=result.get("inchikey"),
        molecular_formula=result.get("formula"),
        molecular_weight=result.get("mass"),
        name=result.get("label"),
        url=result.get("url"),
    )


def _build_resolved_entry(smiles: Optional[str]) -> DatabaseEntry:
    """Build a 'Resolved' reference entry from the user's input SMILES using RDKit.

    This entry represents the locally-computed canonical identifiers (SMILES,
    InChI, InChIKey) so users can compare what RDKit produces against what each
    external database stores.
    """
    if not smiles:
        return DatabaseEntry(database="Resolved", found=False)
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return DatabaseEntry(database="Resolved", found=False)
        canonical = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)
        Chem.Kekulize(mol)
        kekulized = Chem.MolToSmiles(mol, kekuleSmiles=True)
        return DatabaseEntry(
            database="Resolved",
            found=True,
            canonical_smiles=canonical,
            kekulized_smiles=kekulized,
            inchi=inchi,
            inchikey=inchikey,
        )
    except Exception:
        return DatabaseEntry(database="Resolved", found=False)


async def compare_across_databases(
    smiles: Optional[str] = None,
    inchikey: Optional[str] = None,
) -> ConsistencyResult:
    """Compare compound representation across PubChem, ChEMBL, COCONUT, and Wikidata.

    Args:
        smiles: SMILES string of the compound.
        inchikey: InChIKey of the compound.

    Returns:
        ConsistencyResult with per-database entries, field comparisons, and verdict.
    """
    if not smiles and not inchikey:
        return ConsistencyResult(overall_verdict="no_data", summary="No input provided")

    # Build a "Resolved" reference entry from the input using RDKit
    resolved_entry = _build_resolved_entry(smiles)

    # Generate InChIKey from SMILES if not provided
    if smiles and not inchikey:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                inchikey = Chem.MolToInchiKey(mol)
        except Exception:
            pass

    # Query all databases in parallel
    pc_result, chembl_result, coconut_result, wikidata_result = await asyncio.gather(
        get_compound_info(PubChemRequest(smiles=smiles, inchikey=inchikey)),
        get_bioactivity(ChEMBLRequest(smiles=smiles, inchikey=inchikey)),
        lookup_natural_product(COCONUTRequest(smiles=smiles, inchikey=inchikey)),
        _lookup_wikidata(inchikey=inchikey),
        return_exceptions=True,
    )

    # Handle exceptions gracefully
    if isinstance(pc_result, Exception):
        pc_result = PubChemResult(found=False)
    if isinstance(chembl_result, Exception):
        chembl_result = ChEMBLResult(found=False)
    if isinstance(coconut_result, Exception):
        coconut_result = COCONUTResult(found=False)
    if isinstance(wikidata_result, Exception):
        wikidata_result = None

    # Build database entries — external databases first, Resolved last
    entries: list[DatabaseEntry] = [
        _build_pubchem_entry(pc_result),
        _build_chembl_entry(chembl_result),
        _build_coconut_entry(coconut_result),
        _build_wikidata_entry(wikidata_result),
        resolved_entry,
    ]

    found_entries = [e for e in entries if e.found]
    # External database entries only (exclude Resolved for verdict logic)
    external_found = [e for e in found_entries if e.database != "Resolved"]

    if len(external_found) == 0:
        return ConsistencyResult(
            entries=entries,
            overall_verdict="no_data",
            summary="Compound not found in any database",
        )

    if len(external_found) == 1 and not resolved_entry.found:
        db_name = external_found[0].database
        return ConsistencyResult(
            entries=entries,
            overall_verdict="consistent",
            summary=f"Found only in {db_name} — no cross-database comparison possible",
        )

    # Build comparisons for structural identifiers only (SMILES, InChIKey, InChI).
    # Formula and MW are excluded: databases use different formats (Wikidata
    # stores Unicode subscripts, monoisotopic mass) which create false mismatches,
    # and these fields are not displayed in the comparison panel.
    comparisons: list[PropertyComparison] = []

    # SMILES — specialized comparator with RDKit canonicalization detail
    smiles_values = {e.database: e.canonical_smiles for e in found_entries}
    if any(v is not None for v in smiles_values.values()):
        comparisons.append(_compare_smiles(smiles_values))

    # InChIKey — specialized comparator with layer analysis
    ik_values = {e.database: e.inchikey for e in found_entries}
    if any(v is not None for v in ik_values.values()):
        comparisons.append(_compare_inchikey(ik_values))

    # InChI — specialized comparator with layer analysis
    inchi_values = {e.database: e.inchi for e in found_entries}
    if any(v is not None for v in inchi_values.values()):
        comparisons.append(_compare_inchi(inchi_values))

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
