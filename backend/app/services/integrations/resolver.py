"""
Universal Identifier Resolution orchestrator.

Accepts any chemical identifier (SMILES, InChI, InChIKey, PubChem CID, ChEMBL ID,
CAS number, DrugBank ID, ChEBI ID, UNII, Wikipedia URL, or compound name) and
resolves it to a canonical structure with cross-database references.
"""

import logging
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from app.schemas.integrations import CrossReferences, ResolvedCompound
from app.services.integrations.chebi import ChEBIClient
from app.services.integrations.chembl import ChEMBLClient
from app.services.integrations.identifier_detect import detect_identifier_type
from app.services.integrations.pubchem import PubChemClient
from app.services.integrations.unichem import SRC_DRUGBANK, UniChemClient
from app.services.integrations.wikidata import WikidataClient

logger = logging.getLogger(__name__)


def _canonicalize(smiles: str) -> Optional[dict]:
    """Canonicalize SMILES via RDKit and compute basic properties.

    Returns dict with canonical_smiles, inchi, inchikey, molecular_formula,
    molecular_weight, or None if invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        canonical = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol) if inchi else None
        return {
            "canonical_smiles": canonical,
            "inchi": inchi,
            "inchikey": inchikey,
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.ExactMolWt(mol), 4),
        }
    except Exception:
        return None


async def _get_xrefs(inchikey: Optional[str]) -> dict:
    """Fetch cross-references from UniChem by InChIKey, returning empty dict on failure."""
    if not inchikey:
        return {}
    unichem = UniChemClient()
    return await unichem.get_cross_references(inchikey)


def _build_xrefs(xrefs_data: dict, **overrides: object) -> CrossReferences:
    """Build CrossReferences from UniChem data with optional field overrides."""
    pubchem_cid_raw = overrides.pop("pubchem_cid", xrefs_data.get("pubchem_cid"))
    pubchem_cid = int(pubchem_cid_raw) if pubchem_cid_raw else None
    return CrossReferences(
        pubchem_cid=pubchem_cid,
        chembl_id=overrides.pop("chembl_id", xrefs_data.get("chembl_id")),
        drugbank_id=overrides.pop("drugbank_id", xrefs_data.get("drugbank_id")),
        chebi_id=overrides.pop("chebi_id", xrefs_data.get("chebi_id")),
        kegg_id=overrides.pop("kegg_id", xrefs_data.get("kegg_id")),
        **overrides,
    )


async def _resolve_from_pubchem_properties(
    properties: dict,
    *,
    id_type: str,
    cid: int,
    chain_label: str,
    confidence: str = "high",
) -> ResolvedCompound:
    """Build a ResolvedCompound from PubChem properties dict.

    Shared by _resolve_inchikey, _resolve_pubchem_cid, and _resolve_via_pubchem_name.
    """
    smiles = properties.get("CanonicalSMILES")
    props = _canonicalize(smiles) if smiles else None
    inchikey = properties.get("InChIKey") or (props["inchikey"] if props else None)
    xrefs_data = await _get_xrefs(inchikey)

    return ResolvedCompound(
        resolved=True,
        identifier_type_detected=id_type,
        canonical_smiles=props["canonical_smiles"] if props else smiles,
        inchi=properties.get("InChI") or (props["inchi"] if props else None),
        inchikey=inchikey,
        molecular_formula=properties.get("MolecularFormula")
        or (props["molecular_formula"] if props else None),
        molecular_weight=properties.get("MolecularWeight"),
        iupac_name=properties.get("IUPACName"),
        resolution_source="pubchem",
        resolution_chain=[chain_label],
        cross_references=_build_xrefs(xrefs_data, pubchem_cid=cid),
        confidence=confidence,
    )


async def _resolve_smiles(smiles: str) -> Optional[ResolvedCompound]:
    """Resolve a SMILES string locally with RDKit."""
    props = _canonicalize(smiles)
    if not props:
        return None
    return ResolvedCompound(
        resolved=True,
        identifier_type_detected="smiles",
        canonical_smiles=props["canonical_smiles"],
        inchi=props["inchi"],
        inchikey=props["inchikey"],
        molecular_formula=props["molecular_formula"],
        molecular_weight=props["molecular_weight"],
        resolution_source="rdkit",
        resolution_chain=["SMILES → RDKit canonicalization"],
        confidence="high",
    )


async def _resolve_inchi(inchi: str) -> Optional[ResolvedCompound]:
    """Resolve an InChI string locally with RDKit."""
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol)
        props = _canonicalize(smiles)
        if not props:
            return None
        return ResolvedCompound(
            resolved=True,
            identifier_type_detected="inchi",
            canonical_smiles=props["canonical_smiles"],
            inchi=inchi,
            inchikey=props["inchikey"],
            molecular_formula=props["molecular_formula"],
            molecular_weight=props["molecular_weight"],
            resolution_source="rdkit",
            resolution_chain=["InChI → RDKit → SMILES"],
            confidence="high",
        )
    except Exception:
        return None


async def _resolve_inchikey(inchikey: str) -> Optional[ResolvedCompound]:
    """Resolve an InChIKey via PubChem + UniChem cross-references."""
    pc = PubChemClient()
    cid = await pc.search_by_inchikey(inchikey)
    if cid is None:
        return None

    properties = await pc.get_compound_properties(cid)
    if not properties:
        return None

    return await _resolve_from_pubchem_properties(
        properties,
        id_type="inchikey",
        cid=cid,
        chain_label=f"InChIKey → PubChem CID:{cid} → SMILES",
    )


async def _resolve_chembl_id(chembl_id: str) -> Optional[ResolvedCompound]:
    """Resolve a ChEMBL ID via ChEMBL API."""
    client = ChEMBLClient()
    mol_data = await client.get_molecule(chembl_id.upper())
    if not mol_data:
        return None

    structures = mol_data.get("molecule_structures") or {}
    smiles = structures.get("canonical_smiles")
    inchi = structures.get("standard_inchi")
    inchikey = structures.get("standard_inchi_key")

    props = _canonicalize(smiles) if smiles else None
    xrefs_data = await _get_xrefs(inchikey)
    mol_props = mol_data.get("molecule_properties") or {}

    mw_raw = mol_props.get("molecular_weight")
    molecular_weight = float(mw_raw) if mw_raw else (props["molecular_weight"] if props else None)

    return ResolvedCompound(
        resolved=True,
        identifier_type_detected="chembl_id",
        canonical_smiles=props["canonical_smiles"] if props else smiles,
        inchi=inchi,
        inchikey=inchikey,
        molecular_formula=mol_props.get("full_molecular_formula")
        or (props["molecular_formula"] if props else None),
        molecular_weight=molecular_weight,
        iupac_name=mol_data.get("pref_name"),
        resolution_source="chembl",
        resolution_chain=[f"{chembl_id.upper()} → ChEMBL API → SMILES"],
        cross_references=_build_xrefs(xrefs_data, chembl_id=chembl_id.upper()),
        confidence="high",
    )


async def _resolve_pubchem_cid(cid_str: str) -> Optional[ResolvedCompound]:
    """Resolve a PubChem CID."""
    clean = cid_str.upper().replace("CID:", "").strip()
    try:
        cid = int(clean)
    except ValueError:
        return None

    pc = PubChemClient()
    properties = await pc.get_compound_properties(cid)
    if not properties:
        return None

    return await _resolve_from_pubchem_properties(
        properties,
        id_type="pubchem_cid",
        cid=cid,
        chain_label=f"CID:{cid} → PubChem API → SMILES",
    )


async def _resolve_via_pubchem_name(identifier: str, id_type: str) -> Optional[ResolvedCompound]:
    """Resolve CAS, UNII, or similar via PubChem name search."""
    pc = PubChemClient()
    cid = await pc.search_by_name(identifier)
    if cid is None:
        return None

    properties = await pc.get_compound_properties(cid)
    if not properties:
        return None

    return await _resolve_from_pubchem_properties(
        properties,
        id_type=id_type,
        cid=cid,
        chain_label=f"{identifier} → PubChem name search → CID:{cid} → SMILES",
        confidence="medium",
    )


async def _resolve_drugbank(drugbank_id: str) -> Optional[ResolvedCompound]:
    """Resolve a DrugBank ID via UniChem → InChIKey → PubChem."""
    unichem = UniChemClient()
    inchikey = await unichem.resolve_to_inchikey(drugbank_id, src_id=SRC_DRUGBANK)
    if not inchikey:
        return None

    result = await _resolve_inchikey(inchikey)
    if result:
        result.identifier_type_detected = "drugbank_id"
        result.resolution_chain = [
            f"{drugbank_id} → UniChem → InChIKey:{inchikey} → PubChem → SMILES"
        ]
        result.cross_references.drugbank_id = drugbank_id
    return result


async def _resolve_chebi(chebi_id: str) -> Optional[ResolvedCompound]:
    """Resolve a ChEBI ID via ChEBI REST API."""
    clean_id = chebi_id.replace("CHEBI:", "").strip()
    client = ChEBIClient()
    data = await client.get_compound(clean_id)
    if not data or not data.get("smiles"):
        return None

    props = _canonicalize(data["smiles"])
    if not props:
        return None

    inchikey = data.get("inchikey") or props["inchikey"]
    xrefs_data = await _get_xrefs(inchikey)

    return ResolvedCompound(
        resolved=True,
        identifier_type_detected="chebi_id",
        canonical_smiles=props["canonical_smiles"],
        inchi=data.get("inchi") or props["inchi"],
        inchikey=inchikey,
        molecular_formula=props["molecular_formula"],
        molecular_weight=props["molecular_weight"],
        iupac_name=data.get("name"),
        resolution_source="chebi",
        resolution_chain=[f"CHEBI:{clean_id} → ChEBI API → SMILES"],
        cross_references=_build_xrefs(xrefs_data, chebi_id=f"CHEBI:{clean_id}"),
        confidence="high",
    )


async def _resolve_wikipedia(url: str) -> Optional[ResolvedCompound]:
    """Resolve a Wikipedia URL via Wikidata SPARQL."""
    client = WikidataClient()
    data = await client.resolve_from_wikipedia(url)
    if not data or not data.get("smiles"):
        return None

    props = _canonicalize(data["smiles"])
    if not props:
        return None

    return ResolvedCompound(
        resolved=True,
        identifier_type_detected="wikipedia",
        canonical_smiles=props["canonical_smiles"],
        inchi=data.get("inchi") or props["inchi"],
        inchikey=data.get("inchikey") or props["inchikey"],
        molecular_formula=props["molecular_formula"],
        molecular_weight=props["molecular_weight"],
        iupac_name=data.get("label"),
        resolution_source="wikidata",
        resolution_chain=["Wikipedia → Wikidata SPARQL → SMILES"],
        cross_references=CrossReferences(
            wikipedia_url=url,
            cas=data.get("cas"),
        ),
        confidence="medium",
    )


async def _resolve_name(name: str) -> Optional[ResolvedCompound]:
    """Resolve a compound name via OPSIN/PubChem."""
    try:
        from app.services.iupac.converter import name_to_smiles

        smiles, source = name_to_smiles(name)
        if smiles:
            props = _canonicalize(smiles)
            if props:
                return ResolvedCompound(
                    resolved=True,
                    identifier_type_detected="name",
                    canonical_smiles=props["canonical_smiles"],
                    inchi=props["inchi"],
                    inchikey=props["inchikey"],
                    molecular_formula=props["molecular_formula"],
                    molecular_weight=props["molecular_weight"],
                    iupac_name=name,
                    resolution_source=source or "pubchem",
                    resolution_chain=[f'"{name}" → {source or "pubchem"} → SMILES'],
                    confidence="medium",
                )
    except Exception:
        pass

    # Fallback: try PubChem name search directly
    return await _resolve_via_pubchem_name(name, "name")


# Dispatcher mapping identifier types to resolver functions.
# Functions that accept only one arg can be referenced directly; those needing
# a fixed second arg use a lambda for the partial application.
_RESOLVERS = {
    "smiles": _resolve_smiles,
    "inchi": _resolve_inchi,
    "inchikey": _resolve_inchikey,
    "chembl_id": _resolve_chembl_id,
    "pubchem_cid": _resolve_pubchem_cid,
    "cas": lambda ident: _resolve_via_pubchem_name(ident, "cas"),
    "unii": lambda ident: _resolve_via_pubchem_name(ident, "unii"),
    "drugbank_id": _resolve_drugbank,
    "chebi_id": _resolve_chebi,
    "wikipedia": _resolve_wikipedia,
    "name": _resolve_name,
}


async def resolve_identifier(
    identifier: str, identifier_type: str = "auto"
) -> ResolvedCompound:
    """Resolve any chemical identifier to a canonical structure.

    Args:
        identifier: The raw input string (SMILES, name, CAS, ChEMBL ID, etc.)
        identifier_type: Explicit type hint, or "auto" for auto-detection.

    Returns:
        ResolvedCompound with canonical SMILES, properties, and cross-references.
        If resolution fails, resolved=False.
    """
    detected_type = detect_identifier_type(identifier, identifier_type)
    resolver = _RESOLVERS.get(detected_type)

    if not resolver:
        return ResolvedCompound(
            resolved=False, identifier_type_detected=detected_type
        )

    try:
        result = await resolver(identifier)
        if result:
            return result
    except Exception:
        logger.exception("Resolution failed for identifier: %s (type: %s)", identifier, detected_type)

    return ResolvedCompound(resolved=False, identifier_type_detected=detected_type)
