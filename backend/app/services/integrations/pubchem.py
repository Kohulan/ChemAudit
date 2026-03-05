"""
PubChem integration client.

PubChem is a public database of chemical compounds and their properties,
maintained by the NIH National Library of Medicine.
https://pubchem.ncbi.nlm.nih.gov/

This client provides compound lookup and cross-reference functionality.
"""

from typing import Optional

import httpx
from rdkit import Chem

from app.core.config import settings
from app.schemas.integrations import PubChemRequest, PubChemResult


class PubChemClient:
    """Client for PubChem PUG REST API."""

    def __init__(self):
        self.base_url = settings.PUBCHEM_API_URL
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    def _extract_first_cid(self, data: dict) -> Optional[int]:
        """Extract the first CID from a PubChem identifier response."""
        cids = data.get("IdentifierList", {}).get("CID", [])
        # PubChem returns CID 0 to indicate "not found" — treat as None
        return cids[0] if cids and cids[0] != 0 else None

    async def search_by_smiles(self, smiles: str) -> Optional[int]:
        """Search PubChem by SMILES string, returning the first CID or None.

        Tries POST first, then falls back to canonical SMILES via GET.
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # POST with raw SMILES
                response = await client.post(
                    f"{self.base_url}/compound/smiles/cids/JSON",
                    data={"smiles": smiles},
                )
                response.raise_for_status()
                cid = self._extract_first_cid(response.json())
                if cid is not None:
                    return cid

                # Fallback: canonicalize with RDKit then GET
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canon = Chem.MolToSmiles(mol)
                    resp2 = await client.get(
                        f"{self.base_url}/compound/smiles/{canon}/cids/JSON",
                    )
                    resp2.raise_for_status()
                    return self._extract_first_cid(resp2.json())
        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            pass
        return None

    async def search_by_name(self, name: str) -> Optional[int]:
        """Search PubChem by compound name (common, trade, IUPAC, CAS, UNII, etc.)."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/compound/name/{name}/cids/JSON",
                )
                response.raise_for_status()
                return self._extract_first_cid(response.json())
        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            return None

    async def search_by_inchikey(self, inchikey: str) -> Optional[int]:
        """Search PubChem by InChIKey, returning the first CID or None."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/compound/inchikey/{inchikey}/cids/JSON",
                )
                response.raise_for_status()
                return self._extract_first_cid(response.json())
        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            return None

    async def get_compound_properties(self, cid: int) -> Optional[dict]:
        """Get compound properties by CID."""
        try:
            prop_names = (
                "MolecularFormula,MolecularWeight,"
                "IsomericSMILES,CanonicalSMILES,"
                "InChI,InChIKey,IUPACName"
            )
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/compound/cid/{cid}/property/{prop_names}/JSON",
                )
                response.raise_for_status()

                props_list = response.json().get("PropertyTable", {}).get("Properties", [])
                return props_list[0] if props_list else None
        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            return None

    async def get_synonyms(self, cid: int, max_synonyms: int = 10) -> list[str]:
        """Get compound synonyms by CID (up to max_synonyms)."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/compound/cid/{cid}/synonyms/JSON",
                )
                response.raise_for_status()

                synonyms = (
                    response.json()
                    .get("InformationList", {})
                    .get("Information", [{}])[0]
                    .get("Synonym", [])
                )
                return synonyms[:max_synonyms]
        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            return []


async def get_compound_info(request: PubChemRequest) -> PubChemResult:
    """Get compound information from PubChem.

    Searches by InChIKey first (most specific), falls back to SMILES.
    """
    client = PubChemClient()
    cid = None

    if request.inchikey:
        cid = await client.search_by_inchikey(request.inchikey)

    if cid is None and request.smiles:
        # Generate InChIKey from SMILES for more reliable search
        try:
            mol = Chem.MolFromSmiles(request.smiles)
            if mol:
                inchikey = Chem.MolToInchiKey(mol)
                cid = await client.search_by_inchikey(inchikey)
                if cid is None:
                    cid = await client.search_by_smiles(request.smiles)
        except Exception:
            cid = await client.search_by_smiles(request.smiles)

    if cid is None:
        return PubChemResult(found=False)

    properties = await client.get_compound_properties(cid)
    if properties is None:
        return PubChemResult(found=False)

    synonyms = await client.get_synonyms(cid, max_synonyms=10)

    # Prefer isomeric SMILES (preserves stereochemistry) over canonical/connectivity.
    # PubChem API key renames: IsomericSMILES → SMILES, CanonicalSMILES → ConnectivitySMILES.
    # Accept both old and new key names for resilience.
    smiles_value = (
        properties.get("IsomericSMILES")
        or properties.get("SMILES")
        or properties.get("CanonicalSMILES")
        or properties.get("ConnectivitySMILES")
    )

    return PubChemResult(
        found=True,
        cid=cid,
        iupac_name=properties.get("IUPACName"),
        molecular_formula=properties.get("MolecularFormula"),
        molecular_weight=properties.get("MolecularWeight"),
        canonical_smiles=smiles_value,
        inchi=properties.get("InChI"),
        inchikey=properties.get("InChIKey"),
        synonyms=synonyms or None,
        url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
    )
