"""
COCONUT (COlleCtion of Open Natural ProdUcTs) database integration client.

COCONUT is a comprehensive database of natural products with >400,000 entries.
https://coconut.naturalproducts.net/

This client provides lookup functionality for natural product information.
"""
import httpx
from typing import Optional
from rdkit import Chem

from app.core.config import settings
from app.schemas.integrations import COCONUTRequest, COCONUTResult


class COCONUTClient:
    """
    Client for COCONUT natural products database API.

    Provides search and lookup functionality for natural products.
    """

    def __init__(self):
        self.base_url = settings.COCONUT_API_URL
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    async def search_by_smiles(self, smiles: str) -> Optional[dict]:
        """
        Search COCONUT by SMILES string.

        Args:
            smiles: SMILES string to search

        Returns:
            Compound data dict if found, None otherwise
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/search/smiles",
                    params={"smiles": smiles},
                )
                response.raise_for_status()

                data = response.json()
                if data and len(data) > 0:
                    return data[0]  # Return first match
                return None

        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            # External API failure - return None gracefully
            return None

    async def search_by_inchikey(self, inchikey: str) -> Optional[dict]:
        """
        Search COCONUT by InChIKey.

        Args:
            inchikey: InChIKey to search

        Returns:
            Compound data dict if found, None otherwise
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/search/inchikey/{inchikey}",
                )
                response.raise_for_status()

                data = response.json()
                return data if data else None

        except (httpx.HTTPError, KeyError, ValueError):
            # External API failure - return None gracefully
            return None

    async def get_compound(self, coconut_id: str) -> Optional[dict]:
        """
        Get compound by COCONUT ID.

        Args:
            coconut_id: COCONUT compound ID

        Returns:
            Compound data dict if found, None otherwise
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/compound/{coconut_id}",
                )
                response.raise_for_status()

                return response.json()

        except (httpx.HTTPError, KeyError, ValueError):
            # External API failure - return None gracefully
            return None


async def lookup_natural_product(request: COCONUTRequest) -> COCONUTResult:
    """
    Look up molecule in COCONUT natural products database.

    Searches by InChIKey first (more specific), falls back to SMILES.

    Args:
        request: COCONUT lookup request with SMILES or InChIKey

    Returns:
        Natural product information if found
    """
    client = COCONUTClient()

    # Try InChIKey first (most specific)
    if request.inchikey:
        data = await client.search_by_inchikey(request.inchikey)
        if data:
            return _parse_coconut_result(data, found=True)

    # Try SMILES
    if request.smiles:
        # Generate InChIKey from SMILES for more reliable search
        try:
            mol = Chem.MolFromSmiles(request.smiles)
            if mol:
                inchikey = Chem.MolToInchiKey(mol)
                data = await client.search_by_inchikey(inchikey)
                if data:
                    return _parse_coconut_result(data, found=True)

                # Fallback to SMILES search
                data = await client.search_by_smiles(request.smiles)
                if data:
                    return _parse_coconut_result(data, found=True)
        except Exception:
            # RDKit error - try direct SMILES search anyway
            data = await client.search_by_smiles(request.smiles)
            if data:
                return _parse_coconut_result(data, found=True)

    # Not found
    return COCONUTResult(found=False)


def _parse_coconut_result(data: dict, found: bool) -> COCONUTResult:
    """Parse COCONUT API response into result schema."""
    coconut_id = data.get("coconut_id")

    return COCONUTResult(
        found=found,
        coconut_id=coconut_id,
        name=data.get("name") or data.get("iupac_name"),
        smiles=data.get("smiles") or data.get("canonical_smiles"),
        inchikey=data.get("inchikey"),
        molecular_formula=data.get("molecular_formula"),
        molecular_weight=data.get("molecular_weight"),
        organism=data.get("organism") or data.get("biological_source"),
        organism_type=data.get("organism_type"),
        nplikeness=data.get("nplikeness") or data.get("np_likeness"),
        url=f"https://coconut.naturalproducts.net/compound/{coconut_id}" if coconut_id else None,
    )
