"""
SureChEMBL patent presence lookup client.

SureChEMBL is a freely available database of chemical compounds extracted
from patent literature. This client uses a UniChem-first strategy to resolve
InChIKeys to SCHEMBL identifiers, then optionally enriches from the
SureChEMBL direct API for patent count data.

https://www.surechembl.org/
"""

import logging
from typing import Optional

import httpx
from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi

from app.core.config import settings
from app.schemas.integrations import SureChEMBLRequest, SureChEMBLResult

logger = logging.getLogger(__name__)

SURECHEMBL_BASE_URL = "https://www.surechembl.org/api"
UNICHEM_BASE_URL = "https://www.ebi.ac.uk/unichem/rest"
SRC_SURECHEMBL = 15  # UniChem source ID for SureChEMBL
SURECHEMBL_ENTRY_URL = "https://www.surechembl.org/chemical"


class SureChEMBLClient:
    """
    Client for SureChEMBL patent presence lookups.

    Uses a two-stage strategy:
    1. UniChem cross-reference to get SCHEMBL compound ID from InChIKey
    2. SureChEMBL direct API for enrichment (patent count, properties)

    Falls back to UniChem-only result when the SureChEMBL API is unavailable.
    """

    def __init__(self):
        self.surechembl_url = SURECHEMBL_BASE_URL
        self.unichem_url = UNICHEM_BASE_URL
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    async def _get_schembl_id_via_unichem(self, inchikey: str) -> Optional[str]:
        """
        Resolve InChIKey to SCHEMBL compound ID via UniChem cross-reference.

        Args:
            inchikey: Standard InChIKey (27 characters).

        Returns:
            SCHEMBL compound ID if found in UniChem, None otherwise.
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(f"{self.unichem_url}/inchikey/{inchikey}")
                response.raise_for_status()
                data = response.json()

            if not isinstance(data, list):
                return None

            for entry in data:
                if str(entry.get("src_id")) == str(SRC_SURECHEMBL):
                    return entry.get("src_compound_id")

            return None

        except (httpx.HTTPError, KeyError, ValueError):
            logger.debug("UniChem lookup failed for %s", inchikey)
            return None

    async def _get_compound_data(self, schembl_id: str) -> Optional[dict]:
        """
        Fetch compound data from SureChEMBL direct API.

        Args:
            schembl_id: SCHEMBL compound identifier (e.g. "SCHEMBL12345").

        Returns:
            Compound data dict if API returns status OK, None otherwise.
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.surechembl_url}/chemical/id/{schembl_id}"
                )
                response.raise_for_status()
                data = response.json()

            if data.get("status") == "OK":
                return data.get("data")
            return None

        except (httpx.HTTPError, KeyError, ValueError):
            logger.debug("SureChEMBL API unavailable for %s", schembl_id)
            return None

    async def _build_rich_result(self, schembl_id: str, compound_data: dict) -> dict:
        """
        Build enriched result from SureChEMBL compound data.

        Args:
            schembl_id: SCHEMBL compound identifier.
            compound_data: Raw compound data from SureChEMBL API.

        Returns:
            Dict with found, schembl_id, url, patent_count, source, and
            optional structural fields.
        """
        patent_count = int(compound_data.get("global_frequency", 0))

        result: dict = {
            "found": True,
            "schembl_id": schembl_id,
            "url": f"{SURECHEMBL_ENTRY_URL}/{schembl_id}",
            "patent_count": patent_count,
            "source": "surechembl_api",
        }

        # Include optional structural fields if present
        for field in ("smiles", "inchi", "inchikey", "molecular_weight"):
            value = compound_data.get(field)
            if value is not None:
                result[field] = value

        return result

    async def lookup_by_inchikey(self, inchikey: str) -> dict:
        """
        Look up a compound in SureChEMBL by InChIKey.

        Strategy:
        1. Resolve InChIKey to SCHEMBL ID via UniChem
        2. Normalize SCHEMBL ID prefix
        3. Enrich from SureChEMBL direct API if available
        4. Fall back to UniChem-only result if API unavailable

        Args:
            inchikey: Standard InChIKey (27 characters).

        Returns:
            Dict with at minimum 'found' key. When found, includes
            schembl_id, url, and source. May include patent_count
            and structural fields when SureChEMBL API is available.
        """
        # Step 1: Resolve via UniChem
        schembl_id = await self._get_schembl_id_via_unichem(inchikey)
        if schembl_id is None:
            return {"found": False}

        # Step 2: Normalize SCHEMBL ID prefix
        if schembl_id.isdigit():
            schembl_id = f"SCHEMBL{schembl_id}"

        # Step 3: Try SureChEMBL direct API for enrichment
        compound_data = await self._get_compound_data(schembl_id)
        if compound_data is not None:
            return await self._build_rich_result(schembl_id, compound_data)

        # Step 4: Fallback to UniChem-only result
        return {
            "found": True,
            "schembl_id": schembl_id,
            "url": f"{SURECHEMBL_ENTRY_URL}/{schembl_id}",
            "source": "unichem_only",
        }


async def lookup_surechembl(request: SureChEMBLRequest) -> SureChEMBLResult:
    """
    Look up molecule in SureChEMBL patent database.

    Convenience function wrapping SureChEMBLClient. Accepts SMILES or InChIKey
    input, deriving InChIKey from SMILES if needed.

    Args:
        request: SureChEMBL lookup request with SMILES or InChIKey.

    Returns:
        SureChEMBLResult with patent presence information.
    """
    inchikey = request.inchikey

    # Derive InChIKey from SMILES if not provided directly
    if not inchikey and request.smiles:
        try:
            mol = Chem.MolFromSmiles(request.smiles)
            if mol is None:
                return SureChEMBLResult(found=False)
            inchikey = rdkit_inchi.MolToInchiKey(mol)
        except Exception:
            return SureChEMBLResult(found=False)

    if not inchikey:
        return SureChEMBLResult(found=False)

    client = SureChEMBLClient()
    result = await client.lookup_by_inchikey(inchikey)
    return SureChEMBLResult(**result)
