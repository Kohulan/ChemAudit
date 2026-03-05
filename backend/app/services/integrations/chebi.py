"""
ChEBI integration client.

ChEBI (Chemical Entities of Biological Interest) is a dictionary of molecular
entities focused on small chemical compounds. Maintained by EMBL-EBI.
https://www.ebi.ac.uk/chebi/
"""

import re
from typing import Optional

import httpx

from app.core.config import settings

CHEBI_REST_URL = "https://www.ebi.ac.uk/chebi/searchId.do"


class ChEBIClient:
    """Client for ChEBI REST API."""

    def __init__(self):
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    async def get_compound(self, chebi_id: str) -> Optional[dict]:
        """Get compound data by ChEBI ID.

        Args:
            chebi_id: Numeric ChEBI ID (without "CHEBI:" prefix).

        Returns:
            Dict with keys: smiles, inchikey, inchi, name. None if not found.
        """
        clean_id = chebi_id.replace("CHEBI:", "").strip()

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    CHEBI_REST_URL,
                    params={"chebiId": f"CHEBI:{clean_id}"},
                    headers={"Accept": "application/xml"},
                )
                response.raise_for_status()
                xml_text = response.text

            # Simple XML extraction (avoid heavy lxml dependency)
            smiles = self._extract_xml_value(xml_text, "smiles")
            inchikey = self._extract_xml_value(xml_text, "inchiKey")
            inchi = self._extract_xml_value(xml_text, "inchi")
            name = self._extract_xml_value(xml_text, "chebiAsciiName")

            if not smiles and not inchikey:
                return None

            return {
                "smiles": smiles,
                "inchikey": inchikey,
                "inchi": inchi,
                "name": name,
            }

        except (httpx.HTTPError, KeyError, ValueError):
            return None

    @staticmethod
    def _extract_xml_value(xml_text: str, tag: str) -> Optional[str]:
        """Extract a value from a simple XML tag."""
        match = re.search(rf"<{tag}>(.*?)</{tag}>", xml_text)
        return match.group(1) if match else None
