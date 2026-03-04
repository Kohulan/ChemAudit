"""
UniChem integration client.

UniChem is a large-scale non-redundant mapping of pointers between chemical
structures and EMBL-EBI chemistry resources. It provides cross-references
between PubChem, ChEMBL, DrugBank, ChEBI, KEGG, and many more databases.
https://www.ebi.ac.uk/unichem/
"""

from typing import Optional

import httpx

from app.core.config import settings

# UniChem source IDs for databases we care about
SRC_CHEMBL = 1
SRC_DRUGBANK = 2
SRC_PUBCHEM = 22
SRC_CHEBI = 7
SRC_KEGG = 6

_SRC_TO_KEY = {
    SRC_CHEMBL: "chembl_id",
    SRC_DRUGBANK: "drugbank_id",
    SRC_PUBCHEM: "pubchem_cid",
    SRC_CHEBI: "chebi_id",
    SRC_KEGG: "kegg_id",
}

UNICHEM_BASE_URL = "https://www.ebi.ac.uk/unichem/rest"


class UniChemClient:
    """Client for UniChem REST API for cross-database ID mapping."""

    def __init__(self):
        self.base_url = UNICHEM_BASE_URL
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    async def get_cross_references(self, inchikey: str) -> dict[str, Optional[str]]:
        """Get cross-database IDs for a compound by InChIKey.

        Args:
            inchikey: Standard InChIKey (27 chars, e.g. BSYNRYMUTXBXSQ-UHFFFAOYSA-N)

        Returns:
            Dict with keys: chembl_id, drugbank_id, pubchem_cid, chebi_id, kegg_id.
            Values are None if not found in that database.
        """
        result: dict[str, Optional[str]] = {v: None for v in _SRC_TO_KEY.values()}

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/inchikey/{inchikey}",
                )
                response.raise_for_status()
                data = response.json()

            if not isinstance(data, list):
                return result

            for entry in data:
                src_id_raw = entry.get("src_id")
                if src_id_raw is None:
                    continue
                src_id = int(src_id_raw)
                key = _SRC_TO_KEY.get(src_id)
                if key:
                    compound_id = entry.get("src_compound_id")
                    if key == "chembl_id" and compound_id:
                        if not compound_id.startswith("CHEMBL"):
                            compound_id = f"CHEMBL{compound_id}"
                    if key == "chebi_id" and compound_id:
                        if not compound_id.startswith("CHEBI:"):
                            compound_id = f"CHEBI:{compound_id}"
                    result[key] = compound_id

        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            pass

        return result

    async def resolve_to_inchikey(self, compound_id: str, src_id: int) -> Optional[str]:
        """Resolve a database-specific ID to an InChIKey via UniChem.

        Args:
            compound_id: The compound ID in the source database.
            src_id: UniChem source ID (e.g. 2 for DrugBank).

        Returns:
            InChIKey string if resolved, None otherwise.
        """
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    f"{self.base_url}/src_compound_id/{compound_id}/{src_id}",
                )
                response.raise_for_status()
                data = response.json()

            if isinstance(data, list) and len(data) > 0:
                inchikey = data[0].get("src_compound_id")
                if inchikey and len(inchikey) == 27 and "-" in inchikey:
                    return inchikey

        except (httpx.HTTPError, KeyError, ValueError, IndexError):
            pass

        return None
