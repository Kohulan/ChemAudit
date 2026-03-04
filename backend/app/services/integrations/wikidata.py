"""
Wikidata integration client.

Resolves Wikipedia URLs to chemical structure data via Wikidata SPARQL.
Wikidata stores structured chemical data (SMILES, InChIKey, CAS) linked
to Wikipedia articles.
"""

import re
import urllib.parse
from typing import Optional

import httpx

from app.core.config import settings

WIKIDATA_SPARQL_URL = "https://query.wikidata.org/sparql"


class WikidataClient:
    """Client for Wikidata SPARQL endpoint for chemical data lookup."""

    def __init__(self):
        self.timeout = settings.EXTERNAL_API_TIMEOUT

    def _extract_title(self, url: str) -> Optional[str]:
        """Extract article title from a Wikipedia URL."""
        match = re.search(r"wikipedia\.org/wiki/(.+?)(?:#.*)?$", url)
        if not match:
            return None
        return urllib.parse.unquote(match.group(1)).replace("_", " ")

    async def resolve_from_wikipedia(self, url: str) -> Optional[dict]:
        """Resolve a Wikipedia URL to chemical structure data.

        Args:
            url: Wikipedia article URL.

        Returns:
            Dict with keys: smiles, inchikey, cas, inchi, label. None if not chemical.
        """
        title = self._extract_title(url)
        if not title:
            return None

        # SPARQL: find Wikidata item linked to this Wikipedia article,
        # then get its chemical identifiers
        query = f"""
        SELECT ?smiles ?inchikey ?cas ?inchi ?label WHERE {{
          ?article schema:about ?item ;
                   schema:isPartOf <https://en.wikipedia.org/> ;
                   schema:name "{title}"@en .
          OPTIONAL {{ ?item wdt:P233 ?smiles . }}
          OPTIONAL {{ ?item wdt:P235 ?inchikey . }}
          OPTIONAL {{ ?item wdt:P231 ?cas . }}
          OPTIONAL {{ ?item wdt:P234 ?inchi . }}
          OPTIONAL {{ ?item rdfs:label ?label . FILTER(LANG(?label) = "en") }}
        }} LIMIT 1
        """

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(
                    WIKIDATA_SPARQL_URL,
                    params={"query": query, "format": "json"},
                    headers={"User-Agent": "ChemAudit/2.0 (chemical validation tool)"},
                )
                response.raise_for_status()
                data = response.json()

            bindings = data.get("results", {}).get("bindings", [])
            if not bindings:
                return None

            entry = bindings[0]
            # Must have at least SMILES or InChIKey
            smiles = entry.get("smiles", {}).get("value")
            inchikey = entry.get("inchikey", {}).get("value")
            if not smiles and not inchikey:
                return None

            return {
                "smiles": smiles,
                "inchikey": inchikey,
                "cas": entry.get("cas", {}).get("value"),
                "inchi": entry.get("inchi", {}).get("value"),
                "label": entry.get("label", {}).get("value"),
            }

        except (httpx.HTTPError, KeyError, ValueError):
            return None
