"""
External integrations service.

Provides clients for COCONUT, PubChem, ChEMBL, UniChem, Wikidata, ChEBI,
identifier detection, universal resolver, and cross-database comparator.
"""

from app.services.integrations.chebi import ChEBIClient
from app.services.integrations.chembl import ChEMBLClient, get_bioactivity
from app.services.integrations.coconut import COCONUTClient, lookup_natural_product
from app.services.integrations.comparator import compare_across_databases
from app.services.integrations.identifier_detect import detect_identifier_type
from app.services.integrations.pubchem import PubChemClient, get_compound_info
from app.services.integrations.resolver import resolve_identifier
from app.services.integrations.unichem import UniChemClient
from app.services.integrations.wikidata import WikidataClient

__all__ = [
    "COCONUTClient",
    "lookup_natural_product",
    "PubChemClient",
    "get_compound_info",
    "ChEMBLClient",
    "get_bioactivity",
    "UniChemClient",
    "WikidataClient",
    "ChEBIClient",
    "detect_identifier_type",
    "resolve_identifier",
    "compare_across_databases",
]
