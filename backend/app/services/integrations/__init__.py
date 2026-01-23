"""
External integrations service.

Provides clients for DECIMER, COCONUT, PubChem, and ChEMBL.
"""
from app.services.integrations.decimer import DECIMERClient, validate_ocsr_result
from app.services.integrations.coconut import COCONUTClient, lookup_natural_product
from app.services.integrations.pubchem import PubChemClient, get_compound_info
from app.services.integrations.chembl import ChEMBLClient, get_bioactivity

__all__ = [
    "DECIMERClient",
    "validate_ocsr_result",
    "COCONUTClient",
    "lookup_natural_product",
    "PubChemClient",
    "get_compound_info",
    "ChEMBLClient",
    "get_bioactivity",
]
