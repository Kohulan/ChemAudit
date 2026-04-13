"""
External Integrations API Routes

Endpoints for COCONUT, PubChem, ChEMBL, Wikidata, and SureChEMBL integrations,
plus cross-database comparison.
"""

from typing import Optional

from fastapi import APIRouter, Depends, Request
from pydantic import BaseModel, model_validator

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.integrations import (
    ChEMBLRequest,
    ChEMBLResult,
    COCONUTRequest,
    COCONUTResult,
    ConsistencyResult,
    PubChemRequest,
    PubChemResult,
    SureChEMBLRequest,
    SureChEMBLResult,
    WikidataRequest,
    WikidataResult,
)
from app.services.integrations import (
    get_bioactivity,
    get_compound_info,
    lookup_natural_product,
    lookup_wikidata,
)
from app.services.integrations.comparator import compare_across_databases
from app.services.integrations.surechembl import lookup_surechembl

router = APIRouter()


@router.post("/integrations/coconut/lookup", response_model=COCONUTResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_coconut(
    request: Request,
    body: COCONUTRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Look up molecule in COCONUT natural products database.

    COCONUT (COlleCtion of Open Natural ProdUcTs) contains >400,000
    natural product structures with organism and literature information.

    Rate limits:
    - Anonymous: 10 requests/minute (shared with other endpoints)
    - API key: 300 requests/minute (shared with other endpoints)
    - This endpoint: 30 requests/minute (specific limit)

    Args:
        request: FastAPI request (required for rate limiting)
        body: COCONUT lookup request with SMILES or InChIKey
        api_key: Optional API key for higher rate limits

    Returns:
        Natural product information if found in COCONUT

    Example:
        ```json
        {
            "found": true,
            "coconut_id": "CNP0123456",
            "name": "Caffeine",
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "molecular_weight": 194.19,
            "organism": "Coffea arabica",
            "url": "https://coconut.naturalproducts.net/compounds/CNP0123456"
        }
        ```
    """
    return await lookup_natural_product(body)


@router.post("/integrations/pubchem/lookup", response_model=PubChemResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_pubchem(
    request: Request,
    body: PubChemRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Cross-reference molecule with PubChem database.

    PubChem is the NIH's public repository of chemical structures,
    properties, and biological activities.

    Rate limits:
    - Anonymous: 10 requests/minute (shared with other endpoints)
    - API key: 300 requests/minute (shared with other endpoints)
    - This endpoint: 30 requests/minute (specific limit)

    Args:
        request: FastAPI request (required for rate limiting)
        body: PubChem lookup request with SMILES or InChIKey
        api_key: Optional API key for higher rate limits

    Returns:
        Compound information if found in PubChem

    Example:
        ```json
        {
            "found": true,
            "cid": 702,
            "iupac_name": "ethanol",
            "molecular_formula": "C2H6O",
            "molecular_weight": 46.07,
            "synonyms": ["ethanol", "ethyl alcohol", "alcohol"],
            "url": "https://pubchem.ncbi.nlm.nih.gov/compound/702"
        }
        ```
    """
    return await get_compound_info(body)


@router.post("/integrations/chembl/bioactivity", response_model=ChEMBLResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_chembl_bioactivity(
    request: Request, body: ChEMBLRequest, api_key: Optional[str] = Depends(get_api_key)
):
    """
    Look up ChEMBL bioactivity data for molecule.

    ChEMBL is a manually curated database of bioactive molecules with
    drug-like properties, including bioactivity data from scientific literature.

    Rate limits:
    - Anonymous: 10 requests/minute (shared with other endpoints)
    - API key: 300 requests/minute (shared with other endpoints)
    - This endpoint: 30 requests/minute (specific limit)

    Args:
        request: FastAPI request (required for rate limiting)
        body: ChEMBL lookup request with SMILES or InChIKey
        api_key: Optional API key for higher rate limits

    Returns:
        Molecule information and bioactivity data if found in ChEMBL

    Example:
        ```json
        {
            "found": true,
            "chembl_id": "CHEMBL25",
            "pref_name": "ASPIRIN",
            "max_phase": 4,
            "bioactivity_count": 1234,
            "bioactivities": [
                {
                    "target_chembl_id": "CHEMBL240",
                    "target_name": "Cyclooxygenase-1",
                    "activity_type": "IC50",
                    "activity_value": 100.0,
                    "activity_unit": "nM"
                }
            ],
            "url": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25"
        }
        ```
    """
    return await get_bioactivity(body)


@router.post("/integrations/wikidata/lookup", response_model=WikidataResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_wikidata_endpoint(
    request: Request,
    body: WikidataRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Look up molecule in Wikidata via SPARQL.

    Wikidata stores structured chemical data (SMILES, InChIKey, CAS,
    molecular formula) linked to Wikipedia articles and other knowledge bases.

    Rate limits:
    - This endpoint: 30 requests/minute (specific limit)

    Args:
        request: FastAPI request (required for rate limiting)
        body: Wikidata lookup request with SMILES or InChIKey
        api_key: Optional API key for higher rate limits

    Returns:
        Compound information if found in Wikidata
    """
    return await lookup_wikidata(body)


class CompareRequest(BaseModel):
    """Request body for cross-database comparison."""

    smiles: Optional[str] = None
    inchikey: Optional[str] = None

    @model_validator(mode="after")
    def require_at_least_one(self):
        if not self.smiles and not self.inchikey:
            raise ValueError("At least one of 'smiles' or 'inchikey' is required")
        return self


@router.post("/integrations/compare", response_model=ConsistencyResult)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def compare_databases(
    request: Request,
    body: CompareRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """Compare molecule representation across PubChem, ChEMBL, and COCONUT."""
    return await compare_across_databases(smiles=body.smiles, inchikey=body.inchikey)


@router.post("/integrations/surechembl/lookup", response_model=SureChEMBLResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_surechembl_patents(
    request: Request,
    body: SureChEMBLRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Look up molecule in SureChEMBL patent database.

    Checks patent literature presence via UniChem cross-reference to SureChEMBL,
    with optional enrichment from SureChEMBL direct API for patent count data.

    Returns binary present/absent result with a direct link to the SureChEMBL
    entry when found.
    """
    return await lookup_surechembl(body)
