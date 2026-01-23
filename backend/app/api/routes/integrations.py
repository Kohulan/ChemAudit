"""
External Integrations API Routes

Endpoints for DECIMER, COCONUT, PubChem, and ChEMBL integrations.
"""
from fastapi import APIRouter, Request, Depends
from typing import Optional

from app.schemas.integrations import (
    DECIMERRequest,
    DECIMERValidation,
    COCONUTRequest,
    COCONUTResult,
    PubChemRequest,
    PubChemResult,
    ChEMBLRequest,
    ChEMBLResult,
)
from app.services.integrations import (
    validate_ocsr_result,
    lookup_natural_product,
    get_compound_info,
    get_bioactivity,
)
from app.core.rate_limit import limiter, get_rate_limit_key
from app.core.security import get_api_key


router = APIRouter()


@router.post("/integrations/decimer/validate", response_model=DECIMERValidation)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def validate_decimer_ocsr(
    req: Request,
    request: DECIMERRequest,
    api_key: Optional[str] = Depends(get_api_key)
):
    """
    Validate DECIMER OCSR output.

    DECIMER (Deep Learning for Chemical Image Recognition) converts
    chemical structure images to SMILES. This endpoint validates the
    OCSR result and provides confidence-adjusted validation.

    Rate limits:
    - Anonymous: 10 requests/minute (shared with other endpoints)
    - API key: 300 requests/minute (shared with other endpoints)
    - This endpoint: 30 requests/minute (specific limit)

    Args:
        req: FastAPI request (required for rate limiting)
        request: DECIMER validation request with SMILES and optional confidence
        api_key: Optional API key for higher rate limits

    Returns:
        Validation result with canonical SMILES and identifiers

    Example:
        ```json
        {
            "smiles": "CCO",
            "confidence": 0.95,
            "is_valid": true,
            "validation_message": "High confidence OCSR result - valid structure",
            "canonical_smiles": "CCO",
            "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        }
        ```
    """
    return validate_ocsr_result(request)


@router.post("/integrations/coconut/lookup", response_model=COCONUTResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_coconut(
    req: Request,
    request: COCONUTRequest,
    api_key: Optional[str] = Depends(get_api_key)
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
        req: FastAPI request (required for rate limiting)
        request: COCONUT lookup request with SMILES or InChIKey
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
            "url": "https://coconut.naturalproducts.net/compound/CNP0123456"
        }
        ```
    """
    return await lookup_natural_product(request)


@router.post("/integrations/pubchem/lookup", response_model=PubChemResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_pubchem(
    req: Request,
    request: PubChemRequest,
    api_key: Optional[str] = Depends(get_api_key)
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
        req: FastAPI request (required for rate limiting)
        request: PubChem lookup request with SMILES or InChIKey
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
    return await get_compound_info(request)


@router.post("/integrations/chembl/bioactivity", response_model=ChEMBLResult)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def lookup_chembl_bioactivity(
    req: Request,
    request: ChEMBLRequest,
    api_key: Optional[str] = Depends(get_api_key)
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
        req: FastAPI request (required for rate limiting)
        request: ChEMBL lookup request with SMILES or InChIKey
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
    return await get_bioactivity(request)
