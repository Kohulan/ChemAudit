"""
Validation API Routes

Endpoints for single molecule validation with Redis caching support.
"""
from fastapi import APIRouter, HTTPException, Request, Depends
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi as rdkit_inchi, rdMolDescriptors
from redis.asyncio import Redis
import time
from typing import Optional

from app.schemas.validation import ValidationRequest, ValidationResponse, MoleculeInfo, CheckResultSchema
from app.services.parser.molecule_parser import parse_molecule, MoleculeFormat
from app.services.validation.engine import validation_engine
from app.core.rate_limit import limiter, get_rate_limit_key
from app.core.security import get_api_key
from app.core.config import settings
from app.core.cache import (
    validation_cache_key,
    get_cached_validation,
    set_cached_validation,
)


router = APIRouter()


async def get_redis(request: Request) -> Optional[Redis]:
    """
    Get Redis client from app state.

    Returns None if Redis is not available (caching will be skipped).
    """
    if hasattr(request.app.state, "redis") and request.app.state.redis:
        return request.app.state.redis
    return None


@router.post("/validate", response_model=ValidationResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def validate_molecule(
    req: Request,
    request: ValidationRequest,
    api_key: Optional[str] = Depends(get_api_key)
):
    """
    Validate a single molecule.

    Results are cached by InChIKey for 1 hour (configurable).
    Repeated validation of the same molecule returns cached result.

    Rate limits:
    - Anonymous: 10 requests/minute
    - API key: 300 requests/minute

    Args:
        req: FastAPI request (required for rate limiting)
        request: Validation request with molecule and options
        api_key: Optional API key for higher rate limits

    Returns:
        ValidationResponse with validation results and score

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()

    # Parse molecule
    format_map = {
        "smiles": MoleculeFormat.SMILES,
        "inchi": MoleculeFormat.INCHI,
        "mol": MoleculeFormat.MOL,
        "auto": None
    }
    input_format = format_map.get(request.format)

    parse_result = parse_molecule(request.molecule, input_format)

    if not parse_result.success or parse_result.mol is None:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "Failed to parse molecule",
                "errors": parse_result.errors,
                "warnings": parse_result.warnings,
                "format_detected": parse_result.format_detected.value
            }
        )

    mol = parse_result.mol

    # Extract molecule info (needed for InChIKey)
    mol_info = extract_molecule_info(mol, request.molecule)

    # Check cache if enabled
    redis = await get_redis(req)
    cache_key = None
    cached_from_redis = False

    if settings.VALIDATION_CACHE_ENABLED and redis and mol_info.inchikey:
        cache_key = validation_cache_key(mol_info.inchikey, request.checks)
        cached = await get_cached_validation(redis, cache_key)

        if cached:
            # Return cached result with updated execution time
            execution_time = int((time.time() - start_time) * 1000)
            cached["execution_time_ms"] = execution_time
            cached["cached"] = True
            return ValidationResponse(**cached)

    # Run validation (cache miss or caching disabled)
    results, score = validation_engine.validate(mol, request.checks)

    # Convert to schema
    check_results = [
        CheckResultSchema(
            check_name=r.check_name,
            passed=r.passed,
            severity=r.severity,
            message=r.message,
            affected_atoms=r.affected_atoms,
            details=r.details
        )
        for r in results
    ]

    execution_time = int((time.time() - start_time) * 1000)

    response = ValidationResponse(
        molecule_info=mol_info,
        overall_score=score,
        issues=[c for c in check_results if not c.passed],
        all_checks=check_results,
        execution_time_ms=execution_time
    )

    # Cache the result if enabled
    if settings.VALIDATION_CACHE_ENABLED and redis and cache_key:
        # Convert to dict for caching (exclude execution_time as it varies)
        cache_data = response.model_dump()
        cache_data.pop("execution_time_ms", None)
        cache_data.pop("cached", None)
        await set_cached_validation(redis, cache_key, cache_data)

    return response


@router.get("/checks")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def list_checks(
    request: Request,
    api_key: Optional[str] = Depends(get_api_key)
):
    """
    List available validation checks.

    Rate limits:
    - Anonymous: 10 requests/minute
    - API key: 300 requests/minute

    Returns:
        Dictionary mapping check categories to check names
    """
    return validation_engine.list_checks()


def extract_molecule_info(mol: Chem.Mol, input_smiles: str) -> MoleculeInfo:
    """
    Extract molecule properties.

    Args:
        mol: RDKit molecule object
        input_smiles: Original input string

    Returns:
        MoleculeInfo with molecular properties
    """
    try:
        canonical = Chem.MolToSmiles(mol)
        mol_inchi = rdkit_inchi.MolToInchi(mol)
        mol_inchikey = rdkit_inchi.MolToInchiKey(mol) if mol_inchi else None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
    except Exception:
        canonical = None
        mol_inchi = None
        mol_inchikey = None
        formula = None
        mw = None
        num_atoms = None

    return MoleculeInfo(
        input_smiles=input_smiles,
        canonical_smiles=canonical,
        inchi=mol_inchi,
        inchikey=mol_inchikey,
        molecular_formula=formula,
        molecular_weight=mw,
        num_atoms=num_atoms
    )
