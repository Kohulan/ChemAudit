"""
Validation API Routes

Endpoints for single molecule validation with Redis caching support.
Supports both synchronous and async (Celery priority queue) validation.
"""

import time
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Query, Request
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi as rdkit_inchi
from redis.asyncio import Redis

from app.core.cache import (
    get_cached_validation,
    set_cached_validation,
    validation_cache_key,
)
from app.core.config import settings
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.validation import (
    CheckResultSchema,
    InputInterpretation,
    MoleculeInfo,
    SimilarityRequest,
    SimilarityResponse,
    ValidationRequest,
    ValidationResponse,
)
from app.services.batch.tasks import validate_single_molecule
from app.services.parser.molecule_parser import MoleculeFormat, parse_molecule
from app.services.validation.engine import validation_engine

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
    request: Request,
    body: ValidationRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Validate a single molecule.

    Results are cached by InChIKey for 1 hour (configurable).
    Repeated validation of the same molecule returns cached result.

    Rate limits:
    - Anonymous: 10 requests/minute
    - API key: 300 requests/minute

    Args:
        request: FastAPI request (required for rate limiting)
        body: Validation request with molecule and options
        api_key: Optional API key for higher rate limits

    Returns:
        ValidationResponse with validation results and score

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()

    # IUPAC name detection and conversion
    input_interpretation = None
    molecule_input = body.molecule
    effective_input_type = body.input_type or "auto"

    if effective_input_type in ("auto", "iupac"):
        from app.services.iupac.converter import (
            detect_input_type,
            is_opsin_available,
            iupac_to_smiles,
        )

        if effective_input_type == "iupac":
            detected = "iupac"
        else:
            detected = detect_input_type(molecule_input)

        if detected == "iupac" or (detected == "ambiguous" and effective_input_type == "auto"):
            # Attempt IUPAC conversion if OPSIN available
            if is_opsin_available():
                converted = iupac_to_smiles(molecule_input)
                if converted:
                    input_interpretation = InputInterpretation(
                        detected_as="iupac",
                        original_input=molecule_input,
                        converted_smiles=converted,
                        conversion_source="opsin",
                    )
                    molecule_input = converted
                elif detected == "iupac" and effective_input_type == "iupac":
                    raise HTTPException(
                        status_code=400,
                        detail={
                            "error": "IUPAC name conversion failed",
                            "message": f"Could not convert '{molecule_input}' to SMILES via OPSIN",
                        },
                    )
                # If ambiguous and OPSIN fails, fall through to SMILES parsing
            elif detected == "iupac" and effective_input_type == "iupac":
                raise HTTPException(
                    status_code=503,
                    detail={
                        "error": "OPSIN unavailable",
                        "message": "IUPAC name input requires OPSIN which is not initialized",
                    },
                )

    # Parse molecule
    format_map = {
        "smiles": MoleculeFormat.SMILES,
        "inchi": MoleculeFormat.INCHI,
        "mol": MoleculeFormat.MOL,
        "iupac": MoleculeFormat.SMILES,  # Already converted above
        "auto": None,
    }
    input_format = format_map.get(body.format)

    parse_result = parse_molecule(molecule_input, input_format)

    if not parse_result.success or parse_result.mol is None:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "Failed to parse molecule",
                "errors": parse_result.errors,
                "warnings": parse_result.warnings,
                "format_detected": parse_result.format_detected.value,
            },
        )

    mol = parse_result.mol

    # Extract molecule info (needed for InChIKey)
    mol_info = extract_molecule_info(mol, body.molecule, body.preserve_aromatic)

    # Check cache if enabled
    redis = await get_redis(request)
    cache_key = None

    if settings.VALIDATION_CACHE_ENABLED and redis and mol_info.inchikey:
        cache_key = validation_cache_key(mol_info.inchikey, body.checks)
        cached = await get_cached_validation(redis, cache_key)

        if cached:
            # Return cached result with updated execution time
            execution_time = int((time.time() - start_time) * 1000)
            cached["execution_time_ms"] = execution_time
            cached["cached"] = True
            return ValidationResponse(**cached)

    # Run validation (cache miss or caching disabled)
    results, score = validation_engine.validate(mol, body.checks)

    # Convert to schema
    check_results = [
        CheckResultSchema(
            check_name=r.check_name,
            passed=r.passed,
            severity=r.severity,
            message=r.message,
            affected_atoms=r.affected_atoms,
            details=r.details,
        )
        for r in results
    ]

    execution_time = int((time.time() - start_time) * 1000)

    response = ValidationResponse(
        molecule_info=mol_info,
        overall_score=score,
        issues=[c for c in check_results if not c.passed],
        all_checks=check_results,
        execution_time_ms=execution_time,
        input_interpretation=input_interpretation,
    )

    # Cache the result if enabled
    if settings.VALIDATION_CACHE_ENABLED and redis and cache_key:
        # Convert to dict for caching (exclude execution_time as it varies)
        cache_data = response.model_dump()
        cache_data.pop("execution_time_ms", None)
        cache_data.pop("cached", None)
        await set_cached_validation(redis, cache_key, cache_data)

    return response


@router.post("/validate/async")
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def validate_molecule_async(
    request: Request,
    body: ValidationRequest,
    timeout: int = Query(default=30, ge=1, le=60, description="Timeout in seconds"),
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Validate a molecule using the Celery priority queue.

    This endpoint submits validation to the high_priority Celery queue,
    ensuring it won't be blocked by large batch processing jobs.
    Use this when batch jobs are running and you need responsive single validation.

    Args:
        request: FastAPI request
        body: Validation request with molecule and options
        timeout: Maximum time to wait for result (1-60 seconds)
        api_key: Optional API key for higher rate limits

    Returns:
        Validation results from the priority queue worker

    Raises:
        HTTPException: If validation times out or fails
    """
    start_time = time.time()

    # Parse molecule to get canonical SMILES
    format_map = {
        "smiles": MoleculeFormat.SMILES,
        "inchi": MoleculeFormat.INCHI,
        "mol": MoleculeFormat.MOL,
        "auto": None,
    }
    input_format = format_map.get(body.format)

    parse_result = parse_molecule(body.molecule, input_format)

    if not parse_result.success or parse_result.mol is None:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "Failed to parse molecule",
                "errors": parse_result.errors,
            },
        )

    # Get canonical SMILES for the task
    canonical_smiles = Chem.MolToSmiles(parse_result.mol)

    # Submit to priority queue
    task = validate_single_molecule.delay(canonical_smiles, body.checks)

    try:
        # Wait for result with timeout
        result = task.get(timeout=timeout)

        # Extract molecule info from the parsed molecule
        mol_info = extract_molecule_info(
            parse_result.mol, body.molecule, body.preserve_aromatic
        )

        execution_time = int((time.time() - start_time) * 1000)

        return {
            "molecule_info": (
                mol_info.model_dump()
                if hasattr(mol_info, "model_dump")
                else mol_info.__dict__
            ),
            "overall_score": result.get("validation", {}).get("overall_score", 0),
            "issues": result.get("validation", {}).get("issues", []),
            "alerts": result.get("alerts"),
            "scoring": result.get("scoring"),
            "execution_time_ms": execution_time,
            "queue": "high_priority",
        }

    except Exception as e:
        # Revoke the task if it's still running
        task.revoke(terminate=True)
        raise HTTPException(
            status_code=504,
            detail=f"Validation timed out or failed: {str(e)}",
        )


@router.get("/checks")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def list_checks(request: Request, api_key: Optional[str] = Depends(get_api_key)):
    """
    List available validation checks.

    Rate limits:
    - Anonymous: 10 requests/minute
    - API key: 300 requests/minute

    Returns:
        Dictionary mapping check categories to check names
    """
    return validation_engine.list_checks()


@router.post("/validate/similarity", response_model=SimilarityResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def compute_similarity(
    request: Request,
    body: SimilarityRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Compute ECFP4 Tanimoto similarity between two molecules.

    Uses Morgan fingerprints (radius=2, 2048 bits) and Tanimoto coefficient.
    """
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    mol_a = Chem.MolFromSmiles(body.smiles_a)
    if mol_a is None:
        raise HTTPException(status_code=400, detail=f"Cannot parse smiles_a: {body.smiles_a}")

    mol_b = Chem.MolFromSmiles(body.smiles_b)
    if mol_b is None:
        raise HTTPException(status_code=400, detail=f"Cannot parse smiles_b: {body.smiles_b}")

    fp_a = morgan_gen.GetFingerprint(mol_a)
    fp_b = morgan_gen.GetFingerprint(mol_b)

    similarity = DataStructs.TanimotoSimilarity(fp_a, fp_b)
    bits_a = fp_a.GetNumOnBits()
    bits_b = fp_b.GetNumOnBits()

    # Common bits = |A ∩ B| derived from Tanimoto: T = |A∩B| / |A∪B|
    # |A∪B| = |A| + |B| - |A∩B|, so |A∩B| = T * (|A| + |B|) / (1 + T)
    if similarity > 0:
        common_bits = int(round(similarity * (bits_a + bits_b) / (1 + similarity)))
    else:
        common_bits = 0

    return SimilarityResponse(
        tanimoto_similarity=round(similarity, 6),
        common_bits=common_bits,
        bits_a=bits_a,
        bits_b=bits_b,
    )


def extract_molecule_info(
    mol: Chem.Mol, input_smiles: str, preserve_aromatic: bool = False
) -> MoleculeInfo:
    """
    Extract molecule properties.

    Args:
        mol: RDKit molecule object
        input_smiles: Original input string
        preserve_aromatic: If True, output aromatic SMILES notation (lowercase atoms)

    Returns:
        MoleculeInfo with molecular properties
    """
    try:
        # When preserve_aromatic is True, we don't kekulize (default RDKit behavior)
        # When False (default), we kekulize to get explicit double bonds
        if preserve_aromatic:
            # Use default RDKit behavior which preserves aromatic notation
            canonical = Chem.MolToSmiles(mol)
        else:
            # Kekulize to get explicit double bonds (C1=CC=CC=C1 instead of c1ccccc1)
            # Note: RDKit's MolToSmiles uses aromatic notation by default
            # We need to use the kekuleSmiles parameter to get kekulized form
            try:
                Chem.Kekulize(mol, clearAromaticFlags=False)
                canonical = Chem.MolToSmiles(mol, kekuleSmiles=True)
            except Exception:
                # If kekulization fails, fall back to aromatic form
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
        num_atoms=num_atoms,
    )
