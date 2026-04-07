"""
Structure Quality Diagnostics API Routes

Five endpoints for SMILES error diagnostics, InChI layer diff,
format round-trip lossiness, cross-pipeline comparison, and file pre-validation.
"""

from typing import Optional

from fastapi import APIRouter, Depends, File, HTTPException, Request, UploadFile

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.diagnostics import (
    CrossPipelineRequest,
    CrossPipelineResponse,
    FilePreValidationResponse,
    InChIDiffRequest,
    InChIDiffResponse,
    RoundTripRequest,
    RoundTripResponse,
    SMILESDiagnosticsRequest,
    SMILESDiagnosticsResponse,
)
from app.services.diagnostics import (
    check_roundtrip,
    compare_pipelines,
    diagnose_smiles,
    diff_inchi_layers,
    prevalidate_csv,
    prevalidate_sdf,
)

router = APIRouter()

# Maximum allowed upload file size (50 MB default)
_MAX_FILE_SIZE_BYTES = 50 * 1024 * 1024


@router.post("/diagnostics/smiles", response_model=SMILESDiagnosticsResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def smiles_diagnostics(
    request: Request,
    body: SMILESDiagnosticsRequest,
    api_key: Optional[str] = Depends(get_api_key),
) -> SMILESDiagnosticsResponse:
    """
    Diagnose a SMILES string for position-specific errors and fix suggestions.

    Uses dual strategy: RDKit DetectChemistryProblems for parseable-but-problematic
    SMILES, and log capture for unparseable SMILES. Returns error types, character
    positions, and ranked fix suggestions.

    Args:
        body: Request containing the raw SMILES string.

    Returns:
        SMILESDiagnosticsResponse with valid flag, canonical SMILES, warnings, and errors.

    Raises:
        HTTPException: 500 if an unexpected server error occurs.
    """
    try:
        result = diagnose_smiles(body.smiles)
        return SMILESDiagnosticsResponse(**result)
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail={"error": "Diagnostics failed", "detail": str(exc)},
        ) from exc


@router.post("/diagnostics/inchi-diff", response_model=InChIDiffResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def inchi_diff(
    request: Request,
    body: InChIDiffRequest,
    api_key: Optional[str] = Depends(get_api_key),
) -> InChIDiffResponse:
    """
    Compare two InChI strings layer by layer.

    Parses each InChI into its constituent layers (formula, connections, hydrogens,
    stereo, isotope, etc.) and produces a per-layer diff table showing matches and
    differences. Pure string comparison — no RDKit required.

    Args:
        body: Request with two InChI strings (inchi_a and inchi_b).

    Returns:
        InChIDiffResponse with identical flag, layer_rows diff table, and parsed layers.

    Raises:
        HTTPException: 500 if an unexpected server error occurs.
    """
    try:
        result = diff_inchi_layers(body.inchi_a, body.inchi_b)
        return InChIDiffResponse(**result)
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail={"error": "Diagnostics failed", "detail": str(exc)},
        ) from exc


@router.post("/diagnostics/roundtrip", response_model=RoundTripResponse)
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def roundtrip_check(
    request: Request,
    body: RoundTripRequest,
    api_key: Optional[str] = Depends(get_api_key),
) -> RoundTripResponse:
    """
    Check whether a SMILES string round-trips losslessly through an intermediate format.

    Supported routes:
    - smiles_inchi_smiles: SMILES -> InChI -> SMILES (detects stereo/isotope loss)
    - smiles_mol_smiles: SMILES -> MOL block -> SMILES (detects stereo/charge loss)

    Args:
        body: Request with SMILES string and conversion route.

    Returns:
        RoundTripResponse with lossy flag, detected losses, and intermediate representation.

    Raises:
        HTTPException: 500 if an unexpected server error occurs.
    """
    try:
        result = check_roundtrip(body.smiles, body.route)
        return RoundTripResponse(**result)
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail={"error": "Diagnostics failed", "detail": str(exc)},
        ) from exc


@router.post("/diagnostics/cross-pipeline", response_model=CrossPipelineResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def cross_pipeline(
    request: Request,
    body: CrossPipelineRequest,
    api_key: Optional[str] = Depends(get_api_key),
) -> CrossPipelineResponse:
    """
    Compare standardization output across three pipelines.

    Runs the input molecule through:
    1. RDKit MolStandardize (Cleanup + LargestFragment + Uncharger + TautomerEnumerator)
    2. ChEMBL-style pipeline (existing StandardizationPipeline)
    3. Minimal pipeline (parse + sanitize only)

    Compares SMILES, InChIKey, MW, formula, charge, and stereo count across all three,
    reporting disagreements and structural conflicts.

    Args:
        body: Request with molecule string (SMILES/InChI/MOL) and format hint.

    Returns:
        CrossPipelineResponse with per-pipeline results and property comparison table.

    Raises:
        HTTPException: 400 if molecule cannot be parsed; 500 on unexpected error.
    """
    try:
        result = compare_pipelines(body.molecule)
        return CrossPipelineResponse(**result)
    except ValueError as exc:
        raise HTTPException(
            status_code=400,
            detail={"error": "Invalid molecule input", "detail": str(exc)},
        ) from exc
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail={"error": "Diagnostics failed", "detail": str(exc)},
        ) from exc


@router.post("/diagnostics/file-prevalidate", response_model=FilePreValidationResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def file_prevalidate(
    request: Request,
    file: UploadFile = File(...),
    api_key: Optional[str] = Depends(get_api_key),
) -> FilePreValidationResponse:
    """
    Pre-validate an SDF or CSV file for structural integrity issues.

    For SDF files: checks M END terminators, counts line format, and encoding.
    For CSV files: checks encoding (UTF-8/Latin-1 fallback), SMILES column detection,
    and empty row detection. Both formats scan for security-suspicious patterns.

    Args:
        file: Uploaded file (.sdf, .sd, or .csv extension).

    Returns:
        FilePreValidationResponse with issue list, counts, and validity verdict.

    Raises:
        HTTPException: 400 if file type is unsupported or file is too large; 500 on unexpected error.
    """
    try:
        content = await file.read()

        # Enforce file size limit
        if len(content) > _MAX_FILE_SIZE_BYTES:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": "File too large",
                    "detail": f"Maximum allowed size is {_MAX_FILE_SIZE_BYTES // (1024 * 1024)} MB",
                },
            )

        # Determine file type from filename extension
        filename = file.filename or ""
        filename_lower = filename.lower()
        if filename_lower.endswith(".sdf") or filename_lower.endswith(".sd"):
            result = prevalidate_sdf(content)
        elif filename_lower.endswith(".csv"):
            result = prevalidate_csv(content)
        else:
            raise HTTPException(
                status_code=400,
                detail={"error": "Unsupported file type. Use .sdf or .csv"},
            )

        return FilePreValidationResponse(**result)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail={"error": "Diagnostics failed", "detail": str(exc)},
        ) from exc
