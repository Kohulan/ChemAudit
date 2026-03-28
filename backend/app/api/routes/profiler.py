"""Compound Profiler API routes.

Provides compound profiling endpoints for PFI, #stars, Abbott bioavailability,
consensus LogP, skin permeation, 3D shape, custom MPO, ligand efficiency,
and SA comparison. Per D-25: combined /full endpoint + granular per-metric.
Per D-26: /full does NOT call compute_3d_shape (lazy — use /shape-3d separately).
"""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Request
from rdkit import Chem

from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.profiler import (
    LERequest,
    MPORequest,
    ProfileRequest,
    SAComparisonRequest,
    Shape3DRequest,
)
from app.services.profiler.compound_profile import compute_3d_shape, compute_full_profile
from app.services.profiler.ligand_efficiency import compute_ligand_efficiency
from app.services.profiler.mpo_scoring import compute_cns_mpo, compute_custom_mpo
from app.services.profiler.sa_comparison import compute_sa_comparison

router = APIRouter(prefix="/profiler", tags=["profiler"])


@router.post("/full")
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def full_profile(
    request: Request,
    body: ProfileRequest,
    api_key: str | None = Depends(get_api_key),
):
    """
    Compute a full compound profile for a molecule.

    Returns PFI, #stars, Abbott bioavailability score, consensus LogP,
    skin permeation, SA comparison (synthesizability), and CNS MPO.

    Per D-26: 3D shape descriptors are NOT included. Use /shape-3d for those.

    Args:
        body: ProfileRequest with SMILES string.

    Returns:
        Composite profile dict with all core metrics.

    Raises:
        HTTPException: 400 if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(body.smiles)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail={"error": "Invalid SMILES", "smiles": body.smiles},
        )

    profile = compute_full_profile(mol)
    sa_comparison = compute_sa_comparison(mol, body.smiles)
    cns_mpo = compute_cns_mpo(mol)

    return {
        **profile,
        "sa_comparison": sa_comparison,
        "cns_mpo": cns_mpo,
    }


@router.post("/shape-3d")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def shape_3d(
    request: Request,
    body: Shape3DRequest,
    api_key: str | None = Depends(get_api_key),
):
    """
    Compute 3D shape descriptors for a molecule.

    Uses ETKDGv3 conformer generation with MMFF94 optimization. Returns
    PMI ratios (npr1, npr2), Plane of Best Fit (pbf), and shape classification.
    This is a separate endpoint from /full (lazy computation per D-26).

    Args:
        body: Shape3DRequest with SMILES string.

    Returns:
        3D shape descriptor dict, or failure flag if conformer generation fails.

    Raises:
        HTTPException: 400 if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(body.smiles)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail={"error": "Invalid SMILES", "smiles": body.smiles},
        )

    result = compute_3d_shape(mol)
    if result is None:
        return {
            "3d_conformer_failed": True,
            "message": "Conformer generation failed after 3 attempts",
        }

    return {
        "3d_conformer_failed": False,
        **result,
    }


@router.post("/sa-comparison")
@limiter.limit("20/minute", key_func=get_rate_limit_key)
async def sa_comparison(
    request: Request,
    body: SAComparisonRequest,
    api_key: str | None = Depends(get_api_key),
):
    """
    Compare SA Score, SCScore, and SYBA synthesizability scores side-by-side.

    Returns all three scores with graceful fallbacks:
    - SA Score: always available (RDKit Contrib bundled)
    - SCScore: available if scscore is pip-installed
    - SYBA: available if syba is installed (runs via subprocess for GPL-3.0 isolation)

    Args:
        body: SAComparisonRequest with SMILES string.

    Returns:
        Comparison dict with sa_score, scscore, syba, and rascore (placeholder).

    Raises:
        HTTPException: 400 if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(body.smiles)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail={"error": "Invalid SMILES", "smiles": body.smiles},
        )

    return compute_sa_comparison(mol, body.smiles)


@router.post("/efficiency")
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def ligand_efficiency(
    request: Request,
    body: LERequest,
    api_key: str | None = Depends(get_api_key),
):
    """
    Compute extended ligand efficiency metrics from measured activity data.

    Returns LE, LLE, LELP, BEI, and SEI metrics.

    Activity types supported: IC50_nM, IC50_uM, Ki_nM, pIC50, pKd.

    Args:
        body: LERequest with SMILES, activity value, and activity type.

    Returns:
        Ligand efficiency dict with pIC50, LE, LLE, LELP, BEI, SEI.

    Raises:
        HTTPException: 400 if SMILES is invalid or activity_type is unknown.
    """
    mol = Chem.MolFromSmiles(body.smiles)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail={"error": "Invalid SMILES", "smiles": body.smiles},
        )

    result = compute_ligand_efficiency(mol, body.activity_value, body.activity_type)
    if "error" in result:
        raise HTTPException(
            status_code=400,
            detail={"error": result["error"]},
        )

    return result


@router.post("/mpo")
@limiter.limit("30/minute", key_func=get_rate_limit_key)
async def custom_mpo(
    request: Request,
    body: MPORequest,
    api_key: str | None = Depends(get_api_key),
):
    """
    Compute a custom MPO (Multi-Parameter Optimization) score.

    Uses user-defined desirability functions with sigmoid, ramp, or step shapes.
    Each property entry defines a range and weight for scoring.

    Args:
        body: MPORequest with SMILES and profile list.

    Returns:
        Custom MPO result with score, max_score, normalized score, and components.

    Raises:
        HTTPException: 400 if SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(body.smiles)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail={"error": "Invalid SMILES", "smiles": body.smiles},
        )

    profile_dicts = [p.model_dump() for p in body.profile]
    return compute_custom_mpo(mol, profile_dicts)
