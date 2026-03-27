"""
Safety Assessment API Routes

Endpoints for CYP metabolism soft-spot prediction, hERG liability assessment,
beyond-Rule-of-5 (bRo5), REOS filter, and complexity percentile evaluation.
"""

import time
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Request

from app.api.routes.validation import extract_molecule_info
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.safety import (
    Bro5ResultSchema,
    Bro5ViolationSchema,
    ComplexityPropertySchema,
    ComplexityResultSchema,
    CypResultSchema,
    CypSiteSchema,
    HergResultSchema,
    ReosResultSchema,
    ReosViolationSchema,
    SafetyAssessRequest,
    SafetyAssessResponse,
    SafetySummaryResponse,
)
from app.services.alerts.complexity_filter import compute_complexity_percentile
from app.services.alerts.unified_screen import unified_screen
from app.services.parser.molecule_parser import MoleculeFormat, parse_molecule
from app.services.safety.bro5 import compute_bro5
from app.services.safety.cyp_softspots import screen_cyp_softspots
from app.services.safety.herg_risk import compute_herg_risk
from app.services.safety.reos import compute_reos

router = APIRouter()

# ---------------------------------------------------------------------------
# Molecule parse helper (shared by both endpoints)
# ---------------------------------------------------------------------------

_FORMAT_MAP = {
    "smiles": MoleculeFormat.SMILES,
    "inchi": MoleculeFormat.INCHI,
    "mol": MoleculeFormat.MOL,
    "auto": None,
}


def _parse_or_raise(molecule: str, fmt: str):
    """Parse molecule string or raise HTTP 400 with error details."""
    input_format = _FORMAT_MAP.get(fmt)
    parse_result = parse_molecule(molecule, input_format)
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
    return parse_result.mol


# ---------------------------------------------------------------------------
# Schema conversion helpers
# ---------------------------------------------------------------------------


def _build_cyp_result(raw: list[dict]) -> CypResultSchema:
    """Convert raw CYP service output to CypResultSchema."""
    sites = [
        CypSiteSchema(
            site_name=s["site_name"],
            reaction_type=s["reaction_type"],
            matched_atoms=s["matched_atoms"],
        )
        for s in raw
    ]
    return CypResultSchema(sites=sites, n_sites=len(sites))


def _build_herg_result(raw: dict) -> HergResultSchema:
    """Convert raw hERG service output to HergResultSchema."""
    return HergResultSchema(
        herg_risk=raw["herg_risk"],
        risk_score=raw["risk_score"],
        max_score=raw.get("max_score", 4),
        flags=raw.get("flags", []),
        descriptors=raw.get("descriptors", {}),
    )


def _build_bro5_result(raw: dict) -> Bro5ResultSchema:
    """Convert raw bRo5 service output to Bro5ResultSchema."""
    violations = [
        Bro5ViolationSchema(
            property=v["property"],
            value=float(v["value"]),
            threshold=float(v["threshold"]),
            direction=v["direction"],
        )
        for v in raw.get("violations", [])
    ]
    return Bro5ResultSchema(
        applicable=raw["applicable"],
        passed=raw["passed"],
        message=raw.get("message"),
        violations=violations,
        values={k: float(v) for k, v in raw.get("values", {}).items()},
    )


def _build_reos_result(raw: dict) -> ReosResultSchema:
    """Convert raw REOS service output to ReosResultSchema."""
    violations = [
        ReosViolationSchema(
            property=v["property"],
            value=float(v["value"]),
            range=[float(r) for r in v["range"]],
            exceeded=v["exceeded"],
        )
        for v in raw.get("violations", [])
    ]
    return ReosResultSchema(
        passed=raw["passed"],
        violations=violations,
        n_violations=raw.get("n_violations", len(violations)),
        descriptors={k: float(v) for k, v in raw.get("descriptors", {}).items()},
    )


def _build_complexity_result(raw: dict) -> ComplexityResultSchema:
    """Convert raw complexity service output to ComplexityResultSchema."""
    properties = {
        prop: ComplexityPropertySchema(
            value=float(pdata["value"]),
            p5=float(pdata["p5"]),
            p95=float(pdata["p95"]),
            outlier=pdata["outlier"],
            direction=pdata.get("direction"),
        )
        for prop, pdata in raw.get("properties", {}).items()
    }
    return ComplexityResultSchema(
        properties=properties,
        n_outliers=raw.get("n_outliers", 0),
        outlier_properties=raw.get("outlier_properties", []),
        within_range=raw.get("within_range", True),
    )


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------


@router.post("/safety/assess", response_model=SafetyAssessResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def safety_assess(
    request: Request,
    body: SafetyAssessRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Full safety assessment for a molecule.

    Evaluates CYP metabolism soft-spots, hERG liability (amphiphile rule),
    beyond-Rule-of-5 (bRo5) thresholds, REOS 7-property filter, and
    complexity percentile against commercial compound distributions.

    Args:
        body: Request with molecule string and format hint

    Returns:
        SafetyAssessResponse with results for all 5 assessment modules

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()
    mol = _parse_or_raise(body.molecule, body.format)
    mol_info = extract_molecule_info(mol, body.molecule)

    cyp_raw = screen_cyp_softspots(mol)
    herg_raw = compute_herg_risk(mol)
    bro5_raw = compute_bro5(mol)
    reos_raw = compute_reos(mol)
    complexity_raw = compute_complexity_percentile(mol)

    execution_time = int((time.time() - start_time) * 1000)

    return SafetyAssessResponse(
        status="completed",
        molecule_info=mol_info,
        cyp_softspots=_build_cyp_result(cyp_raw),
        herg=_build_herg_result(herg_raw),
        bro5=_build_bro5_result(bro5_raw),
        reos=_build_reos_result(reos_raw),
        complexity=_build_complexity_result(complexity_raw),
        execution_time_ms=execution_time,
    )


@router.post("/safety/summary", response_model=SafetySummaryResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def safety_summary(
    request: Request,
    body: SafetyAssessRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Lightweight safety summary for badge display.

    Returns traffic-light status strings for hERG, bRo5, and REOS, plus
    total alert counts for use in SingleValidation/Profiler badge display.
    CYP is always reported as "default" (informational, not pass/fail).

    Traffic-light mappings:
    - cyp_status: always "default"
    - herg_status: "success" (score<=1), "warning" (score==2), "error" (score>=3)
    - bro5_status: "default" (not applicable), "success" (passed), "error" (failed)
    - reos_status: "success" (0 violations), "warning" (1 violation), "error" (2+ violations)

    Args:
        body: Request with molecule string and format hint

    Returns:
        SafetySummaryResponse with status strings and alert counts

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    mol = _parse_or_raise(body.molecule, body.format)

    screen_result = unified_screen(mol)
    herg_raw = compute_herg_risk(mol)
    bro5_raw = compute_bro5(mol)
    reos_raw = compute_reos(mol)
    complexity_raw = compute_complexity_percentile(mol)

    # Map hERG risk score to traffic-light status
    herg_score = herg_raw["risk_score"]
    if herg_score <= 1:
        herg_status = "success"
    elif herg_score == 2:
        herg_status = "warning"
    else:
        herg_status = "error"

    # Map bRo5 result to traffic-light status
    if not bro5_raw["applicable"]:
        bro5_status = "default"
    elif bro5_raw["passed"]:
        bro5_status = "success"
    else:
        bro5_status = "error"

    # Map REOS violations to traffic-light status
    n_reos = reos_raw["n_violations"]
    if n_reos == 0:
        reos_status = "success"
    elif n_reos == 1:
        reos_status = "warning"
    else:
        reos_status = "error"

    return SafetySummaryResponse(
        status="completed",
        total_alerts=screen_result["total_raw"],
        has_critical=screen_result["has_critical"],
        cyp_status="default",
        herg_status=herg_status,
        bro5_status=bro5_status,
        reos_status=reos_status,
        complexity_outliers=complexity_raw["n_outliers"],
    )
