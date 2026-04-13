"""
Structural Alert Screening API Routes

Endpoints for screening molecules against PAINS, BRENK, and other
structural alert pattern catalogs.

IMPORTANT: Alerts are warnings for investigation, not automatic rejections.
87 FDA-approved drugs contain PAINS patterns.
"""

import time
from typing import Optional

from fastapi import APIRouter, Depends, HTTPException, Request
from rdkit import Chem

from app.api.routes.validation import extract_molecule_info
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.alerts import (
    AlertResultSchema,
    AlertScreenRequest,
    AlertScreenResponse,
    AlertSeverity,
    CatalogInfoSchema,
    CatalogListResponse,
    ConcernGroupSchema,
    UnifiedScreenRequest,
    UnifiedScreenResponse,
)
from app.services.alerts.alert_manager import AlertResult, alert_manager
from app.services.alerts.filter_catalog import AVAILABLE_CATALOGS
from app.services.alerts.unified_screen import unified_screen
from app.services.parser.molecule_parser import MoleculeFormat, parse_molecule

router = APIRouter()

_FORMAT_MAP = {
    "smiles": MoleculeFormat.SMILES,
    "inchi": MoleculeFormat.INCHI,
    "mol": MoleculeFormat.MOL,
    "auto": None,
}


def _parse_or_raise(molecule: str, fmt: str) -> Chem.Mol:
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


def _alert_to_schema(alert: AlertResult) -> AlertResultSchema:
    """Convert an AlertResult dataclass to AlertResultSchema."""
    return AlertResultSchema(
        pattern_name=alert.pattern_name,
        description=alert.description,
        severity=AlertSeverity(alert.severity.value),
        matched_atoms=alert.matched_atoms,
        catalog_source=alert.catalog_source,
        smarts=alert.smarts,
        reference=alert.reference,
        scope=alert.scope,
        filter_set=alert.filter_set,
        catalog_description=alert.catalog_description,
        category=alert.category,
        concern_group=getattr(alert, "concern_group", None),
    )


@router.post("/alerts", response_model=AlertScreenResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def screen_alerts(
    request: Request,
    body: AlertScreenRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Screen a molecule for structural alerts.

    Screens the input molecule against specified pattern catalogs
    (default: PAINS) and returns matched alerts with atom indices
    for highlighting.

    Args:
        body: Alert screening request with molecule and catalog selection

    Returns:
        AlertScreenResponse with matched alerts and screening metadata

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()

    mol = _parse_or_raise(body.molecule, body.format)
    mol_info = extract_molecule_info(mol, body.molecule)

    screening_result = alert_manager.screen(mol, catalogs=body.catalogs)
    alerts = [_alert_to_schema(a) for a in screening_result.alerts]

    execution_time = int((time.time() - start_time) * 1000)

    return AlertScreenResponse(
        status="completed",
        molecule_info=mol_info,
        alerts=alerts,
        total_alerts=screening_result.total_alerts,
        screened_catalogs=screening_result.screened_catalogs,
        has_critical=screening_result.has_critical,
        has_warning=screening_result.has_warning,
        execution_time_ms=execution_time,
    )


@router.get("/alerts/catalogs", response_model=CatalogListResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def list_catalogs(
    request: Request, api_key: Optional[str] = Depends(get_api_key)
):
    """
    List available structural alert catalogs.

    Returns information about all available pattern catalogs
    including PAINS (A/B/C), BRENK, NIH, and ZINC.

    Returns:
        CatalogListResponse with available catalogs and their descriptions
    """
    catalogs = {}
    for cat_type, cat_info in AVAILABLE_CATALOGS.items():
        catalogs[cat_type] = CatalogInfoSchema(
            name=cat_info["name"],
            description=cat_info["description"],
            pattern_count=cat_info["pattern_count"],
            severity=cat_info["severity"],
            note=cat_info.get("note"),
            reference=cat_info.get("reference"),
            scope=cat_info.get("scope"),
            doi=cat_info.get("doi"),
            pmid=cat_info.get("pmid"),
        )

    return CatalogListResponse(
        catalogs=catalogs,
        default_catalogs=["PAINS"],
    )


@router.post("/alerts/quick-check")
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def quick_check_alerts(
    request: Request,
    body: AlertScreenRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Quick check if molecule has any structural alerts.

    Faster than full screening - only checks for presence of alerts,
    not specific patterns or atom indices.

    Args:
        body: Alert screening request with molecule and catalog selection

    Returns:
        Dictionary with has_alerts boolean and checked catalogs
    """
    mol = _parse_or_raise(body.molecule, body.format)
    has_alerts = alert_manager.has_alerts(mol, catalogs=body.catalogs)

    return {
        "has_alerts": has_alerts,
        "checked_catalogs": body.catalogs,
    }


@router.post("/alerts/screen", response_model=UnifiedScreenResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def unified_screen_endpoint(
    request: Request,
    body: UnifiedScreenRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Unified structural alert screening across ALL catalogs.

    Screens the input molecule against all available alert sources:
    PAINS, BRENK, NIH, ZINC, ChEMBL sub-catalogs, custom SMARTS (21 patterns),
    Kazius toxicophores (29 patterns), and NIBR Novartis filters.

    Returns both a flat raw-alert list and a deduplicated concern-group view.

    Args:
        body: Request with molecule string and format hint

    Returns:
        UnifiedScreenResponse with alerts, concern_groups, and summary flags

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()

    mol = _parse_or_raise(body.molecule, body.format)
    mol_info = extract_molecule_info(mol, body.molecule)

    result = unified_screen(mol)

    alerts_schema = [_alert_to_schema(a) for a in result["alerts"]]

    concern_groups_schema: dict[str, ConcernGroupSchema] = {}
    for group_name, group_data in result["concern_groups"].items():
        concern_groups_schema[group_name] = ConcernGroupSchema(
            name=group_name,
            count=group_data["count"],
            severity=group_data["severity"],
            alerts=[_alert_to_schema(a) for a in group_data["alerts"]],
        )

    execution_time = int((time.time() - start_time) * 1000)

    return UnifiedScreenResponse(
        status="completed",
        molecule_info=mol_info,
        alerts=alerts_schema,
        concern_groups=concern_groups_schema,
        total_raw=result["total_raw"],
        total_deduped=result["total_deduped"],
        screened_catalogs=result["screened_catalogs"],
        has_critical=result["has_critical"],
        has_warning=result["has_warning"],
        execution_time_ms=execution_time,
    )
