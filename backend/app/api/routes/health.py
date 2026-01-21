"""
Health check endpoint.

Provides system status, version information, and RDKit availability.
"""
from fastapi import APIRouter
from app.schemas.common import HealthResponse
from app.core.config import settings

router = APIRouter()


@router.get("/health", response_model=HealthResponse)
async def health_check() -> HealthResponse:
    """
    Health check endpoint.

    Returns:
        HealthResponse with status, app info, and RDKit version
    """
    # Get RDKit version
    rdkit_version = None
    try:
        from rdkit import rdBase
        rdkit_version = rdBase.rdkitVersion
    except ImportError:
        pass

    return HealthResponse(
        status="healthy",
        app_name=settings.APP_NAME,
        app_version=settings.APP_VERSION,
        rdkit_version=rdkit_version,
    )
