"""
Structural Alert Screening Services

Provides PAINS, BRENK, and other structural alert pattern screening
using RDKit FilterCatalog.
"""

from .filter_catalog import get_filter_catalog, AVAILABLE_CATALOGS
from .alert_manager import AlertManager, alert_manager

__all__ = [
    "get_filter_catalog",
    "AVAILABLE_CATALOGS",
    "AlertManager",
    "alert_manager",
]
