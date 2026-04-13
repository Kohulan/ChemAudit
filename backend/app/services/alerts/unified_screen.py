"""
Unified Structural Alert Screening Orchestrator

Combines all alert sources (FilterCatalog-based: PAINS/Brenk/NIH/ZINC/ChEMBL,
custom SMARTS, Kazius toxicophores, and NIBR) into a single orchestrated screen.

Provides two views of results:
  1. ``alerts`` — flat list of all raw AlertResult objects (undeduplicated),
     preserving every catalog-source attribution for the catalog-source view.
  2. ``concern_groups`` — alerts grouped and deduplicated by functional concern.
     Within each group, alerts with the same (concern_group, pattern_name_normalized)
     key are merged; multiple catalog_source attributions are collected.

Note: complexity_filter is intentionally NOT called here — complexity percentile
is a separate assessment served by the safety/assess endpoint (Plan 02).
"""

from __future__ import annotations

from collections.abc import Iterable

from rdkit import Chem

from .alert_manager import AlertSeverity, alert_manager
from .custom_smarts import screen_custom_smarts
from .kazius_rules import screen_kazius
from .nibr_filters import screen_nibr
from .pattern_categories import classify_pattern

# Severity ordering for "worst severity" calculation
_SEVERITY_ORDER = {
    AlertSeverity.CRITICAL: 2,
    AlertSeverity.WARNING: 1,
    AlertSeverity.INFO: 0,
}


def _severity_value(severity: AlertSeverity) -> int:
    """Return numeric value for severity comparison."""
    return _SEVERITY_ORDER.get(severity, 0)


def _worst_severity(severities: Iterable[AlertSeverity]) -> str:
    """Return the worst (highest) severity from a collection."""
    best = max(severities, key=_severity_value, default=AlertSeverity.INFO)
    return best.value


def unified_screen(mol: Chem.Mol) -> dict:
    """
    Perform unified structural alert screening across all available sources.

    Sources screened:
      1. All FilterCatalog catalogs (PAINS, Brenk, NIH, ZINC, ChEMBL sub-catalogs)
         via AlertManager with ``catalogs=["ALL"]``.
      2. 21 custom SMARTS patterns (custom_smarts.py).
      3. 29 Kazius mutagenicity toxicophores (kazius_rules.py).
      4. NIBR Novartis screening deck filters (nibr_filters.py).

    Note: complexity_filter is NOT called — it belongs in the safety/assess
    endpoint (Plan 02) and is a separate concern from substructure alerts.

    Deduplication strategy:
      Within each concern group, alerts are deduplicated by the key
      ``(concern_group, pattern_name_normalized)`` where
      ``pattern_name_normalized = pattern_name.lower().replace(" ", "_")``.
      The first occurrence is kept; additional catalog_source values from
      duplicate entries are merged into the kept alert's sources list.

    Args:
        mol: RDKit molecule to screen.

    Returns:
        Dictionary with the following structure::

            {
                "alerts": list[AlertResult],  # all raw matches (undeduplicated)
                "concern_groups": {
                    "<group_name>": {
                        "alerts": list[AlertResult],   # deduplicated
                        "count": int,
                        "severity": str,  # worst severity in group
                    },
                    ...
                },
                "total_raw": int,
                "total_deduped": int,
                "screened_catalogs": list[str],
                "has_critical": bool,
                "has_warning": bool,
            }
    """
    if mol is None:
        return {
            "alerts": [],
            "concern_groups": {},
            "total_raw": 0,
            "total_deduped": 0,
            "screened_catalogs": [],
            "has_critical": False,
            "has_warning": False,
        }

    # ------------------------------------------------------------------ #
    # 1. FilterCatalog-based alerts (PAINS, Brenk, NIH, ZINC, ChEMBL)    #
    # ------------------------------------------------------------------ #
    screening_result = alert_manager.screen(mol, catalogs=["ALL"])
    for alert in screening_result.alerts:
        if alert.concern_group is None:
            alert.concern_group = classify_pattern(alert.pattern_name, alert.catalog_source)

    # ------------------------------------------------------------------ #
    # 2–4. New alert sources                                              #
    # ------------------------------------------------------------------ #
    custom_alerts = screen_custom_smarts(mol)
    kazius_alerts = screen_kazius(mol)
    try:
        nibr_alerts = screen_nibr(mol)
    except (ImportError, ModuleNotFoundError):
        nibr_alerts = []

    # Flat raw list (catalog-source view preserves all entries)
    raw_alerts = (
        list(screening_result.alerts)
        + custom_alerts
        + kazius_alerts
        + nibr_alerts
    )

    # Track all screened catalog sources
    screened_catalogs: list[str] = list(screening_result.screened_catalogs)
    for src in ["CUSTOM", "KAZIUS", "NIBR"]:
        if src not in screened_catalogs:
            screened_catalogs.append(src)

    # ------------------------------------------------------------------ #
    # 5–6. Build concern_groups with deduplication                        #
    # ------------------------------------------------------------------ #
    # Intermediate: group all alerts by concern_group
    groups: dict[str, list] = {}
    for alert in raw_alerts:
        group_key = alert.concern_group or "Unwanted Functionality"
        groups.setdefault(group_key, []).append(alert)

    concern_groups: dict[str, dict] = {}
    total_deduped = 0

    for group_name, group_alerts in groups.items():
        # Deduplicate within the group
        seen: dict[str, int] = {}  # normalized_key -> index in deduped list
        deduped: list = []

        for alert in group_alerts:
            norm_key = alert.pattern_name.lower().replace(" ", "_")
            dedup_key = f"{group_name}|{norm_key}"

            if dedup_key not in seen:
                seen[dedup_key] = len(deduped)
                # Make a copy-like view — we mutate sources list in place
                deduped.append(alert)
            else:
                # Merge: collect additional catalog_source into first occurrence
                existing = deduped[seen[dedup_key]]
                existing_src = existing.catalog_source
                new_src = alert.catalog_source
                if new_src != existing_src:
                    # Store merged sources in a supplementary attribute
                    if not hasattr(existing, "_extra_sources"):
                        existing._extra_sources = [existing_src]
                    if new_src not in existing._extra_sources:
                        existing._extra_sources.append(new_src)

        # Worst severity in this group
        severities = [a.severity for a in deduped]
        worst_sev = _worst_severity(severities)

        concern_groups[group_name] = {
            "alerts": deduped,
            "count": len(deduped),
            "severity": worst_sev,
        }
        total_deduped += len(deduped)

    # ------------------------------------------------------------------ #
    # 7. Summary flags                                                    #
    # ------------------------------------------------------------------ #
    has_critical = any(
        _severity_value(a.severity) >= _severity_value(AlertSeverity.CRITICAL)
        for a in raw_alerts
    )
    has_warning = any(
        _severity_value(a.severity) >= _severity_value(AlertSeverity.WARNING)
        for a in raw_alerts
    )

    return {
        "alerts": raw_alerts,
        "concern_groups": concern_groups,
        "total_raw": len(raw_alerts),
        "total_deduped": total_deduped,
        "screened_catalogs": screened_catalogs,
        "has_critical": has_critical,
        "has_warning": has_warning,
    }
