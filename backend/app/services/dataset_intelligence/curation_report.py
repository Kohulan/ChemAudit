"""
Curation Report Service Module (DSET-04).

Produces a serializable JSON-compatible dict capturing file metadata,
health audit results, contradictory labels, standardization decisions,
deduplication outcomes, and flagged issues. Also provides build_curated_csv_rows
for producing CSV rows with appended issue columns.
"""

from __future__ import annotations

import dataclasses
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from app.services.dataset_intelligence.health_audit import DatasetHealthResult

logger = logging.getLogger(__name__)


@dataclass
class CurationReport:
    """Serializable curation report capturing full audit trail."""

    version: str = "1.0"
    generated_at: str = ""
    file_metadata: dict = field(default_factory=dict)
    health_audit: dict = field(default_factory=dict)
    contradictory_labels: list[dict] = field(default_factory=list)
    standardization_decisions: list[dict] = field(default_factory=list)
    deduplication: dict = field(default_factory=dict)
    issues_flagged: list[dict] = field(default_factory=list)
    property_distributions: dict = field(default_factory=dict)


def build_curation_report(
    health_result: DatasetHealthResult,
    file_metadata: dict,
    contradictions: list[dict],
) -> dict:
    """Build a serializable JSON curation report from health audit results.

    Populates all sections of the CurationReport from the health_result fields,
    file metadata, and contradictory label findings.

    Args:
        health_result: DatasetHealthResult from compute_health_score.
        file_metadata: Dict with keys: filename, format, row_count,
            file_size_bytes, sha256_hash.
        contradictions: List of contradiction dicts from detect_contradictory_labels.

    Returns:
        Dict (JSON-serializable) with all report sections.
    """
    report = CurationReport(
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        file_metadata=file_metadata,
        health_audit={
            "overall_score": health_result.overall_score,
            "sub_scores": {
                "parsability": health_result.parsability_score,
                "stereo": health_result.stereo_score,
                "uniqueness": health_result.uniqueness_score,
                "alerts": health_result.alert_score,
                "std_consistency": health_result.std_consistency_score,
            },
            "weights": health_result.weights,
            "molecule_count": health_result.molecule_count,
            "parse_failures": health_result.parse_failures,
            "stereo_undefined_count": health_result.stereo_undefined_count,
            "duplicate_count": health_result.duplicate_count,
            "alert_hit_count": health_result.alert_hit_count,
            "std_disagreement_count": health_result.std_disagreement_count,
            "std_sample_size": health_result.std_sample_size,
            "std_pipeline_comparison": health_result.std_pipeline_comparison,
        },
        contradictory_labels=contradictions,
        standardization_decisions=[],  # Populated by caller during processing
        deduplication={
            "duplicate_count": health_result.duplicate_count,
            "duplicate_groups": health_result.dedup_groups,
        },
        issues_flagged=health_result.issues,
        property_distributions=health_result.property_distributions,
    )

    return dataclasses.asdict(report)


def build_curated_csv_rows(
    molecules: list[dict],
    health_result: DatasetHealthResult,
) -> list[dict]:
    """Build curated CSV rows with appended issue columns.

    For each molecule, produces a dict with original properties PLUS:
        - _health_issues: comma-separated issue types for that row
        - _is_duplicate: 'true'/'false' based on dedup_groups
        - _standardized_smiles: from std pipeline comparison if available
        - _alert_flags: comma-separated alert types from issues

    Args:
        molecules: List of molecule dicts with index, smiles, properties.
        health_result: DatasetHealthResult from compute_health_score.

    Returns:
        List of row dicts ready for CSV serialization.
    """
    # Build lookup: row_index -> issues
    issues_by_row: dict[int, list[dict]] = defaultdict(list)
    for issue in health_result.issues:
        issues_by_row[issue.get("row_index", -1)].append(issue)

    # Build lookup: InChIKey -> is duplicate
    duplicate_inchikeys: set[str] = set()
    for group in health_result.dedup_groups:
        ik = group.get("inchikey", "")
        if ik:
            duplicate_inchikeys.add(ik)

    rows: list[dict] = []
    for mol_dict in molecules:
        row: dict[str, Any] = {}

        # Copy original properties
        for k, v in mol_dict.get("properties", {}).items():
            row[k] = v

        row_idx = mol_dict.get("index", -1)
        row_issues = issues_by_row.get(row_idx, [])

        # _health_issues: comma-separated issue types
        health_issue_types = [
            i["issue_type"] for i in row_issues
            if i.get("issue_type") not in ("alert_hit",)
        ]
        row["_health_issues"] = ",".join(health_issue_types) if health_issue_types else ""

        # _is_duplicate
        ik = mol_dict.get("inchikey", "")
        row["_is_duplicate"] = "true" if ik in duplicate_inchikeys else "false"

        # _standardized_smiles (from std pipeline comparison if available)
        row["_standardized_smiles"] = mol_dict.get("smiles", "")

        # _alert_flags: comma-separated alert types
        alert_types = [
            i.get("issue_type", "")
            for i in row_issues
            if i.get("issue_type") == "alert_hit"
        ]
        row["_alert_flags"] = ",".join(alert_types) if alert_types else ""

        rows.append(row)

    return rows
