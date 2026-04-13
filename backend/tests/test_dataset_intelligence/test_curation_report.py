"""Tests for the curation report service module."""

import json

from app.services.dataset_intelligence.curation_report import (
    CurationReport,
    build_curated_csv_rows,
    build_curation_report,
)
from app.services.dataset_intelligence.health_audit import DatasetHealthResult


def _make_health_result(**overrides) -> DatasetHealthResult:
    """Helper to create a DatasetHealthResult with sensible defaults."""
    defaults = {
        "overall_score": 75.0,
        "parsability_score": 0.95,
        "stereo_score": 0.90,
        "uniqueness_score": 0.85,
        "alert_score": 0.80,
        "std_consistency_score": 0.70,
        "weights": {
            "parsability": 0.25,
            "stereo": 0.15,
            "uniqueness": 0.20,
            "alerts": 0.20,
            "std_consistency": 0.20,
        },
        "molecule_count": 100,
        "parse_failures": 5,
        "stereo_undefined_count": 10,
        "duplicate_count": 15,
        "alert_hit_count": 20,
        "std_disagreement_count": 30,
        "std_sample_size": 100,
        "issues": [
            {
                "row_index": 0,
                "smiles": "CCO",
                "issue_type": "parse_failure",
                "severity": "critical",
                "description": "Cannot parse SMILES",
            },
            {
                "row_index": 5,
                "smiles": "c1ccccc1",
                "issue_type": "alert_hit",
                "severity": "warning",
                "description": "Triggered 2 alerts",
            },
        ],
        "property_distributions": {
            "mw": {"bins": [100.0, 200.0, 300.0], "counts": [5, 10]},
            "logp": {"bins": [0.0, 2.0, 4.0], "counts": [8, 12]},
            "tpsa": {"bins": [0.0, 50.0, 100.0], "counts": [7, 13]},
        },
        "std_pipeline_comparison": {
            "sample_size": 100,
            "agree_count": 70,
            "disagree_count": 30,
        },
        "contradictions": [],
        "dedup_groups": [{"inchikey": "IK_A", "rows": [0, 5]}],
        "numeric_columns": ["activity"],
    }
    defaults.update(overrides)
    return DatasetHealthResult(**defaults)


class TestBuildCurationReport:
    """Tests for build_curation_report function."""

    def test_required_keys(self):
        """Report contains all required top-level keys."""
        health_result = _make_health_result()
        file_metadata = {
            "filename": "test.csv",
            "format": "csv",
            "row_count": 100,
            "file_size_bytes": 5000,
            "sha256_hash": "abc123",
        }
        report = build_curation_report(health_result, file_metadata, [])
        assert "version" in report
        assert "generated_at" in report
        assert "file_metadata" in report
        assert "health_audit" in report
        assert "contradictory_labels" in report
        assert "standardization_decisions" in report
        assert "deduplication" in report
        assert "issues_flagged" in report

    def test_json_serializable(self):
        """Report is JSON serializable — json.dumps does not raise."""
        health_result = _make_health_result()
        file_metadata = {
            "filename": "test.csv",
            "format": "csv",
            "row_count": 100,
            "file_size_bytes": 5000,
            "sha256_hash": "abc123",
        }
        report = build_curation_report(health_result, file_metadata, [])
        # This MUST not raise
        serialized = json.dumps(report)
        assert isinstance(serialized, str)

    def test_generated_at_iso_format(self):
        """generated_at is in ISO 8601 format."""
        health_result = _make_health_result()
        file_metadata = {
            "filename": "test.csv",
            "format": "csv",
            "row_count": 100,
            "file_size_bytes": 5000,
            "sha256_hash": "abc123",
        }
        report = build_curation_report(health_result, file_metadata, [])
        generated_at = report["generated_at"]
        assert generated_at.endswith("Z")
        # Should contain date separator
        assert "T" in generated_at

    def test_version_is_1_0(self):
        """Report version is 1.0."""
        health_result = _make_health_result()
        file_metadata = {
            "filename": "test.csv",
            "format": "csv",
            "row_count": 100,
            "file_size_bytes": 5000,
            "sha256_hash": "abc123",
        }
        report = build_curation_report(health_result, file_metadata, [])
        assert report["version"] == "1.0"

    def test_health_audit_sub_scores(self):
        """Health audit section has overall score and sub-scores."""
        health_result = _make_health_result()
        file_metadata = {
            "filename": "test.csv",
            "format": "csv",
            "row_count": 100,
            "file_size_bytes": 5000,
            "sha256_hash": "abc123",
        }
        report = build_curation_report(health_result, file_metadata, [])
        ha = report["health_audit"]
        assert "overall_score" in ha
        assert "sub_scores" in ha
        assert "parsability" in ha["sub_scores"]
        assert "weights" in ha


class TestBuildCuratedCsvRows:
    """Tests for build_curated_csv_rows function."""

    def test_appends_issue_columns(self):
        """Curated rows have _health_issues, _is_duplicate, _standardized_smiles, _alert_flags."""
        health_result = _make_health_result()
        molecules = [
            {
                "index": 0,
                "smiles": "CCO",
                "mol": None,
                "inchikey": "IK_A",
                "properties": {"activity": "5.0"},
            },
            {
                "index": 5,
                "smiles": "c1ccccc1",
                "mol": None,
                "inchikey": "IK_B",
                "properties": {"activity": "10.0"},
            },
        ]
        rows = build_curated_csv_rows(molecules, health_result)
        assert len(rows) == 2
        for row in rows:
            assert "_health_issues" in row
            assert "_is_duplicate" in row
            assert "_standardized_smiles" in row
            assert "_alert_flags" in row

    def test_duplicate_detection_in_rows(self):
        """Molecule in dedup_groups is marked as duplicate."""
        health_result = _make_health_result(
            dedup_groups=[{"inchikey": "IK_A", "rows": [0, 5]}]
        )
        molecules = [
            {
                "index": 0,
                "smiles": "CCO",
                "mol": None,
                "inchikey": "IK_A",
                "properties": {},
            },
            {
                "index": 1,
                "smiles": "c1ccccc1",
                "mol": None,
                "inchikey": "IK_B",
                "properties": {},
            },
        ]
        rows = build_curated_csv_rows(molecules, health_result)
        # First molecule (IK_A) is in dedup group
        assert rows[0]["_is_duplicate"] == "true"
        # Second molecule (IK_B) is NOT in dedup group
        assert rows[1]["_is_duplicate"] == "false"

    def test_preserves_original_properties(self):
        """Original properties are preserved in curated rows."""
        health_result = _make_health_result()
        molecules = [
            {
                "index": 0,
                "smiles": "CCO",
                "mol": None,
                "inchikey": "IK_A",
                "properties": {"activity": "5.0", "name": "ethanol"},
            },
        ]
        rows = build_curated_csv_rows(molecules, health_result)
        assert rows[0]["activity"] == "5.0"
        assert rows[0]["name"] == "ethanol"


class TestCurationReportDataclass:
    """Tests for CurationReport dataclass."""

    def test_default_values(self):
        """CurationReport has sensible defaults."""
        report = CurationReport()
        assert report.version == "1.0"
        assert report.generated_at == ""
        assert isinstance(report.file_metadata, dict)
        assert isinstance(report.contradictory_labels, list)
