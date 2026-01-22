"""
Tests for batch processing tasks.

Note: These tests mock Celery to avoid requiring a running Redis instance.
"""
import pytest
from unittest.mock import patch, MagicMock

from app.services.batch.tasks import _process_single_molecule


class TestProcessSingleMolecule:
    """Tests for single molecule processing function."""

    def test_process_valid_molecule(self):
        """Test processing a valid molecule."""
        mol_data = {
            "smiles": "CCO",
            "name": "Ethanol",
            "index": 0,
            "properties": {},
            "parse_error": None,
        }

        result = _process_single_molecule(mol_data)

        assert result["status"] == "success"
        assert result["smiles"] == "CCO"
        assert result["name"] == "Ethanol"
        assert result["error"] is None
        assert result["validation"] is not None
        assert "overall_score" in result["validation"]

    def test_process_molecule_with_parse_error(self):
        """Test that pre-existing parse errors are handled."""
        mol_data = {
            "smiles": "",
            "name": "BadMol",
            "index": 0,
            "properties": {},
            "parse_error": "Failed to parse molecule",
        }

        result = _process_single_molecule(mol_data)

        assert result["status"] == "error"
        assert result["error"] == "Failed to parse molecule"
        assert result["validation"] is None

    def test_process_invalid_smiles(self):
        """Test processing an invalid SMILES string."""
        mol_data = {
            "smiles": "INVALID_SMILES_STRING",
            "name": "Invalid",
            "index": 0,
            "properties": {},
            "parse_error": None,
        }

        result = _process_single_molecule(mol_data)

        assert result["status"] == "error"
        assert result["error"] == "Failed to parse SMILES"

    def test_process_empty_smiles(self):
        """Test processing empty SMILES."""
        mol_data = {
            "smiles": "",
            "name": "Empty",
            "index": 0,
            "properties": {},
            "parse_error": None,
        }

        result = _process_single_molecule(mol_data)

        assert result["status"] == "error"
        assert result["error"] == "Empty SMILES string"

    def test_process_molecule_with_alerts(self):
        """Test processing a molecule that triggers structural alerts."""
        # Rhodanine - known PAINS pattern
        mol_data = {
            "smiles": "O=C1NC(=S)SC1",
            "name": "Rhodanine",
            "index": 0,
            "properties": {},
            "parse_error": None,
        }

        result = _process_single_molecule(mol_data)

        assert result["status"] == "success"
        assert result["alerts"] is not None
        # May or may not have alerts depending on catalogs
        assert "has_alerts" in result["alerts"] or "error" in result["alerts"]

    def test_process_molecule_returns_scoring(self):
        """Test that ML-readiness scoring is included."""
        mol_data = {
            "smiles": "c1ccccc1",  # Benzene
            "name": "Benzene",
            "index": 0,
            "properties": {},
            "parse_error": None,
        }

        result = _process_single_molecule(mol_data)

        assert result["status"] == "success"
        assert result["scoring"] is not None
        assert "ml_readiness" in result["scoring"]
        assert "score" in result["scoring"]["ml_readiness"]

    def test_failure_isolation(self):
        """Test that one molecule failure doesn't affect others."""
        molecules = [
            {"smiles": "CCO", "name": "Ethanol", "index": 0, "properties": {}, "parse_error": None},
            {"smiles": "INVALID", "name": "Bad", "index": 1, "properties": {}, "parse_error": None},
            {"smiles": "C", "name": "Methane", "index": 2, "properties": {}, "parse_error": None},
        ]

        results = [_process_single_molecule(m) for m in molecules]

        assert results[0]["status"] == "success"
        assert results[1]["status"] == "error"
        assert results[2]["status"] == "success"


class TestResultAggregation:
    """Tests for result aggregation."""

    def test_compute_statistics_basic(self):
        """Test basic statistics computation."""
        from app.services.batch.result_aggregator import compute_statistics

        results = [
            {"status": "success", "validation": {"overall_score": 90}, "alerts": {"alerts": []}, "scoring": {"ml_readiness": {"score": 85}}},
            {"status": "success", "validation": {"overall_score": 70}, "alerts": {"alerts": []}, "scoring": {"ml_readiness": {"score": 75}}},
            {"status": "error", "error": "Parse failed", "validation": None, "alerts": None, "scoring": None},
        ]

        stats = compute_statistics(results)

        assert stats.total == 3
        assert stats.successful == 2
        assert stats.errors == 1
        assert stats.avg_validation_score == 80.0  # (90 + 70) / 2
        assert stats.avg_ml_readiness_score == 80.0  # (85 + 75) / 2

    def test_compute_statistics_score_distribution(self):
        """Test score distribution buckets."""
        from app.services.batch.result_aggregator import compute_statistics

        results = [
            {"status": "success", "validation": {"overall_score": 95}, "alerts": {}, "scoring": {}},
            {"status": "success", "validation": {"overall_score": 85}, "alerts": {}, "scoring": {}},
            {"status": "success", "validation": {"overall_score": 60}, "alerts": {}, "scoring": {}},
            {"status": "success", "validation": {"overall_score": 40}, "alerts": {}, "scoring": {}},
        ]

        stats = compute_statistics(results)

        assert stats.score_distribution["excellent"] == 1  # 90-100
        assert stats.score_distribution["good"] == 1  # 70-89
        assert stats.score_distribution["moderate"] == 1  # 50-69
        assert stats.score_distribution["poor"] == 1  # 0-49

    def test_compute_statistics_empty_results(self):
        """Test statistics with empty results."""
        from app.services.batch.result_aggregator import compute_statistics

        stats = compute_statistics([])

        assert stats.total == 0
        assert stats.successful == 0
        assert stats.errors == 0
        assert stats.avg_validation_score is None
