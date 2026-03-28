"""
Integration tests for the QSAR-Ready API endpoints.

Tests the single-molecule endpoint with valid SMILES, invalid SMILES,
and custom config. Batch upload/status/results tests are skipped if
Redis is not available.
"""

import os

import pytest

# Set required env vars before importing app modules
os.environ.setdefault("SECRET_KEY", "test-secret-key-for-ci-validation-only")
os.environ.setdefault("DEBUG", "True")


def _is_redis_available() -> bool:
    """Check if Redis is available for batch endpoint tests."""
    try:
        import redis

        r = redis.from_url("redis://localhost:6379/0")
        r.ping()
        return True
    except Exception:
        return False


REDIS_AVAILABLE = _is_redis_available()


@pytest.fixture(scope="module")
def client():
    """Create a TestClient for the FastAPI app."""
    from fastapi.testclient import TestClient

    from app.main import app

    with TestClient(app) as c:
        yield c


# =============================================================================
# Single molecule endpoint tests
# =============================================================================


class TestQSARSingleEndpoint:
    """Tests for POST /api/v1/qsar-ready/single."""

    def test_single_pipeline_endpoint(self, client):
        """POST with aspirin SMILES returns 200 with full pipeline result."""
        # Aspirin SMILES
        aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
        payload = {
            "smiles": aspirin,
            "config": {
                "enable_metals": True,
                "enable_desalt": True,
                "enable_normalize": True,
                "enable_neutralize": True,
                "enable_tautomer": True,
                "enable_stereo_strip": False,
                "enable_isotope_strip": True,
                "min_heavy_atoms": 3,
                "max_heavy_atoms": 100,
                "max_mw": 1500.0,
                "remove_inorganics": True,
            },
        }
        response = client.post("/api/v1/qsar-ready/single", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        assert "original_smiles" in data
        assert "curated_smiles" in data
        assert "steps" in data
        assert "inchikey_changed" in data
        assert "status" in data
        assert data["status"] == "ok"
        assert data["original_smiles"] == aspirin
        assert data["curated_smiles"] is not None
        # Pipeline should have 10 steps
        assert len(data["steps"]) == 10, f"Expected 10 steps, got {len(data['steps'])}"

    def test_single_invalid_smiles(self, client):
        """POST with invalid SMILES returns 200 with status='rejected'."""
        payload = {
            "smiles": "not_a_valid_smiles_!!!",
            "config": {
                "enable_metals": True,
                "enable_desalt": True,
                "enable_normalize": True,
                "enable_neutralize": True,
                "enable_tautomer": True,
                "enable_stereo_strip": False,
                "enable_isotope_strip": True,
                "min_heavy_atoms": 3,
                "max_heavy_atoms": 100,
                "max_mw": 1500.0,
                "remove_inorganics": True,
            },
        }
        response = client.post("/api/v1/qsar-ready/single", json=payload)
        # Pipeline handles errors gracefully — returns 200 with error/rejected status
        assert response.status_code == 200, response.text

        data = response.json()
        assert data["status"] in ("rejected", "error"), (
            f"Expected rejected or error status for invalid SMILES, got {data['status']}"
        )
        assert data["curated_smiles"] is None

    def test_single_custom_config(self, client):
        """POST with enable_stereo_strip=True strips stereo from L-leucine."""
        # L-leucine with stereo center
        l_leucine = "CC(C)C[C@H](N)C(=O)O"
        payload = {
            "smiles": l_leucine,
            "config": {
                "enable_metals": True,
                "enable_desalt": True,
                "enable_normalize": True,
                "enable_neutralize": True,
                "enable_tautomer": True,
                "enable_stereo_strip": True,  # Should strip stereo
                "enable_isotope_strip": True,
                "min_heavy_atoms": 3,
                "max_heavy_atoms": 100,
                "max_mw": 1500.0,
                "remove_inorganics": True,
            },
        }
        response = client.post("/api/v1/qsar-ready/single", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        assert data["status"] == "ok"
        # Stereo step should have status "applied" or "no_change"
        stereo_step = next(
            (s for s in data["steps"] if s["step_name"] == "stereo"), None
        )
        assert stereo_step is not None, "Stereo step not found in results"
        assert stereo_step["enabled"] is True
        # Curated SMILES should not contain stereo markers if stereo was stripped
        if stereo_step["status"] == "applied":
            assert "@" not in (data["curated_smiles"] or ""), (
                "Stereo marker still present after stripping"
            )

    def test_single_inchikey_tracking(self, client):
        """POST returns both original_inchikey and standardized_inchikey."""
        payload = {
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
            "config": {
                "enable_metals": True,
                "enable_desalt": True,
                "enable_normalize": True,
                "enable_neutralize": True,
                "enable_tautomer": True,
                "enable_stereo_strip": False,
                "enable_isotope_strip": True,
                "min_heavy_atoms": 3,
                "max_heavy_atoms": 100,
                "max_mw": 1500.0,
                "remove_inorganics": True,
            },
        }
        response = client.post("/api/v1/qsar-ready/single", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        # Both InChIKey fields must be present (D-14 locked decision)
        assert "original_inchikey" in data
        assert "standardized_inchikey" in data
        assert "inchikey_changed" in data
        # Aspirin has no stereo/salts — InChIKeys should be equal
        assert data["inchikey_changed"] is False

    def test_single_step_structure(self, client):
        """Each step result has the required fields."""
        payload = {
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "config": {
                "enable_metals": True,
                "enable_desalt": True,
                "enable_normalize": True,
                "enable_neutralize": True,
                "enable_tautomer": True,
                "enable_stereo_strip": False,
                "enable_isotope_strip": True,
                "min_heavy_atoms": 3,
                "max_heavy_atoms": 100,
                "max_mw": 1500.0,
                "remove_inorganics": True,
            },
        }
        response = client.post("/api/v1/qsar-ready/single", json=payload)
        assert response.status_code == 200

        data = response.json()
        for step in data["steps"]:
            assert "step_name" in step
            assert "step_index" in step
            assert "enabled" in step
            assert "status" in step
            assert step["status"] in ("applied", "no_change", "skipped", "error")


# =============================================================================
# Schema validation tests
# =============================================================================


class TestQSARBatchSchemas:
    """Tests for QSAR batch Pydantic schemas."""

    def test_batch_status_schema_has_eta_seconds(self):
        """QSARBatchStatusResponse must accept eta_seconds field."""
        from app.schemas.qsar_ready import QSARBatchStatusResponse
        response = QSARBatchStatusResponse(
            job_id="test-uuid",
            status="processing",
            progress=50,
            processed=25,
            total=50,
            eta_seconds=120,
        )
        assert response.eta_seconds == 120

    def test_batch_status_schema_eta_seconds_optional(self):
        """eta_seconds should default to None when not provided."""
        from app.schemas.qsar_ready import QSARBatchStatusResponse
        response = QSARBatchStatusResponse(
            job_id="test-uuid",
            status="pending",
        )
        assert response.eta_seconds is None

    def test_batch_results_schema_has_total_results(self):
        """QSARBatchResultsResponse must accept total_results field."""
        from app.schemas.qsar_ready import (
            QSARBatchResultsResponse,
            QSARBatchSummary,
            QSARReadyConfigSchema,
        )
        config = QSARReadyConfigSchema()
        summary = QSARBatchSummary(total=10, ok=8, rejected=1, duplicate=1, error=0)
        response = QSARBatchResultsResponse(
            job_id="test-uuid",
            status="complete",
            config=config,
            summary=summary,
            total_results=10,
        )
        assert response.total_results == 10


# =============================================================================
# Batch endpoint tests (skipped if Redis not available)
# =============================================================================


@pytest.mark.skipif(not REDIS_AVAILABLE, reason="Redis not available")
class TestQSARBatchEndpoints:
    """Tests for batch upload, status, and results endpoints (requires Redis)."""

    def test_batch_upload_smiles_text(self, client):
        """POST /qsar-ready/batch/upload with smiles_text returns job_id."""
        smiles_text = "CC(=O)OC1=CC=CC=C1C(=O)O\nCC(C)C[C@H](N)C(=O)O"
        config_json = '{"enable_metals":true,"enable_desalt":true,"enable_normalize":true,"enable_neutralize":true,"enable_tautomer":true,"enable_stereo_strip":false,"enable_isotope_strip":true,"min_heavy_atoms":3,"max_heavy_atoms":100,"max_mw":1500.0,"remove_inorganics":true}'

        response = client.post(
            "/api/v1/qsar-ready/batch/upload",
            data={"config": config_json, "smiles_text": smiles_text},
        )
        assert response.status_code == 200, response.text
        data = response.json()
        assert "job_id" in data
        assert data["total_molecules"] == 2
        assert data["status"] == "pending"

    def test_batch_status_not_found(self, client):
        """GET /qsar-ready/batch/{job_id}/status returns 404 for unknown job."""
        fake_job_id = "00000000-0000-0000-0000-000000000000"
        response = client.get(f"/api/v1/qsar-ready/batch/{fake_job_id}/status")
        assert response.status_code == 404

    def test_batch_invalid_job_id(self, client):
        """GET with invalid job_id format returns 400."""
        response = client.get("/api/v1/qsar-ready/batch/not-a-uuid/status")
        assert response.status_code == 400
