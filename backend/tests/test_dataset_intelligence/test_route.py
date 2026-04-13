"""
Integration tests for the Dataset Intelligence API routes (Phase 12, Plan 02).

Tests:
- POST /api/v1/dataset/upload with CSV file
- POST /api/v1/dataset/upload with invalid file type
- GET  /api/v1/dataset/{job_id}/status for nonexistent job
- GET  /api/v1/dataset/{job_id}/results before completion (202)
- GET  /api/v1/dataset/{job_id}/download/report for nonexistent job
- GET  /api/v1/dataset/{job_id}/download/csv for nonexistent job
"""

import io
import json
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def client():
    """Create a TestClient with mocked settings (DEBUG=True for tests)."""
    with patch.dict("os.environ", {"DEBUG": "true"}):
        from app.main import app

        return TestClient(app)


@pytest.fixture
def sample_csv_bytes():
    """Small CSV file with 3 molecule rows."""
    content = "SMILES,Name,Activity\n"
    content += "CCO,ethanol,1.5\n"
    content += "CC(=O)O,acetic_acid,2.3\n"
    content += "c1ccccc1,benzene,0.8\n"
    return content.encode("utf-8")


class TestUploadEndpoint:
    """Tests for POST /api/v1/dataset/upload."""

    def test_upload_csv_returns_job_id(self, client, sample_csv_bytes):
        """Upload a valid CSV file and verify response contains job_id."""
        mock_task = MagicMock()
        mock_task.delay = MagicMock()

        with (
            patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis,
            patch(
                "app.services.dataset_intelligence.batch_processor.process_dataset_audit",
                mock_task,
            ),
        ):
            mock_r = MagicMock()
            mock_redis.return_value = mock_r

            response = client.post(
                "/api/v1/dataset/upload",
                files={"file": ("test.csv", io.BytesIO(sample_csv_bytes), "text/csv")},
            )

        assert response.status_code == 200
        data = response.json()
        assert "job_id" in data
        assert data["filename"] == "test.csv"
        assert data["file_type"] == "csv"
        assert data["status"] == "pending"
        assert data["message"] == "Dataset audit job submitted"

        # Verify Celery task was dispatched
        mock_task.delay.assert_called_once()
        call_args = mock_task.delay.call_args
        assert call_args[0][1] == "test.csv"  # filename
        assert call_args[0][2] == "csv"  # file_type

    def test_upload_invalid_file_type(self, client):
        """Upload a .txt file returns 400."""
        response = client.post(
            "/api/v1/dataset/upload",
            files={
                "file": ("test.txt", io.BytesIO(b"not a csv"), "text/plain")
            },
        )
        assert response.status_code == 400
        data = response.json()
        assert "Invalid file type" in data["detail"]

    def test_upload_sdf_returns_job_id(self, client):
        """Upload a valid SDF file returns 200 with job_id."""
        mock_task = MagicMock()
        mock_task.delay = MagicMock()

        sdf_content = b"""
  Mrv2211 08092307142D

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
"""
        with (
            patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis,
            patch(
                "app.services.dataset_intelligence.batch_processor.process_dataset_audit",
                mock_task,
            ),
        ):
            mock_r = MagicMock()
            mock_redis.return_value = mock_r

            response = client.post(
                "/api/v1/dataset/upload",
                files={
                    "file": (
                        "test.sdf",
                        io.BytesIO(sdf_content),
                        "chemical/x-mdl-sdfile",
                    )
                },
            )

        assert response.status_code == 200
        data = response.json()
        assert data["file_type"] == "sdf"
        assert data["status"] == "pending"


class TestStatusEndpoint:
    """Tests for GET /api/v1/dataset/{job_id}/status."""

    def test_nonexistent_job_returns_404(self, client):
        """Querying status for a nonexistent job returns 404."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.hgetall.return_value = {}
            mock_redis.return_value = mock_r

            response = client.get(f"/api/v1/dataset/{fake_job_id}/status")

        assert response.status_code == 404

    def test_invalid_job_id_format_returns_400(self, client):
        """Invalid UUID format returns 400."""
        response = client.get("/api/v1/dataset/not-a-uuid/status")
        assert response.status_code == 400

    def test_processing_job_returns_status(self, client):
        """Querying status for a processing job returns correct progress."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.hgetall.return_value = {
                "status": "processing",
                "progress": "45.5",
                "current_stage": "health_alerts",
            }
            mock_redis.return_value = mock_r

            response = client.get(f"/api/v1/dataset/{fake_job_id}/status")

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "processing"
        assert data["progress"] == 45.5
        assert data["current_stage"] == "health_alerts"


class TestResultsEndpoint:
    """Tests for GET /api/v1/dataset/{job_id}/results."""

    def test_processing_job_returns_202(self, client):
        """Requesting results while job is processing returns 202."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.hgetall.return_value = {
                "status": "processing",
                "progress": "50",
            }
            mock_r.get.return_value = None
            mock_redis.return_value = mock_r

            response = client.get(f"/api/v1/dataset/{fake_job_id}/results")

        assert response.status_code == 202
        data = response.json()
        assert data["status"] == "processing"

    def test_complete_job_returns_results(self, client):
        """Completed job returns full audit results."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        results = {
            "health_audit": {
                "overall_score": 85.0,
                "sub_scores": [
                    {"name": "parsability", "score": 1.0, "weight": 0.25,
                     "count": 0, "total": 3},
                ],
                "weights": {"parsability": 0.25},
                "molecule_count": 3,
                "issues": [],
                "property_distributions": {},
                "std_pipeline_comparison": {},
                "std_sample_size": 3,
                "dedup_groups": [],
            },
            "contradictions": [],
            "numeric_columns": [{"name": "Activity", "priority": 1}],
            "curation_report": {"version": "1.0"},
            "curated_csv_available": True,
        }

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.hgetall.return_value = {"status": "complete", "progress": "100"}
            mock_r.get.return_value = json.dumps(results)
            mock_redis.return_value = mock_r

            response = client.get(f"/api/v1/dataset/{fake_job_id}/results")

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "complete"
        assert data["health_audit"]["overall_score"] == 85.0
        assert data["curated_csv_available"] is True


class TestDownloadEndpoints:
    """Tests for download endpoints."""

    def test_download_report_nonexistent_returns_404(self, client):
        """Downloading report for nonexistent job returns 404."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.get.return_value = None
            mock_redis.return_value = mock_r

            response = client.get(
                f"/api/v1/dataset/{fake_job_id}/download/report"
            )

        assert response.status_code == 404

    def test_download_csv_nonexistent_returns_404(self, client):
        """Downloading CSV for nonexistent job returns 404."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.get.return_value = None
            mock_redis.return_value = mock_r

            response = client.get(
                f"/api/v1/dataset/{fake_job_id}/download/csv"
            )

        assert response.status_code == 404

    def test_download_report_success(self, client):
        """Downloading report for completed job returns JSON attachment."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        results = {
            "curation_report": {
                "version": "1.0",
                "generated_at": "2026-03-28T10:00:00Z",
            },
        }

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.get.return_value = json.dumps(results)
            mock_redis.return_value = mock_r

            response = client.get(
                f"/api/v1/dataset/{fake_job_id}/download/report"
            )

        assert response.status_code == 200
        assert response.headers["content-type"] == "application/json"
        assert "attachment" in response.headers["content-disposition"]
        data = response.json()
        assert data["version"] == "1.0"

    def test_download_csv_success(self, client):
        """Downloading curated CSV for completed job returns CSV attachment."""
        fake_job_id = "12345678-1234-1234-1234-123456789abc"

        results = {
            "curated_csv_rows": [
                {"Activity": "1.5", "_health_issues": "", "_is_duplicate": "false",
                 "_standardized_smiles": "CCO", "_alert_flags": ""},
                {"Activity": "2.3", "_health_issues": "", "_is_duplicate": "false",
                 "_standardized_smiles": "CC(=O)O", "_alert_flags": ""},
            ],
        }

        with patch("app.api.routes.dataset_intelligence._get_redis") as mock_redis:
            mock_r = MagicMock()
            mock_r.get.return_value = json.dumps(results)
            mock_redis.return_value = mock_r

            response = client.get(
                f"/api/v1/dataset/{fake_job_id}/download/csv"
            )

        assert response.status_code == 200
        assert "text/csv" in response.headers["content-type"]
        assert "attachment" in response.headers["content-disposition"]
        # Verify CSV content has header + 2 data rows
        lines = response.text.strip().split("\n")
        assert len(lines) == 3  # header + 2 rows
