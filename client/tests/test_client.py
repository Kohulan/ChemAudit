"""
Tests for ChemAudit API client.

Uses respx for HTTP mocking to test client behavior without a live server.
"""

import json

import pytest
import httpx
import respx
from pathlib import Path
from chemaudit import ChemAuditClient
from chemaudit.exceptions import (
    APIError,
    RateLimitError,
    AuthenticationError,
    ValidationError,
    BatchJobNotFoundError,
    TimeoutError,
)


@pytest.fixture
def base_url():
    """Base URL for testing."""
    return "http://test-api.example.com"


@pytest.fixture
def client(base_url):
    """Create test client."""
    return ChemAuditClient(base_url=base_url, api_key="test-key")


@pytest.fixture
def mock_molecule_info():
    """Mock molecule_info block shared by most responses."""
    return {"input_smiles": "CCO", "canonical_smiles": "CCO"}


@pytest.fixture
def mock_validation_response():
    """Mock validation response data."""
    return {
        "status": "completed",
        "molecule_info": {
            "input_smiles": "CCO",
            "canonical_smiles": "CCO",
            "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "molecular_formula": "C2H6O",
            "molecular_weight": 46.069,
            "num_atoms": 9,
        },
        "overall_score": 95,
        "issues": [],
        "all_checks": [
            {
                "check_name": "parsability",
                "passed": True,
                "severity": "info",
                "message": "Molecule parsed successfully",
                "affected_atoms": [],
                "details": {},
            }
        ],
        "execution_time_ms": 15,
    }


@pytest.fixture
def mock_batch_upload_response():
    """Mock batch upload response."""
    return {
        "job_id": "test-job-123",
        "status": "pending",
        "total_molecules": 5,
        "message": "Job submitted successfully",
    }


@pytest.fixture
def mock_batch_status_response():
    """Mock batch status response."""
    return {
        "job_id": "test-job-123",
        "status": "processing",
        "progress": 60,
        "processed": 3,
        "total": 5,
        "eta_seconds": 10,
        "error_message": None,
    }


@pytest.fixture
def mock_batch_result_response():
    """Mock batch result response."""
    return {
        "job_id": "test-job-123",
        "status": "complete",
        "statistics": {
            "total": 5,
            "successful": 4,
            "errors": 1,
            "avg_validation_score": 85.5,
            "avg_ml_readiness_score": 78.2,
            "score_distribution": {
                "excellent": 2,
                "good": 2,
                "moderate": 0,
                "poor": 1,
            },
            "alert_summary": {
                "PAINS": 1,
                "BRENK": 0,
            },
            "processing_time_seconds": 12.5,
        },
        "results": [
            {
                "smiles": "CCO",
                "name": "Ethanol",
                "index": 0,
                "status": "success",
                "error": None,
                "validation": {"overall_score": 95},
            },
            {
                "smiles": "c1ccccc1",
                "name": "Benzene",
                "index": 1,
                "status": "success",
                "error": None,
                "validation": {"overall_score": 80},
            },
        ],
        "page": 1,
        "page_size": 100,
        "total_results": 5,
        "total_pages": 1,
    }


@respx.mock
def test_validate_success(client, base_url, mock_validation_response):
    """Test successful validation."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(200, json=mock_validation_response)
    )

    result = client.validate("CCO")

    assert route.called
    assert result.overall_score == 95
    assert result.molecule_info.canonical_smiles == "CCO"
    assert len(result.issues) == 0


@respx.mock
def test_validate_sends_molecule_and_format(client, base_url, mock_validation_response):
    """Test validation request body."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(200, json=mock_validation_response)
    )

    client.validate(molecule="CCO", format="smiles")

    request = route.calls.last.request
    assert json.loads(request.content) == {"molecule": "CCO", "format": "smiles"}
    assert request.headers["Content-Type"] == "application/json"


@respx.mock
def test_validate_parses_all_server_severities(client, base_url, mock_validation_response):
    """Server severities are lowercase and include "pass"."""
    checks = [
        {**mock_validation_response["all_checks"][0], "severity": sev}
        for sev in ("critical", "error", "warning", "info", "pass")
    ]
    respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(200, json={**mock_validation_response, "all_checks": checks})
    )

    result = client.validate("CCO")

    assert [c.severity.value for c in result.all_checks] == [
        "critical", "error", "warning", "info", "pass"
    ]


@respx.mock
def test_screen_alerts(client, base_url, mock_molecule_info):
    """Alert screening sends `molecule` and parses the server response."""
    route = respx.post(f"{base_url}/api/v1/alerts").mock(
        return_value=httpx.Response(
            200,
            json={
                "status": "completed",
                "molecule_info": mock_molecule_info,
                "alerts": [
                    {
                        "pattern_name": "quinone_A",
                        "description": "Quinone",
                        "severity": "warning",
                        "matched_atoms": [0, 1],
                        "catalog_source": "PAINS",
                        "smarts": "O=C1C=CC(=O)C=C1",
                    }
                ],
                "total_alerts": 1,
                "screened_catalogs": ["PAINS", "BRENK"],
                "has_critical": False,
                "has_warning": True,
                "execution_time_ms": 4,
            },
        )
    )

    result = client.screen_alerts("O=C1C=CC(=O)C=C1", catalogs=["PAINS", "BRENK"])

    body = json.loads(route.calls.last.request.content)
    assert body == {"molecule": "O=C1C=CC(=O)C=C1", "catalogs": ["PAINS", "BRENK"]}
    assert result.alerts[0].catalog_source == "PAINS"
    assert result.screened_catalogs == ["PAINS", "BRENK"]


@respx.mock
def test_score(client, base_url, mock_molecule_info):
    """Scoring sends `molecule` and exposes score blocks as attributes."""
    route = respx.post(f"{base_url}/api/v1/score").mock(
        return_value=httpx.Response(
            200,
            json={
                "status": "completed",
                "molecule_info": mock_molecule_info,
                "ml_readiness": {"score": 70, "label": "Good"},
                "execution_time_ms": 9,
            },
        )
    )

    result = client.score("CCO", include=["ml_readiness"])

    body = json.loads(route.calls.last.request.content)
    assert body == {"molecule": "CCO", "include": ["ml_readiness"]}
    assert result.ml_readiness["score"] == 70
    assert result.np_likeness is None


@respx.mock
def test_standardize(client, base_url, mock_molecule_info):
    """Standardization sends `molecule` + options and parses the nested result."""
    route = respx.post(f"{base_url}/api/v1/standardize").mock(
        return_value=httpx.Response(
            200,
            json={
                "molecule_info": mock_molecule_info,
                "result": {
                    "original_smiles": "CCN.Cl",
                    "standardized_smiles": "CCN",
                    "success": True,
                    "steps_applied": [{"step_name": "get_parent", "applied": True}],
                },
                "execution_time_ms": 12,
            },
        )
    )

    result = client.standardize("CCN.Cl", tautomer=True)

    body = json.loads(route.calls.last.request.content)
    assert body == {"molecule": "CCN.Cl", "options": {"include_tautomer": True}}
    assert result.result.standardized_smiles == "CCN"


@respx.mock
def test_session_cookie_is_not_sent_back(base_url, mock_validation_response):
    """A server session cookie would make CSRF protection reject later POSTs."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(
            200,
            json=mock_validation_response,
            headers={"Set-Cookie": "chemaudit_sid=abc; Path=/"},
        )
    )
    anonymous = ChemAuditClient(base_url=base_url)

    anonymous.validate("CCO")
    anonymous.validate("CCO")

    assert "cookie" not in route.calls.last.request.headers


@respx.mock
def test_error_message_from_error_body(client, base_url):
    """Errors with {"error", "message"} bodies (CSRF, rate limit) keep the message."""
    respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(
            403, json={"error": "csrf_validation_failed", "message": "Invalid CSRF token"}
        )
    )

    with pytest.raises(APIError, match="Invalid CSRF token"):
        client.validate("CCO")


@respx.mock
def test_error_with_non_json_body(client, base_url):
    """An HTML error page (e.g. from nginx) still raises APIError."""
    respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(413, text="<html>Request Entity Too Large</html>")
    )

    with pytest.raises(APIError) as exc_info:
        client.validate("CCO")

    assert exc_info.value.status_code == 413


@respx.mock
def test_rate_limit_error_no_retry(client, base_url):
    """Test rate limit error when max retries exhausted."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(
            429, headers={"Retry-After": "1"}, json={"detail": "Rate limit exceeded"}
        )
    )

    with pytest.raises(RateLimitError) as exc_info:
        client.validate("CCO")

    # Should retry 3 times (max_retries=3)
    assert route.call_count == 4  # 1 initial + 3 retries
    assert exc_info.value.retry_after == 1


@respx.mock
def test_rate_limit_error_with_retry(client, base_url, mock_validation_response):
    """Test rate limit error with successful retry."""
    # First call returns 429, second call succeeds
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        side_effect=[
            httpx.Response(429, headers={"Retry-After": "0"}),
            httpx.Response(200, json=mock_validation_response),
        ]
    )

    result = client.validate("CCO")

    assert route.call_count == 2
    assert result.overall_score == 95


@respx.mock
def test_authentication_error(client, base_url):
    """Test authentication error."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(401, json={"detail": "Invalid API key"})
    )

    with pytest.raises(AuthenticationError):
        client.validate("CCO")

    assert route.called


@respx.mock
def test_validation_error(client, base_url):
    """Test validation error (422)."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(422, json={"detail": "Invalid SMILES string"})
    )

    with pytest.raises(ValidationError):
        client.validate("INVALID")

    assert route.called


@respx.mock
def test_api_error_404(client, base_url):
    """Test 404 not found error."""
    route = respx.get(f"{base_url}/api/v1/batch/nonexistent/status").mock(
        return_value=httpx.Response(404, json={"detail": "Job not found"})
    )

    with pytest.raises(BatchJobNotFoundError) as exc_info:
        client.get_batch_status("nonexistent")

    assert route.called
    assert exc_info.value.status_code == 404


@respx.mock
def test_api_key_header(client, base_url, mock_validation_response):
    """Test that API key is sent in header."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(200, json=mock_validation_response)
    )

    client.validate("CCO")

    request = route.calls.last.request
    assert request.headers["X-API-Key"] == "test-key"


@respx.mock
def test_batch_submit(client, base_url, mock_batch_upload_response, tmp_path):
    """Test batch file submission."""
    # Create temporary CSV file
    csv_file = tmp_path / "test.csv"
    csv_file.write_text("smiles,name\nCCO,Ethanol\n")

    route = respx.post(f"{base_url}/api/v1/batch/upload").mock(
        return_value=httpx.Response(200, json=mock_batch_upload_response)
    )

    response = client.submit_batch(
        str(csv_file), smiles_column="smiles", include_standardization=True
    )

    assert route.called
    request = route.calls.last.request
    assert request.headers["Content-Type"].startswith("multipart/form-data")
    body = request.read()
    assert b'name="file"; filename="test.csv"' in body
    assert b'name="include_standardization"\r\n\r\ntrue' in body
    assert response.job_id == "test-job-123"
    assert response.total_molecules == 5


@respx.mock
def test_batch_status(client, base_url, mock_batch_status_response):
    """Test getting batch status."""
    route = respx.get(f"{base_url}/api/v1/batch/test-job-123/status").mock(
        return_value=httpx.Response(200, json=mock_batch_status_response)
    )

    status = client.get_batch_status("test-job-123")

    assert route.called
    assert status.job_id == "test-job-123"
    assert status.progress == 60
    assert status.processed == 3
    assert status.total == 5


@respx.mock
def test_batch_results(client, base_url, mock_batch_result_response):
    """Test getting batch results."""
    route = respx.get(f"{base_url}/api/v1/batch/test-job-123").mock(
        return_value=httpx.Response(200, json=mock_batch_result_response)
    )

    results = client.get_batch_results("test-job-123", page=1, page_size=100)

    assert route.called
    assert results.job_id == "test-job-123"
    assert len(results.results) == 2
    assert results.total_results == 5
    assert results.statistics.avg_validation_score == 85.5


@respx.mock
def test_batch_results_with_filters(client, base_url, mock_batch_result_response):
    """Test getting batch results with filters."""
    route = respx.get(f"{base_url}/api/v1/batch/test-job-123").mock(
        return_value=httpx.Response(200, json=mock_batch_result_response)
    )

    results = client.get_batch_results(
        "test-job-123",
        page=1,
        page_size=50,
        score_min=70,
        score_max=90,
        status="success",
    )

    assert route.called
    request = route.calls.last.request
    params = dict(request.url.params)
    assert params["min_score"] == "70"
    assert params["max_score"] == "90"
    assert params["status_filter"] == "success"


@respx.mock
def test_wait_for_batch_complete(client, base_url):
    """Test waiting for batch completion."""
    # Mock status calls: processing -> processing -> complete
    route = respx.get(f"{base_url}/api/v1/batch/test-job-123/status").mock(
        side_effect=[
            httpx.Response(
                200,
                json={
                    "job_id": "test-job-123",
                    "status": "processing",
                    "progress": 30,
                    "processed": 1,
                    "total": 3,
                },
            ),
            httpx.Response(
                200,
                json={
                    "job_id": "test-job-123",
                    "status": "processing",
                    "progress": 60,
                    "processed": 2,
                    "total": 3,
                },
            ),
            httpx.Response(
                200,
                json={
                    "job_id": "test-job-123",
                    "status": "complete",
                    "progress": 100,
                    "processed": 3,
                    "total": 3,
                },
            ),
        ]
    )

    job = client.wait_for_batch("test-job-123", poll_interval=0.1, timeout=10.0)

    assert route.call_count == 3
    assert job.status.value == "complete"
    assert job.progress == 100


@respx.mock
def test_wait_for_batch_failed(client, base_url):
    """Test waiting for batch that fails."""
    route = respx.get(f"{base_url}/api/v1/batch/test-job-123/status").mock(
        return_value=httpx.Response(
            200,
            json={
                "job_id": "test-job-123",
                "status": "failed",
                "progress": 50,
                "processed": 1,
                "total": 2,
                "error_message": "Processing error",
            },
        )
    )

    with pytest.raises(APIError) as exc_info:
        client.wait_for_batch("test-job-123", poll_interval=0.1)

    assert "failed" in str(exc_info.value).lower()


@respx.mock
def test_iter_batch_results(client, base_url):
    """Test iterating through all batch results."""
    # Mock two pages of results
    page1 = {
        "job_id": "test-job-123",
        "status": "complete",
        "statistics": None,
        "results": [
            {
                "smiles": "CCO",
                "name": None,
                "index": 0,
                "status": "success",
                "error": None,
                "validation": {"overall_score": 95},
            }
        ],
        "page": 1,
        "page_size": 1,
        "total_results": 2,
        "total_pages": 2,
    }

    page2 = {
        "job_id": "test-job-123",
        "status": "complete",
        "statistics": None,
        "results": [
            {
                "smiles": "c1ccccc1",
                "name": None,
                "index": 1,
                "status": "success",
                "error": None,
                "validation": {"overall_score": 80},
            }
        ],
        "page": 2,
        "page_size": 1,
        "total_results": 2,
        "total_pages": 2,
    }

    route = respx.get(f"{base_url}/api/v1/batch/test-job-123").mock(
        side_effect=[
            httpx.Response(200, json=page1),
            httpx.Response(200, json=page2),
        ]
    )

    results = list(client.iter_batch_results("test-job-123", page_size=1))

    assert route.call_count == 2
    assert len(results) == 2
    assert results[0].smiles == "CCO"
    assert results[1].smiles == "c1ccccc1"


@respx.mock
def test_context_manager(base_url, mock_validation_response):
    """Test client as context manager."""
    route = respx.post(f"{base_url}/api/v1/validate").mock(
        return_value=httpx.Response(200, json=mock_validation_response)
    )

    with ChemAuditClient(base_url=base_url) as client:
        result = client.validate("CCO")
        assert result.overall_score == 95

    # Client should be closed after context exit
    assert route.called


@respx.mock
def test_export_batch(client, base_url, tmp_path):
    """Test exporting batch results."""
    export_content = b"smiles,name,score\nCCO,Ethanol,95\n"

    route = respx.get(f"{base_url}/api/v1/batch/test-job-123/export").mock(
        return_value=httpx.Response(200, content=export_content)
    )

    output_path = tmp_path / "export.csv"
    path = client.export_batch(
        "test-job-123", format="csv", output_path=str(output_path)
    )

    assert route.called
    assert Path(path).exists()
    assert Path(path).read_bytes() == export_content

    # Check query parameters
    request = route.calls.last.request
    params = dict(request.url.params)
    assert params["format"] == "csv"


@respx.mock
def test_export_batch_uses_server_filename(client, base_url, tmp_path, monkeypatch):
    """Without output_path, the file gets the server's name and extension."""
    monkeypatch.chdir(tmp_path)
    respx.get(f"{base_url}/api/v1/batch/test-job-123/export").mock(
        return_value=httpx.Response(
            200,
            content=b"%PDF-1.7",
            headers={"Content-Disposition": 'attachment; filename="report_test.pdf"'},
        )
    )

    path = client.export_batch("test-job-123", format="pdf")

    assert path == "report_test.pdf"
    assert (tmp_path / "report_test.pdf").read_bytes() == b"%PDF-1.7"


@respx.mock
def test_get_batch_stats(client, base_url):
    """Statistics come from the dedicated /stats endpoint."""
    route = respx.get(f"{base_url}/api/v1/batch/test-job-123/stats").mock(
        return_value=httpx.Response(
            200, json={"total": 3, "successful": 3, "errors": 0, "avg_validation_score": 86.7}
        )
    )

    stats = client.get_batch_stats("test-job-123")

    assert route.called
    assert stats.avg_validation_score == 86.7


def test_client_properties(client, base_url):
    """Test client properties."""
    assert client.base_url == base_url
    assert client.api_key == "test-key"
    assert client.timeout == 30.0
    assert client.max_retries == 3
