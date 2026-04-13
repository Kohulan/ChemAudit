"""
Integration tests for the GenChem Filter API endpoints (Phase 11).

Tests:
- POST /api/v1/genchem/filter       — funnel pipeline (sync path ≤1000)
- POST /api/v1/genchem/score        — composite 0-1 scorer
- POST /api/v1/genchem/reinvent-score — REINVENT 4 contract

Ibuprofen (CC(C)Cc1ccc(cc1)C(C)C(=O)O, MW~206) is used as the reference
drug-like SMILES because it reliably passes the drug_like property thresholds
(min_mw=200). Ethanol (CCO, MW=46) would be rejected by the property stage.
"""

import os

# Set required env vars before importing app modules
os.environ.setdefault("SECRET_KEY", "test-secret-key-for-ci-validation-only")
os.environ.setdefault("DEBUG", "True")

import pytest

# ---------------------------------------------------------------------------
# SMILES fixtures
# ---------------------------------------------------------------------------

# Ibuprofen: MW~206, drug-like, no PAINS/Brenk/Kazius alerts expected
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
# Aspirin: MW~180, below drug_like min_mw=200 — rejected at property stage
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
# Totally invalid SMILES
INVALID_SMILES = "not_valid###"


@pytest.fixture(scope="module")
def client():
    """Create a TestClient for the FastAPI app."""
    from fastapi.testclient import TestClient

    from app.main import app

    with TestClient(app) as c:
        yield c


# =============================================================================
# Filter endpoint tests
# =============================================================================


class TestGenChemFilterEndpoint:
    """Tests for POST /api/v1/genchem/filter."""

    def test_filter_endpoint(self, client):
        """POST with valid SMILES and drug_like preset returns funnel result."""
        payload = {
            "smiles_list": [IBUPROFEN, ASPIRIN],
            "preset": "drug_like",
        }
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        assert "input_count" in data
        assert "output_count" in data
        assert "stages" in data
        assert "molecules" in data
        assert data["input_count"] == 2
        # At least 6 stages expected (parse, valence, alerts, property, sa, dedup)
        assert len(data["stages"]) >= 6

    def test_filter_stage_structure(self, client):
        """Each stage result has the required fields."""
        payload = {"smiles_list": [IBUPROFEN], "preset": "drug_like"}
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 200

        data = response.json()
        for stage in data["stages"]:
            assert "stage_name" in stage
            assert "stage_index" in stage
            assert "input_count" in stage
            assert "passed_count" in stage
            assert "rejected_count" in stage
            assert "enabled" in stage

    def test_filter_molecule_structure(self, client):
        """Each molecule result has the required fields."""
        payload = {"smiles_list": [IBUPROFEN, INVALID_SMILES], "preset": "drug_like"}
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 200

        data = response.json()
        for mol in data["molecules"]:
            assert "smiles" in mol
            assert "status" in mol
            assert mol["status"] in ("passed", "rejected", "duplicate", "error")

    def test_filter_invalid_smiles_rejected(self, client):
        """Invalid SMILES is rejected at parse stage."""
        payload = {"smiles_list": [INVALID_SMILES], "preset": "drug_like"}
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 200

        data = response.json()
        assert data["output_count"] == 0
        assert data["molecules"][0]["status"] == "rejected"
        assert data["molecules"][0]["failed_at"] == "parse"

    def test_filter_invalid_preset(self, client):
        """POST with unknown preset returns 400."""
        payload = {"smiles_list": [IBUPROFEN], "preset": "nonexistent_preset"}
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 400

    def test_filter_with_config(self, client):
        """POST with explicit config dict (no preset) returns 200."""
        config = {
            "min_mw": 100.0,
            "max_mw": 800.0,
            "min_logp": -5.0,
            "max_logp": 8.0,
            "max_tpsa": 200.0,
            "max_rot_bonds": 15,
            "max_rings": None,
            "max_sa_score": 7.0,
            "use_pains": False,
            "use_brenk": False,
            "use_kazius": False,
            "use_nibr": False,
            "enable_novelty": False,
            "novelty_threshold": 0.85,
            "weight_validity": 0.4,
            "weight_qed": 0.3,
            "weight_alert_free": 0.1,
            "weight_sa": 0.2,
        }
        payload = {"smiles_list": [IBUPROFEN, ASPIRIN], "config": config}
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        assert data["input_count"] == 2
        assert "stages" in data
        assert "molecules" in data

    def test_filter_dedup(self, client):
        """Duplicate SMILES is tracked as 'duplicate' status."""
        payload = {
            "smiles_list": [IBUPROFEN, IBUPROFEN],
            "preset": "permissive",
        }
        response = client.post("/api/v1/genchem/filter", json=payload)
        assert response.status_code == 200

        data = response.json()
        statuses = [m["status"] for m in data["molecules"]]
        assert "duplicate" in statuses, f"Expected duplicate in {statuses}"

    def test_filter_all_presets(self, client):
        """All 4 valid presets return 200."""
        for preset in ("drug_like", "lead_like", "fragment_like", "permissive"):
            payload = {"smiles_list": [IBUPROFEN], "preset": preset}
            response = client.post("/api/v1/genchem/filter", json=payload)
            assert response.status_code == 200, f"Preset '{preset}' failed: {response.text}"


# =============================================================================
# Score endpoint tests
# =============================================================================


class TestGenChemScoreEndpoint:
    """Tests for POST /api/v1/genchem/score."""

    def test_score_endpoint(self, client):
        """POST with valid and invalid SMILES returns scores list."""
        payload = {"smiles_list": [IBUPROFEN, INVALID_SMILES]}
        response = client.post("/api/v1/genchem/score", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        assert "scores" in data
        assert len(data["scores"]) == 2

    def test_score_float_for_valid(self, client):
        """Score for valid SMILES is a float in [0, 1]."""
        payload = {"smiles_list": [IBUPROFEN]}
        response = client.post("/api/v1/genchem/score", json=payload)
        assert response.status_code == 200

        data = response.json()
        score = data["scores"][0]
        assert score is not None, "Expected float score for valid SMILES"
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_score_null_for_invalid(self, client):
        """Score for invalid SMILES is null (None), not 0.0 — per D-14/Pitfall 4."""
        payload = {"smiles_list": [INVALID_SMILES]}
        response = client.post("/api/v1/genchem/score", json=payload)
        assert response.status_code == 200

        data = response.json()
        score = data["scores"][0]
        assert score is None, (
            f"Expected null for invalid SMILES (D-14 Pitfall 4), got {score!r}"
        )

    def test_score_preserves_order(self, client):
        """Scores list preserves input order — null at correct positions."""
        payload = {
            "smiles_list": [IBUPROFEN, INVALID_SMILES, IBUPROFEN],
            "preset": "drug_like",
        }
        response = client.post("/api/v1/genchem/score", json=payload)
        assert response.status_code == 200

        data = response.json()
        scores = data["scores"]
        assert len(scores) == 3
        assert scores[0] is not None  # valid
        assert scores[1] is None  # invalid
        assert scores[2] is not None  # valid

    def test_score_invalid_preset_returns_400(self, client):
        """POST with unknown preset returns 400."""
        payload = {"smiles_list": [IBUPROFEN], "preset": "nonexistent"}
        response = client.post("/api/v1/genchem/score", json=payload)
        assert response.status_code == 400


# =============================================================================
# REINVENT score endpoint tests
# =============================================================================


class TestGenChemReinventScoreEndpoint:
    """Tests for POST /api/v1/genchem/reinvent-score."""

    def test_reinvent_score_endpoint(self, client):
        """POST raw list returns {output: {successes_list: [...]}} contract."""
        payload = [
            {"input_string": IBUPROFEN, "query_id": "0"},
            {"input_string": INVALID_SMILES, "query_id": "1"},
        ]
        response = client.post("/api/v1/genchem/reinvent-score", json=payload)
        assert response.status_code == 200, response.text

        data = response.json()
        assert "output" in data
        assert "successes_list" in data["output"]

    def test_reinvent_omits_invalid(self, client):
        """Invalid SMILES (query_id='1') is omitted from successes_list — not scored 0.0."""
        payload = [
            {"input_string": IBUPROFEN, "query_id": "0"},
            {"input_string": INVALID_SMILES, "query_id": "1"},
        ]
        response = client.post("/api/v1/genchem/reinvent-score", json=payload)
        assert response.status_code == 200

        data = response.json()
        successes = data["output"]["successes_list"]
        query_ids = [s["query_id"] for s in successes]

        # query_id "1" (invalid SMILES) must be absent
        assert "1" not in query_ids, (
            f"Invalid SMILES query_id '1' should be omitted, got successes: {successes}"
        )
        # query_id "0" (valid SMILES) must be present
        assert "0" in query_ids, (
            f"Valid SMILES query_id '0' should be in successes, got: {query_ids}"
        )

    def test_reinvent_success_item_structure(self, client):
        """Each success item has query_id and output_value fields."""
        payload = [{"input_string": IBUPROFEN, "query_id": "test-id"}]
        response = client.post("/api/v1/genchem/reinvent-score", json=payload)
        assert response.status_code == 200

        data = response.json()
        successes = data["output"]["successes_list"]
        assert len(successes) == 1
        item = successes[0]
        assert "query_id" in item
        assert "output_value" in item
        assert item["query_id"] == "test-id"
        assert isinstance(item["output_value"], float)
        assert 0.0 <= item["output_value"] <= 1.0

    def test_reinvent_all_invalid_returns_empty_list(self, client):
        """All-invalid input returns empty successes_list (not null/error)."""
        payload = [
            {"input_string": INVALID_SMILES, "query_id": "0"},
            {"input_string": "another###invalid", "query_id": "1"},
        ]
        response = client.post("/api/v1/genchem/reinvent-score", json=payload)
        assert response.status_code == 200

        data = response.json()
        successes = data["output"]["successes_list"]
        assert successes == [], (
            f"Expected empty successes_list for all-invalid input, got: {successes}"
        )

    def test_reinvent_preset_query_param(self, client):
        """preset query parameter is accepted (default 'drug_like')."""
        payload = [{"input_string": IBUPROFEN, "query_id": "0"}]
        for preset in ("drug_like", "lead_like", "fragment_like", "permissive"):
            response = client.post(
                f"/api/v1/genchem/reinvent-score?preset={preset}", json=payload
            )
            assert response.status_code == 200, (
                f"Preset '{preset}' failed: {response.text}"
            )

    def test_reinvent_invalid_preset_returns_400(self, client):
        """Unknown preset query param returns 400."""
        payload = [{"input_string": IBUPROFEN, "query_id": "0"}]
        response = client.post(
            "/api/v1/genchem/reinvent-score?preset=nonexistent", json=payload
        )
        assert response.status_code == 400


# =============================================================================
# Batch status/results — minimal tests without Redis
# =============================================================================


class TestGenChemBatchEndpoints:
    """Tests for batch status and results endpoints (404/400 paths only)."""

    def test_batch_status_not_found(self, client):
        """GET status for unknown job returns 404."""
        fake_id = "00000000-0000-0000-0000-000000000000"
        response = client.get(f"/api/v1/genchem/batch/{fake_id}/status")
        # 404 or 503 (if Redis not available) — both are acceptable non-200 responses
        assert response.status_code in (404, 503), response.text

    def test_batch_invalid_uuid_returns_400(self, client):
        """GET status with non-UUID job_id returns 400."""
        response = client.get("/api/v1/genchem/batch/not-a-uuid/status")
        assert response.status_code == 400

    def test_batch_results_not_found(self, client):
        """GET results for unknown job returns 404."""
        fake_id = "00000000-0000-0000-0000-000000000001"
        response = client.get(f"/api/v1/genchem/batch/{fake_id}/results")
        assert response.status_code in (404, 503), response.text

    def test_batch_download_invalid_format(self, client):
        """GET download with unsupported format returns 400."""
        fake_id = "00000000-0000-0000-0000-000000000002"
        response = client.get(f"/api/v1/genchem/batch/{fake_id}/download/unsupported_format")
        assert response.status_code == 400
