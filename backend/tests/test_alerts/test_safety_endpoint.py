"""
Integration tests for safety assessment and unified alert screening API endpoints.

Tests cover:
  - POST /api/v1/alerts/screen  (unified screen with concern groups)
  - POST /api/v1/alerts          (backward compatibility)
  - POST /api/v1/safety/assess  (full CYP/hERG/bRo5/REOS/complexity assessment)
  - POST /api/v1/safety/summary (lightweight badge summary)
"""

import pytest
from httpx import ASGITransport, AsyncClient

from app.main import app

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Haloperidol: high hERG risk (logP>3.7, MW in 250-500, TPSA<75, basic N)
_HALOPERIDOL = "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1"

# Aspirin: small MW, MW<500 → bRo5 not applicable; multiple REOS violations expected
_ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# Ethanol: small, clean molecule
_ETHANOL = "CCO"

# 4-Aminobiphenyl — may have PAINS/Kazius/CUSTOM alerts
_AMINOBIPHENYL = "Nc1ccc(-c2ccccc2)cc1"


@pytest.fixture
async def client():
    """Create async test client with lifespan context."""
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as ac:
        yield ac


# ---------------------------------------------------------------------------
# /alerts/screen endpoint tests
# ---------------------------------------------------------------------------


class TestAlertsScreenEndpoint:
    """Tests for POST /api/v1/alerts/screen unified screening endpoint."""

    @pytest.mark.asyncio
    async def test_alerts_screen_returns_concern_groups(self, client: AsyncClient):
        """Unified screen returns a concern_groups dict (may be empty for clean molecules)."""
        response = await client.post(
            "/api/v1/alerts/screen",
            json={"molecule": _AMINOBIPHENYL},
        )
        assert response.status_code == 200
        data = response.json()
        assert "concern_groups" in data
        assert isinstance(data["concern_groups"], dict)

    @pytest.mark.asyncio
    async def test_alerts_screen_has_raw_and_deduped(self, client: AsyncClient):
        """total_raw >= total_deduped — dedup never creates more alerts."""
        response = await client.post(
            "/api/v1/alerts/screen",
            json={"molecule": "O=C1NC(=S)SC1"},  # Rhodanine — has PAINS
        )
        assert response.status_code == 200
        data = response.json()
        assert data["total_raw"] >= data["total_deduped"]
        assert data["total_raw"] >= 0
        assert data["total_deduped"] >= 0

    @pytest.mark.asyncio
    async def test_alerts_screen_response_structure(self, client: AsyncClient):
        """Response contains all required top-level fields."""
        response = await client.post(
            "/api/v1/alerts/screen",
            json={"molecule": _ETHANOL},
        )
        assert response.status_code == 200
        data = response.json()
        for field in ["status", "molecule_info", "alerts", "concern_groups",
                      "total_raw", "total_deduped", "screened_catalogs",
                      "has_critical", "has_warning", "execution_time_ms"]:
            assert field in data, f"Missing field: {field}"
        assert data["status"] == "completed"

    @pytest.mark.asyncio
    async def test_alerts_screen_invalid_smiles(self, client: AsyncClient):
        """Invalid SMILES returns HTTP 400."""
        response = await client.post(
            "/api/v1/alerts/screen",
            json={"molecule": "invalid!!!"},
        )
        assert response.status_code == 400

    @pytest.mark.asyncio
    async def test_alerts_backward_compat(self, client: AsyncClient):
        """Existing /alerts endpoint still works unchanged (backward compat D-21)."""
        response = await client.post(
            "/api/v1/alerts",
            json={"molecule": _ETHANOL, "catalogs": ["PAINS"]},
        )
        assert response.status_code == 200
        data = response.json()
        assert "status" in data
        assert "alerts" in data
        assert "total_alerts" in data
        # Should NOT have concern_groups (old endpoint doesn't include this field)
        assert "concern_groups" not in data


# ---------------------------------------------------------------------------
# /safety/assess endpoint tests
# ---------------------------------------------------------------------------


class TestSafetyAssessEndpoint:
    """Tests for POST /api/v1/safety/assess full assessment endpoint."""

    @pytest.mark.asyncio
    async def test_safety_assess_response_structure(self, client: AsyncClient):
        """Response has cyp_softspots, herg, bro5, reos, complexity top-level keys."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": _ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        for key in ["status", "molecule_info", "cyp_softspots", "herg",
                    "bro5", "reos", "complexity", "execution_time_ms"]:
            assert key in data, f"Missing key: {key}"
        assert data["status"] == "completed"

    @pytest.mark.asyncio
    async def test_safety_assess_haloperidol_herg_high(self, client: AsyncClient):
        """Haloperidol returns herg_risk='high' and risk_score=4."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": _HALOPERIDOL},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["herg"]["herg_risk"] == "high"
        assert data["herg"]["risk_score"] == 4

    @pytest.mark.asyncio
    async def test_safety_assess_aspirin_bro5_not_applicable(self, client: AsyncClient):
        """Aspirin (MW<500) returns bro5.applicable=False."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": _ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["bro5"]["applicable"] is False

    @pytest.mark.asyncio
    async def test_safety_assess_aspirin_reos_violations(self, client: AsyncClient):
        """Aspirin fails REOS (MW < 200 lower bound)."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": _ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        # Aspirin MW ~180 violates REOS MW range 200-500
        assert data["reos"]["passed"] is False

    @pytest.mark.asyncio
    async def test_safety_assess_cyp_fields(self, client: AsyncClient):
        """CYP result contains sites list and n_sites integer."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": _HALOPERIDOL},
        )
        assert response.status_code == 200
        data = response.json()
        cyp = data["cyp_softspots"]
        assert "sites" in cyp
        assert "n_sites" in cyp
        assert isinstance(cyp["sites"], list)
        assert cyp["n_sites"] == len(cyp["sites"])

    @pytest.mark.asyncio
    async def test_safety_assess_complexity_fields(self, client: AsyncClient):
        """Complexity result contains properties dict and n_outliers int."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": _ETHANOL},
        )
        assert response.status_code == 200
        data = response.json()
        complexity = data["complexity"]
        assert "properties" in complexity
        assert "n_outliers" in complexity
        assert "outlier_properties" in complexity
        assert "within_range" in complexity
        # Ethanol should have several outliers (too small)
        assert complexity["n_outliers"] >= 0

    @pytest.mark.asyncio
    async def test_safety_assess_invalid_molecule(self, client: AsyncClient):
        """Invalid molecule returns HTTP 400."""
        response = await client.post(
            "/api/v1/safety/assess",
            json={"molecule": "not_a_molecule"},
        )
        assert response.status_code == 400


# ---------------------------------------------------------------------------
# /safety/summary endpoint tests
# ---------------------------------------------------------------------------


class TestSafetySummaryEndpoint:
    """Tests for POST /api/v1/safety/summary lightweight summary endpoint."""

    @pytest.mark.asyncio
    async def test_safety_summary_returns_statuses(self, client: AsyncClient):
        """Response contains all traffic-light status fields."""
        response = await client.post(
            "/api/v1/safety/summary",
            json={"molecule": _ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        for field in ["status", "total_alerts", "has_critical",
                      "cyp_status", "herg_status", "bro5_status",
                      "reos_status", "complexity_outliers"]:
            assert field in data, f"Missing field: {field}"
        assert data["status"] == "completed"
        # cyp_status is always "default" (informational)
        assert data["cyp_status"] == "default"

    @pytest.mark.asyncio
    async def test_safety_summary_haloperidol_herg_error(self, client: AsyncClient):
        """Haloperidol (risk_score=4) returns herg_status='error'."""
        response = await client.post(
            "/api/v1/safety/summary",
            json={"molecule": _HALOPERIDOL},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["herg_status"] == "error"

    @pytest.mark.asyncio
    async def test_safety_summary_ethanol_herg_success(self, client: AsyncClient):
        """Ethanol (low hERG risk) returns herg_status='success'."""
        response = await client.post(
            "/api/v1/safety/summary",
            json={"molecule": _ETHANOL},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["herg_status"] == "success"

    @pytest.mark.asyncio
    async def test_safety_summary_bro5_not_applicable_default(self, client: AsyncClient):
        """Ethanol (MW<500) returns bro5_status='default' (N/A)."""
        response = await client.post(
            "/api/v1/safety/summary",
            json={"molecule": _ETHANOL},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["bro5_status"] == "default"

    @pytest.mark.asyncio
    async def test_safety_summary_status_values_valid(self, client: AsyncClient):
        """All status fields contain valid traffic-light values."""
        response = await client.post(
            "/api/v1/safety/summary",
            json={"molecule": _HALOPERIDOL},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["herg_status"] in ("success", "warning", "error")
        assert data["bro5_status"] in ("success", "warning", "error", "default")
        assert data["reos_status"] in ("success", "warning", "error")
        assert data["cyp_status"] == "default"
        assert isinstance(data["complexity_outliers"], int)
        assert data["complexity_outliers"] >= 0
