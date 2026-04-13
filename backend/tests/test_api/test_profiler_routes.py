"""
Integration tests for /api/v1/profiler endpoints.

Tests all 5 profiler endpoints:
- POST /api/v1/profiler/full — composite profile (PFI, stars, Abbott, consensus LogP, etc.)
- POST /api/v1/profiler/shape-3d — 3D shape descriptors (lazy per D-26)
- POST /api/v1/profiler/sa-comparison — SA Score / SCScore / SYBA comparison
- POST /api/v1/profiler/efficiency — extended ligand efficiency (LE, LLE, LELP, BEI, SEI)
- POST /api/v1/profiler/mpo — custom MPO scoring

Reference molecule: Aspirin (CC(=O)Oc1ccccc1C(=O)O)
  MW=180.16, heavy atoms=13, LogP≈1.31, TPSA=63.6, HBD=1
  Expected PFI≈2.31 (cLogP + 1 aromatic ring), stars=0, Abbott=85%
  Expected SA Score≈1.58 (easy to synthesize)
  Expected LE≈0.754 (1.4 * pIC50(100nM) / 13 HA = 1.4 * 7.0 / 13)
"""

import pytest

# All tests use the shared async 'client' fixture from conftest.py

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


class TestFullProfile:
    """Tests for POST /api/v1/profiler/full"""

    @pytest.mark.asyncio
    async def test_full_profile_aspirin(self, client):
        """Full profile for aspirin should return all core metric keys."""
        response = await client.post(
            "/api/v1/profiler/full",
            json={"smiles": ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()

        # All required top-level keys must be present
        assert "pfi" in data
        assert "stars" in data
        assert "abbott" in data
        assert "consensus_logp" in data
        assert "skin_permeation" in data
        assert "sa_comparison" in data
        assert "cns_mpo" in data

        # PFI verification: cLogP≈1.31, 1 aromatic ring → PFI≈2.31
        pfi = data["pfi"]
        assert "pfi" in pfi
        assert "risk" in pfi
        assert abs(pfi["pfi"] - 2.31) < 0.1, f"Expected PFI≈2.31, got {pfi['pfi']}"
        assert pfi["risk"] == "low"

        # Stars: aspirin is well within 95th-percentile → 0 violations
        stars = data["stars"]
        assert "stars" in stars
        assert "details" in stars
        assert stars["stars"] == 0
        assert len(stars["details"]) == 8

        # Abbott: aspirin has 0 Lipinski violations, charge=0, TPSA≈63.6 → 85%
        abbott = data["abbott"]
        assert "probability_pct" in abbott
        assert abbott["probability_pct"] == 85

    @pytest.mark.asyncio
    async def test_full_profile_no_3d(self, client):
        """Full profile must NOT include shape_3d key (lazy per D-26)."""
        response = await client.post(
            "/api/v1/profiler/full",
            json={"smiles": ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        assert "shape_3d" not in data, "shape_3d should NOT be in /full response (D-26)"

    @pytest.mark.asyncio
    async def test_full_profile_invalid_smiles(self, client):
        """Invalid SMILES should return 400."""
        response = await client.post(
            "/api/v1/profiler/full",
            json={"smiles": "invalid_smiles_xyz"},
        )
        assert response.status_code == 400
        data = response.json()
        assert "detail" in data
        assert "error" in data["detail"]

    @pytest.mark.asyncio
    async def test_full_profile_cns_mpo_structure(self, client):
        """CNS MPO result should have score, max_score, and components."""
        response = await client.post(
            "/api/v1/profiler/full",
            json={"smiles": ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        cns_mpo = data["cns_mpo"]
        assert "score" in cns_mpo
        assert "max_score" in cns_mpo
        assert cns_mpo["max_score"] == 4.0
        assert "components" in cns_mpo
        assert 0 <= cns_mpo["score"] <= 4.0

    @pytest.mark.asyncio
    async def test_full_profile_consensus_logp_approx_flag(self, client):
        """Consensus LogP must always disclose xlogp3_is_approximation=True."""
        response = await client.post(
            "/api/v1/profiler/full",
            json={"smiles": ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()
        cl = data["consensus_logp"]
        assert cl["xlogp3_is_approximation"] is True


class TestShape3D:
    """Tests for POST /api/v1/profiler/shape-3d"""

    @pytest.mark.asyncio
    async def test_shape3d_aspirin(self, client):
        """3D shape descriptors for aspirin should succeed with valid values."""
        response = await client.post(
            "/api/v1/profiler/shape-3d",
            json={"smiles": ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()

        # Conformer generation for aspirin should succeed
        assert data["3d_conformer_failed"] is False

        # All shape descriptor keys must be present
        assert "npr1" in data
        assert "npr2" in data
        assert "pbf" in data
        assert "pmi1" in data
        assert "pmi2" in data
        assert "pmi3" in data
        assert "shape_class" in data

        # NPR values must be in [0, 1]
        assert 0 <= data["npr1"] <= 1
        assert 0 <= data["npr2"] <= 1

        # Shape class must be one of the 3 categories
        assert data["shape_class"] in ("sphere", "disc", "rod")

    @pytest.mark.asyncio
    async def test_shape3d_invalid_smiles(self, client):
        """Invalid SMILES should return 400."""
        response = await client.post(
            "/api/v1/profiler/shape-3d",
            json={"smiles": "not_valid"},
        )
        assert response.status_code == 400


class TestSAComparison:
    """Tests for POST /api/v1/profiler/sa-comparison"""

    @pytest.mark.asyncio
    async def test_sa_comparison_aspirin(self, client):
        """SA comparison for aspirin: SA Score must be available and easy."""
        response = await client.post(
            "/api/v1/profiler/sa-comparison",
            json={"smiles": ASPIRIN},
        )
        assert response.status_code == 200
        data = response.json()

        # SA Score is always available
        assert "sa_score" in data
        sa = data["sa_score"]
        assert sa["available"] is True
        assert "score" in sa
        # Aspirin SA Score ≈ 1.58 (very easy to synthesize)
        assert abs(sa["score"] - 1.58) < 0.15, f"Expected SA Score≈1.58, got {sa['score']}"
        assert sa["classification"] == "easy"

        # SCScore and SYBA keys must be present (may not be available)
        assert "scscore" in data
        assert "syba" in data

        # All must have "available" key
        assert "available" in data["scscore"]
        assert "available" in data["syba"]

    @pytest.mark.asyncio
    async def test_sa_comparison_invalid_smiles(self, client):
        """Invalid SMILES should return 400."""
        response = await client.post(
            "/api/v1/profiler/sa-comparison",
            json={"smiles": "INVALID@@"},
        )
        assert response.status_code == 400


class TestLigandEfficiency:
    """Tests for POST /api/v1/profiler/efficiency"""

    @pytest.mark.asyncio
    async def test_efficiency_aspirin(self, client):
        """Ligand efficiency for aspirin with IC50=100nM: LE≈0.754."""
        response = await client.post(
            "/api/v1/profiler/efficiency",
            json={
                "smiles": ASPIRIN,
                "activity_value": 100,
                "activity_type": "IC50_nM",
            },
        )
        assert response.status_code == 200
        data = response.json()

        # All expected keys
        assert "pIC50" in data
        assert "LE" in data
        assert "LLE" in data
        assert "LELP" in data
        assert "BEI" in data
        assert "SEI" in data

        # pIC50 = -log10(100e-9) = 7.0
        assert abs(data["pIC50"] - 7.0) < 0.01

        # LE = 1.4 * pIC50 / HA = 1.4 * 7.0 / 13 = 0.754
        assert abs(data["LE"] - 0.754) < 0.01, f"Expected LE≈0.754, got {data['LE']}"

    @pytest.mark.asyncio
    async def test_efficiency_invalid_type(self, client):
        """Invalid activity_type should return 422 (Pydantic validation)."""
        response = await client.post(
            "/api/v1/profiler/efficiency",
            json={
                "smiles": ASPIRIN,
                "activity_value": 100,
                "activity_type": "INVALID_TYPE",
            },
        )
        assert response.status_code == 422

    @pytest.mark.asyncio
    async def test_efficiency_invalid_smiles(self, client):
        """Invalid SMILES should return 400."""
        response = await client.post(
            "/api/v1/profiler/efficiency",
            json={
                "smiles": "not_a_molecule",
                "activity_value": 100,
                "activity_type": "IC50_nM",
            },
        )
        assert response.status_code == 400

    @pytest.mark.asyncio
    async def test_efficiency_pki_type(self, client):
        """pKd activity type should be accepted and treated like pIC50."""
        response = await client.post(
            "/api/v1/profiler/efficiency",
            json={
                "smiles": ASPIRIN,
                "activity_value": 7.0,
                "activity_type": "pKd",
            },
        )
        assert response.status_code == 200
        data = response.json()
        assert abs(data["pIC50"] - 7.0) < 0.01


class TestCustomMPO:
    """Tests for POST /api/v1/profiler/mpo"""

    @pytest.mark.asyncio
    async def test_mpo_custom(self, client):
        """Custom MPO with MW ramp profile should return a numeric score."""
        response = await client.post(
            "/api/v1/profiler/mpo",
            json={
                "smiles": ASPIRIN,
                "profile": [
                    {
                        "property": "MW",
                        "low": 100,
                        "high": 500,
                        "weight": 1.0,
                        "shape": "ramp",
                    }
                ],
            },
        )
        assert response.status_code == 200
        data = response.json()

        assert "score" in data
        assert "max_score" in data
        assert "normalized" in data
        assert "components" in data
        assert isinstance(data["score"], (int, float))
        assert 0 <= data["normalized"] <= 1
        assert len(data["components"]) == 1
        assert data["components"][0]["property"] == "MW"

    @pytest.mark.asyncio
    async def test_mpo_multi_property(self, client):
        """Custom MPO with multiple properties should sum contributions."""
        response = await client.post(
            "/api/v1/profiler/mpo",
            json={
                "smiles": ASPIRIN,
                "profile": [
                    {"property": "MW", "low": 100, "high": 500, "weight": 1.0, "shape": "ramp"},
                    {"property": "LogP", "low": -1, "high": 5, "weight": 1.0, "shape": "sigmoid"},
                    {"property": "TPSA", "low": 20, "high": 130, "weight": 1.0, "shape": "ramp"},
                ],
            },
        )
        assert response.status_code == 200
        data = response.json()
        assert len(data["components"]) == 3
        assert data["max_score"] == 3.0

    @pytest.mark.asyncio
    async def test_mpo_invalid_smiles(self, client):
        """Invalid SMILES should return 400."""
        response = await client.post(
            "/api/v1/profiler/mpo",
            json={
                "smiles": "BAD_SMILES",
                "profile": [
                    {"property": "MW", "low": 100, "high": 500, "weight": 1.0, "shape": "ramp"}
                ],
            },
        )
        assert response.status_code == 400
