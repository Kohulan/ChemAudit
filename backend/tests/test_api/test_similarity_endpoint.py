"""
Tests for POST /api/v1/validate/similarity endpoint.

Tests ECFP4 Tanimoto similarity computation between two molecules.
"""

import pytest


class TestSimilarityEndpoint:
    """Test POST /api/v1/validate/similarity endpoint"""

    @pytest.mark.asyncio
    async def test_similar_molecules(self, client):
        """Should return high similarity for structurally similar molecules"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "CCO", "smiles_b": "CCCO"},
        )
        assert response.status_code == 200
        data = response.json()

        assert "tanimoto_similarity" in data
        assert 0.0 <= data["tanimoto_similarity"] <= 1.0
        assert data["fingerprint_type"] == "ECFP4"
        assert data["radius"] == 2
        assert data["n_bits"] == 2048
        assert isinstance(data["common_bits"], int)
        assert isinstance(data["bits_a"], int)
        assert isinstance(data["bits_b"], int)
        assert data["common_bits"] >= 0
        assert data["bits_a"] > 0
        assert data["bits_b"] > 0

    @pytest.mark.asyncio
    async def test_identical_molecules(self, client):
        """Should return 1.0 similarity for identical molecules"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "c1ccccc1", "smiles_b": "C1=CC=CC=C1"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["tanimoto_similarity"] == 1.0
        assert data["common_bits"] == data["bits_a"] == data["bits_b"]

    @pytest.mark.asyncio
    async def test_dissimilar_molecules(self, client):
        """Should return low similarity for structurally different molecules"""
        # Benzene vs a long aliphatic chain
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "c1ccccc1", "smiles_b": "CCCCCCCCCCCCCCCC"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["tanimoto_similarity"] < 0.5

    @pytest.mark.asyncio
    async def test_invalid_smiles_a(self, client):
        """Should return 400 for invalid smiles_a"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "not_a_smiles", "smiles_b": "CCO"},
        )
        assert response.status_code == 400
        data = response.json()
        assert "smiles_a" in data["detail"]

    @pytest.mark.asyncio
    async def test_invalid_smiles_b(self, client):
        """Should return 400 for invalid smiles_b"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "CCO", "smiles_b": "not_a_smiles"},
        )
        assert response.status_code == 400
        data = response.json()
        assert "smiles_b" in data["detail"]

    @pytest.mark.asyncio
    async def test_empty_smiles_rejected(self, client):
        """Should return 422 for empty SMILES strings"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "", "smiles_b": "CCO"},
        )
        assert response.status_code == 422

    @pytest.mark.asyncio
    async def test_missing_field(self, client):
        """Should return 422 when a required field is missing"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "CCO"},
        )
        assert response.status_code == 422

    @pytest.mark.asyncio
    async def test_common_bits_consistency(self, client):
        """Common bits should not exceed either fingerprint's on-bit count"""
        response = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "c1ccc(O)cc1", "smiles_b": "c1ccc(N)cc1"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["common_bits"] <= data["bits_a"]
        assert data["common_bits"] <= data["bits_b"]

    @pytest.mark.asyncio
    async def test_symmetry(self, client):
        """Tanimoto(A, B) should equal Tanimoto(B, A)"""
        resp1 = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "c1ccccc1O", "smiles_b": "CCO"},
        )
        resp2 = await client.post(
            "/api/v1/validate/similarity",
            json={"smiles_a": "CCO", "smiles_b": "c1ccccc1O"},
        )
        assert resp1.status_code == 200
        assert resp2.status_code == 200
        assert resp1.json()["tanimoto_similarity"] == resp2.json()["tanimoto_similarity"]
