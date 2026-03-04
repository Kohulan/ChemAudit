"""Tests for the /resolve API endpoint."""

from unittest.mock import AsyncMock, patch

import pytest

from app.schemas.integrations import CrossReferences, ResolvedCompound


@pytest.mark.asyncio
async def test_resolve_smiles(client):
    """Test resolving a SMILES string via API."""
    response = await client.post("/api/v1/resolve", json={"identifier": "CCO"})
    assert response.status_code == 200
    data = response.json()
    assert data["resolved"] is True
    assert data["identifier_type_detected"] == "smiles"
    assert data["canonical_smiles"] is not None


@pytest.mark.asyncio
async def test_resolve_with_explicit_type(client):
    """Test resolving with explicit identifier type."""
    response = await client.post(
        "/api/v1/resolve",
        json={"identifier": "CCO", "identifier_type": "smiles"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["resolved"] is True


@pytest.mark.asyncio
async def test_resolve_empty_identifier(client):
    """Test that empty identifier returns 422."""
    response = await client.post("/api/v1/resolve", json={"identifier": ""})
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_resolve_chembl_id(client):
    """Test resolving ChEMBL ID via API."""
    mock_result = ResolvedCompound(
        resolved=True,
        identifier_type_detected="chembl_id",
        canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        resolution_source="chembl",
        resolution_chain=["CHEMBL25 → ChEMBL API → SMILES"],
        cross_references=CrossReferences(),
        confidence="high",
    )

    with patch(
        "app.api.routes.resolve.resolve_identifier",
        new_callable=AsyncMock,
        return_value=mock_result,
    ):
        response = await client.post(
            "/api/v1/resolve", json={"identifier": "CHEMBL25"}
        )
    assert response.status_code == 200
    assert response.json()["canonical_smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
