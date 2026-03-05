"""Tests for the /integrations/compare API endpoint."""

from unittest.mock import AsyncMock, patch

import pytest

from app.schemas.integrations import ConsistencyResult


@pytest.mark.asyncio
async def test_compare_with_smiles(client):
    """Test comparison endpoint with SMILES input."""
    mock_result = ConsistencyResult(
        entries=[], comparisons=[], overall_verdict="consistent", summary="All match"
    )
    with patch(
        "app.api.routes.integrations.compare_across_databases",
        new_callable=AsyncMock,
        return_value=mock_result,
    ):
        response = await client.post(
            "/api/v1/integrations/compare",
            json={"smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        )
    assert response.status_code == 200
    assert response.json()["overall_verdict"] == "consistent"


@pytest.mark.asyncio
async def test_compare_no_input(client):
    """Test that missing both smiles and inchikey returns 422."""
    response = await client.post("/api/v1/integrations/compare", json={})
    assert response.status_code == 422
