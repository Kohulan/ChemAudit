"""
Test health check endpoint
"""
import pytest


@pytest.mark.asyncio
async def test_health_check(client):
    """Test GET /health returns 200 and correct status"""
    response = await client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "healthy"


@pytest.mark.asyncio
async def test_root_endpoint(client):
    """Test GET / returns API info"""
    response = await client.get("/")
    assert response.status_code == 200
    data = response.json()
    assert "name" in data
    assert "version" in data
    assert data["name"] == "ChemStructVal API"
