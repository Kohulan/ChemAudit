"""
Tests for COCONUT integration client.
"""
import pytest
from unittest.mock import AsyncMock, patch

from app.services.integrations.coconut import COCONUTClient, lookup_natural_product
from app.schemas.integrations import COCONUTRequest


class TestCOCONUTClient:
    """Tests for COCONUTClient."""

    @pytest.mark.asyncio
    async def test_search_by_smiles_success(self):
        """Test successful search by SMILES."""
        client = COCONUTClient()
        mock_response = AsyncMock()
        mock_response.json.return_value = [
            {
                "coconut_id": "CNP0123456",
                "name": "Caffeine",
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            }
        ]
        mock_response.raise_for_status = AsyncMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(return_value=mock_response)

            result = await client.search_by_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

            assert result is not None
            assert result["coconut_id"] == "CNP0123456"

    @pytest.mark.asyncio
    async def test_search_by_smiles_not_found(self):
        """Test search by SMILES not found."""
        client = COCONUTClient()
        mock_response = AsyncMock()
        mock_response.json.return_value = []
        mock_response.raise_for_status = AsyncMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(return_value=mock_response)

            result = await client.search_by_smiles("CCO")

            assert result is None

    @pytest.mark.asyncio
    async def test_search_by_inchikey_success(self):
        """Test successful search by InChIKey."""
        client = COCONUTClient()
        mock_response = AsyncMock()
        mock_response.json.return_value = {
            "coconut_id": "CNP0123456",
            "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        }
        mock_response.raise_for_status = AsyncMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(return_value=mock_response)

            result = await client.search_by_inchikey("RYYVLZVUVIJVGH-UHFFFAOYSA-N")

            assert result is not None
            assert result["coconut_id"] == "CNP0123456"

    @pytest.mark.asyncio
    async def test_search_api_failure(self):
        """Test API failure returns None gracefully."""
        client = COCONUTClient()
        mock_response = AsyncMock()
        mock_response.raise_for_status.side_effect = Exception("API error")

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(return_value=mock_response)

            result = await client.search_by_smiles("CCO")

            assert result is None


class TestLookupNaturalProduct:
    """Tests for lookup_natural_product function."""

    @pytest.mark.asyncio
    async def test_lookup_by_inchikey(self):
        """Test lookup by InChIKey."""
        mock_data = {
            "coconut_id": "CNP0123456",
            "name": "Caffeine",
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
            "molecular_weight": 194.19,
        }

        with patch("app.services.integrations.coconut.COCONUTClient") as mock_client_class:
            mock_client = AsyncMock()
            mock_client.search_by_inchikey = AsyncMock(return_value=mock_data)
            mock_client_class.return_value = mock_client

            request = COCONUTRequest(inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N")
            result = await lookup_natural_product(request)

            assert result.found is True
            assert result.coconut_id == "CNP0123456"
            assert result.name == "Caffeine"
            assert result.url is not None

    @pytest.mark.asyncio
    async def test_lookup_by_smiles(self):
        """Test lookup by SMILES."""
        mock_data = {
            "coconut_id": "CNP0123456",
            "name": "Ethanol",
            "smiles": "CCO",
        }

        with patch("app.services.integrations.coconut.COCONUTClient") as mock_client_class:
            mock_client = AsyncMock()
            mock_client.search_by_inchikey = AsyncMock(return_value=mock_data)
            mock_client_class.return_value = mock_client

            request = COCONUTRequest(smiles="CCO")
            result = await lookup_natural_product(request)

            assert result.found is True

    @pytest.mark.asyncio
    async def test_lookup_not_found(self):
        """Test lookup not found."""
        with patch("app.services.integrations.coconut.COCONUTClient") as mock_client_class:
            mock_client = AsyncMock()
            mock_client.search_by_inchikey = AsyncMock(return_value=None)
            mock_client.search_by_smiles = AsyncMock(return_value=None)
            mock_client_class.return_value = mock_client

            request = COCONUTRequest(smiles="CCO")
            result = await lookup_natural_product(request)

            assert result.found is False
            assert result.coconut_id is None
