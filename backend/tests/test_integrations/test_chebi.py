"""Tests for ChEBI integration client."""

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from app.services.integrations.chebi import ChEBIClient


class TestChEBIClient:
    """Tests for ChEBIClient."""

    @pytest.mark.asyncio
    async def test_get_compound_success(self):
        """Test successful compound lookup."""
        client = ChEBIClient()
        mock_response = MagicMock()
        mock_response.text = (
            "<return><smiles>CC(=O)Oc1ccccc1C(=O)O</smiles>"
            "<inchiKey>BSYNRYMUTXBXSQ-UHFFFAOYSA-N</inchiKey>"
            "<chebiAsciiName>aspirin</chebiAsciiName>"
            "<inchi>InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12"
            "/h2-5H,1H3,(H,11,12)</inchi></return>"
        )
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.get_compound("15365")

        assert result is not None
        assert result["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"

    @pytest.mark.asyncio
    async def test_get_compound_not_found(self):
        """Test lookup when compound not found."""
        client = ChEBIClient()
        mock_response = MagicMock()
        mock_response.text = "<return></return>"
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.get_compound("999999999")

        assert result is None
