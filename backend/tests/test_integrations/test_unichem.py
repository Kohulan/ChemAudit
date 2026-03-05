"""Tests for UniChem integration client."""

from unittest.mock import AsyncMock, MagicMock, patch

import httpx
import pytest

from app.services.integrations.unichem import UniChemClient


class TestUniChemClient:
    """Tests for UniChemClient."""

    @pytest.mark.asyncio
    async def test_get_cross_references_success(self):
        """Test successful cross-reference lookup by InChIKey."""
        client = UniChemClient()
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {"src_id": "1", "src_compound_id": "2244"},     # ChEMBL
            {"src_id": "22", "src_compound_id": "2244"},     # PubChem
            {"src_id": "2", "src_compound_id": "DB00945"},   # DrugBank
            {"src_id": "7", "src_compound_id": "15365"},     # ChEBI
        ]
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.get_cross_references("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")

        assert result is not None
        assert result.get("chembl_id") is not None
        assert result.get("drugbank_id") == "DB00945"

    @pytest.mark.asyncio
    async def test_get_cross_references_not_found(self):
        """Test lookup when InChIKey not found."""
        client = UniChemClient()
        mock_response = MagicMock()
        mock_response.json.return_value = []
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.get_cross_references("INVALID-INCHIKEY-X")

        assert result is not None
        assert result.get("chembl_id") is None

    @pytest.mark.asyncio
    async def test_api_failure_returns_empty(self):
        """Test graceful failure."""
        client = UniChemClient()
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = httpx.HTTPStatusError(
            "Not found", request=MagicMock(), response=MagicMock()
        )

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.get_cross_references("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")

        assert result is not None  # Returns empty dict, never None
        assert result.get("chembl_id") is None

    @pytest.mark.asyncio
    async def test_resolve_drugbank_to_inchikey(self):
        """Test resolving DrugBank ID to InChIKey via UniChem."""
        client = UniChemClient()
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {"src_compound_id": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}
        ]
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.resolve_to_inchikey("DB00945", src_id=2)

        assert result == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
