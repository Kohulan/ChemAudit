"""Tests for Wikidata integration client."""

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from app.services.integrations.wikidata import WikidataClient


class TestWikidataClient:
    """Tests for WikidataClient."""

    @pytest.mark.asyncio
    async def test_resolve_wikipedia_url(self):
        """Test resolving Wikipedia URL to chemical data."""
        client = WikidataClient()
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": {
                "bindings": [
                    {
                        "smiles": {"value": "CC(=O)Oc1ccccc1C(=O)O"},
                        "inchikey": {"value": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"},
                        "cas": {"value": "50-78-2"},
                        "label": {"value": "aspirin"},
                    }
                ]
            }
        }
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.resolve_from_wikipedia(
                "https://en.wikipedia.org/wiki/Aspirin"
            )

        assert result is not None
        assert result["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
        assert result["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    @pytest.mark.asyncio
    async def test_resolve_not_found(self):
        """Test resolution when article has no chemical data."""
        client = WikidataClient()
        mock_response = MagicMock()
        mock_response.json.return_value = {"results": {"bindings": []}}
        mock_response.raise_for_status = MagicMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                return_value=mock_response
            )
            result = await client.resolve_from_wikipedia(
                "https://en.wikipedia.org/wiki/Dog"
            )

        assert result is None

    @pytest.mark.asyncio
    async def test_extract_article_title(self):
        """Test extraction of article title from URL."""
        client = WikidataClient()
        assert client._extract_title("https://en.wikipedia.org/wiki/Aspirin") == "Aspirin"
        assert (
            client._extract_title("https://en.wikipedia.org/wiki/Acetylsalicylic_acid")
            == "Acetylsalicylic acid"
        )
