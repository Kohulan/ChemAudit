"""
Tests for DECIMER integration client.
"""
import pytest
from unittest.mock import AsyncMock, patch

from app.services.integrations.decimer import DECIMERClient, validate_ocsr_result
from app.schemas.integrations import DECIMERRequest


class TestDECIMERClient:
    """Tests for DECIMERClient."""

    @pytest.mark.asyncio
    async def test_ocsr_from_image_success(self):
        """Test successful OCSR from image."""
        client = DECIMERClient()
        mock_response = AsyncMock()
        mock_response.json.return_value = {"smiles": "CCO"}
        mock_response.raise_for_status = AsyncMock()

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.post = AsyncMock(return_value=mock_response)

            result = await client.ocsr_from_image(b"fake_image_data")

            assert result == "CCO"

    @pytest.mark.asyncio
    async def test_ocsr_from_image_failure(self):
        """Test OCSR failure returns None gracefully."""
        client = DECIMERClient()
        mock_response = AsyncMock()
        mock_response.raise_for_status.side_effect = Exception("API error")

        with patch("httpx.AsyncClient") as mock_client:
            mock_client.return_value.__aenter__.return_value.post = AsyncMock(return_value=mock_response)

            result = await client.ocsr_from_image(b"fake_image_data")

            assert result is None


class TestValidateOCSRResult:
    """Tests for validate_ocsr_result function."""

    def test_validate_valid_smiles_high_confidence(self):
        """Test validation of valid SMILES with high confidence."""
        request = DECIMERRequest(smiles="CCO", confidence=0.95)
        result = validate_ocsr_result(request)

        assert result.is_valid is True
        assert result.confidence == 0.95
        assert "High confidence" in result.validation_message
        assert result.canonical_smiles is not None
        assert result.inchikey is not None

    def test_validate_valid_smiles_moderate_confidence(self):
        """Test validation of valid SMILES with moderate confidence."""
        request = DECIMERRequest(smiles="c1ccccc1", confidence=0.75)
        result = validate_ocsr_result(request)

        assert result.is_valid is True
        assert result.confidence == 0.75
        assert "Moderate confidence" in result.validation_message

    def test_validate_valid_smiles_low_confidence(self):
        """Test validation of valid SMILES with low confidence."""
        request = DECIMERRequest(smiles="CCO", confidence=0.5)
        result = validate_ocsr_result(request)

        assert result.is_valid is True
        assert result.confidence == 0.5
        assert "Low confidence" in result.validation_message
        assert "strongly recommended" in result.validation_message

    def test_validate_valid_smiles_no_confidence(self):
        """Test validation of valid SMILES without confidence score."""
        request = DECIMERRequest(smiles="CCO")
        result = validate_ocsr_result(request)

        assert result.is_valid is True
        assert result.confidence is None
        assert "no confidence score" in result.validation_message

    def test_validate_invalid_smiles(self):
        """Test validation of invalid SMILES."""
        request = DECIMERRequest(smiles="INVALID_SMILES_123")
        result = validate_ocsr_result(request)

        assert result.is_valid is False
        assert "cannot parse" in result.validation_message
        assert result.canonical_smiles is None

    def test_validate_invalid_chemistry(self):
        """Test validation of parseable but invalid chemistry."""
        # SMILES that parses but fails sanitization
        request = DECIMERRequest(smiles="C1CCC", confidence=0.8)
        result = validate_ocsr_result(request)

        assert result.is_valid is False
        assert "Invalid chemistry" in result.validation_message
