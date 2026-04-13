"""
Tests for SureChEMBL patent presence lookup integration.

Tests cover:
- SureChEMBLClient UniChem resolution and SureChEMBL API enrichment
- Fallback behavior when SureChEMBL API is unavailable
- lookup_surechembl convenience function (SMILES/InChIKey input)
- FastAPI endpoint integration
"""

from unittest.mock import AsyncMock, MagicMock, patch

import httpx
import pytest

from app.schemas.integrations import SureChEMBLRequest, SureChEMBLResult
from app.services.integrations.surechembl import SureChEMBLClient, lookup_surechembl


# =============================================================================
# SureChEMBLClient Tests
# =============================================================================


class TestSureChEMBLClientUniChem:
    """Tests for _get_schembl_id_via_unichem method."""

    @pytest.mark.asyncio
    async def test_get_schembl_id_via_unichem_found(self):
        """UniChem returns a SCHEMBL ID for src_id=15."""
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {"src_id": "1", "src_compound_id": "CHEMBL25"},
            {"src_id": "15", "src_compound_id": "SCHEMBL12345"},
        ]
        mock_response.raise_for_status = MagicMock()

        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            with patch("app.services.integrations.surechembl.httpx.AsyncClient") as mock_client:
                mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                    return_value=mock_response
                )

                client = SureChEMBLClient()
                result = await client._get_schembl_id_via_unichem(
                    "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
                )

                assert result == "SCHEMBL12345"

    @pytest.mark.asyncio
    async def test_get_schembl_id_via_unichem_not_found(self):
        """UniChem returns entries but none with src_id=15."""
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {"src_id": "1", "src_compound_id": "CHEMBL25"},
            {"src_id": "22", "src_compound_id": "702"},
        ]
        mock_response.raise_for_status = MagicMock()

        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            with patch("app.services.integrations.surechembl.httpx.AsyncClient") as mock_client:
                mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                    return_value=mock_response
                )

                client = SureChEMBLClient()
                result = await client._get_schembl_id_via_unichem(
                    "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
                )

                assert result is None

    @pytest.mark.asyncio
    async def test_get_schembl_id_via_unichem_numeric_id(self):
        """UniChem returns a numeric-only ID without SCHEMBL prefix."""
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {"src_id": "15", "src_compound_id": "12345"},
        ]
        mock_response.raise_for_status = MagicMock()

        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            with patch("app.services.integrations.surechembl.httpx.AsyncClient") as mock_client:
                mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                    return_value=mock_response
                )

                client = SureChEMBLClient()
                result = await client._get_schembl_id_via_unichem(
                    "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
                )

                # Returns raw "12345" — prefix normalization happens in lookup_by_inchikey
                assert result == "12345"

    @pytest.mark.asyncio
    async def test_get_schembl_id_via_unichem_error(self):
        """UniChem request raises HTTPError — returns None gracefully."""
        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            with patch("app.services.integrations.surechembl.httpx.AsyncClient") as mock_client:
                mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                    side_effect=httpx.ConnectError("Connection refused")
                )

                client = SureChEMBLClient()
                result = await client._get_schembl_id_via_unichem(
                    "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
                )

                assert result is None


class TestSureChEMBLClientCompoundData:
    """Tests for _get_compound_data method."""

    @pytest.mark.asyncio
    async def test_get_compound_data_success(self):
        """SureChEMBL API returns status OK with compound data."""
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "status": "OK",
            "data": {
                "chemical_id": "SCHEMBL12345",
                "global_frequency": 42,
                "smiles": "CCO",
            },
        }
        mock_response.raise_for_status = MagicMock()

        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            with patch("app.services.integrations.surechembl.httpx.AsyncClient") as mock_client:
                mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                    return_value=mock_response
                )

                client = SureChEMBLClient()
                result = await client._get_compound_data("SCHEMBL12345")

                assert result is not None
                assert result["chemical_id"] == "SCHEMBL12345"
                assert result["global_frequency"] == 42

    @pytest.mark.asyncio
    async def test_get_compound_data_api_unavailable(self):
        """SureChEMBL API raises HTTPError — returns None."""
        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            with patch("app.services.integrations.surechembl.httpx.AsyncClient") as mock_client:
                mock_client.return_value.__aenter__.return_value.get = AsyncMock(
                    side_effect=httpx.ConnectError("Service unavailable")
                )

                client = SureChEMBLClient()
                result = await client._get_compound_data("SCHEMBL12345")

                assert result is None


class TestSureChEMBLClientLookup:
    """Tests for lookup_by_inchikey method."""

    @pytest.mark.asyncio
    async def test_lookup_by_inchikey_full_success(self):
        """Both UniChem and SureChEMBL API succeed — rich result."""
        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            client = SureChEMBLClient()

            # Mock UniChem returning SCHEMBL ID
            client._get_schembl_id_via_unichem = AsyncMock(return_value="SCHEMBL12345")
            # Mock SureChEMBL API returning compound data
            client._get_compound_data = AsyncMock(
                return_value={
                    "chemical_id": "SCHEMBL12345",
                    "global_frequency": 42,
                    "smiles": "CCO",
                    "molecular_weight": 46.07,
                }
            )

            result = await client.lookup_by_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

            assert result["found"] is True
            assert result["schembl_id"] == "SCHEMBL12345"
            assert result["patent_count"] == 42
            assert result["source"] == "surechembl_api"
            assert "SCHEMBL12345" in result["url"]
            assert result["smiles"] == "CCO"
            assert result["molecular_weight"] == 46.07

    @pytest.mark.asyncio
    async def test_lookup_by_inchikey_unichem_only_fallback(self):
        """UniChem succeeds but SureChEMBL API is unavailable — fallback."""
        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            client = SureChEMBLClient()
            client._get_schembl_id_via_unichem = AsyncMock(return_value="SCHEMBL12345")
            client._get_compound_data = AsyncMock(return_value=None)

            result = await client.lookup_by_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

            assert result["found"] is True
            assert result["source"] == "unichem_only"
            assert result["schembl_id"] == "SCHEMBL12345"
            assert "SCHEMBL12345" in result["url"]
            assert "patent_count" not in result

    @pytest.mark.asyncio
    async def test_lookup_by_inchikey_not_found(self):
        """UniChem returns no SCHEMBL entry — compound not in SureChEMBL."""
        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            client = SureChEMBLClient()
            client._get_schembl_id_via_unichem = AsyncMock(return_value=None)

            result = await client.lookup_by_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

            assert result["found"] is False

    @pytest.mark.asyncio
    async def test_lookup_by_inchikey_numeric_id_normalization(self):
        """Numeric-only UniChem ID gets SCHEMBL prefix added."""
        with patch("app.services.integrations.surechembl.settings") as mock_settings:
            mock_settings.EXTERNAL_API_TIMEOUT = 30

            client = SureChEMBLClient()
            client._get_schembl_id_via_unichem = AsyncMock(return_value="12345")
            client._get_compound_data = AsyncMock(return_value=None)

            result = await client.lookup_by_inchikey("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

            assert result["found"] is True
            assert result["schembl_id"] == "SCHEMBL12345"


# =============================================================================
# lookup_surechembl convenience function tests
# =============================================================================


class TestLookupSurechembl:
    """Tests for the module-level lookup_surechembl convenience function."""

    @pytest.mark.asyncio
    async def test_lookup_surechembl_with_inchikey(self):
        """Lookup with InChIKey directly."""
        mock_result = {
            "found": True,
            "schembl_id": "SCHEMBL12345",
            "url": "https://www.surechembl.org/chemical/SCHEMBL12345",
            "source": "unichem_only",
        }

        with patch(
            "app.services.integrations.surechembl.SureChEMBLClient"
        ) as mock_client_class:
            mock_client = AsyncMock()
            mock_client.lookup_by_inchikey = AsyncMock(return_value=mock_result)
            mock_client_class.return_value = mock_client

            request = SureChEMBLRequest(inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
            result = await lookup_surechembl(request)

            assert isinstance(result, SureChEMBLResult)
            assert result.found is True
            assert result.schembl_id == "SCHEMBL12345"
            mock_client.lookup_by_inchikey.assert_called_once_with(
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
            )

    @pytest.mark.asyncio
    async def test_lookup_surechembl_with_smiles(self):
        """Lookup with SMILES — derives InChIKey internally."""
        mock_result = {
            "found": True,
            "schembl_id": "SCHEMBL99999",
            "url": "https://www.surechembl.org/chemical/SCHEMBL99999",
            "source": "unichem_only",
        }

        with patch(
            "app.services.integrations.surechembl.SureChEMBLClient"
        ) as mock_client_class:
            mock_client = AsyncMock()
            mock_client.lookup_by_inchikey = AsyncMock(return_value=mock_result)
            mock_client_class.return_value = mock_client

            request = SureChEMBLRequest(smiles="CCO")
            result = await lookup_surechembl(request)

            assert isinstance(result, SureChEMBLResult)
            assert result.found is True
            # Verify lookup_by_inchikey was called with a valid InChIKey
            call_args = mock_client.lookup_by_inchikey.call_args[0]
            assert len(call_args[0]) == 27  # Standard InChIKey length
            assert "-" in call_args[0]

    @pytest.mark.asyncio
    async def test_lookup_surechembl_invalid_smiles(self):
        """Invalid SMILES returns found=False immediately."""
        request = SureChEMBLRequest(smiles="INVALID_NOT_A_SMILES")
        result = await lookup_surechembl(request)

        assert isinstance(result, SureChEMBLResult)
        assert result.found is False

    @pytest.mark.asyncio
    async def test_lookup_surechembl_no_input(self):
        """No SMILES or InChIKey returns found=False."""
        request = SureChEMBLRequest()
        result = await lookup_surechembl(request)

        assert isinstance(result, SureChEMBLResult)
        assert result.found is False


# =============================================================================
# Endpoint integration tests
# =============================================================================


class TestSureChEMBLEndpoint:
    """Tests for the /integrations/surechembl/lookup endpoint."""

    @pytest.mark.asyncio
    async def test_surechembl_endpoint_returns_schema(self):
        """POST endpoint returns SureChEMBLResult-shaped response."""
        mock_result = SureChEMBLResult(
            found=True,
            schembl_id="SCHEMBL12345",
            url="https://www.surechembl.org/chemical/SCHEMBL12345",
            source="unichem_only",
        )

        with patch(
            "app.api.routes.integrations.lookup_surechembl",
            new_callable=AsyncMock,
            return_value=mock_result,
        ):
            from app.main import app

            async with httpx.AsyncClient(
                transport=httpx.ASGITransport(app=app),
                base_url="http://testserver",
            ) as ac:
                response = await ac.post(
                    "/api/v1/integrations/surechembl/lookup",
                    json={"inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"},
                )

            assert response.status_code == 200
            data = response.json()
            assert "found" in data
            assert data["found"] is True
            assert data["schembl_id"] == "SCHEMBL12345"
            assert data["source"] == "unichem_only"
