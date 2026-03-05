"""Tests for the cross-database comparator service."""

from unittest.mock import AsyncMock, patch

import pytest

from app.schemas.integrations import (
    ChEMBLResult,
    COCONUTResult,
    PubChemResult,
)
from app.services.integrations.comparator import compare_across_databases


class TestCompareAcrossDatabases:
    """Tests for compare_across_databases."""

    @pytest.mark.asyncio
    async def test_consistent_results(self):
        """Test that identical SMILES across databases yields 'consistent'."""
        with patch(
            "app.services.integrations.comparator.get_compound_info",
            new_callable=AsyncMock,
        ) as mock_pc, patch(
            "app.services.integrations.comparator.get_bioactivity",
            new_callable=AsyncMock,
        ) as mock_chembl, patch(
            "app.services.integrations.comparator.lookup_natural_product",
            new_callable=AsyncMock,
        ) as mock_coconut, patch(
            "app.services.integrations.comparator._lookup_wikidata",
            new_callable=AsyncMock,
        ) as mock_wikidata:
            mock_pc.return_value = PubChemResult(
                found=True, cid=2244,
                canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
                inchi="InChI=1S/C9H8O4", inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                molecular_formula="C9H8O4", molecular_weight=180.16,
                iupac_name="aspirin",
            )
            mock_chembl.return_value = ChEMBLResult(
                found=True, chembl_id="CHEMBL25",
                molecular_formula="C9H8O4", molecular_weight=180.16,
            )
            mock_coconut.return_value = COCONUTResult(found=False)
            mock_wikidata.return_value = None

            result = await compare_across_databases(
                smiles="CC(=O)Oc1ccccc1C(=O)O"
            )

        assert result.overall_verdict in ("consistent", "minor_differences")
        assert len(result.entries) >= 2

    @pytest.mark.asyncio
    async def test_no_databases_found(self):
        """Test when no database has the compound."""
        with patch(
            "app.services.integrations.comparator.get_compound_info",
            new_callable=AsyncMock,
        ) as mock_pc, patch(
            "app.services.integrations.comparator.get_bioactivity",
            new_callable=AsyncMock,
        ) as mock_chembl, patch(
            "app.services.integrations.comparator.lookup_natural_product",
            new_callable=AsyncMock,
        ) as mock_coconut, patch(
            "app.services.integrations.comparator._lookup_wikidata",
            new_callable=AsyncMock,
        ) as mock_wikidata:
            mock_pc.return_value = PubChemResult(found=False)
            mock_chembl.return_value = ChEMBLResult(found=False)
            mock_coconut.return_value = COCONUTResult(found=False)
            mock_wikidata.return_value = None

            result = await compare_across_databases(smiles="C#N")

        assert result.overall_verdict == "no_data"
        # Resolved entry is present but external databases are not found
        resolved = [e for e in result.entries if e.database == "Resolved"]
        assert len(resolved) == 1

    @pytest.mark.asyncio
    async def test_consistent_with_wikidata(self):
        """Test consistency when all 4 databases return matching data."""
        with patch(
            "app.services.integrations.comparator.get_compound_info",
            new_callable=AsyncMock,
        ) as mock_pc, patch(
            "app.services.integrations.comparator.get_bioactivity",
            new_callable=AsyncMock,
        ) as mock_chembl, patch(
            "app.services.integrations.comparator.lookup_natural_product",
            new_callable=AsyncMock,
        ) as mock_coconut, patch(
            "app.services.integrations.comparator._lookup_wikidata",
            new_callable=AsyncMock,
        ) as mock_wikidata:
            mock_pc.return_value = PubChemResult(
                found=True, cid=2244,
                canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
                inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                molecular_formula="C9H8O4", molecular_weight=180.16,
            )
            mock_chembl.return_value = ChEMBLResult(
                found=True, chembl_id="CHEMBL25",
                canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
                inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                molecular_formula="C9H8O4", molecular_weight=180.16,
            )
            mock_coconut.return_value = COCONUTResult(
                found=True,
                smiles="CC(=O)Oc1ccccc1C(=O)O",
                inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                molecular_formula="C9H8O4", molecular_weight=180.16,
            )
            mock_wikidata.return_value = {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "formula": "C9H8O4",
                "mass": 180.16,
                "label": "aspirin",
            }

            result = await compare_across_databases(
                smiles="CC(=O)Oc1ccccc1C(=O)O"
            )

        assert result.overall_verdict == "consistent"
        # 5 entries: Resolved + 4 external databases
        assert len(result.entries) == 5
        found = [e for e in result.entries if e.found]
        assert len(found) == 5
        db_names = {e.database for e in found}
        assert "Wikidata" in db_names
        assert "Resolved" in db_names

    @pytest.mark.asyncio
    async def test_no_input(self):
        """Test with no input returns no_data."""
        result = await compare_across_databases()
        assert result.overall_verdict == "no_data"
