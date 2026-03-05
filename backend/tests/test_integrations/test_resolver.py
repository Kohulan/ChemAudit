"""Tests for the Universal Identifier Resolver orchestrator."""

from unittest.mock import AsyncMock, patch

import pytest

from app.services.integrations.resolver import resolve_identifier


class TestResolveIdentifier:
    """Tests for resolve_identifier."""

    @pytest.mark.asyncio
    async def test_resolve_smiles(self):
        """Test resolving a SMILES string."""
        result = await resolve_identifier("CCO")
        assert result.resolved is True
        assert result.identifier_type_detected == "smiles"
        assert result.canonical_smiles is not None
        assert result.inchikey is not None
        assert result.confidence == "high"

    @pytest.mark.asyncio
    async def test_resolve_inchikey(self):
        """Test resolving an InChIKey — needs cross-ref lookup."""
        with patch(
            "app.services.integrations.resolver.UniChemClient"
        ) as mock_unichem_cls:
            mock_unichem = AsyncMock()
            mock_unichem.get_cross_references = AsyncMock(return_value={
                "chembl_id": "CHEMBL545", "drugbank_id": None,
                "pubchem_cid": "702", "chebi_id": None, "kegg_id": None,
            })
            mock_unichem_cls.return_value = mock_unichem

            with patch(
                "app.services.integrations.resolver.PubChemClient"
            ) as mock_pc_cls:
                mock_pc = AsyncMock()
                mock_pc.search_by_inchikey = AsyncMock(return_value=702)
                mock_pc.get_compound_properties = AsyncMock(return_value={
                    "CanonicalSMILES": "CCO",
                    "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    "MolecularFormula": "C2H6O",
                    "MolecularWeight": 46.07,
                    "IUPACName": "ethanol",
                })
                mock_pc.get_synonyms = AsyncMock(return_value=[])
                mock_pc_cls.return_value = mock_pc

                result = await resolve_identifier("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

        assert result.resolved is True
        assert result.identifier_type_detected == "inchikey"
        assert result.canonical_smiles == "CCO"

    @pytest.mark.asyncio
    async def test_resolve_chembl_id(self):
        """Test resolving a ChEMBL ID."""
        with patch(
            "app.services.integrations.resolver.ChEMBLClient"
        ) as mock_chembl_cls:
            mock_chembl = AsyncMock()
            mock_chembl.get_molecule = AsyncMock(return_value={
                "molecule_chembl_id": "CHEMBL25",
                "pref_name": "ASPIRIN",
                "molecule_structures": {
                    "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "standard_inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                    "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                },
                "molecule_properties": {
                    "full_molecular_formula": "C9H8O4",
                    "molecular_weight": "180.16",
                },
            })
            mock_chembl_cls.return_value = mock_chembl

            with patch(
                "app.services.integrations.resolver.UniChemClient"
            ) as mock_unichem_cls:
                mock_unichem = AsyncMock()
                mock_unichem.get_cross_references = AsyncMock(return_value={
                    "chembl_id": "CHEMBL25", "drugbank_id": "DB00945",
                    "pubchem_cid": "2244", "chebi_id": "CHEBI:15365", "kegg_id": None,
                })
                mock_unichem_cls.return_value = mock_unichem

                result = await resolve_identifier("CHEMBL25")

        assert result.resolved is True
        assert result.identifier_type_detected == "chembl_id"
        assert result.canonical_smiles is not None
        assert result.cross_references.drugbank_id == "DB00945"

    @pytest.mark.asyncio
    async def test_resolve_invalid_returns_not_resolved(self):
        """Test that unresolvable input returns resolved=False."""
        with patch(
            "app.services.integrations.resolver._resolve_name",
            new_callable=AsyncMock,
            return_value=None,
        ):
            result = await resolve_identifier("xyzzy_not_a_molecule_12345")
        assert result.resolved is False

    @pytest.mark.asyncio
    async def test_resolve_cas_number(self):
        """Test resolving a CAS number via PubChem name search."""
        with patch(
            "app.services.integrations.resolver.PubChemClient"
        ) as mock_pc_cls:
            mock_pc = AsyncMock()
            mock_pc.search_by_name = AsyncMock(return_value=2244)
            mock_pc.get_compound_properties = AsyncMock(return_value={
                "CanonicalSMILES": "CC(=O)Oc1ccccc1C(=O)O",
                "InChI": "InChI=1S/C9H8O4",
                "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                "MolecularFormula": "C9H8O4",
                "MolecularWeight": 180.16,
                "IUPACName": "aspirin",
            })
            mock_pc.get_synonyms = AsyncMock(return_value=["aspirin"])
            mock_pc_cls.return_value = mock_pc

            with patch(
                "app.services.integrations.resolver.UniChemClient"
            ) as mock_unichem_cls:
                mock_unichem = AsyncMock()
                mock_unichem.get_cross_references = AsyncMock(return_value={
                    "chembl_id": None, "drugbank_id": None,
                    "pubchem_cid": None, "chebi_id": None, "kegg_id": None,
                })
                mock_unichem_cls.return_value = mock_unichem

                result = await resolve_identifier("50-78-2")

        assert result.resolved is True
        assert result.identifier_type_detected == "cas"
