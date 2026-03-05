"""Tests for identifier type detection."""


from app.services.integrations.identifier_detect import detect_identifier_type


class TestDetectIdentifierType:
    """Tests for detect_identifier_type."""

    def test_pubchem_cid_numeric(self):
        assert detect_identifier_type("2244") == "pubchem_cid"

    def test_pubchem_cid_prefixed(self):
        assert detect_identifier_type("CID:2244") == "pubchem_cid"

    def test_chembl_id(self):
        assert detect_identifier_type("CHEMBL25") == "chembl_id"

    def test_chembl_id_lowercase(self):
        assert detect_identifier_type("chembl25") == "chembl_id"

    def test_cas_number(self):
        assert detect_identifier_type("50-78-2") == "cas"

    def test_cas_number_long(self):
        assert detect_identifier_type("9002-93-1") == "cas"

    def test_drugbank_id(self):
        assert detect_identifier_type("DB00945") == "drugbank_id"

    def test_chebi_id_prefixed(self):
        assert detect_identifier_type("CHEBI:15365") == "chebi_id"

    def test_unii(self):
        assert detect_identifier_type("R16CO5Y76E") == "unii"

    def test_inchikey(self):
        assert detect_identifier_type("BSYNRYMUTXBXSQ-UHFFFAOYSA-N") == "inchikey"

    def test_inchi_string(self):
        assert detect_identifier_type(
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        ) == "inchi"

    def test_wikipedia_url(self):
        assert detect_identifier_type("https://en.wikipedia.org/wiki/Aspirin") == "wikipedia"

    def test_smiles_with_rings(self):
        assert detect_identifier_type("CC(=O)Oc1ccccc1C(=O)O") == "smiles"

    def test_smiles_simple(self):
        assert detect_identifier_type("CCO") == "smiles"

    def test_compound_name(self):
        assert detect_identifier_type("aspirin") == "name"

    def test_compound_name_with_spaces(self):
        assert detect_identifier_type("acetylsalicylic acid") == "name"

    def test_explicit_type_override(self):
        assert detect_identifier_type("2244", identifier_type="pubchem_cid") == "pubchem_cid"

    def test_explicit_type_name(self):
        assert detect_identifier_type("CCO", identifier_type="name") == "name"
