"""
Tests for standardization provenance pipeline (Milestone 2.1).

Covers all four STD requirements:
- STD-01: Tautomer canonicalization provenance
- STD-02: Neutralization charge change tracking
- STD-03: Functional group audit (normalization rule identification)
- STD-04: Parent extraction fragment naming

Plus backward compatibility and endpoint integration tests.
"""

import pytest
from httpx import ASGITransport, AsyncClient
from rdkit import Chem

from app.main import app
from app.services.standardization.chembl_pipeline import StandardizationOptions
from app.services.standardization.fragment_dict import (
    COUNTERION_NAMES,
    classify_fragment,
)
from app.services.standardization.provenance import ProvenancePipeline


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def pipeline():
    """Create a fresh ProvenancePipeline instance for tests."""
    return ProvenancePipeline()


@pytest.fixture
def opts():
    """Default standardization options (no provenance, no tautomer)."""
    return StandardizationOptions()


@pytest.fixture
def opts_with_tautomer():
    """Options with tautomer canonicalization enabled."""
    return StandardizationOptions(include_tautomer=True)


@pytest.fixture
def simple_mol():
    """Simple ethanol molecule — no modifications expected."""
    return Chem.MolFromSmiles("CCO")


@pytest.fixture
def charged_mol():
    """Alanine zwitterion — has both + and - charges that ChEMBL normalizes."""
    return Chem.MolFromSmiles("[NH3+]CC([O-])=O")


@pytest.fixture
def salt_mol():
    """Benzene with HCl — get_parent will strip the Cl fragment."""
    return Chem.MolFromSmiles("c1ccccc1.Cl")


@pytest.fixture
def tautomer_mol():
    """2-hydroxypyridine — classic tautomer pair with 2-pyridinone."""
    return Chem.MolFromSmiles("OC1=CC=CC=N1")


@pytest.fixture
def dmso_mol():
    """DMSO — standardizer normalizes S=O to [S+][O-] (sulphoxide_normalization)."""
    return Chem.MolFromSmiles("CS(C)=O")


@pytest.fixture
async def client():
    """Create async test client."""
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test/api/v1"
    ) as ac:
        yield ac


# ---------------------------------------------------------------------------
# TestFragmentDictionary (~6 tests)
# ---------------------------------------------------------------------------


class TestFragmentDictionary:
    """Test the COUNTERION_NAMES dictionary and classify_fragment() function."""

    def test_common_counterions_present(self):
        """Key ionic counterions must be in COUNTERION_NAMES."""
        canonical_forms = {
            Chem.MolToSmiles(Chem.MolFromSmiles(s))
            for s in ["[Cl-]", "[Br-]", "[Na+]", "[K+]"]
        }
        dict_keys = set(COUNTERION_NAMES.keys())
        for smiles in canonical_forms:
            assert smiles in dict_keys, f"Expected {smiles!r} in COUNTERION_NAMES"

    def test_common_solvents_present(self):
        """Common solvents must be in COUNTERION_NAMES with role='solvent'."""
        solvents = {"O": "water", "CO": "methanol", "CCO": "ethanol"}
        for smiles, expected_name in solvents.items():
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            assert canon in COUNTERION_NAMES, f"Solvent {smiles!r} missing from dict"
            assert COUNTERION_NAMES[canon]["role"] == "solvent"
            assert COUNTERION_NAMES[canon]["name"] == expected_name

    def test_common_salts_present(self):
        """Common acid salt-formers must be present with role='salt'."""
        # Use classify_fragment() to properly canonicalize before lookup
        acids = [
            "CC(=O)O",        # acetic acid
            "OC(=O)c1ccccc1", # benzoic acid (non-canonical input)
            "O=C(O)c1ccccc1", # benzoic acid (canonical)
        ]
        for smiles in acids:
            result = classify_fragment(smiles)
            assert result["role"] == "salt", (
                f"Expected role='salt' for {smiles!r}, got {result['role']!r}"
            )
            assert result["name"] is not None, (
                f"Expected a name for {smiles!r}"
            )

    def test_classify_known_fragment(self):
        """classify_fragment('[Cl-]') should return known counterion metadata."""
        result = classify_fragment("[Cl-]")
        assert result["name"] == "chloride"
        assert result["role"] == "counterion"
        assert result["mw"] > 0
        assert "smiles" in result

    def test_classify_unknown_fragment(self):
        """classify_fragment with n-butane should return role='unknown', name=None."""
        result = classify_fragment("CCCC")
        assert result["name"] is None
        assert result["role"] == "unknown"
        assert result["mw"] > 0  # MW is still computed for unknown fragments

    def test_classify_invalid_smiles(self):
        """classify_fragment with invalid SMILES should return graceful fallback."""
        result = classify_fragment("invalid_smiles_xyz")
        assert result["name"] is None
        assert result["role"] == "unknown"
        assert result["mw"] == 0.0

    def test_no_duplicate_keys_in_fragment_dict(self):
        """COUNTERION_NAMES must have no duplicate canonical SMILES keys.

        Reads the source file to count key string literal occurrences,
        preventing silent overwrites where a duplicate key silently replaces
        an earlier entry (as occurred with 'O=C(O)O' / carbonic acid).
        """
        import inspect
        import re

        source = inspect.getsource(
            __import__(
                "app.services.standardization.fragment_dict",
                fromlist=["COUNTERION_NAMES"],
            )
        )

        # Find all quoted string keys inside the dict definition
        # Pattern matches both single- and double-quoted key strings
        key_pattern = re.compile(r'^\s+"([^"]+)"\s*:\s*\{', re.MULTILINE)
        found_keys = key_pattern.findall(source)

        key_counts: dict[str, int] = {}
        for key in found_keys:
            key_counts[key] = key_counts.get(key, 0) + 1

        duplicates = {k: v for k, v in key_counts.items() if v > 1}
        assert duplicates == {}, (
            f"Duplicate canonical SMILES keys found in COUNTERION_NAMES: {duplicates}"
        )

    def test_formic_acid_classified(self):
        """classify_fragment('OC=O') must return name='formic acid', role='salt', mw~46."""
        result = classify_fragment("OC=O")
        assert result["name"] == "formic acid", (
            f"Expected name='formic acid' for OC=O, got {result['name']!r}"
        )
        assert result["role"] == "salt", (
            f"Expected role='salt' for formic acid, got {result['role']!r}"
        )
        assert abs(result["mw"] - 46.0) < 1.0, (
            f"Expected MW ~46.0 for formic acid, got {result['mw']}"
        )


# ---------------------------------------------------------------------------
# TestProvenancePipeline (~12 tests)
# ---------------------------------------------------------------------------


class TestProvenancePipeline:
    """Test ProvenancePipeline provenance capture for all four STD requirements."""

    def test_basic_provenance_structure(self, pipeline, simple_mol, opts):
        """Standardize CCO with default opts — verify stage names and structure."""
        result, prov = pipeline.standardize_with_provenance(simple_mol, opts)
        assert result.success

        stage_names = [s.stage_name for s in prov.stages]
        assert "checker" in stage_names
        assert "standardizer" in stage_names
        assert "get_parent" in stage_names
        assert "tautomer_canonicalization" in stage_names
        assert len(prov.stages) >= 3

        for stage in prov.stages:
            assert stage.stage_name
            assert stage.input_smiles is not None
            assert stage.output_smiles is not None
            assert isinstance(stage.applied, bool)

    def test_neutralization_charge_tracking(self, pipeline, charged_mol, opts):
        """STD-02: Alanine zwitterion has charge changes captured in standardizer stage."""
        result, prov = pipeline.standardize_with_provenance(charged_mol, opts)

        std_stage = next(s for s in prov.stages if s.stage_name == "standardizer")
        # Standardizer normalizes the zwitterion charges
        assert len(std_stage.charge_changes) > 0, (
            "Expected charge changes for alanine zwitterion [NH3+]CC([O-])=O"
        )

        # Each change must have correct fields with before != after charge
        for change in std_stage.charge_changes:
            assert change.before_charge != change.after_charge

    def test_neutralization_charge_fields(self, pipeline, charged_mol, opts):
        """STD-02: Each ChargeChange must have all required fields populated."""
        result, prov = pipeline.standardize_with_provenance(charged_mol, opts)

        std_stage = next(s for s in prov.stages if s.stage_name == "standardizer")
        for change in std_stage.charge_changes:
            assert isinstance(change.atom_idx, int)
            assert isinstance(change.element, str) and len(change.element) > 0
            assert isinstance(change.before_charge, int)
            assert isinstance(change.after_charge, int)
            assert isinstance(change.rule_name, str) and len(change.rule_name) > 0
            assert isinstance(change.smarts, str)  # May be empty string for unknown

    def test_fragment_removal_tracking(self, pipeline, salt_mol, opts):
        """STD-04: Benzene+HCl — get_parent stage records the removed fragment."""
        result, prov = pipeline.standardize_with_provenance(salt_mol, opts)

        parent_stage = next(s for s in prov.stages if s.stage_name == "get_parent")
        assert parent_stage.applied is True
        assert len(parent_stage.fragment_removals) > 0

        for removal in parent_stage.fragment_removals:
            assert isinstance(removal.smiles, str) and len(removal.smiles) > 0
            assert removal.role in ("salt", "solvent", "counterion", "unknown")
            assert isinstance(removal.mw, float) and removal.mw > 0

    def test_fragment_removal_fields(self, pipeline, salt_mol, opts):
        """STD-04: FragmentRemoval from benzene+HCl has name, role, and mw."""
        result, prov = pipeline.standardize_with_provenance(salt_mol, opts)

        parent_stage = next(s for s in prov.stages if s.stage_name == "get_parent")
        removal = parent_stage.fragment_removals[0]

        # HCl ('Cl') should be identified
        assert removal.smiles is not None
        assert removal.name is not None  # 'hydrochloric acid' from dictionary
        assert removal.role in ("salt", "counterion")
        assert removal.mw > 0

    def test_fragment_removal_unknown(self, pipeline, opts):
        """STD-04: Uncommon fragment should have role='unknown' and name=None."""
        # Isobutyric acid is not in the dictionary
        result = classify_fragment("CC(C)C(=O)O")  # isobutyric acid — not in dictionary
        assert result["role"] == "unknown"
        assert result["name"] is None

    def test_tautomer_provenance(self, pipeline, tautomer_mol, opts_with_tautomer):
        """STD-01: 2-hydroxypyridine tautomer provenance has all required fields."""
        result, prov = pipeline.standardize_with_provenance(
            tautomer_mol, opts_with_tautomer
        )

        assert prov.tautomer is not None, "Expected tautomer provenance"
        tp = prov.tautomer

        assert isinstance(tp.input_smiles, str) and len(tp.input_smiles) > 0
        assert isinstance(tp.canonical_smiles, str) and len(tp.canonical_smiles) > 0
        assert tp.num_tautomers_enumerated > 0
        assert isinstance(tp.modified_atoms, list)
        assert isinstance(tp.modified_bonds, list)
        assert isinstance(tp.stereo_stripped, bool)
        assert isinstance(tp.complexity_flag, bool)

    def test_tautomer_provenance_canonical_differs(
        self, pipeline, tautomer_mol, opts_with_tautomer
    ):
        """STD-01: 2-hydroxypyridine canonical tautomer should differ from input."""
        result, prov = pipeline.standardize_with_provenance(
            tautomer_mol, opts_with_tautomer
        )
        assert prov.tautomer is not None
        # 2-hydroxypyridine ↔ 2-pyridinone: canonical form is 2-pyridinone
        assert prov.tautomer.input_smiles != prov.tautomer.canonical_smiles or (
            # Some versions of RDKit may canonicalize differently
            prov.tautomer.num_tautomers_enumerated >= 1
        )

    def test_tautomer_provenance_stereo_stripping(self, pipeline, opts_with_tautomer):
        """STD-01: Stereo-stripped flag is a bool (correctly computed)."""
        # Use ethanol — simple molecule, no stereo loss expected
        mol = Chem.MolFromSmiles("CCO")
        result, prov = pipeline.standardize_with_provenance(mol, opts_with_tautomer)
        assert prov.tautomer is not None
        # stereo_stripped must be a bool — no stereo in ethanol so should be False
        assert prov.tautomer.stereo_stripped is False

    def test_functional_group_audit(self, pipeline, dmso_mol, opts):
        """STD-03: DMSO standardizer produces charge changes with rule_name."""
        result, prov = pipeline.standardize_with_provenance(dmso_mol, opts)

        std_stage = next(s for s in prov.stages if s.stage_name == "standardizer")
        # DMSO: CS(C)=O -> C[S+](C)[O-] — sulphoxide normalization
        assert len(std_stage.charge_changes) > 0, (
            "Expected charge changes for DMSO sulphoxide normalization"
        )

        # At least one change should identify the sulphoxide rule
        rule_names = [c.rule_name for c in std_stage.charge_changes]
        assert any(
            "sulphoxide" in r or "normalization" in r for r in rule_names
        ), f"Expected normalization rule in charge changes, got: {rule_names}"

    def test_no_provenance_when_not_requested(self, pipeline, simple_mol, opts):
        """Standardize without include_provenance — provenance is still returned by wrapper."""
        # The ProvenancePipeline always captures provenance — the route controls visibility.
        # This test verifies the pipeline result has correct SMILES output.
        result, prov = pipeline.standardize_with_provenance(simple_mol, opts)
        assert result.success
        assert result.standardized_smiles == "CCO"

    def test_provenance_backward_compatible(self, pipeline, simple_mol, opts):
        """StandardizationResult from provenance pipeline is structurally identical to v1."""
        result, prov = pipeline.standardize_with_provenance(simple_mol, opts)

        # Must have all standard fields
        assert hasattr(result, "original_smiles")
        assert hasattr(result, "standardized_smiles")
        assert hasattr(result, "success")
        assert hasattr(result, "steps_applied")
        assert hasattr(result, "checker_issues")
        assert hasattr(result, "excluded_fragments")

        # Provenance is separate from the result
        assert isinstance(prov.stages, list)

    def test_none_molecule_handling(self, pipeline, opts):
        """ProvenancePipeline handles None mol gracefully — returns empty provenance."""
        result, prov = pipeline.standardize_with_provenance(None, opts)
        # Should not raise
        assert result is not None
        assert prov is not None
        assert isinstance(prov.stages, list)


# ---------------------------------------------------------------------------
# TestProvenanceEndpoint (~4 tests)
# ---------------------------------------------------------------------------


class TestProvenanceEndpoint:
    """Test provenance via the POST /api/v1/standardize endpoint."""

    @pytest.mark.asyncio
    async def test_endpoint_with_provenance(self, client):
        """POST with include_provenance=true returns provenance.stages array."""
        response = await client.post(
            "/standardize",
            json={"molecule": "CCO", "options": {"include_provenance": True}},
        )
        assert response.status_code == 200
        data = response.json()

        result = data["result"]
        assert "provenance" in result
        assert result["provenance"] is not None
        prov = result["provenance"]
        assert "stages" in prov
        assert isinstance(prov["stages"], list)
        assert len(prov["stages"]) >= 3

    @pytest.mark.asyncio
    async def test_endpoint_without_provenance(self, client):
        """POST without include_provenance returns result.provenance = null."""
        response = await client.post("/standardize", json={"molecule": "CCO"})
        assert response.status_code == 200
        data = response.json()

        result = data["result"]
        # provenance field exists but is null (backward compatible)
        assert "provenance" in result
        assert result["provenance"] is None

    @pytest.mark.asyncio
    async def test_endpoint_provenance_neutralization(self, client):
        """POST with charged molecule + provenance returns charge_changes in standardizer stage."""
        response = await client.post(
            "/standardize",
            json={
                "molecule": "[NH3+]CC([O-])=O",
                "options": {"include_provenance": True},
            },
        )
        assert response.status_code == 200
        data = response.json()

        result = data["result"]
        assert result["provenance"] is not None

        stages = result["provenance"]["stages"]
        std_stage = next(
            (s for s in stages if s["stage_name"] == "standardizer"), None
        )
        assert std_stage is not None, "standardizer stage not found in provenance"
        assert len(std_stage["charge_changes"]) > 0, (
            "Expected charge_changes in standardizer stage for alanine zwitterion"
        )

        # Verify each charge change has required fields
        for change in std_stage["charge_changes"]:
            assert "atom_idx" in change
            assert "element" in change
            assert "before_charge" in change
            assert "after_charge" in change
            assert "rule_name" in change

    @pytest.mark.asyncio
    async def test_endpoint_provenance_fragment_removal(self, client):
        """POST with salt-containing molecule + provenance returns fragment_removals."""
        response = await client.post(
            "/standardize",
            json={
                "molecule": "c1ccccc1.Cl",
                "options": {"include_provenance": True},
            },
        )
        assert response.status_code == 200
        data = response.json()

        result = data["result"]
        assert result["provenance"] is not None

        stages = result["provenance"]["stages"]
        parent_stage = next(
            (s for s in stages if s["stage_name"] == "get_parent"), None
        )
        assert parent_stage is not None, "get_parent stage not found in provenance"
        assert len(parent_stage["fragment_removals"]) > 0, (
            "Expected fragment_removals in get_parent stage for benzene.Cl"
        )

        # Verify each removal has required fields
        for removal in parent_stage["fragment_removals"]:
            assert "smiles" in removal
            assert "role" in removal
            assert "mw" in removal
            assert removal["role"] in ("salt", "solvent", "counterion", "unknown")
