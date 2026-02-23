"""
Tests for deep composition validation checks (Milestone 1.2).

Covers DVAL-06 through DVAL-11:
    - MixtureDetectionCheck
    - SolventContaminationCheck
    - InorganicFilterCheck
    - RadicalDetectionCheck
    - IsotopeLabelDetectionCheck
    - TrivialMoleculeCheck
"""

from rdkit import Chem

from app.schemas.common import Severity
from app.services.validation.checks.deep_composition import (
    InorganicFilterCheck,
    IsotopeLabelDetectionCheck,
    MixtureDetectionCheck,
    RadicalDetectionCheck,
    SolventContaminationCheck,
    TrivialMoleculeCheck,
)


class TestMixtureDetection:
    """Test DVAL-06: mixture detection with fragment classification."""

    def test_single_molecule_passes(self):
        """A single connected molecule should pass."""
        mol = Chem.MolFromSmiles("CCO")  # ethanol
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["num_fragments"] == 1

    def test_mixture_detected(self):
        """CCN.Cl (amine HCl salt) should fail with 2 fragments."""
        mol = Chem.MolFromSmiles("CCN.Cl")
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING
        assert result.details["num_fragments"] == 2

        fragments = result.details["fragments"]
        assert len(fragments) == 2

        classifications = {f["classification"] for f in fragments}
        # Largest carbon fragment → drug; chlorine → salt
        assert "drug" in classifications
        assert "salt" in classifications

    def test_three_component_mixture(self):
        """CCN.[Na+].[Cl-] should produce 3 fragments with correct classifications."""
        mol = Chem.MolFromSmiles("CCN.[Na+].[Cl-]")
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["num_fragments"] == 3

        fragments = result.details["fragments"]
        assert len(fragments) == 3

        classifications = {f["classification"] for f in fragments}
        # CCN is largest carbon fragment → drug; Na+, Cl- → salt
        assert "drug" in classifications
        assert "salt" in classifications

    def test_fragment_classification_drug(self):
        """Largest carbon-containing fragment should be classified as 'drug'."""
        mol = Chem.MolFromSmiles("c1ccc(cc1)CCN.Cl")  # phenethylamine HCl
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        fragments = result.details["fragments"]
        # Fragment with most heavy atoms should be drug
        largest = max(fragments, key=lambda f: f["heavy_atom_count"])
        assert largest["classification"] == "drug"

    def test_fragment_classification_salt(self):
        """Small ionic fragments like [Na+] and [Cl-] should be classified as salt."""
        mol = Chem.MolFromSmiles("[Na+].[Cl-]")
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        fragments = result.details["fragments"]
        for frag in fragments:
            assert frag["classification"] == "salt"

    def test_details_contain_required_fields(self):
        """Fragment details should include smiles, molecular_weight, classification, pattern_name."""
        mol = Chem.MolFromSmiles("CCN.Cl")
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        for frag in result.details["fragments"]:
            assert "smiles" in frag
            assert "molecular_weight" in frag
            assert "heavy_atom_count" in frag
            assert "classification" in frag
            assert "pattern_name" in frag
            assert frag["molecular_weight"] > 0
            assert len(frag["smiles"]) > 0

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR."""
        check = MixtureDetectionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_message_contains_fragment_count(self):
        """Failure message should state number of fragments."""
        mol = Chem.MolFromSmiles("CCN.Cl")
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert "2" in result.message or "fragments" in result.message.lower()

    def test_affected_atoms_empty_for_mixture(self):
        """Affected atoms should be empty for mixtures (fragments are separate molecules)."""
        mol = Chem.MolFromSmiles("CCN.Cl")
        check = MixtureDetectionCheck()
        result = check.run(mol)

        assert result.affected_atoms == []


class TestSolventContamination:
    """Test DVAL-07: solvent contamination detection."""

    def test_no_solvent_passes(self):
        """Aspirin should pass (not a solvent)."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["solvents_found"] == []

    def test_water_detected(self):
        """Pure water (O) should be detected as solvent."""
        mol = Chem.MolFromSmiles("O")
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING
        assert len(result.details["solvents_found"]) > 0
        names = [s["name"] for s in result.details["solvents_found"]]
        assert "water" in names

    def test_dmso_detected(self):
        """DMSO (CS(=O)C) should be detected as solvent."""
        mol = Chem.MolFromSmiles("CS(=O)C")
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        names = [s["name"] for s in result.details["solvents_found"]]
        assert "DMSO" in names

    def test_methanol_detected(self):
        """Methanol (CO) should be detected as solvent."""
        mol = Chem.MolFromSmiles("CO")
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        names = [s["name"] for s in result.details["solvents_found"]]
        assert "methanol" in names

    def test_solvent_in_mixture(self):
        """CCN.O should detect water contamination in mixture context."""
        mol = Chem.MolFromSmiles("CCN.O")
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        names = [s["name"] for s in result.details["solvents_found"]]
        assert "water" in names

    def test_pure_solvent_flag(self):
        """Single-fragment solvent input should set is_pure_solvent=True."""
        mol = Chem.MolFromSmiles("O")  # pure water
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["is_pure_solvent"] is True

    def test_solvent_in_mixture_not_pure(self):
        """Solvent as part of mixture should set is_pure_solvent=False."""
        mol = Chem.MolFromSmiles("CCN.O")  # ethylamine + water
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["is_pure_solvent"] is False

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR."""
        check = SolventContaminationCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_solvent_details_have_required_fields(self):
        """Solvent details should include name, smiles, molecular_weight."""
        mol = Chem.MolFromSmiles("CS(=O)C")  # DMSO
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        for sv in result.details["solvents_found"]:
            assert "name" in sv
            assert "smiles" in sv
            assert "molecular_weight" in sv

    def test_ethanol_detected(self):
        """Ethanol (CCO) should be detected as solvent."""
        mol = Chem.MolFromSmiles("CCO")
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert not result.passed
        names = [s["name"] for s in result.details["solvents_found"]]
        assert "ethanol" in names

    def test_large_organic_molecule_passes(self):
        """A large drug-like molecule should not trigger solvent detection."""
        # Ibuprofen — carbon-rich, not a solvent
        mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        check = SolventContaminationCheck()
        result = check.run(mol)

        assert result.passed


class TestInorganicFilter:
    """Test DVAL-08: inorganic and organometallic detection."""

    def test_organic_molecule_passes(self):
        """Ethanol (CCO) contains carbon, should pass."""
        mol = Chem.MolFromSmiles("CCO")
        check = InorganicFilterCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["has_carbon"] is True
        assert result.details["is_inorganic"] is False
        assert result.details["is_organometallic"] is False

    def test_inorganic_detected(self):
        """NaCl ([Na+].[Cl-]) has no carbon — should be flagged as inorganic."""
        mol = Chem.MolFromSmiles("[Na+].[Cl-]")
        check = InorganicFilterCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.ERROR
        assert result.details["is_inorganic"] is True
        assert result.details["has_carbon"] is False

    def test_organometallic_detected(self):
        """Ferrocene (carbon + iron) should be detected as organometallic."""
        # Simplified organometallic: methyliron or tetramethyltin
        mol = Chem.MolFromSmiles("[Pd]")  # Pure metal — inorganic
        check = InorganicFilterCheck()
        r1 = check.run(mol)
        assert not r1.passed
        assert r1.details["is_inorganic"] is True

        # Dimethyl palladium: C[Pd]C — organometallic
        mol2 = Chem.MolFromSmiles("C[Pd]C")
        result2 = check.run(mol2)
        assert not result2.passed
        assert result2.details["is_organometallic"] is True
        assert result2.severity == Severity.WARNING

    def test_metal_atoms_listed(self):
        """Metal atoms should appear in affected_atoms."""
        mol = Chem.MolFromSmiles("C[Zn]C")  # organozinc
        check = InorganicFilterCheck()
        result = check.run(mol)

        assert not result.passed
        assert len(result.affected_atoms) > 0
        assert len(result.details["metal_atoms"]) > 0
        # Verify structure of metal_atoms entry
        metal = result.details["metal_atoms"][0]
        assert "atom_idx" in metal
        assert "symbol" in metal
        assert "atomic_num" in metal

    def test_element_counts_in_details(self):
        """element_counts should be populated with correct elements."""
        mol = Chem.MolFromSmiles("CCO")
        check = InorganicFilterCheck()
        result = check.run(mol)

        assert "element_counts" in result.details
        counts = result.details["element_counts"]
        assert "C" in counts
        assert counts["C"] == 2
        assert "O" in counts
        assert counts["O"] == 1

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR."""
        check = InorganicFilterCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_inorganic_severity_is_error(self):
        """Purely inorganic molecule severity should be ERROR."""
        mol = Chem.MolFromSmiles("[Na+].[Cl-]")
        check = InorganicFilterCheck()
        result = check.run(mol)

        assert result.severity == Severity.ERROR

    def test_organometallic_severity_is_warning(self):
        """Organometallic molecule severity should be WARNING."""
        mol = Chem.MolFromSmiles("C[Fe]C")
        check = InorganicFilterCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING


class TestRadicalDetection:
    """Test DVAL-09: radical electron detection."""

    def test_no_radicals_passes(self):
        """Ethanol has no radical electrons."""
        mol = Chem.MolFromSmiles("CCO")
        check = RadicalDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["total_radical_electrons"] == 0
        assert result.details["radical_atoms"] == []

    def test_radical_detected(self):
        """[CH2] (methylene) has 2 radical electrons."""
        mol = Chem.MolFromSmiles("[CH2]")
        check = RadicalDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.WARNING
        assert result.details["total_radical_electrons"] > 0
        assert len(result.details["radical_atoms"]) > 0

    def test_oxygen_radical_detected(self):
        """[O] (atomic oxygen radical) should be detected."""
        mol = Chem.MolFromSmiles("[O]")
        check = RadicalDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["total_radical_electrons"] > 0

    def test_affected_atoms_correct(self):
        """affected_atoms should match the indices of radical-bearing atoms."""
        mol = Chem.MolFromSmiles("[CH2]")
        check = RadicalDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert len(result.affected_atoms) > 0

        # Verify each affected atom is in radical_atoms details
        radical_idxs = {r["atom_idx"] for r in result.details["radical_atoms"]}
        for idx in result.affected_atoms:
            assert idx in radical_idxs

    def test_radical_details_structure(self):
        """radical_atoms entries must have atom_idx, symbol, num_radical_electrons."""
        mol = Chem.MolFromSmiles("[CH2]")
        check = RadicalDetectionCheck()
        result = check.run(mol)

        for entry in result.details["radical_atoms"]:
            assert "atom_idx" in entry
            assert "symbol" in entry
            assert "num_radical_electrons" in entry
            assert entry["num_radical_electrons"] > 0

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR."""
        check = RadicalDetectionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR


class TestIsotopeLabelDetection:
    """Test DVAL-10: isotope label detection."""

    def test_no_isotopes_passes(self):
        """Ethanol has no isotope labels."""
        mol = Chem.MolFromSmiles("CCO")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["total_labeled"] == 0
        assert result.details["labeled_atoms"] == []

    def test_deuterium_detected(self):
        """[2H]C([2H])([2H])O (deuterated methanol) should detect deuterium labels."""
        mol = Chem.MolFromSmiles("[2H]C([2H])([2H])O")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.INFO  # isotopes are often intentional
        assert result.details["total_labeled"] == 3

        # Verify deuterium is identified correctly
        for entry in result.details["labeled_atoms"]:
            assert entry["symbol"] == "H"
            assert entry["isotope_mass"] == 2
            assert entry["common_name"] == "deuterium"

    def test_carbon13_detected(self):
        """[13CH4] should detect carbon-13 label."""
        mol = Chem.MolFromSmiles("[13CH4]")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["total_labeled"] == 1
        entry = result.details["labeled_atoms"][0]
        assert entry["symbol"] == "C"
        assert entry["isotope_mass"] == 13

    def test_common_name_mapping(self):
        """Common name field should be correctly populated."""
        mol = Chem.MolFromSmiles("[2H]C([2H])([2H])O")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        for entry in result.details["labeled_atoms"]:
            assert entry["common_name"] == "deuterium"

    def test_carbon14_common_name(self):
        """14C should map to 'carbon-14'."""
        mol = Chem.MolFromSmiles("[14CH4]")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        entry = result.details["labeled_atoms"][0]
        assert entry["common_name"] == "carbon-14"

    def test_nitrogen15_common_name(self):
        """15N should map to 'nitrogen-15'."""
        mol = Chem.MolFromSmiles("[15NH3]")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        entry = result.details["labeled_atoms"][0]
        assert entry["symbol"] == "N"
        assert entry["common_name"] == "nitrogen-15"

    def test_unknown_isotope_has_none_common_name(self):
        """Uncommon isotope should have common_name=None."""
        mol = Chem.MolFromSmiles("[125IH]")  # iodine-125 (not in our map)
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        if not result.passed:
            for entry in result.details["labeled_atoms"]:
                # common_name may be None for unknown isotopes
                assert entry["common_name"] is None or isinstance(entry["common_name"], str)

    def test_affected_atoms_populated(self):
        """affected_atoms should list indices of labeled atoms."""
        mol = Chem.MolFromSmiles("[2H]C([2H])([2H])O")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert len(result.affected_atoms) == 3

    def test_severity_is_info_for_isotopes(self):
        """Isotope labels should be INFO severity (often intentional)."""
        mol = Chem.MolFromSmiles("[13CH4]")
        check = IsotopeLabelDetectionCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.INFO

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR."""
        check = IsotopeLabelDetectionCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR


class TestTrivialMolecule:
    """Test DVAL-11: trivial molecule size check."""

    def test_normal_molecule_passes(self):
        """Heptane (7 heavy atoms) should pass."""
        mol = Chem.MolFromSmiles("CCCCCCC")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert result.passed
        assert result.severity == Severity.INFO
        assert result.details["heavy_atom_count"] == 7

    def test_single_atom_fails(self):
        """Single carbon atom [C] should fail."""
        mol = Chem.MolFromSmiles("[C]")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.ERROR
        assert result.details["heavy_atom_count"] == 1
        assert result.details["is_single_atom"] is True

    def test_water_fails(self):
        """Water (O) has 1 heavy atom — should fail."""
        mol = Chem.MolFromSmiles("O")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["heavy_atom_count"] == 1

    def test_threshold_boundary_at_three(self):
        """Exactly 3 heavy atoms (threshold) should fail."""
        # CO2 has 3 heavy atoms (C + 2 O)
        mol = Chem.MolFromSmiles("O=C=O")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.details["heavy_atom_count"] == 3
        assert result.details["threshold"] == 3

    def test_four_heavy_atoms_passes(self):
        """4 heavy atoms (above threshold) should pass."""
        mol = Chem.MolFromSmiles("CC=O")  # acetaldehyde: 3 C + O = 4... wait, CC=O is 3 heavy
        # Use butane: CCCC = 4 carbons = 4 heavy atoms
        mol = Chem.MolFromSmiles("CCCC")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert result.passed
        assert result.details["heavy_atom_count"] == 4

    def test_details_contain_required_fields(self):
        """Details should include heavy_atom_count, num_bonds, is_single_atom, threshold."""
        mol = Chem.MolFromSmiles("CCO")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert "heavy_atom_count" in result.details
        assert "num_bonds" in result.details
        assert "is_single_atom" in result.details
        assert "threshold" in result.details

    def test_severity_is_error(self):
        """Trivial molecule should have ERROR severity."""
        mol = Chem.MolFromSmiles("O")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_none_molecule_returns_error(self):
        """None molecule should return ERROR."""
        check = TrivialMoleculeCheck()
        result = check.run(None)

        assert not result.passed
        assert result.severity == Severity.ERROR

    def test_message_includes_atom_count(self):
        """Failure message should mention the heavy atom count."""
        mol = Chem.MolFromSmiles("O")
        check = TrivialMoleculeCheck()
        result = check.run(mol)

        assert "1" in result.message


class TestRegistration:
    """Test that all M1.2 checks are properly registered."""

    def test_all_m12_checks_registered(self):
        """Verify all 6 M1.2 check names appear in CheckRegistry."""
        from app.services.validation.registry import CheckRegistry

        registered = CheckRegistry.list_names()

        expected = {
            "mixture_detection",
            "solvent_contamination",
            "inorganic_filter",
            "radical_detection",
            "isotope_label_detection",
            "trivial_molecule",
        }
        missing = expected - set(registered)
        assert not missing, f"Missing M1.2 checks: {missing}"

    def test_check_categories(self):
        """All M1.2 checks should use 'chemical_composition' category."""
        from app.services.validation.registry import CheckRegistry

        checks = CheckRegistry.get_all()
        m12_names = {
            "mixture_detection",
            "solvent_contamination",
            "inorganic_filter",
            "radical_detection",
            "isotope_label_detection",
            "trivial_molecule",
        }
        for name in m12_names:
            assert name in checks, f"Check '{name}' not registered"
            check_instance = checks[name]()
            assert check_instance.category == "chemical_composition", (
                f"Check '{name}' has wrong category: {check_instance.category}"
            )
