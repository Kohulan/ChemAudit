"""
Provenance Pipeline for Standardization Intelligence.

Wraps the existing ChEMBL StandardizationPipeline to capture atom-level diffs
at each pipeline stage. Does NOT modify StandardizationPipeline internals.

Implements Milestone 2.1 provenance tracking for:
- STD-01: Tautomer canonicalization provenance
- STD-02: Neutralization charge change tracking
- STD-03: Functional group audit (normalization rule identification)
- STD-04: Parent extraction fragment naming
"""

from typing import Optional

from chembl_structure_pipeline import checker, get_parent_mol, standardizer
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from app.schemas.standardization import (
    BondChange,
    ChargeChange,
    FragmentRemoval,
    ProvStageRecord,
    RadicalChange,
    RingChange,
    StereoCenterDetailSchema,
    StereoProvenance,
    StandardizationProvenance,
    TautomerProvenance,
)
from app.services.standardization.chembl_pipeline import (
    StandardizationOptions,
    StandardizationPipeline,
    StandardizationResult,
)
from app.services.standardization.fragment_dict import classify_fragment
from app.services.standardization.stereo_tracker import StereoTracker

# Mapping from (element, change_type) patterns to normalization rule information.
# Used to identify which ChEMBL normalization rule produced a given change.
# Format: key is a string pattern, value has rule_name and smarts.
NORMALIZE_RULE_NAMES: dict[str, dict[str, str]] = {
    # Nitro group normalization: N(=O)=O -> [N+](=O)[O-]
    "N_charge_increase": {
        "rule_name": "nitro_normalization",
        "smarts": "[N;X3:1](=[O:2])=[O:3]>>[N+:1](=[O:2])[O-:3]",
    },
    # Sulphoxide normalization: S(=O) -> [S+][O-]
    "S_charge_increase": {
        "rule_name": "sulphoxide_normalization",
        "smarts": "[S;X3:1](=[O:2])>>[S+:1][O-:2]",
    },
    # Oxide anion from nitrogen: nitroxide
    "N_charge_increase_O_charge_decrease": {
        "rule_name": "nitroxide_normalization",
        "smarts": "[N:1]=[O:2]>>[N+:1][O-:2]",
    },
    # Diazo/diazonium normalization
    "N_diazo_normalization": {
        "rule_name": "diazo_normalization",
        "smarts": "[C:1]=[N:2]=[N:3]>>[C:1][N+:2]#[N:3]",
    },
    # Azide normalization
    "N_azide_normalization": {
        "rule_name": "azide_normalization",
        "smarts": "[N:1]=[N:2]=[N:3]>>[N-:1]=[N+:2]=[N:3]",
    },
    # Phosphate normalization
    "P_charge_increase": {
        "rule_name": "phosphate_normalization",
        "smarts": "[P;X4:1](=[O:2])>>[P+:1][O-:2]",
    },
    # Carboxylate normalization (deprotonation)
    "O_charge_decrease_acidic": {
        "rule_name": "carboxylate_normalization",
        "smarts": "[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])[O-:3]",
    },
    # Amine protonation/deprotonation
    "N_charge_decrease_amine": {
        "rule_name": "amine_neutralization",
        "smarts": "[N+:1]>>[N:1]",
    },
    # Metal salt normalization
    "metal_disconnection": {
        "rule_name": "metal_disconnection",
        "smarts": "[M:1]-[A:2]>>[M:1].[A:2]",
    },
}


def _identify_normalization_rule(
    atom: Chem.Atom,
    before_charge: int,
    after_charge: int,
    before_radicals: int,
    after_radicals: int,
) -> tuple[str, str]:
    """
    Attempt to identify the normalization rule that caused a change.

    Args:
        atom: The atom that changed.
        before_charge: Formal charge before standardization.
        after_charge: Formal charge after standardization.
        before_radicals: Radical electrons before.
        after_radicals: Radical electrons after.

    Returns:
        Tuple of (rule_name, smarts). Falls back to "unknown_normalization".
    """
    element = atom.GetSymbol()
    charge_delta = after_charge - before_charge

    # Nitrogen increasing charge (nitro, nitroxide, etc.)
    if element == "N" and charge_delta > 0:
        return (
            NORMALIZE_RULE_NAMES["N_charge_increase"]["rule_name"],
            NORMALIZE_RULE_NAMES["N_charge_increase"]["smarts"],
        )
    # Sulfur increasing charge (sulphoxide)
    if element == "S" and charge_delta > 0:
        return (
            NORMALIZE_RULE_NAMES["S_charge_increase"]["rule_name"],
            NORMALIZE_RULE_NAMES["S_charge_increase"]["smarts"],
        )
    # Phosphorus increasing charge (phosphate)
    if element == "P" and charge_delta > 0:
        return (
            NORMALIZE_RULE_NAMES["P_charge_increase"]["rule_name"],
            NORMALIZE_RULE_NAMES["P_charge_increase"]["smarts"],
        )
    # Oxygen decreasing charge (carboxylate)
    if element == "O" and charge_delta < 0:
        return (
            NORMALIZE_RULE_NAMES["O_charge_decrease_acidic"]["rule_name"],
            NORMALIZE_RULE_NAMES["O_charge_decrease_acidic"]["smarts"],
        )
    # Nitrogen decreasing charge (amine neutralization)
    if element == "N" and charge_delta < 0:
        return (
            NORMALIZE_RULE_NAMES["N_charge_decrease_amine"]["rule_name"],
            NORMALIZE_RULE_NAMES["N_charge_decrease_amine"]["smarts"],
        )

    return "unknown_normalization", ""


class ProvenancePipeline:
    """
    Provenance-capturing wrapper around StandardizationPipeline.

    Runs the same ChEMBL pipeline stages but captures atom-level diffs
    at each boundary. Does NOT modify StandardizationPipeline internals.

    Usage:
        pipeline = ProvenancePipeline()
        result, provenance = pipeline.standardize_with_provenance(mol, options)
    """

    def __init__(self) -> None:
        """Initialize the provenance pipeline."""
        self._pipeline = StandardizationPipeline()
        self._tautomer_enumerator = rdMolStandardize.TautomerEnumerator()

    def standardize_with_provenance(
        self,
        mol: Optional[Chem.Mol],
        options: Optional[StandardizationOptions] = None,
        dval_results: Optional[dict] = None,
    ) -> tuple[StandardizationResult, StandardizationProvenance]:
        """
        Standardize a molecule and capture per-stage provenance records.

        Runs all four pipeline stages (checker, standardizer, get_parent,
        tautomer) and records atom-level diffs at each boundary.

        Args:
            mol: RDKit molecule. If None, returns error result.
            options: Standardization options.
            dval_results: Optional prior deep-validation results for cross-referencing.
                Accepted keys: 'undefined_stereo' (DVAL-01), 'tautomer_detection' (DVAL-03).
                Each value should contain a 'count' field. When None, dval_cross_refs remain
                empty (backward compatible).

        Returns:
            Tuple of (StandardizationResult, StandardizationProvenance).
            The StandardizationResult is identical to what StandardizationPipeline
            would return (backward compatible). The provenance captures per-stage diffs.
        """
        if options is None:
            options = StandardizationOptions()

        # Run the existing pipeline first to get the canonical result
        pipeline_result = self._pipeline.standardize(mol, options)

        # If molecule is None, return empty provenance
        if mol is None:
            return pipeline_result, StandardizationProvenance(stages=[])

        # Now run each stage individually to capture provenance
        stages: list[ProvStageRecord] = []
        tautomer_provenance: Optional[TautomerProvenance] = None
        stereo_summary: Optional[StereoProvenance] = None

        try:
            # Stage 1: Checker provenance
            checker_stage = self._capture_checker_provenance(mol)
            stages.append(checker_stage)

            # Stage 2: Standardizer provenance (charge/radical/bond diffs)
            std_mol, standardizer_stage = self._capture_standardizer_provenance(mol)
            stages.append(standardizer_stage)

            if std_mol is None:
                std_mol = mol

            # Stage 3: GetParent provenance (fragment removal with classification)
            parent_mol, parent_stage = self._capture_get_parent_provenance(std_mol)
            stages.append(parent_stage)

            if parent_mol is None:
                parent_mol = std_mol

            # Stage 4: Tautomer provenance (optional)
            # Populate DVAL-03 cross-ref on tautomer stage when dval_results provided.
            taut_dval_cross_refs: list[str] = []
            if dval_results and "tautomer_detection" in dval_results:
                count = dval_results["tautomer_detection"].get("count", 0)
                taut_dval_cross_refs.append(f"DVAL-03: {count} tautomers enumerated")

            if options.include_tautomer:
                final_mol, taut_stage, tautomer_provenance = (
                    self._capture_tautomer_provenance(parent_mol)
                )
                taut_stage.dval_cross_refs = taut_dval_cross_refs
                stages.append(taut_stage)
            else:
                stages.append(
                    ProvStageRecord(
                        stage_name="tautomer_canonicalization",
                        input_smiles=self._safe_smiles(parent_mol),
                        output_smiles=self._safe_smiles(parent_mol),
                        applied=False,
                        dval_cross_refs=taut_dval_cross_refs,
                    )
                )
                final_mol = parent_mol

            # Stereo summary with per-center detail (STD-06)
            if pipeline_result.stereo_comparison:
                sc = pipeline_result.stereo_comparison
                # Build per-center detail by comparing overall before/after stereo
                stereo_before_overall = StereoTracker.get_stereo_info(mol)
                stereo_after_overall = StereoTracker.get_stereo_info(final_mol)
                stereo_comparison_detailed = StereoTracker.compare(
                    stereo_before_overall, stereo_after_overall, reason="standardization"
                )
                # Convert StereoCenterDetail dataclasses to Pydantic schemas
                per_center_schemas = [
                    StereoCenterDetailSchema(
                        atom_idx=cd.atom_idx,
                        type=cd.type,
                        before_config=cd.before_config,
                        after_config=cd.after_config,
                        reason=cd.reason,
                    )
                    for cd in stereo_comparison_detailed.per_center_detail
                ]
                # Populate DVAL cross-references when prior validation results are provided.
                # This links standardization stereo provenance to deep-validation findings,
                # allowing the UI to show "DVAL-01 found N undefined stereocenters" alongside
                # the stereo summary. When dval_results is None, the list remains empty
                # (backward compatible — same behavior as before this feature).
                stereo_dval_cross_refs: list[str] = []
                if dval_results and "undefined_stereo" in dval_results:
                    count = dval_results["undefined_stereo"].get("count", 0)
                    stereo_dval_cross_refs.append(
                        f"DVAL-01: {count} undefined stereocenters detected"
                    )
                stereo_summary = StereoProvenance(
                    stereo_stripped=sc.has_stereo_loss,
                    centers_lost=sc.stereocenters_lost,
                    bonds_lost=sc.double_bond_stereo_lost,
                    per_center=per_center_schemas,
                    dval_cross_refs=stereo_dval_cross_refs,
                )

        except Exception:
            # If provenance capture fails, return empty provenance
            # (the main result is still valid from the pipeline run above)
            pass

        provenance = StandardizationProvenance(
            stages=stages,
            tautomer=tautomer_provenance,
            stereo_summary=stereo_summary,
        )
        return pipeline_result, provenance

    def _safe_smiles(self, mol: Optional[Chem.Mol]) -> str:
        """Safely convert mol to SMILES, returning empty string on failure."""
        if mol is None:
            return ""
        try:
            return Chem.MolToSmiles(mol)
        except Exception:
            return ""

    def _capture_checker_provenance(self, mol: Chem.Mol) -> ProvStageRecord:
        """
        Capture provenance for the checker stage.

        The checker only detects issues — it does not modify the molecule.
        Issues are reported as checker_issues on the main result.

        Args:
            mol: Input molecule.

        Returns:
            ProvStageRecord for the checker stage.
        """
        input_smiles = self._safe_smiles(mol)
        try:
            molblock = Chem.MolToMolBlock(mol)
            issues = checker.check_molblock(molblock)
            applied = True
            # Checker does not change the molecule — output = input
            return ProvStageRecord(
                stage_name="checker",
                input_smiles=input_smiles,
                output_smiles=input_smiles,
                applied=applied,
                # checker_issues captured on main result, not in provenance changes
            )
        except Exception:
            return ProvStageRecord(
                stage_name="checker",
                input_smiles=input_smiles,
                output_smiles=input_smiles,
                applied=False,
            )

    def _capture_ring_aromaticity_changes(
        self, before_mol: Chem.Mol, after_mol: Chem.Mol
    ) -> list[RingChange]:
        """
        Detect ring-level aromaticity changes between before and after molecules (STD-05).

        A ring is "aromatic" if ALL its atoms are aromatic, "kekulized" otherwise.
        Only valid when atom count is preserved (same guard as charge tracking).

        Args:
            before_mol: Molecule before standardization.
            after_mol: Molecule after standardization.

        Returns:
            List of RingChange records for rings that changed aromaticity type.
        """
        if before_mol is None or after_mol is None:
            return []

        # Only valid when atom count is preserved
        if before_mol.GetNumAtoms() != after_mol.GetNumAtoms():
            return []

        ring_changes: list[RingChange] = []

        try:
            before_rings = list(before_mol.GetRingInfo().AtomRings())
            after_rings = list(after_mol.GetRingInfo().AtomRings())

            # Match rings by their sorted atom index sets
            # After standardization with atom count preserved, rings correspond 1:1
            for ring_atoms_tuple in before_rings:
                ring_atoms = list(ring_atoms_tuple)

                # Check aromaticity before: a ring is aromatic if ALL atoms are aromatic
                before_aromatic = all(
                    before_mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms
                )
                before_type = "aromatic" if before_aromatic else "kekulized"

                # Check aromaticity after (same atom indices — count preserved)
                after_aromatic = all(
                    after_mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms
                )
                after_type = "aromatic" if after_aromatic else "kekulized"

                if before_type != after_type:
                    ring_changes.append(
                        RingChange(
                            ring_atoms=sorted(ring_atoms),
                            ring_size=len(ring_atoms),
                            before_type=before_type,
                            after_type=after_type,
                        )
                    )
        except Exception:
            # Ring detection may fail on unusual molecules — return empty list
            pass

        return ring_changes

    def _capture_standardizer_provenance(
        self, mol: Chem.Mol
    ) -> tuple[Optional[Chem.Mol], ProvStageRecord]:
        """
        Capture provenance for the standardizer stage.

        Compares atom-by-atom before/after standardization:
        - Formal charge changes (ChargeChange) — STD-02
        - Radical electron changes (RadicalChange) — STD-03
        - Bond type changes (BondChange) — STD-03
        - Aromaticity changes — STD-03

        Only valid when atom count is preserved after standardization.

        Args:
            mol: Input molecule.

        Returns:
            Tuple of (standardized_mol_or_None, ProvStageRecord).
        """
        input_smiles = self._safe_smiles(mol)
        charge_changes: list[ChargeChange] = []
        radical_changes: list[RadicalChange] = []
        bond_changes: list[BondChange] = []

        try:
            std_mol = standardizer.standardize_mol(mol)
            if std_mol is None:
                return None, ProvStageRecord(
                    stage_name="standardizer",
                    input_smiles=input_smiles,
                    output_smiles=input_smiles,
                    applied=False,
                )

            output_smiles = self._safe_smiles(std_mol)
            applied = input_smiles != output_smiles

            # Atom-level diff — only valid when atom count is preserved
            if mol.GetNumAtoms() == std_mol.GetNumAtoms():
                for idx in range(mol.GetNumAtoms()):
                    before_atom = mol.GetAtomWithIdx(idx)
                    after_atom = std_mol.GetAtomWithIdx(idx)

                    element = before_atom.GetSymbol()
                    before_charge = before_atom.GetFormalCharge()
                    after_charge = after_atom.GetFormalCharge()
                    before_radicals = before_atom.GetNumRadicalElectrons()
                    after_radicals = after_atom.GetNumRadicalElectrons()

                    # Charge changes (STD-02, STD-03)
                    if before_charge != after_charge:
                        rule_name, smarts = _identify_normalization_rule(
                            before_atom,
                            before_charge,
                            after_charge,
                            before_radicals,
                            after_radicals,
                        )
                        charge_changes.append(
                            ChargeChange(
                                atom_idx=idx,
                                element=element,
                                before_charge=before_charge,
                                after_charge=after_charge,
                                rule_name=rule_name,
                                smarts=smarts,
                            )
                        )

                    # Radical electron changes (STD-03)
                    if before_radicals != after_radicals:
                        radical_changes.append(
                            RadicalChange(
                                atom_idx=idx,
                                element=element,
                                before_radicals=before_radicals,
                                after_radicals=after_radicals,
                            )
                        )

                # Bond-level diff for bond type changes (STD-03)
                if mol.GetNumBonds() == std_mol.GetNumBonds():
                    for bond_idx in range(mol.GetNumBonds()):
                        before_bond = mol.GetBondWithIdx(bond_idx)
                        after_bond = std_mol.GetBondWithIdx(bond_idx)

                        before_type = str(before_bond.GetBondType()).split(".")[-1]
                        after_type = str(after_bond.GetBondType()).split(".")[-1]

                        if before_type != after_type:
                            rule_name = "unknown_normalization"
                            # Try to identify rule from connected atoms
                            atom1 = before_bond.GetBeginAtom()
                            atom2 = before_bond.GetEndAtom()
                            connected_elements = {atom1.GetSymbol(), atom2.GetSymbol()}
                            if "N" in connected_elements and "O" in connected_elements:
                                rule_name = "nitro_normalization"
                            elif "S" in connected_elements and "O" in connected_elements:
                                rule_name = "sulphoxide_normalization"

                            bond_changes.append(
                                BondChange(
                                    bond_idx=bond_idx,
                                    atom1_idx=before_bond.GetBeginAtomIdx(),
                                    atom2_idx=before_bond.GetEndAtomIdx(),
                                    before_type=before_type,
                                    after_type=after_type,
                                    rule_name=rule_name,
                                )
                            )

            # Ring aromaticity tracking (STD-05)
            ring_changes = self._capture_ring_aromaticity_changes(mol, std_mol)

            return std_mol, ProvStageRecord(
                stage_name="standardizer",
                input_smiles=input_smiles,
                output_smiles=output_smiles,
                applied=applied,
                charge_changes=charge_changes,
                radical_changes=radical_changes,
                bond_changes=bond_changes,
                ring_changes=ring_changes,
            )

        except Exception:
            return None, ProvStageRecord(
                stage_name="standardizer",
                input_smiles=input_smiles,
                output_smiles=input_smiles,
                applied=False,
            )

    def _capture_get_parent_provenance(
        self, mol: Chem.Mol
    ) -> tuple[Optional[Chem.Mol], ProvStageRecord]:
        """
        Capture provenance for the get_parent stage.

        Uses SMILES fragment diffing (not atom-idx comparison) to find
        removed fragments and classify each via classify_fragment().
        This avoids atom-idx comparison pitfall when atom count changes.

        Args:
            mol: Input molecule (post-standardizer).

        Returns:
            Tuple of (parent_mol_or_None, ProvStageRecord).
        """
        input_smiles = self._safe_smiles(mol)
        fragment_removals: list[FragmentRemoval] = []

        try:
            parent_mol, exclude = get_parent_mol(mol)
            if parent_mol is None:
                return None, ProvStageRecord(
                    stage_name="get_parent",
                    input_smiles=input_smiles,
                    output_smiles=input_smiles,
                    applied=False,
                )

            output_smiles = self._safe_smiles(parent_mol)
            applied = input_smiles != output_smiles

            # Fragment diff using SMILES splitting — per research recommendation
            # Avoids atom-idx comparison when atom count changes after removal
            if applied and input_smiles:
                before_parts = set(input_smiles.split("."))
                after_parts = set(output_smiles.split("."))
                removed_smiles = before_parts - after_parts

                for frag_smiles in removed_smiles:
                    classification = classify_fragment(frag_smiles)
                    fragment_removals.append(
                        FragmentRemoval(
                            smiles=classification["smiles"],
                            name=classification["name"],
                            role=classification["role"],
                            mw=classification["mw"],
                        )
                    )

            return parent_mol, ProvStageRecord(
                stage_name="get_parent",
                input_smiles=input_smiles,
                output_smiles=output_smiles,
                applied=applied,
                fragment_removals=fragment_removals,
            )

        except Exception:
            return None, ProvStageRecord(
                stage_name="get_parent",
                input_smiles=input_smiles,
                output_smiles=input_smiles,
                applied=False,
            )

    def _capture_tautomer_provenance(
        self, mol: Chem.Mol
    ) -> tuple[Optional[Chem.Mol], ProvStageRecord, TautomerProvenance]:
        """
        Capture provenance for the tautomer canonicalization stage (STD-01).

        Enumerates tautomers first to get modifiedAtoms/modifiedBonds and count.
        Then canonicalizes. Detects stereo stripping by comparing chiral centers.

        Per research pitfall: uses result.tautomers (Python-iterable Mol objects),
        NOT result.smiles (C++ vector). Does NOT include all tautomer SMILES.
        Sets complexity_flag=True if > 100 tautomers.

        Args:
            mol: Input molecule (post-get_parent).

        Returns:
            Tuple of (canonical_mol_or_None, ProvStageRecord, TautomerProvenance).
        """
        input_smiles = self._safe_smiles(mol)

        try:
            # Enumerate first to get modified atoms/bonds and tautomer count
            enum_result = self._tautomer_enumerator.Enumerate(mol)
            # Use Python-iterable tautomers (NOT .smiles C++ vector)
            tautomers = list(enum_result.tautomers)
            num_tautomers = len(tautomers)

            # Get modified atoms and bonds from enumeration
            modified_atoms: list[int] = list(enum_result.modifiedAtoms)
            modified_bonds: list[int] = list(enum_result.modifiedBonds)

            # Get stereo info before canonicalization
            stereo_before = StereoTracker.get_stereo_info(mol)

            # Canonicalize
            canon_mol = self._tautomer_enumerator.Canonicalize(mol)
            if canon_mol is None:
                canon_mol = mol

            output_smiles = self._safe_smiles(canon_mol)

            # Detect stereo stripping
            stereo_after = StereoTracker.get_stereo_info(canon_mol)
            stereo_comparison = StereoTracker.compare(stereo_before, stereo_after)
            stereo_stripped = stereo_comparison.has_stereo_loss

            complexity_flag = num_tautomers > 100
            applied = input_smiles != output_smiles

            taut_prov = TautomerProvenance(
                input_smiles=input_smiles,
                canonical_smiles=output_smiles,
                num_tautomers_enumerated=num_tautomers,
                modified_atoms=modified_atoms,
                modified_bonds=modified_bonds,
                stereo_stripped=stereo_stripped,
                complexity_flag=complexity_flag,
            )

            return canon_mol, ProvStageRecord(
                stage_name="tautomer_canonicalization",
                input_smiles=input_smiles,
                output_smiles=output_smiles,
                applied=applied,
            ), taut_prov

        except Exception:
            # Tautomer provenance failure — return minimal data
            taut_prov = TautomerProvenance(
                input_smiles=input_smiles,
                canonical_smiles=input_smiles,
                num_tautomers_enumerated=0,
                modified_atoms=[],
                modified_bonds=[],
                stereo_stripped=False,
            )
            return mol, ProvStageRecord(
                stage_name="tautomer_canonicalization",
                input_smiles=input_smiles,
                output_smiles=input_smiles,
                applied=False,
            ), taut_prov
