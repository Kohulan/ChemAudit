"""
QSAR-Ready Curation Pipeline.

A 10-step configurable curation pipeline for preparing molecules for QSAR modeling.

Steps (in order):
  1. parse      - SMILES → RDKit Mol, validate parsability
  2. metals     - Disconnect metal-organic bonds (MetalDisconnector)
  3. desalt     - Keep largest fragment (LargestFragmentChooser)
  4. normalize  - Normalize functional groups (Normalizer)
  5. neutralize - Neutralize charges (Uncharger)
  6. tautomer   - Canonical tautomer (TautomerEnumerator)
  7. stereo     - Strip stereochemistry (RemoveStereochemistry) [opt]
  8. isotope    - Strip isotope labels (SetIsotope(0)) [opt]
  9. filter     - Composition filter (HA count, MW, carbon check)
  10. canonical  - Generate canonical SMILES (MolToSmiles)

References:
  - Mansouri et al. 2024 — QSAR-2D/3D preset configurations
  - RDKit MolStandardize API (rdMolStandardize module)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import inchi as rdkit_inchi
from rdkit.Chem.MolStandardize import rdMolStandardize


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class QSARStepResult:
    """Result of a single pipeline step."""

    step_name: str
    """One of: parse, metals, desalt, normalize, neutralize, tautomer, stereo, isotope, filter, canonical."""

    step_index: int
    """1-based step index (1 through 10)."""

    enabled: bool
    """Was this step enabled in the config?"""

    status: str
    """One of: 'applied', 'no_change', 'skipped', 'error'."""

    before_smiles: Optional[str] = None
    """SMILES before step (None for skipped/parse steps)."""

    after_smiles: Optional[str] = None
    """SMILES after step (None for skipped/error steps)."""

    detail: Optional[str] = None
    """Rejection reason, error message, or change description."""


@dataclass
class QSARReadyConfig:
    """Configuration for the QSAR-Ready curation pipeline.

    Default configuration enables all curation steps except stereo stripping.
    Use classmethods qsar_2d(), qsar_3d(), or minimal() for preset configurations.
    """

    enable_metals: bool = True
    """Disconnect metal-organic bonds (step 2)."""

    enable_desalt: bool = True
    """Keep largest fragment / remove salts (step 3)."""

    enable_normalize: bool = True
    """Normalize functional groups (step 4)."""

    enable_neutralize: bool = True
    """Neutralize charges (step 5)."""

    enable_tautomer: bool = True
    """Canonicalize tautomers (step 6)."""

    enable_stereo_strip: bool = False
    """Strip stereochemistry (step 7). True for QSAR-2D, False for QSAR-3D."""

    enable_isotope_strip: bool = True
    """Strip isotope labels (step 8)."""

    min_heavy_atoms: int = 3
    """Minimum heavy atom count for composition filter (step 9)."""

    max_heavy_atoms: int = 100
    """Maximum heavy atom count for composition filter (step 9)."""

    max_mw: float = 1500.0
    """Maximum molecular weight for composition filter (step 9)."""

    remove_inorganics: bool = True
    """Reject molecules with no carbon atoms in composition filter (step 9)."""

    @classmethod
    def qsar_2d(cls) -> "QSARReadyConfig":
        """QSAR-2D preset: all steps on, strip stereo + isotopes.

        Industry standard (Mansouri et al. 2024) for 2D QSAR model preparation.
        Stereo and isotope information removed for consistent 2D representations.
        """
        return cls(
            enable_metals=True,
            enable_desalt=True,
            enable_normalize=True,
            enable_neutralize=True,
            enable_tautomer=True,
            enable_stereo_strip=True,
            enable_isotope_strip=True,
        )

    @classmethod
    def qsar_3d(cls) -> "QSARReadyConfig":
        """QSAR-3D preset: all steps on, preserve stereo, strip isotopes.

        Industry standard (Mansouri et al. 2024) for 3D QSAR model preparation.
        Stereochemistry preserved; isotope information removed.
        """
        return cls(
            enable_metals=True,
            enable_desalt=True,
            enable_normalize=True,
            enable_neutralize=True,
            enable_tautomer=True,
            enable_stereo_strip=False,
            enable_isotope_strip=True,
        )

    @classmethod
    def minimal(cls) -> "QSARReadyConfig":
        """Minimal preset: parse + desalt + normalize + canonical only.

        Runs only the core structure-normalizing steps. Skips metals, neutralize,
        tautomer, stereo, and isotope for lightest-touch curation.
        """
        return cls(
            enable_metals=False,
            enable_desalt=True,
            enable_normalize=True,
            enable_neutralize=False,
            enable_tautomer=False,
            enable_stereo_strip=False,
            enable_isotope_strip=False,
            remove_inorganics=False,
        )


@dataclass
class QSARReadyResult:
    """Result of processing a molecule through the QSAR-Ready pipeline."""

    original_smiles: str
    """The input SMILES string."""

    original_inchikey: Optional[str] = None
    """InChIKey computed from parsed mol BEFORE any pipeline steps."""

    curated_smiles: Optional[str] = None
    """Canonical SMILES output after all pipeline steps (None if rejected/error)."""

    standardized_inchikey: Optional[str] = None
    """InChIKey computed from final mol AFTER all pipeline steps."""

    inchikey_changed: bool = False
    """True if original_inchikey != standardized_inchikey (D-14, locked decision)."""

    status: str = "error"
    """One of: 'ok', 'rejected', 'duplicate', 'error'."""

    rejection_reason: Optional[str] = None
    """Human-readable reason for rejection (if status='rejected' or 'duplicate')."""

    steps: list = field(default_factory=list)
    """List of QSARStepResult, one per pipeline step (always 10 entries when status != 'error')."""


# ---------------------------------------------------------------------------
# Internal step helpers
# ---------------------------------------------------------------------------


def _compute_inchikey(mol: Chem.Mol) -> Optional[str]:
    """Compute InChIKey from an RDKit molecule. Returns None on failure.

    Uses MolToInchiKey(mol) which takes a mol object directly (not the InChI string).
    Confirmed API: rdkit_inchi.MolToInchiKey(mol) -> str
    """
    try:
        return rdkit_inchi.MolToInchiKey(mol)
    except Exception:
        return None


def _run_metals_step(mol: Chem.Mol) -> tuple[Optional[Chem.Mol], str]:
    """Disconnect metal-organic bonds using MetalDisconnector."""
    try:
        result_mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        return result_mol, ""
    except Exception as e:
        return None, str(e)


def _run_desalt_step(mol: Chem.Mol) -> tuple[Optional[Chem.Mol], str]:
    """Keep largest fragment using LargestFragmentChooser."""
    try:
        result_mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
        return result_mol, ""
    except Exception as e:
        return None, str(e)


def _run_normalize_step(mol: Chem.Mol) -> tuple[Optional[Chem.Mol], str]:
    """Normalize functional groups using Normalizer."""
    try:
        result_mol = rdMolStandardize.Normalizer().normalize(mol)
        return result_mol, ""
    except Exception as e:
        return None, str(e)


def _run_neutralize_step(mol: Chem.Mol) -> tuple[Optional[Chem.Mol], str]:
    """Neutralize charges using Uncharger."""
    try:
        result_mol = rdMolStandardize.Uncharger().uncharge(mol)
        return result_mol, ""
    except Exception as e:
        return None, str(e)


def _run_tautomer_step(mol: Chem.Mol, preserve_stereo: bool = True) -> tuple[Optional[Chem.Mol], str]:
    """Canonicalize tautomers using TautomerEnumerator.

    When preserve_stereo=True (QSAR-3D), SetRemoveSp3Stereo(False) prevents
    the enumerator from stripping sp3 stereochemistry during tautomer selection.
    When preserve_stereo=False (QSAR-2D), the stereo step (step 7) handles removal.
    """
    try:
        te = rdMolStandardize.TautomerEnumerator()
        if preserve_stereo:
            te.SetRemoveSp3Stereo(False)
        result_mol = te.Canonicalize(mol)
        return result_mol, ""
    except Exception as e:
        return None, str(e)


def _run_stereo_step(mol: Chem.Mol) -> tuple[Optional[Chem.Mol], str]:
    """Strip stereochemistry using RemoveStereochemistry (in-place).

    Pitfall: RemoveStereochemistry mutates the molecule in-place.
    We must capture before_smiles BEFORE calling it.
    """
    try:
        Chem.RemoveStereochemistry(mol)
        return mol, ""
    except Exception as e:
        return None, str(e)


def _run_isotope_step(mol: Chem.Mol) -> tuple[Optional[Chem.Mol], str]:
    """Strip all isotope labels by setting isotope=0 on every atom.

    Pitfall: This mutates the molecule in-place via atom.SetIsotope(0).
    We must capture before_smiles BEFORE the loop.
    """
    try:
        for atom in mol.GetAtoms():
            if atom.GetIsotope() != 0:
                atom.SetIsotope(0)
        return mol, ""
    except Exception as e:
        return None, str(e)


def _run_filter_step(
    mol: Chem.Mol, config: QSARReadyConfig
) -> tuple[Optional[Chem.Mol], str, Optional[str]]:
    """Run composition filter checks.

    Returns: (mol_or_None, error_detail, rejection_reason)
    - If passes: (mol, "", None)
    - If fails: (None, "", rejection_reason)
    """
    try:
        heavy_atoms = mol.GetNumHeavyAtoms()
        if heavy_atoms < config.min_heavy_atoms:
            return (
                None,
                "",
                f"Too few heavy atoms: {heavy_atoms} < {config.min_heavy_atoms} (minimum)",
            )
        if heavy_atoms > config.max_heavy_atoms:
            return (
                None,
                "",
                f"Too many heavy atoms: {heavy_atoms} > {config.max_heavy_atoms} (maximum)",
            )

        mw = Descriptors.ExactMolWt(mol)
        if mw > config.max_mw:
            return (
                None,
                "",
                f"Molecular weight too high: {mw:.1f} > {config.max_mw:.1f} (maximum)",
            )

        if config.remove_inorganics:
            has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
            if not has_carbon:
                return None, "", "No carbon atoms found (inorganic molecule rejected)"

        return mol, "", None
    except Exception as e:
        return None, str(e), None


# ---------------------------------------------------------------------------
# Step order definition
# ---------------------------------------------------------------------------

# Each tuple: (step_name, step_index, config_attr_name, step_runner_or_None)
# step_runner_or_None is None for special-cased steps (parse, filter, canonical)
_STEP_ORDER = [
    ("parse", 1, None, None),           # Special: SMILES → mol
    ("metals", 2, "enable_metals", _run_metals_step),
    ("desalt", 3, "enable_desalt", _run_desalt_step),
    ("normalize", 4, "enable_normalize", _run_normalize_step),
    ("neutralize", 5, "enable_neutralize", _run_neutralize_step),
    ("tautomer", 6, "enable_tautomer", _run_tautomer_step),
    ("stereo", 7, "enable_stereo_strip", _run_stereo_step),
    ("isotope", 8, "enable_isotope_strip", _run_isotope_step),
    ("filter", 9, None, None),          # Special: composition filter
    ("canonical", 10, None, None),      # Special: final canonical SMILES
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def qsar_ready_single(smiles: str, config: QSARReadyConfig) -> QSARReadyResult:
    """Process a single molecule through the QSAR-Ready 10-step curation pipeline.

    Args:
        smiles: Input SMILES string to curate.
        config: QSARReadyConfig controlling which steps are enabled and filter thresholds.

    Returns:
        QSARReadyResult with per-step provenance, InChIKey tracking, and final
        curated SMILES. Status is 'ok' on success, 'rejected' on filter/parse
        failure, 'error' on unexpected errors.
    """
    result = QSARReadyResult(original_smiles=smiles)

    # -------------------------------------------------------------------------
    # Step 1: Parse
    # -------------------------------------------------------------------------
    parse_step = QSARStepResult(
        step_name="parse",
        step_index=1,
        enabled=True,
        status="error",
    )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        parse_step.status = "error"
        parse_step.detail = "Failed to parse SMILES" if mol is None else "Empty molecule (0 atoms)"
        result.steps.append(parse_step)
        result.status = "rejected"
        result.rejection_reason = parse_step.detail
        return result

    parse_step.status = "no_change"
    parse_step.before_smiles = smiles
    parse_step.after_smiles = Chem.MolToSmiles(mol)
    result.steps.append(parse_step)

    # Compute original_inchikey BEFORE any pipeline transformations (D-14)
    result.original_inchikey = _compute_inchikey(mol)

    # -------------------------------------------------------------------------
    # Steps 2-8: Chemical transformations
    # -------------------------------------------------------------------------
    for step_name, step_index, config_attr, runner in _STEP_ORDER[1:8]:
        enabled = getattr(config, config_attr, True)

        if not enabled:
            result.steps.append(
                QSARStepResult(
                    step_name=step_name,
                    step_index=step_index,
                    enabled=False,
                    status="skipped",
                )
            )
            continue

        # Capture BEFORE state (must be done before in-place mutating steps)
        before_smiles = Chem.MolToSmiles(mol)

        # Run step — tautomer step receives preserve_stereo kwarg
        if step_name == "tautomer":
            new_mol, error_detail = runner(mol, preserve_stereo=not config.enable_stereo_strip)
        else:
            new_mol, error_detail = runner(mol)

        if new_mol is None:
            step_result = QSARStepResult(
                step_name=step_name,
                step_index=step_index,
                enabled=True,
                status="error",
                before_smiles=before_smiles,
                detail=error_detail or f"Step {step_name} returned None",
            )
            result.steps.append(step_result)
            result.status = "error"
            result.rejection_reason = f"Pipeline error at step {step_name}"
            return result

        mol = new_mol

        # Capture AFTER state
        after_smiles = Chem.MolToSmiles(mol)
        status = "applied" if before_smiles != after_smiles else "no_change"

        result.steps.append(
            QSARStepResult(
                step_name=step_name,
                step_index=step_index,
                enabled=True,
                status=status,
                before_smiles=before_smiles,
                after_smiles=after_smiles,
            )
        )

    # -------------------------------------------------------------------------
    # Step 9: Composition filter
    # -------------------------------------------------------------------------
    filter_mol, filter_error, rejection_reason = _run_filter_step(mol, config)

    if filter_error:
        # Unexpected error in filter
        result.steps.append(
            QSARStepResult(
                step_name="filter",
                step_index=9,
                enabled=True,
                status="error",
                detail=filter_error,
            )
        )
        result.status = "error"
        result.rejection_reason = f"Filter step error: {filter_error}"
        return result

    if rejection_reason is not None:
        # Composition filter rejected the molecule
        result.steps.append(
            QSARStepResult(
                step_name="filter",
                step_index=9,
                enabled=True,
                status="error",
                detail=rejection_reason,
            )
        )
        result.status = "rejected"
        result.rejection_reason = rejection_reason
        # Still need to add canonical step as skipped/error
        result.steps.append(
            QSARStepResult(
                step_name="canonical",
                step_index=10,
                enabled=True,
                status="skipped",
                detail="Skipped: molecule rejected by composition filter",
            )
        )
        return result

    result.steps.append(
        QSARStepResult(
            step_name="filter",
            step_index=9,
            enabled=True,
            status="no_change",
            detail="Passed all composition filters",
        )
    )

    # -------------------------------------------------------------------------
    # Step 10: Canonical SMILES
    # -------------------------------------------------------------------------
    try:
        curated_smiles = Chem.MolToSmiles(mol)
    except Exception as e:
        result.steps.append(
            QSARStepResult(
                step_name="canonical",
                step_index=10,
                enabled=True,
                status="error",
                detail=f"Failed to generate canonical SMILES: {e}",
            )
        )
        result.status = "error"
        result.rejection_reason = "Failed to generate canonical SMILES"
        return result

    result.steps.append(
        QSARStepResult(
            step_name="canonical",
            step_index=10,
            enabled=True,
            status="applied",
            after_smiles=curated_smiles,
        )
    )

    # -------------------------------------------------------------------------
    # Finalize: InChIKey comparison (D-14)
    # -------------------------------------------------------------------------
    result.curated_smiles = curated_smiles
    result.standardized_inchikey = _compute_inchikey(mol)
    result.inchikey_changed = result.original_inchikey != result.standardized_inchikey
    result.status = "ok"

    return result


def qsar_ready_batch(
    smiles_list: list[str], config: QSARReadyConfig
) -> list[QSARReadyResult]:
    """Process a list of SMILES through the QSAR-Ready pipeline with deduplication.

    Each molecule is processed independently through qsar_ready_single().
    After processing, molecules with status='ok' are deduplicated by
    standardized_inchikey — the first occurrence is kept as 'ok', subsequent
    occurrences are marked as 'duplicate'.

    Args:
        smiles_list: List of SMILES strings to process.
        config: QSARReadyConfig controlling pipeline steps and filter thresholds.

    Returns:
        List of QSARReadyResult, one per input SMILES (including duplicates with
        status='duplicate').
    """
    results: list[QSARReadyResult] = []
    seen_inchikeys: dict[str, int] = {}  # inchikey -> first-seen index (0-based)

    for idx, smiles in enumerate(smiles_list):
        result = qsar_ready_single(smiles, config)

        if result.status == "ok" and result.standardized_inchikey is not None:
            if result.standardized_inchikey in seen_inchikeys:
                # Mark as duplicate
                first_idx = seen_inchikeys[result.standardized_inchikey]
                result.status = "duplicate"
                result.rejection_reason = (
                    f"Duplicate of molecule {first_idx} (index {first_idx})"
                )
            else:
                seen_inchikeys[result.standardized_inchikey] = idx

        results.append(result)

    return results
