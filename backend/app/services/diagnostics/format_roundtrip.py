"""
Format Round-Trip Lossiness Service (DIAG-03).

Detects information loss when converting SMILES through InChI or MOL block
and back. Checks for stereo center loss, formal charge loss, and isotope loss.
"""

from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters
from rdkit.Chem.inchi import MolFromInchi, MolToInchi
from rdkit.Chem.rdchem import Mol


def _compute_losses(original_mol: Mol, roundtrip_mol: Mol) -> list[dict]:
    """Compute information losses between original and round-tripped molecules.

    Checks for stereo center count, formal charge sum, and isotope atom count.

    Args:
        original_mol: Original RDKit molecule.
        roundtrip_mol: Round-tripped RDKit molecule.

    Returns:
        List of loss dicts with keys: type, description, before, after.
    """
    losses: list[dict] = []

    # Stereo centers
    try:
        stereo_before = len(FindMolChiralCenters(original_mol, includeUnassigned=True))
        stereo_after = len(FindMolChiralCenters(roundtrip_mol, includeUnassigned=True))
        if stereo_before != stereo_after:
            diff = stereo_before - stereo_after
            losses.append({
                "type": "stereo",
                "description": f"{diff} stereo center(s) lost during round-trip conversion",
                "before": stereo_before,
                "after": stereo_after,
            })
    except Exception:
        pass

    # Formal charge
    try:
        charge_before = sum(a.GetFormalCharge() for a in original_mol.GetAtoms())
        charge_after = sum(a.GetFormalCharge() for a in roundtrip_mol.GetAtoms())
        if charge_before != charge_after:
            losses.append({
                "type": "charge",
                "description": f"Formal charge changed from {charge_before:+d} to {charge_after:+d}",
                "before": charge_before,
                "after": charge_after,
            })
    except Exception:
        pass

    # Isotope labels
    try:
        isotope_before = sum(1 for a in original_mol.GetAtoms() if a.GetIsotope() != 0)
        isotope_after = sum(1 for a in roundtrip_mol.GetAtoms() if a.GetIsotope() != 0)
        if isotope_before != isotope_after:
            diff = isotope_before - isotope_after
            losses.append({
                "type": "isotope",
                "description": f"{diff} isotope label(s) lost during round-trip conversion",
                "before": isotope_before,
                "after": isotope_after,
            })
    except Exception:
        pass

    return losses


def _make_error_result(route: str, original_smiles: str, error: str,
                       intermediate: str | None = None) -> dict:
    """Build a standardized error result dict for a failed round-trip route."""
    return {
        "route": route,
        "original_smiles": original_smiles,
        "intermediate": intermediate,
        "roundtrip_smiles": None,
        "lossy": False,
        "losses": [],
        "error": error,
    }


def _roundtrip_via_inchi(mol: Mol) -> tuple[str | None, Mol | None, str | None]:
    """Convert mol to InChI and back. Returns (intermediate, roundtrip_mol, error)."""
    inchi = MolToInchi(mol)
    if inchi is None:
        return None, None, "Failed to generate InChI from molecule"
    roundtrip_mol = MolFromInchi(inchi)
    if roundtrip_mol is None:
        return inchi, None, "Failed to parse molecule from InChI"
    return inchi, roundtrip_mol, None


def _roundtrip_via_molblock(mol: Mol) -> tuple[str | None, Mol | None, str | None]:
    """Convert mol to MOL block and back. Returns (intermediate, roundtrip_mol, error)."""
    mol_block = Chem.MolToMolBlock(mol)
    if mol_block is None:
        return None, None, "Failed to generate MOL block from molecule"
    roundtrip_mol = Chem.MolFromMolBlock(mol_block)
    if roundtrip_mol is None:
        return mol_block, None, "Failed to parse molecule from MOL block"
    return mol_block, roundtrip_mol, None


_ROUTE_FUNCTIONS = {
    "smiles_inchi_smiles": _roundtrip_via_inchi,
    "smiles_mol_smiles": _roundtrip_via_molblock,
}


def check_roundtrip(smiles: str, route: str = "smiles_inchi_smiles") -> dict:
    """Check whether a SMILES string round-trips losslessly through an intermediate format.

    Supported routes:
    - "smiles_inchi_smiles": SMILES → InChI → SMILES
    - "smiles_mol_smiles": SMILES → MOL block → SMILES

    Args:
        smiles: SMILES string to check.
        route: Conversion route identifier. One of "smiles_inchi_smiles" or "smiles_mol_smiles".

    Returns:
        Dict with keys:
            route (str): The route used.
            original_smiles (str): Canonical SMILES of the original molecule.
            intermediate (str | None): InChI string or MOL block, or None on error.
            roundtrip_smiles (str | None): Canonical SMILES after round-trip, or None on error.
            lossy (bool): True if any information was lost.
            losses (list[dict]): List of detected losses.
            error (str | None): Error message if the route failed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return _make_error_result(route, smiles, "Input SMILES could not be parsed by RDKit")

    original_canonical = Chem.MolToSmiles(mol, canonical=True)

    route_fn = _ROUTE_FUNCTIONS.get(route)
    if route_fn is None:
        return _make_error_result(
            route, original_canonical,
            f"Unknown route '{route}'. Supported: smiles_inchi_smiles, smiles_mol_smiles",
        )

    try:
        intermediate, roundtrip_mol, error = route_fn(mol)
    except Exception as exc:
        return _make_error_result(route, original_canonical, str(exc))

    if error is not None:
        return _make_error_result(route, original_canonical, error, intermediate=intermediate)

    roundtrip_smiles = Chem.MolToSmiles(roundtrip_mol, canonical=True)
    losses = _compute_losses(mol, roundtrip_mol)
    return {
        "route": route,
        "original_smiles": original_canonical,
        "intermediate": intermediate,
        "roundtrip_smiles": roundtrip_smiles,
        "lossy": len(losses) > 0,
        "losses": losses,
        "error": None,
    }
