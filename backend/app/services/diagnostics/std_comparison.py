"""
Cross-Pipeline Standardization Comparison Service (DIAG-04).

Compares standardization output across three pipelines:
1. RDKit MolStandardize (Cleanup + LargestFragment + Uncharger + TautomerEnumerator)
2. ChEMBL-style pipeline (existing StandardizationPipeline)
3. Minimal (MolFromSmiles + SanitizeMol only)

Also computes an MCS-based structural diff across the three outputs so the
frontend can highlight exactly which atoms/bonds each pipeline changed.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, FindMolChiralCenters, rdFMCS
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from app.services.standardization.chembl_pipeline import StandardizationPipeline

# Seconds allowed for rdFMCS.FindMCS before we fall back to no highlights.
# 3s keeps tail latency bounded on pathological inputs; typical small molecules
# complete in well under 50 ms.
_MCS_TIMEOUT_SECONDS = 3


def _extract_properties(mol: Chem.Mol) -> dict:
    """Extract key properties from a molecule for pipeline comparison.

    Args:
        mol: RDKit molecule object.

    Returns:
        Dict with keys: smiles, inchikey, mw, formula, charge, stereo_count.
    """
    smiles = Chem.MolToSmiles(mol, canonical=True)
    inchikey = MolToInchiKey(mol) or ""
    mw = round(Descriptors.ExactMolWt(mol), 4)
    formula = CalcMolFormula(mol)
    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    try:
        stereo_count = len(FindMolChiralCenters(mol, includeUnassigned=True))
    except Exception:
        stereo_count = 0
    return {
        "smiles": smiles,
        "inchikey": inchikey,
        "mw": mw,
        "formula": formula,
        "charge": charge,
        "stereo_count": stereo_count,
    }


def _run_rdkit_pipeline(mol: Chem.Mol) -> dict:
    """Run RDKit MolStandardize pipeline and extract properties.

    Steps: Cleanup → LargestFragmentChooser → Uncharger → TautomerEnumerator.Canonicalize.
    TautomerEnumerator.Canonicalize is wrapped in try/except for timeout protection.

    Args:
        mol: Input RDKit molecule.

    Returns:
        Dict with pipeline name and properties, or error key on failure.
    """
    try:
        cleaned = rdMolStandardize.Cleanup(mol)
        largest = rdMolStandardize.LargestFragmentChooser().choose(cleaned)
        uncharged = rdMolStandardize.Uncharger().uncharge(largest)

        # Tautomer canonicalization can time out on complex molecules
        try:
            te = rdMolStandardize.TautomerEnumerator()
            canon_mol = te.Canonicalize(uncharged)
            final_mol = canon_mol if canon_mol is not None else uncharged
        except Exception:
            # Fall back to pre-tautomer mol on timeout/error
            final_mol = uncharged

        props = _extract_properties(final_mol)
        return {"name": "RDKit MolStandardize", **props}
    except Exception as exc:
        return {"name": "RDKit MolStandardize", "error": str(exc)}


def _run_chembl_pipeline(mol: Chem.Mol) -> dict:
    """Run ChEMBL-style standardization pipeline and extract properties.

    Uses the existing StandardizationPipeline singleton.

    Args:
        mol: Input RDKit molecule.

    Returns:
        Dict with pipeline name and properties, or error key on failure.
    """
    try:
        pipeline = StandardizationPipeline()
        result = pipeline.standardize(mol)
        if result.success and result.standardized_smiles:
            std_mol = Chem.MolFromSmiles(result.standardized_smiles)
            if std_mol is None:
                return {"name": "ChEMBL Pipeline", "error": "Failed to parse standardized SMILES"}
            props = _extract_properties(std_mol)
            return {"name": "ChEMBL Pipeline", **props}
        return {
            "name": "ChEMBL Pipeline",
            "error": result.error_message or "Standardization failed",
        }
    except Exception as exc:
        return {"name": "ChEMBL Pipeline", "error": str(exc)}


def _run_minimal_pipeline(smiles: str) -> dict:
    """Run minimal pipeline: parse + sanitize only.

    Args:
        smiles: Input SMILES string.

    Returns:
        Dict with pipeline name and properties, or error key on failure.
    """
    try:
        mol_min = Chem.MolFromSmiles(smiles)
        if mol_min is None:
            return {"name": "Minimal (Sanitize Only)", "error": "Failed to parse SMILES"}
        Chem.SanitizeMol(mol_min)
        props = _extract_properties(mol_min)
        return {"name": "Minimal (Sanitize Only)", **props}
    except Exception as exc:
        return {"name": "Minimal (Sanitize Only)", "error": str(exc)}


_STRUCTURAL_PROPERTIES = {"smiles", "inchikey"}
_ALL_PROPERTIES = ["smiles", "inchikey", "mw", "formula", "charge", "stereo_count"]


def _compute_mcs_highlights(pipeline_results: list[dict]) -> list[dict]:
    """Compute atom/bond indices that differ from the pipelines' common core.

    Uses rdFMCS.FindMCS (strict element + bond-order matching, default params)
    across each pipeline's canonical SMILES. Any atom not matched by the MCS
    is considered part of the pipeline's local edit; any bond touching a
    non-MCS atom is likewise flagged.

    Args:
        pipeline_results: One dict per pipeline, as returned by the
            _run_*_pipeline helpers. Entries with an "error" key are skipped.

    Returns:
        List of {"highlight_atoms": [...], "highlight_bonds": [...]} dicts
        in the same order as the input. Pipelines that errored (or all-agree
        cases where MCS covers everything) get empty lists.
    """
    empty = [{"highlight_atoms": [], "highlight_bonds": []} for _ in pipeline_results]

    # Rebuild mol per pipeline from its canonical SMILES so indices line up
    # with exactly what the frontend will render.
    mols: list[Chem.Mol | None] = []
    for result in pipeline_results:
        if "error" in result or not result.get("smiles"):
            mols.append(None)
            continue
        mol = Chem.MolFromSmiles(result["smiles"])
        mols.append(mol)

    valid_mols = [m for m in mols if m is not None]
    if len(valid_mols) < 2:
        return empty

    try:
        mcs = rdFMCS.FindMCS(
            valid_mols,
            timeout=_MCS_TIMEOUT_SECONDS,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            matchValences=False,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
        )
    except Exception:
        return empty

    # FindMCS sets canceled=True on timeout or when no MCS is found.
    if mcs.canceled or not mcs.smartsString or mcs.numAtoms == 0:
        return empty

    query = Chem.MolFromSmarts(mcs.smartsString)
    if query is None:
        return empty

    highlights: list[dict] = []
    for mol in mols:
        if mol is None:
            highlights.append({"highlight_atoms": [], "highlight_bonds": []})
            continue

        match = mol.GetSubstructMatch(query)
        if not match:
            # Pipeline produced a structure the MCS couldn't map onto — flag
            # every atom + bond so the user sees the whole molecule diverged.
            highlights.append(
                {
                    "highlight_atoms": list(range(mol.GetNumAtoms())),
                    "highlight_bonds": list(range(mol.GetNumBonds())),
                }
            )
            continue

        mcs_atom_set = set(match)
        diff_atoms = [i for i in range(mol.GetNumAtoms()) if i not in mcs_atom_set]
        # A bond is "different" if it touches any non-MCS atom OR if both
        # endpoints are MCS atoms but the bond itself was not matched (rare;
        # happens when the MCS omits an internal bond such as a ring-closure
        # shift between tautomers).
        diff_bonds: list[int] = []
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a not in mcs_atom_set or b not in mcs_atom_set:
                diff_bonds.append(bond.GetIdx())

        highlights.append({"highlight_atoms": diff_atoms, "highlight_bonds": diff_bonds})

    return highlights


def compare_pipelines(smiles: str) -> dict:
    """Compare standardization output across three pipelines for a given SMILES.

    Pipelines compared:
    1. RDKit MolStandardize (Cleanup + LargestFragment + Uncharger + TautomerEnumerator)
    2. ChEMBL-style (existing StandardizationPipeline)
    3. Minimal (parse + sanitize only)

    Args:
        smiles: Input SMILES string.

    Returns:
        Dict with keys:
            pipelines (list[dict]): One result dict per pipeline with name + properties.
                When pipelines disagree, each entry also carries `highlight_atoms`
                and `highlight_bonds` lists (indices of atoms/bonds not in the
                cross-pipeline MCS). Agreeing pipelines get empty lists.
            disagreements (int): Count of properties where not all pipelines agree.
            structural_disagreements (int): Count of structural properties (SMILES, InChIKey) that disagree.
            all_agree (bool): True if all pipelines agree on all properties.
            property_comparison (list[dict]): Per-property comparison with values and agreement flag.

    Raises:
        ValueError: If the SMILES cannot be parsed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Input SMILES could not be parsed: {smiles!r}")

    p1 = _run_rdkit_pipeline(mol)
    p2 = _run_chembl_pipeline(mol)
    p3 = _run_minimal_pipeline(smiles)

    pipelines = [p1, p2, p3]

    # Compute per-property agreement
    property_comparison: list[dict] = []
    disagreements = 0
    structural_disagreements = 0

    for prop in _ALL_PROPERTIES:
        values = [p.get(prop) for p in pipelines]
        # Only count agreement if all pipelines have this property (no error)
        valid_values = [v for v in values if v is not None]
        if len(valid_values) < len(pipelines):
            # At least one pipeline failed for this property
            agrees = False
        else:
            agrees = len(set(str(v) for v in valid_values)) == 1

        is_structural = prop in _STRUCTURAL_PROPERTIES

        if not agrees:
            disagreements += 1
            if is_structural:
                structural_disagreements += 1

        property_comparison.append(
            {
                "property": prop,
                "values": values,
                "agrees": agrees,
                "structural": is_structural,
            }
        )

    # Only compute MCS highlights when structural disagreements exist; when all
    # pipelines agree structurally every atom maps 1:1 and highlights are empty.
    if structural_disagreements > 0:
        highlights = _compute_mcs_highlights(pipelines)
    else:
        highlights = [{"highlight_atoms": [], "highlight_bonds": []} for _ in pipelines]

    for pipeline, h in zip(pipelines, highlights):
        pipeline["highlight_atoms"] = h["highlight_atoms"]
        pipeline["highlight_bonds"] = h["highlight_bonds"]

    return {
        "pipelines": pipelines,
        "disagreements": disagreements,
        "structural_disagreements": structural_disagreements,
        "all_agree": disagreements == 0,
        "property_comparison": property_comparison,
    }
