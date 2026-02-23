"""
Matched Molecular Pair (MMP), Activity Cliff, and Lipophilic Ligand Efficiency Analytics

Background
----------
Matched Molecular Pair (MMP) Analysis:
    Two molecules form a matched molecular pair (MMP) when they differ by exactly one
    structural transformation at a single site. BRICS (Breaking Retrosynthetically
    Interesting Chemical Substructures) fragmentation is used to identify these pairs:
    each molecule is decomposed into BRICS fragments, and two molecules share a common
    core if they share at least one fragment (the core) while differing in another
    (the R-group).

Activity Cliff Detection (SALI index):
    Activity cliffs are pairs of structurally similar molecules with a large difference
    in biological activity. The Structure-Activity Landscape Index (SALI) quantifies this:

        SALI = |delta_activity| / (1 - Tanimoto_similarity)

    High SALI values indicate a dramatic activity change from a small structural change.
    The denominator approaches 0 as molecules become more similar, amplifying the
    activity signal for closely related pairs.

Lipophilic Ligand Efficiency (LLE):
    LLE measures drug-likeness by balancing potency against lipophilicity:

        LLE = pIC50 - LogP

    Higher LLE values indicate better potency per unit of lipophilicity, which correlates
    with lower risk of off-target effects and ADMET liabilities.
"""

import logging
import time
from collections import defaultdict

from rdkit import Chem, DataStructs
from rdkit.Chem import BRICS, Descriptors, rdFingerprintGenerator

logger = logging.getLogger(__name__)

# Morgan fingerprint generator consistent with chemical_space.py
_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Batch size cap: MMP is O(n^2) in the worst case; refuse very large batches
MAX_MMP_BATCH_SIZE = 5000

# Maximum number of MMP pairs returned in response (sorted by Tanimoto descending)
MAX_PAIRS_RETURNED = 1000


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _detect_mmp_pairs(mols: list[Chem.Mol], indices: list[int]) -> list[dict]:
    """
    Detect matched molecular pairs using BRICS single-cut fragmentation.

    For each molecule, compute its BRICS fragments. Two molecules form an MMP
    when they share at least one fragment (common core) while having different
    remaining fragments (different R-groups). A size heuristic filters noise:
    the shared (core) fragment must have more heavy atoms than the variable
    (R-group) fragment.

    Pairs are deduplicated so each (mol_a, mol_b) appears at most once. If
    multiple shared cores exist for the same pair, the entry with the highest
    Tanimoto similarity is retained. Output is capped at MAX_PAIRS_RETURNED
    pairs sorted by Tanimoto descending.

    Args:
        mols: List of RDKit molecule objects (all must be non-None).
        indices: Original batch indices parallel to ``mols``.

    Returns:
        List of dicts matching the MMPPair schema, each with keys:
        mol_a_index, mol_b_index, core_smiles, rgroup_a, rgroup_b, tanimoto.
    """
    t0 = time.perf_counter()

    # Compute BRICS fragments and fingerprints for all molecules
    mol_frags: list[set[str]] = []
    mol_fps = []
    for mol in mols:
        frags = BRICS.BRICSDecompose(mol)
        mol_frags.append(frags)
        mol_fps.append(_morgan_gen.GetFingerprint(mol))

    # Build fragment -> [(position_in_mols_list, remaining_frags), ...] mapping
    # Only index molecules that actually decomposed (len(frags) > 1)
    frag_to_mols: dict[str, list[tuple[int, set[str]]]] = defaultdict(list)
    for pos, (frags, idx) in enumerate(zip(mol_frags, indices)):
        if len(frags) < 2:
            # No BRICS bonds found — molecule cannot be part of an MMP via BRICS
            continue
        for frag in frags:
            frag_mol = Chem.MolFromSmiles(frag)
            if frag_mol is None:
                continue
            remaining = frags - {frag}
            frag_to_mols[frag].append((pos, remaining))

    # For each shared fragment (potential core), build candidate MMP pairs
    # Key: (pos_a, pos_b) with pos_a < pos_b
    # Value: dict with best core and Tanimoto so far
    best_pairs: dict[tuple[int, int], dict] = {}

    for core_smi, mol_list in frag_to_mols.items():
        if len(mol_list) < 2:
            continue

        # Parse core fragment once for heavy-atom count
        core_frag_mol = Chem.MolFromSmiles(core_smi)
        if core_frag_mol is None:
            continue
        core_heavy_atoms = core_frag_mol.GetNumHeavyAtoms()

        # All pairs within this core group
        for i in range(len(mol_list)):
            pos_a, remaining_a = mol_list[i]
            idx_a = indices[pos_a]

            # Check size heuristic for mol_a's R-group(s)
            rgroup_a_heavy = sum(
                Chem.MolFromSmiles(r).GetNumHeavyAtoms()
                for r in remaining_a
                if Chem.MolFromSmiles(r) is not None
            )
            if core_heavy_atoms <= rgroup_a_heavy:
                continue

            for j in range(i + 1, len(mol_list)):
                pos_b, remaining_b = mol_list[j]
                idx_b = indices[pos_b]

                # Check size heuristic for mol_b's R-group(s)
                rgroup_b_heavy = sum(
                    Chem.MolFromSmiles(r).GetNumHeavyAtoms()
                    for r in remaining_b
                    if Chem.MolFromSmiles(r) is not None
                )
                if core_heavy_atoms <= rgroup_b_heavy:
                    continue

                # Tanimoto similarity between the two full molecules
                tanimoto = DataStructs.TanimotoSimilarity(
                    mol_fps[pos_a], mol_fps[pos_b]
                )

                # Normalise pair key so mol_a < mol_b (by original batch index)
                if idx_a < idx_b:
                    pair_key = (idx_a, idx_b)
                    rg_a = ".".join(sorted(remaining_a))
                    rg_b = ".".join(sorted(remaining_b))
                else:
                    pair_key = (idx_b, idx_a)
                    rg_a = ".".join(sorted(remaining_b))
                    rg_b = ".".join(sorted(remaining_a))

                # Keep the entry with the highest Tanimoto for this pair
                if pair_key not in best_pairs or tanimoto > best_pairs[pair_key]["tanimoto"]:
                    best_pairs[pair_key] = {
                        "mol_a_index": pair_key[0],
                        "mol_b_index": pair_key[1],
                        "core_smiles": core_smi,
                        "rgroup_a": rg_a,
                        "rgroup_b": rg_b,
                        "tanimoto": tanimoto,
                    }

    # Sort by Tanimoto descending and cap output
    pairs = sorted(best_pairs.values(), key=lambda p: p["tanimoto"], reverse=True)
    pairs = pairs[:MAX_PAIRS_RETURNED]

    elapsed = time.perf_counter() - t0
    logger.info(
        "MMP pair detection: %d molecules → %d pairs in %.2fs",
        len(mols),
        len(pairs),
        elapsed,
    )
    return pairs


def _compute_activity_cliffs(
    pairs: list[dict], activities: dict[int, float]
) -> list[dict]:
    """
    Compute SALI-based activity cliffs for MMP pairs.

    For each MMP pair where both molecules have known activity values, compute:

        SALI = |delta_activity| / (1 - Tanimoto)

    Pairs where Tanimoto >= 1.0 are skipped (denominator is zero or negative).
    Results are sorted by SALI descending.

    Args:
        pairs: List of MMP pair dicts (from ``_detect_mmp_pairs``).
        activities: Mapping of original batch index → activity value (e.g., pIC50).

    Returns:
        List of dicts matching the ActivityCliff schema, sorted by SALI descending.
    """
    cliffs = []
    for pair in pairs:
        idx_a = pair["mol_a_index"]
        idx_b = pair["mol_b_index"]

        act_a = activities.get(idx_a)
        act_b = activities.get(idx_b)
        if act_a is None or act_b is None:
            continue

        tanimoto = pair["tanimoto"]
        if tanimoto >= 1.0:
            # Perfect similarity — SALI undefined (infinite); skip
            continue

        activity_diff = abs(act_a - act_b)
        sali = activity_diff / (1.0 - tanimoto)

        cliffs.append(
            {
                "mol_a_index": idx_a,
                "mol_b_index": idx_b,
                "sali": sali,
                "tanimoto": tanimoto,
                "activity_diff": activity_diff,
            }
        )

    cliffs.sort(key=lambda c: c["sali"], reverse=True)
    return cliffs


def _compute_lle(results: list[dict], activity_column: str) -> list[dict]:
    """
    Compute Lipophilic Ligand Efficiency (LLE) for each molecule with activity data.

    LLE = activity_value - LogP

    Assumes activity_value is on a -log scale (e.g., pIC50). Higher LLE values
    indicate better potency per unit of lipophilicity.

    Args:
        results: Batch validation result dicts. Each must have a ``smiles`` key
                 for the input SMILES and a ``properties`` dict for activity lookup.
        activity_column: Key in each result's ``properties`` dict containing the
                         activity value.

    Returns:
        List of dicts with keys: molecule_index, activity, logp, lle.
    """
    lle_values = []
    for res in results:
        idx = res.get("molecule_index")
        props = res.get("properties") or {}
        raw_activity = props.get(activity_column)
        if raw_activity is None:
            continue

        try:
            activity_val = float(raw_activity)
        except (ValueError, TypeError):
            continue

        smiles = res.get("standardized_smiles") or res.get("smiles")
        if not smiles:
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        logp = Descriptors.MolLogP(mol)
        lle = activity_val - logp

        lle_values.append(
            {
                "molecule_index": idx,
                "activity": activity_val,
                "logp": logp,
                "lle": lle,
            }
        )

    return lle_values


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_mmp_analysis(
    results: list[dict],
    activity_column: str | None = None,
) -> dict:
    """
    Run full MMP analysis: pair detection, activity cliff scoring, and LLE computation.

    Workflow:
    1.  Filter results to successful molecules only.
    2.  Enforce batch size limit (MAX_MMP_BATCH_SIZE = 5000).
    3.  Detect MMP pairs via BRICS single-cut fragmentation.
    4.  If ``activity_column`` is provided, extract activity values from each
        result's ``properties`` dict, compute activity cliffs (SALI index), and
        compute LLE for each molecule with valid activity data.
    5.  If ``activity_column`` is absent, skip activity cliff and LLE computation
        and return None for those fields.

    Args:
        results: List of batch validation result dicts. Each must have at minimum:
                 - ``molecule_index`` (int): position in the original batch
                 - ``status`` (str): ``"success"`` for processable molecules
                 - ``smiles`` or ``standardized_smiles`` (str): SMILES string
                 - ``properties`` (dict, optional): extra data including activity values
        activity_column: Name of the key in each result's ``properties`` dict that
                         holds the numeric activity value (e.g., ``"pIC50"``). If
                         ``None``, activity cliffs and LLE are skipped.

    Returns:
        Dict matching the MMPResult schema:
        - ``pairs``: list of MMPPair dicts (may be empty)
        - ``activity_cliffs``: list of ActivityCliff dicts or ``None``
        - ``lle_values``: list of LLE dicts or ``None``

        On refusal (batch too large):
        - ``status``: ``"refused"``
        - ``reason``: human-readable message
        - ``molecule_count``: number of successful molecules in the batch
    """
    t0 = time.perf_counter()

    # Filter to successful results with parseable SMILES
    successful = [r for r in results if r.get("status") == "success"]
    n = len(successful)

    if n > MAX_MMP_BATCH_SIZE:
        logger.warning(
            "compute_mmp_analysis: batch of %d molecules exceeds limit of %d — refusing",
            n,
            MAX_MMP_BATCH_SIZE,
        )
        return {
            "status": "refused",
            "reason": f"MMP detection limited to {MAX_MMP_BATCH_SIZE} molecules",
            "molecule_count": n,
        }

    if n == 0:
        logger.info("compute_mmp_analysis: no successful molecules in batch")
        return {
            "pairs": [],
            "activity_cliffs": None,
            "lle_values": None,
        }

    # Parse molecules
    mols: list[Chem.Mol] = []
    indices: list[int] = []
    valid_results: list[dict] = []

    for res in successful:
        smiles = res.get("standardized_smiles") or res.get("smiles")
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        idx = res.get("molecule_index")
        if idx is None:
            continue
        mols.append(mol)
        indices.append(idx)
        valid_results.append(res)

    if not mols:
        return {
            "pairs": [],
            "activity_cliffs": None,
            "lle_values": None,
        }

    # --- MMP pair detection ---
    pairs = _detect_mmp_pairs(mols, indices)

    # --- Activity data (optional) ---
    activity_cliffs = None
    lle_values = None
    note = None

    if activity_column is not None:
        # Extract activity values from properties
        activities: dict[int, float] = {}
        for res in valid_results:
            idx = res.get("molecule_index")
            props = res.get("properties") or {}
            raw = props.get(activity_column)
            if raw is None:
                continue
            try:
                activities[idx] = float(raw)
            except (ValueError, TypeError):
                continue

        if not activities:
            note = f"No valid activity values found in column '{activity_column}'"
            logger.info("compute_mmp_analysis: %s (job has %d molecules)", note, n)
        else:
            activity_cliffs = _compute_activity_cliffs(pairs, activities)
            lle_values = _compute_lle(valid_results, activity_column)

    elapsed = time.perf_counter() - t0
    logger.info(
        "compute_mmp_analysis: %d mols → %d pairs, %s cliffs, %s LLE in %.2fs",
        n,
        len(pairs),
        len(activity_cliffs) if activity_cliffs is not None else "N/A",
        len(lle_values) if lle_values is not None else "N/A",
        elapsed,
    )

    result: dict = {
        "pairs": pairs,
        "activity_cliffs": activity_cliffs,
        "lle_values": lle_values,
    }
    if note:
        result["note"] = note

    return result


# Alias for compatibility with analytics_tasks.py which imports `compute_mmp`
def compute_mmp(
    results: list[dict],
    activity_column: str | None = None,
) -> object:
    """
    Alias for ``compute_mmp_analysis`` that returns an MMPResult-compatible object.

    analytics_tasks.py calls ``result.model_dump()``, so we wrap the output
    dict in a lightweight namespace that exposes ``model_dump()``.

    Args:
        results: Batch validation result dicts.
        activity_column: Optional activity column name.

    Returns:
        Object with a ``model_dump()`` method returning the MMPResult dict.
    """

    class _MMPResultWrapper:
        """Thin wrapper so analytics_tasks.py can call .model_dump()."""

        def __init__(self, data: dict) -> None:
            self._data = data

        def model_dump(self) -> dict:
            """Return the underlying result dict."""
            return self._data

    data = compute_mmp_analysis(results, activity_column=activity_column)
    return _MMPResultWrapper(data)
