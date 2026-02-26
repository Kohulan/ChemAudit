"""
Multi-Level Deduplication Service

Detects duplicate molecules across 4 comparison levels:
  - exact: identical canonical SMILES
  - tautomeric: identical canonical tautomer SMILES
  - stereo_insensitive: identical stereo-stripped InChI
  - salt_form: identical parent-compound SMILES after desalting

Each function is self-contained and can be called independently.
The aggregator ``compute_all_dedup_levels`` calls all four and returns a
dict matching the ``DeduplicationResult`` schema.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from app.schemas.analytics import DeduplicationResult

from rdkit import Chem

logger = logging.getLogger(__name__)

# InChI stereo layers to strip when comparing stereo-insensitive InChIs.
# Standard InChI layers: /t (tetrahedral), /m (config), /s (stereo type)
_STEREO_LAYERS: frozenset[str] = frozenset(["t", "m", "s"])


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _iter_valid(results: list[dict[str, Any]]):
    """
    Yield (original_index, mol) for every result that is a successful parse.

    Skips results where ``status != "success"`` or where the SMILES cannot be
    parsed by RDKit.

    Args:
        results: Raw batch result dicts containing at least ``index``, ``status``,
                 and ``smiles`` keys.

    Yields:
        Tuple of (int index, rdkit.Chem.Mol molecule).
    """
    for result in results:
        if result.get("status") != "success":
            continue
        smiles = result.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        yield result["index"], mol


def _build_groups(level: str, key_index_pairs: list[tuple[str, int]]) -> list[dict]:
    """
    Group (key, index) pairs and return dedup group dicts for groups with >1 member.

    The representative is always ``min(indices)`` (first submitted), not
    ``indices[0]``, because dict insertion order may not match submission
    order when chunks are processed in parallel.

    Args:
        level: Dedup level name (e.g. ``"exact"``).
        key_index_pairs: List of (group_key, original_index) tuples.

    Returns:
        List of dedup-group dicts (only groups with >=2 members).
    """
    buckets: dict[str, list[int]] = defaultdict(list)
    for key, idx in key_index_pairs:
        buckets[key].append(idx)

    groups: list[dict] = []
    for key, indices in buckets.items():
        if len(indices) < 2:
            continue
        rep = min(indices)
        dups = sorted(i for i in indices if i != rep)
        groups.append(
            {
                "level": level,
                "representative_index": rep,
                "duplicate_indices": dups,
                "group_key": key,
                "count": len(indices),
            }
        )
    return groups


# ---------------------------------------------------------------------------
# Public dedup functions
# ---------------------------------------------------------------------------


def compute_exact_dedup(results: list[dict[str, Any]]) -> list[dict]:
    """
    Detect exact structural duplicates by canonical SMILES.

    Two molecules are exact duplicates when RDKit generates the same canonical
    SMILES string for both, regardless of how the SMILES was originally written.

    Args:
        results: Raw batch result dicts (must contain ``index``, ``status``, ``smiles``).

    Returns:
        List of dedup-group dicts with ``level="exact"``.
    """
    pairs: list[tuple[str, int]] = []
    for idx, mol in _iter_valid(results):
        key = Chem.MolToSmiles(mol)
        pairs.append((key, idx))
    return _build_groups("exact", pairs)


def compute_tautomer_dedup(results: list[dict[str, Any]]) -> list[dict]:
    """
    Detect tautomeric duplicates by canonical tautomer SMILES.

    TautomerEnumerator is instantiated locally per call (not at module level)
    because it is not thread-safe (research Pitfall 7).  If tautomer
    canonicalization raises for an individual molecule, the function falls back
    to plain canonical SMILES for that molecule.

    Args:
        results: Raw batch result dicts.

    Returns:
        List of dedup-group dicts with ``level="tautomeric"``.
    """
    from rdkit.Chem.MolStandardize import rdMolStandardize  # noqa: PLC0415

    enumerator = rdMolStandardize.TautomerEnumerator()

    pairs: list[tuple[str, int]] = []
    for idx, mol in _iter_valid(results):
        try:
            canonical_taut = enumerator.Canonicalize(mol)
            key = Chem.MolToSmiles(canonical_taut)
        except Exception:
            logger.debug(
                "compute_tautomer_dedup: tautomer canonicalization failed for index %d;"
                " falling back to canonical SMILES",
                idx,
            )
            key = Chem.MolToSmiles(mol)
        pairs.append((key, idx))
    return _build_groups("tautomeric", pairs)


def compute_stereo_dedup(results: list[dict[str, Any]]) -> list[dict]:
    """
    Detect stereo-insensitive duplicates by stripping stereo layers from InChI.

    Stereo InChI layers /t, /m, /s are removed so that enantiomers and
    diastereomers map to the same key.  Molecules for which InChI generation
    fails are skipped entirely.

    Args:
        results: Raw batch result dicts.

    Returns:
        List of dedup-group dicts with ``level="stereo_insensitive"``.
    """
    from rdkit.Chem.inchi import MolToInchi  # noqa: PLC0415

    pairs: list[tuple[str, int]] = []
    for idx, mol in _iter_valid(results):
        inchi = MolToInchi(mol)
        if not inchi:
            logger.debug(
                "compute_stereo_dedup: InChI generation failed for index %d; skipping",
                idx,
            )
            continue
        # Strip stereo layers: split by '/', remove layers whose first char is
        # a stereo layer identifier, then rejoin.
        parts = inchi.split("/")
        stripped_parts = [
            p for p in parts if not (len(p) > 0 and p[0].lower() in _STEREO_LAYERS)
        ]
        key = "/".join(stripped_parts)
        pairs.append((key, idx))
    return _build_groups("stereo_insensitive", pairs)


def compute_saltform_dedup(results: list[dict[str, Any]]) -> list[dict]:
    """
    Detect salt-form duplicates by comparing canonical SMILES of the parent compound.

    For each molecule the parent compound is obtained via
    ``chembl_structure_pipeline.standardizer.get_parent_mol``.  If pre-computed
    standardization data is available in the result dict
    (``result["standardization"]["standardized_smiles"]``), that SMILES is
    used as a fast path to skip redundant standardization calls.  If parent
    extraction fails for a molecule, the function falls back to the original
    canonical SMILES for that molecule.

    Args:
        results: Raw batch result dicts.

    Returns:
        List of dedup-group dicts with ``level="salt_form"``.
    """
    from chembl_structure_pipeline.standardizer import get_parent_mol  # noqa: PLC0415

    pairs: list[tuple[str, int]] = []
    for idx, mol in _iter_valid(results):
        # Fast path: use already-standardized SMILES if present.
        std_data = (results[idx].get("standardization") or {}) if False else {}
        # Resolve actual result dict by index since _iter_valid yields the
        # original index value, not list position.
        for r in results:
            if r.get("index") == idx and r.get("status") == "success":
                std_data = r.get("standardization") or {}
                break

        pre_smiles = std_data.get("standardized_smiles")
        if pre_smiles:
            std_mol = Chem.MolFromSmiles(pre_smiles)
            if std_mol is not None:
                mol = std_mol  # noqa: PLW2901

        try:
            parent, _exclude = get_parent_mol(mol)
            key = Chem.MolToSmiles(parent)
        except Exception:
            logger.debug(
                "compute_saltform_dedup: parent extraction failed for index %d;"
                " falling back to canonical SMILES",
                idx,
            )
            key = Chem.MolToSmiles(mol)
        pairs.append((key, idx))
    return _build_groups("salt_form", pairs)


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------


def compute_all_dedup_levels(results: list[dict[str, Any]]) -> "DeduplicationResult":
    """
    Run all 4 deduplication levels and return a ``DeduplicationResult`` schema object.

    Total unique count per level is:
        (number of successful molecules) - (sum of duplicate counts across groups)

    Because each duplicate molecule beyond the representative is counted once,
    the formula is:
        unique = total_success - sum(group["count"] - 1 for group in level_groups)

    Args:
        results: Raw batch result dicts for the entire batch job.

    Returns:
        ``DeduplicationResult`` Pydantic model with keys:
        ``exact``, ``tautomeric``, ``stereo_insensitive``, ``salt_form``,
        ``total_unique``.
    """
    from app.schemas.analytics import DeduplicationResult  # noqa: PLC0415

    total_success = sum(1 for r in results if r.get("status") == "success")

    exact_groups = compute_exact_dedup(results)
    tautomeric_groups = compute_tautomer_dedup(results)
    stereo_groups = compute_stereo_dedup(results)
    saltform_groups = compute_saltform_dedup(results)

    def _unique_count(groups: list[dict]) -> int:
        duplicates_removed = sum(g["count"] - 1 for g in groups)
        return max(0, total_success - duplicates_removed)

    return DeduplicationResult(
        exact=exact_groups,
        tautomeric=tautomeric_groups,
        stereo_insensitive=stereo_groups,
        salt_form=saltform_groups,
        total_unique={
            "exact": _unique_count(exact_groups),
            "tautomeric": _unique_count(tautomeric_groups),
            "stereo_insensitive": _unique_count(stereo_groups),
            "salt_form": _unique_count(saltform_groups),
        },
    )
