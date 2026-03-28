"""
Butina Compound Clustering Service

Clusters molecules using the Butina sphere-exclusion algorithm with configurable
Tanimoto distance cutoff. Uses Morgan fingerprints (radius=2, 2048 bits) via the
non-deprecated rdFingerprintGenerator API.

Hard-capped at 1,000 molecules (Phase 13 decision D-03).
"""

from __future__ import annotations

import logging
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

logger = logging.getLogger(__name__)

# Hard cap to prevent combinatorial explosion in distance matrix computation
MAX_MOLECULES = 1000

# Morgan fingerprint generator (module-level singleton for reuse)
_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def compute_butina_clustering(
    results: list[dict[str, Any]],
    distance_cutoff: float = 0.35,
) -> dict:
    """
    Cluster molecules using Butina sphere-exclusion algorithm.

    Extracts SMILES from batch results, generates Morgan fingerprints,
    computes pairwise Tanimoto distances, and clusters using Butina.

    Args:
        results: Batch result dicts, each containing at least ``smiles`` key.
        distance_cutoff: Tanimoto distance threshold for clustering.
            Lower values (e.g. 0.2) produce fewer, larger clusters.
            Higher values (e.g. 0.6) produce more, smaller clusters.
            Default 0.35.

    Returns:
        Dict with keys:
        - ``clusters``: list of cluster dicts with ``cluster_id``,
          ``member_indices`` (original batch indices), ``size``,
          ``representative_index``
        - ``cluster_count``: number of clusters
        - ``singleton_count``: number of clusters with exactly 1 member
        - ``largest_cluster_size``: size of largest cluster
        - ``distance_cutoff``: the cutoff used
    """
    # Step 1: Extract valid molecules and build index mapping
    valid_mols: list[Chem.Mol] = []
    valid_indices: list[int] = []  # valid_indices[i] = original batch index

    for result in results[:MAX_MOLECULES]:
        smiles = result.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        valid_mols.append(mol)
        valid_indices.append(result.get("index", 0))

    n = len(valid_mols)

    # Edge case: no valid molecules
    if n == 0:
        return {
            "clusters": [],
            "cluster_count": 0,
            "singleton_count": 0,
            "largest_cluster_size": 0,
            "distance_cutoff": distance_cutoff,
        }

    # Edge case: single molecule -> 1 cluster, no singletons
    if n == 1:
        return {
            "clusters": [
                {
                    "cluster_id": 0,
                    "member_indices": [valid_indices[0]],
                    "size": 1,
                    "representative_index": valid_indices[0],
                }
            ],
            "cluster_count": 1,
            "singleton_count": 0,
            "largest_cluster_size": 1,
            "distance_cutoff": distance_cutoff,
        }

    # Step 2: Generate Morgan fingerprints
    fps = [_MORGAN_GEN.GetFingerprint(mol) for mol in valid_mols]

    # Step 3: Compute lower-triangle Tanimoto distance list
    dists: list[float] = []
    for i in range(1, n):
        for j in range(0, i):
            dists.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))

    # Step 4: Butina clustering
    cluster_tuples = Butina.ClusterData(dists, n, distThresh=distance_cutoff, isDistData=True)

    # Step 5: Build output, mapping internal indices back to original batch indices
    clusters: list[dict] = []
    largest_size = 0
    singleton_count = 0

    for cid, members in enumerate(cluster_tuples):
        original_indices = [valid_indices[m] for m in members]
        size = len(original_indices)

        if size == 1:
            singleton_count += 1

        if size > largest_size:
            largest_size = size

        clusters.append(
            {
                "cluster_id": cid,
                "member_indices": original_indices,
                "size": size,
                "representative_index": original_indices[0],  # Butina convention: first = centroid
            }
        )

    return {
        "clusters": clusters,
        "cluster_count": len(clusters),
        "singleton_count": singleton_count,
        "largest_cluster_size": largest_size,
        "distance_cutoff": distance_cutoff,
    }
