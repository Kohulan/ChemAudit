"""
Chemical Space Analysis Service

Provides PCA and t-SNE dimensionality reduction, similarity search,
nearest neighbor analysis, and pairwise similarity matrix computation
for batch chemical structure analysis.

All fingerprints use Morgan (ECFP4) via rdFingerprintGenerator (non-deprecated API).
PCA uses randomized SVD (pure numpy, no sklearn).
t-SNE uses openTSNE; batches > 2000 are refused.
"""

import logging
import time
from typing import Optional

import numpy as np
from rdkit import DataStructs
from rdkit.Chem import MolFromSmiles, rdFingerprintGenerator

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Module-level Morgan generator (non-deprecated API)
# ---------------------------------------------------------------------------

_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _parse_and_fingerprint(results: list[dict]) -> tuple[list[int], list]:
    """
    Parse SMILES from batch results and generate Morgan fingerprints.

    Skips results with status != "success" or unparseable SMILES.

    Args:
        results: List of validation result dicts, each expected to have
                 at least 'index', 'status', and 'smiles' keys.

    Returns:
        Tuple of (molecule_indices, fingerprints) as parallel lists.
    """
    indices: list[int] = []
    fps: list = []

    for item in results:
        if item.get("status") != "success":
            continue
        smiles = item.get("smiles") or item.get("original_smiles", "")
        if not smiles:
            continue
        mol = MolFromSmiles(smiles)
        if mol is None:
            continue
        fp = _morgan_gen.GetFingerprint(mol)
        indices.append(item.get("index", len(indices)))
        fps.append(fp)

    return indices, fps


def _fps_to_matrix(fps: list) -> np.ndarray:
    """
    Convert a list of RDKit fingerprints to a float32 numpy matrix (n x 2048).

    Args:
        fps: List of RDKit ExplicitBitVect fingerprints.

    Returns:
        Float32 numpy array of shape (n, 2048).
    """
    n = len(fps)
    X = np.zeros((n, 2048), dtype=np.float32)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, X[i])
    return X


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_pca(results: list[dict], n_components: int = 2) -> dict:
    """
    Compute randomized SVD PCA coordinates for a batch of molecules.

    Uses pure numpy randomized SVD (no sklearn) for reproducibility and
    memory efficiency with high-dimensional fingerprint data.

    Args:
        results: List of batch validation result dicts.
        n_components: Number of PCA components (default 2).

    Returns:
        Dict matching ChemSpaceCoordinates schema with method="pca",
        or error dict if fewer than 3 molecules can be parsed.
    """
    t0 = time.perf_counter()
    indices, fps = _parse_and_fingerprint(results)

    if len(indices) < 3:
        return {
            "error": f"PCA requires at least 3 molecules; got {len(indices)}",
            "molecule_count": len(indices),
        }

    X = _fps_to_matrix(fps)

    # Randomized SVD PCA (per research Pattern 7)
    rng = np.random.default_rng(42)
    X_centered = X - X.mean(axis=0)

    omega = rng.standard_normal((2048, n_components + 10)).astype(np.float32)
    Q, _ = np.linalg.qr(X_centered @ omega)
    B = Q.T @ X_centered
    U, s, _Vt = np.linalg.svd(B, full_matrices=False)

    coords = (Q @ U[:, :n_components]) * s[:n_components]
    variance_explained = (s[:n_components] ** 2 / (s**2).sum()).tolist()

    elapsed = time.perf_counter() - t0
    logger.info(
        "compute_pca: %d molecules, %d components in %.3fs",
        len(indices),
        n_components,
        elapsed,
    )

    return {
        "method": "pca",
        "coordinates": coords.tolist(),
        "molecule_indices": indices,
        "variance_explained": variance_explained,
    }


def compute_tsne(results: list[dict], perplexity: int = 30) -> dict:
    """
    Compute t-SNE 2-D coordinates using openTSNE.

    Enforces a hard limit of 2000 molecules. Auto-adjusts perplexity to
    avoid the perplexity >= n_samples error for small batches.

    Args:
        results: List of batch validation result dicts.
        perplexity: t-SNE perplexity parameter (default 30).

    Returns:
        Dict matching ChemSpaceCoordinates schema with method="tsne",
        or error dict if more than 2000 molecules are in the batch.
    """
    from openTSNE import TSNE

    indices, fps = _parse_and_fingerprint(results)

    if len(indices) > 2000:
        return {
            "error": "t-SNE limited to 2000 molecules",
            "molecule_count": len(indices),
        }

    if len(indices) < 3:
        return {
            "error": f"t-SNE requires at least 3 molecules; got {len(indices)}",
            "molecule_count": len(indices),
        }

    X = _fps_to_matrix(fps)

    # Auto-adjust perplexity: avoids openTSNE crash when perplexity >= n_samples
    perplexity = min(perplexity, max(5, len(indices) // 3))

    t0 = time.perf_counter()
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=42,
        n_jobs=-1,
    )
    coords = tsne.fit(X)

    elapsed = time.perf_counter() - t0
    logger.info(
        "compute_tsne: %d molecules, perplexity=%d in %.3fs",
        len(indices),
        perplexity,
        elapsed,
    )

    return {
        "method": "tsne",
        "coordinates": np.array(coords).tolist(),
        "molecule_indices": indices,
        "variance_explained": None,
    }


def find_similar_molecules(
    results: list[dict],
    query_smiles: Optional[str] = None,
    query_index: Optional[int] = None,
    top_k: int = 10,
) -> dict:
    """
    Find the top-k most similar molecules in the batch to a query molecule.

    Either ``query_smiles`` (external SMILES) or ``query_index`` (index into
    the batch) must be provided; ``query_smiles`` takes precedence.

    Args:
        results: List of batch validation result dicts.
        query_smiles: Optional external query SMILES string.
        query_index: Optional index of the query molecule in the batch.
        top_k: Number of nearest neighbors to return (default 10).

    Returns:
        Dict matching SimilarityResult schema.
    """
    indices, fps = _parse_and_fingerprint(results)

    if not indices:
        return {"query_index": query_index, "query_smiles": query_smiles, "neighbors": []}

    if query_smiles is not None:
        mol = MolFromSmiles(query_smiles)
        if mol is None:
            return {
                "error": f"Could not parse query SMILES: {query_smiles}",
                "query_index": None,
                "query_smiles": query_smiles,
                "neighbors": [],
            }
        query_fp = _morgan_gen.GetFingerprint(mol)
        query_batch_pos = None  # external query — not in batch
    elif query_index is not None:
        # Find position of query_index in the parsed indices list
        if query_index not in indices:
            return {
                "error": f"query_index {query_index} not found in batch or not a successful result",
                "query_index": query_index,
                "query_smiles": None,
                "neighbors": [],
            }
        query_batch_pos = indices.index(query_index)
        query_fp = fps[query_batch_pos]
        # Resolve the SMILES for the response
        for item in results:
            if item.get("index") == query_index:
                query_smiles = item.get("smiles") or item.get("original_smiles", "")
                break
    else:
        return {
            "error": "Either query_smiles or query_index must be provided",
            "query_index": None,
            "query_smiles": None,
            "neighbors": [],
        }

    sims = DataStructs.BulkTanimotoSimilarity(query_fp, fps)

    # Build (batch_pos, similarity) pairs, exclude the query molecule itself
    sim_pairs = [
        (pos, sim)
        for pos, sim in enumerate(sims)
        if pos != query_batch_pos
    ]
    sim_pairs.sort(key=lambda x: x[1], reverse=True)
    top_pairs = sim_pairs[:top_k]

    # Resolve SMILES for each neighbor from results
    index_to_smiles: dict[int, str] = {}
    for item in results:
        idx = item.get("index")
        if idx is not None:
            index_to_smiles[idx] = item.get("smiles") or item.get("original_smiles", "")

    neighbors = [
        {
            "index": indices[pos],
            "similarity": float(sim),
            "smiles": index_to_smiles.get(indices[pos], ""),
        }
        for pos, sim in top_pairs
    ]

    return {
        "query_index": query_index if query_smiles is None else None,
        "query_smiles": query_smiles,
        "neighbors": neighbors,
    }


def compute_nearest_neighbors(results: list[dict]) -> list[dict]:
    """
    Compute the nearest neighbor and isolation score for every molecule.

    Isolation score = 1.0 - max_neighbor_similarity. A high isolation score
    indicates the molecule is structurally distinct from all others in the batch.

    Args:
        results: List of batch validation result dicts.

    Returns:
        List of dicts matching NearestNeighborResult schema, one per molecule.
    """
    t0 = time.perf_counter()
    indices, fps = _parse_and_fingerprint(results)

    if not indices:
        return []

    nn_results: list[dict] = []
    for i, (mol_idx, fp) in enumerate(zip(indices, fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps)
        # Mask self-similarity
        sims_no_self = [s for j, s in enumerate(sims) if j != i]
        if sims_no_self:
            max_sim = float(max(sims_no_self))
            max_pos = max(
                (j for j in range(len(fps)) if j != i),
                key=lambda j: sims[j],
            )
            nearest_idx = indices[max_pos]
        else:
            max_sim = 0.0
            nearest_idx = mol_idx

        nn_results.append(
            {
                "molecule_index": mol_idx,
                "nearest_index": nearest_idx,
                "similarity": max_sim,
                "isolation_score": 1.0 - max_sim,
            }
        )

    elapsed = time.perf_counter() - t0
    logger.info(
        "compute_nearest_neighbors: %d molecules in %.3fs", len(indices), elapsed
    )
    return nn_results


def compute_similarity_matrix(results: list[dict], threshold: float = 0.7) -> dict:
    """
    Compute a pairwise Tanimoto similarity matrix for the batch.

    Size-based representation:
    - ≤ 500 molecules: dense upper-triangle (list of lists of similarity values).
    - 500-2000 molecules: sparse representation (list of {i, j, similarity} dicts
      where similarity > threshold).
    - > 2000 molecules: refused with an error dict.

    Args:
        results: List of batch validation result dicts.
        threshold: Similarity threshold for sparse representation (default 0.7).

    Returns:
        Dict matching SimilarityMatrixResult schema.
    """
    t0 = time.perf_counter()
    indices, fps = _parse_and_fingerprint(results)
    n = len(indices)

    if n > 2000:
        return {
            "error": "Similarity matrix limited to 2000 molecules",
            "molecule_count": n,
        }

    if n == 0:
        return {"size": 0, "representation": "dense", "data": []}

    if n <= 500:
        # Dense upper-triangle: data[i] = list of similarity values for j > i
        data: list = []
        for i in range(n):
            row = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1 :])
            data.append([float(s) for s in row])

        elapsed = time.perf_counter() - t0
        logger.info(
            "compute_similarity_matrix: %d molecules (dense) in %.3fs", n, elapsed
        )
        return {"size": n, "representation": "dense", "data": data}

    else:
        # Sparse: only pairs above threshold
        data = []
        for i in range(n):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1 :])
            for j_offset, sim in enumerate(sims):
                if sim > threshold:
                    data.append({"i": i, "j": i + 1 + j_offset, "similarity": float(sim)})

        elapsed = time.perf_counter() - t0
        logger.info(
            "compute_similarity_matrix: %d molecules (sparse, threshold=%.2f) "
            "in %.3fs — %d pairs above threshold",
            n,
            threshold,
            elapsed,
            len(data),
        )
        return {"size": n, "representation": "sparse", "data": data}


def compute_chemical_space(results: list[dict], method: str = "pca", **kwargs) -> object:
    """
    Dispatcher for chemical space embedding — routes to PCA or t-SNE.

    Used by run_expensive_analytics Celery task.

    Args:
        results: List of batch validation result dicts.
        method: Either "pca" or "tsne" (default "pca").
        **kwargs: Additional keyword arguments forwarded to the chosen method.

    Returns:
        ChemSpaceCoordinates-compatible object or error dict.
    """
    from app.schemas.analytics import ChemSpaceCoordinates

    if method == "tsne":
        raw = compute_tsne(results, **kwargs)
    else:
        raw = compute_pca(results, **kwargs)

    if "error" in raw:
        # Return a minimal object that can be stored; callers should check for errors
        return raw

    return ChemSpaceCoordinates(**raw)
