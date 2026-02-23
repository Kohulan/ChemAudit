"""
Tests for chemical space analysis service.

Covers PCA, t-SNE (with size/perplexity guards), similarity search,
nearest neighbors, and similarity matrix with dense/sparse switching.
"""

import inspect

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_batch(smiles_list: list[str]) -> list[dict]:
    """Create minimal result dicts from a list of SMILES strings."""
    return [
        {
            "index": i,
            "status": "success",
            "smiles": smi,
        }
        for i, smi in enumerate(smiles_list)
    ]


# 10 diverse drug-like molecules
_DIVERSE_10 = [
    "c1ccccc1",            # benzene
    "c1ccc2ccccc2c1",      # naphthalene
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",  # testosterone
    "CN1CCC[C@H]1c2cccnc2",  # nicotine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
    "OC(=O)Cc1ccc(Oc2ccc(CC(O)=O)cc2)cc1",  # fenofibric acid analog
    "c1ccc(cc1)C2CC(=O)c3ccccc3O2",  # flavanone
    "CC(=O)Nc1ccc(O)cc1",  # paracetamol
    "c1cnc2ccccc2n1",      # quinoxaline
]

# 20 molecules (extend diverse 10 with variants)
_DIVERSE_20 = _DIVERSE_10 + [
    "CC1=CC(=O)c2ccccc2C1=O",  # methyl-naphthoquinone
    "O=C(O)c1ccccc1O",         # salicylic acid
    "Nc1ccc(cc1)S(N)(=O)=O",   # sulfanilamide
    "c1ccc2[nH]cccc2c1",       # indole
    "O=C1CCCN1",               # pyrrolidinone
    "CC(O)=O",                 # acetic acid
    "NCC(=O)O",                # glycine
    "c1ccncc1",                # pyridine
    "C1CCNCC1",                # piperidine
    "C1COCCN1",                # morpholine
]


def _make_large_batch(n: int) -> list[dict]:
    """Generate n molecules using varying alkyl/heteroatom SMILES."""
    smiles_list = []
    for i in range(n):
        # Cycle through simple patterns to produce a varied set
        base = i % 20
        if base < 10:
            # Alkyl alcohols: C*n + O
            length = (i // 20) + 1
            smi = "C" * (base + length) + "O"
        else:
            # Alkyl carboxylic acids
            length = (i // 20) + 1
            smi = "C" * (base - 10 + length) + "C(=O)O"
        smiles_list.append(smi)
    return _make_batch(smiles_list)


# ---------------------------------------------------------------------------
# PCA tests
# ---------------------------------------------------------------------------


def test_pca_returns_2d_coordinates():
    """PCA on 10 diverse molecules returns coordinates with correct shape."""
    from app.services.analytics.chemical_space import compute_pca

    result = compute_pca(_make_batch(_DIVERSE_10))

    assert "error" not in result, f"Unexpected error: {result.get('error')}"
    assert result["method"] == "pca"
    assert len(result["coordinates"]) == 10
    assert all(len(coord) == 2 for coord in result["coordinates"])
    assert result["variance_explained"] is not None
    assert len(result["variance_explained"]) == 2
    assert all(isinstance(v, float) for v in result["variance_explained"])
    assert sum(result["variance_explained"]) <= 1.0 + 1e-6


def test_pca_too_few_molecules():
    """PCA on 2 molecules returns an error dict."""
    from app.services.analytics.chemical_space import compute_pca

    result = compute_pca(_make_batch(["c1ccccc1", "CC(=O)O"]))
    assert "error" in result


def test_pca_reproducible():
    """PCA on the same batch twice produces identical coordinates (seeded RNG)."""
    from app.services.analytics.chemical_space import compute_pca

    batch = _make_batch(_DIVERSE_10)
    r1 = compute_pca(batch)
    r2 = compute_pca(batch)

    assert r1["coordinates"] == r2["coordinates"]


# ---------------------------------------------------------------------------
# t-SNE tests
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_tsne_returns_coordinates():
    """t-SNE on 20 diverse molecules returns 2-D coordinates."""
    from app.services.analytics.chemical_space import compute_tsne

    result = compute_tsne(_make_batch(_DIVERSE_20))

    assert "error" not in result, f"Unexpected error: {result.get('error')}"
    assert result["method"] == "tsne"
    assert len(result["coordinates"]) == 20
    assert all(len(coord) == 2 for coord in result["coordinates"])
    # variance_explained is None for t-SNE
    assert result["variance_explained"] is None


def test_tsne_refuses_large_batch():
    """t-SNE refuses batches larger than 2000 molecules."""
    from app.services.analytics.chemical_space import compute_tsne

    large_batch = _make_large_batch(2001)
    result = compute_tsne(large_batch)

    assert "error" in result
    assert "2000" in result["error"] or "limited" in result["error"].lower()


@pytest.mark.slow
def test_tsne_auto_adjusts_perplexity():
    """t-SNE auto-adjusts perplexity for small batches and completes without crashing."""
    from app.services.analytics.chemical_space import compute_tsne

    # 30 molecules with perplexity=30 would crash without auto-adjustment
    batch = _make_batch(_DIVERSE_20 + _DIVERSE_10)  # 30 molecules
    result = compute_tsne(batch, perplexity=30)

    # Should succeed (not crash, no error)
    assert "error" not in result, f"Unexpected error: {result.get('error')}"
    assert result["method"] == "tsne"
    assert len(result["coordinates"]) == 30


# ---------------------------------------------------------------------------
# Similarity search tests
# ---------------------------------------------------------------------------


def test_find_similar_by_smiles():
    """Query by external SMILES returns top_k neighbors sorted by descending similarity."""
    from app.services.analytics.chemical_space import find_similar_molecules

    batch = _make_batch(_DIVERSE_10)
    result = find_similar_molecules(batch, query_smiles="c1ccccc1", top_k=5)

    assert "error" not in result
    assert len(result["neighbors"]) <= 5
    sims = [n["similarity"] for n in result["neighbors"]]
    # Sorted descending
    assert sims == sorted(sims, reverse=True)
    # All similarities are valid floats in [0, 1]
    assert all(0.0 <= s <= 1.0 for s in sims)


def test_find_similar_by_index():
    """Query by batch index excludes the query molecule itself from results."""
    from app.services.analytics.chemical_space import find_similar_molecules

    batch = _make_batch(_DIVERSE_10)
    result = find_similar_molecules(batch, query_index=0, top_k=5)

    assert "error" not in result
    neighbor_indices = [n["index"] for n in result["neighbors"]]
    # Query molecule (index 0) should not appear in its own results
    assert 0 not in neighbor_indices


# ---------------------------------------------------------------------------
# Nearest neighbor tests
# ---------------------------------------------------------------------------


def test_nearest_neighbors():
    """Nearest neighbor analysis returns one result per molecule with valid isolation scores."""
    from app.services.analytics.chemical_space import compute_nearest_neighbors

    batch = _make_batch(_DIVERSE_10[:5])
    results = compute_nearest_neighbors(batch)

    assert len(results) == 5
    for r in results:
        assert "molecule_index" in r
        assert "nearest_index" in r
        assert "similarity" in r
        assert "isolation_score" in r
        assert 0.0 <= r["isolation_score"] <= 1.0
        # isolation_score must equal 1 - similarity
        assert abs(r["isolation_score"] - (1.0 - r["similarity"])) < 1e-6


# ---------------------------------------------------------------------------
# Similarity matrix tests
# ---------------------------------------------------------------------------


def test_similarity_matrix_dense():
    """Batch of 10 molecules produces a dense upper-triangle representation."""
    from app.services.analytics.chemical_space import compute_similarity_matrix

    batch = _make_batch(_DIVERSE_10)
    result = compute_similarity_matrix(batch)

    assert "error" not in result
    assert result["representation"] == "dense"
    assert result["size"] == 10
    # Upper triangle: row i has (n - 1 - i) entries; row 0 has 9, row 9 has 0
    data = result["data"]
    assert len(data) == 10
    for i, row in enumerate(data):
        assert len(row) == 10 - 1 - i
        assert all(isinstance(v, float) for v in row)
        assert all(0.0 <= v <= 1.0 for v in row)


def test_similarity_matrix_sparse():
    """Batch of 600+ molecules produces a sparse representation above threshold."""
    from app.services.analytics.chemical_space import compute_similarity_matrix

    batch = _make_large_batch(600)
    result = compute_similarity_matrix(batch, threshold=0.7)

    assert "error" not in result
    assert result["representation"] == "sparse"
    assert result["size"] == 600
    # Sparse data should only contain pairs above threshold
    for entry in result["data"]:
        assert "i" in entry
        assert "j" in entry
        assert "similarity" in entry
        assert entry["similarity"] > 0.7
        assert entry["i"] < entry["j"]  # upper-triangle only


def test_similarity_matrix_refuses_large():
    """Similarity matrix refuses batches larger than 2000 molecules."""
    from app.services.analytics.chemical_space import compute_similarity_matrix

    large_batch = _make_large_batch(2001)
    result = compute_similarity_matrix(large_batch)

    assert "error" in result
    assert "2000" in result["error"] or "limited" in result["error"].lower()


# ---------------------------------------------------------------------------
# API / implementation correctness tests
# ---------------------------------------------------------------------------


def test_fingerprint_uses_morgan_generator():
    """Verify that rdFingerprintGenerator.GetMorganGenerator is used (not deprecated API)."""
    import app.services.analytics.chemical_space as cspace

    # Check that the module imports rdFingerprintGenerator (non-deprecated)
    source = inspect.getsource(cspace)
    assert "rdFingerprintGenerator" in source
    assert "GetMorganGenerator" in source
    # Ensure the old deprecated call is NOT used
    assert "GetMorganFingerprintAsBitVect" not in source
