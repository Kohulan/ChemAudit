"""
Analytics Pydantic Schemas

Response models for all batch analytics endpoints.
Covers deduplication, scaffold analysis, chemical space, similarity, MMP, and statistics.
"""

from typing import Literal, Optional

from pydantic import BaseModel

# ---------------------------------------------------------------------------
# Shared / status
# ---------------------------------------------------------------------------


class AnalysisStatus(BaseModel):
    """Status of a single analytics computation."""

    status: Literal["pending", "computing", "complete", "failed", "skipped"]
    computed_at: Optional[float] = None
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------


class DeduplicationGroup(BaseModel):
    """A group of structurally equivalent molecules at a given comparison level."""

    level: str
    representative_index: int
    duplicate_indices: list[int]
    group_key: str
    count: int


class DeduplicationResult(BaseModel):
    """Full deduplication result across all comparison levels."""

    exact: list[DeduplicationGroup]
    tautomeric: list[DeduplicationGroup]
    stereo_insensitive: list[DeduplicationGroup]
    salt_form: list[DeduplicationGroup]
    total_unique: dict[str, int]  # counts per level


# ---------------------------------------------------------------------------
# Scaffold analysis
# ---------------------------------------------------------------------------


class ScaffoldGroup(BaseModel):
    """Group of molecules sharing the same Murcko scaffold."""

    scaffold_smiles: str
    generic_scaffold_smiles: str
    molecule_indices: list[int]
    count: int


class ScaffoldResult(BaseModel):
    """Scaffold diversity analysis results."""

    scaffolds: list[ScaffoldGroup]
    unique_scaffold_count: int
    shannon_entropy: float
    frequency_distribution: dict[str, int]


# ---------------------------------------------------------------------------
# R-group decomposition
# ---------------------------------------------------------------------------


class RGroupResult(BaseModel):
    """R-group decomposition results around a common core."""

    core_smarts: str
    decomposition: list[dict]
    unmatched_count: int


# ---------------------------------------------------------------------------
# Chemical space
# ---------------------------------------------------------------------------


class ChemSpaceCoordinates(BaseModel):
    """2-D chemical space embedding coordinates."""

    method: Literal["pca", "tsne"]
    coordinates: list[list[float]]
    molecule_indices: list[int]
    variance_explained: Optional[list[float]] = None  # PCA only


# ---------------------------------------------------------------------------
# Similarity
# ---------------------------------------------------------------------------


class SimilarityHit(BaseModel):
    """Single neighbor from a similarity search."""

    index: int
    similarity: float
    smiles: str


class SimilarityResult(BaseModel):
    """Similarity search results for a query molecule."""

    query_index: Optional[int]
    query_smiles: Optional[str]
    neighbors: list[SimilarityHit]


class NearestNeighborResult(BaseModel):
    """Nearest neighbor and isolation score for a molecule."""

    molecule_index: int
    nearest_index: int
    similarity: float
    isolation_score: float


class SimilarityMatrixResult(BaseModel):
    """Pairwise similarity matrix (dense or sparse)."""

    size: int
    representation: Literal["dense", "sparse"]
    # dense: full upper-triangle list; sparse: list of {i, j, similarity} dicts
    data: list


# ---------------------------------------------------------------------------
# MMP / Activity cliffs
# ---------------------------------------------------------------------------


class MMPPair(BaseModel):
    """Matched molecular pair."""

    mol_a_index: int
    mol_b_index: int
    core_smiles: str
    rgroup_a: str
    rgroup_b: str
    tanimoto: float


class ActivityCliff(BaseModel):
    """Activity cliff between two molecules."""

    mol_a_index: int
    mol_b_index: int
    sali: float
    tanimoto: float
    activity_diff: float


class MMPResult(BaseModel):
    """Matched molecular pair analysis results."""

    pairs: list[MMPPair]
    activity_cliffs: Optional[list[ActivityCliff]] = None
    lle_values: Optional[list[dict]] = None


# ---------------------------------------------------------------------------
# Property statistics
# ---------------------------------------------------------------------------


class PropertyStats(BaseModel):
    """Descriptive statistics for a single molecular property."""

    property_name: str
    mean: float
    median: float
    std: float
    q1: float
    q3: float
    iqr: float
    min: float
    max: float
    count: int


class OutlierInfo(BaseModel):
    """An outlier molecule for a specific property."""

    molecule_index: int
    property_name: str
    value: float
    lower_fence: float
    upper_fence: float


class PropertyCorrelation(BaseModel):
    """Pearson correlation between two molecular properties."""

    property_a: str
    property_b: str
    pearson_r: float


class QualityScore(BaseModel):
    """Composite dataset quality score."""

    score: float
    validity_pct: float
    diversity_pct: float
    druglikeness_pct: float


class StatisticsResult(BaseModel):
    """Comprehensive property statistics for the batch dataset."""

    property_stats: list[PropertyStats]
    correlations: list[PropertyCorrelation]
    outliers: list[OutlierInfo]
    quality_score: QualityScore


# ---------------------------------------------------------------------------
# Top-level response schemas
# ---------------------------------------------------------------------------


class BatchAnalyticsResponse(BaseModel):
    """Response for GET /batch/{job_id}/analytics."""

    job_id: str
    status: dict[str, AnalysisStatus]
    deduplication: Optional[DeduplicationResult] = None
    scaffold: Optional[ScaffoldResult] = None
    chemical_space: Optional[ChemSpaceCoordinates] = None
    similarity_matrix: Optional[SimilarityMatrixResult] = None
    mmp: Optional[MMPResult] = None
    statistics: Optional[StatisticsResult] = None


class AnalyticsTriggerResponse(BaseModel):
    """Response for POST /batch/{job_id}/analytics/{analysis_type}."""

    job_id: str
    analysis_type: str
    status: str
