/**
 * Batch Analytics Types
 *
 * TypeScript interfaces mirroring backend analytics.py schemas.
 * Used by useBatchAnalytics hook and BatchAnalyticsPanel components.
 */

/**
 * Status of a single analytics computation.
 */
export interface AnalysisStatus {
  status: 'pending' | 'computing' | 'complete' | 'failed' | 'skipped';
  computed_at: number | null;
  error: string | null;
}

// ---------------------------------------------------------------------------
// Deduplication
// ---------------------------------------------------------------------------

export interface DeduplicationGroup {
  level: string;
  representative_index: number;
  duplicate_indices: number[];
  group_key: string;
  count: number;
}

export interface DeduplicationResult {
  exact: DeduplicationGroup[];
  tautomeric: DeduplicationGroup[];
  stereo_insensitive: DeduplicationGroup[];
  salt_form: DeduplicationGroup[];
  total_unique: Record<string, number>;
}

// ---------------------------------------------------------------------------
// Scaffold analysis
// ---------------------------------------------------------------------------

export interface ScaffoldGroup {
  scaffold_smiles: string;
  generic_scaffold_smiles: string;
  molecule_indices: number[];
  count: number;
}

export interface ScaffoldResult {
  scaffolds: ScaffoldGroup[];
  unique_scaffold_count: number;
  shannon_entropy: number;
  frequency_distribution: Record<string, number>;
}

// ---------------------------------------------------------------------------
// Chemical space
// ---------------------------------------------------------------------------

export interface ChemSpaceCoordinates {
  method: 'pca' | 'tsne';
  coordinates: number[][];
  molecule_indices: number[];
  variance_explained?: number[];
}

// ---------------------------------------------------------------------------
// Similarity
// ---------------------------------------------------------------------------

export interface SimilarityMatrixResult {
  size: number;
  representation: 'dense' | 'sparse';
  data: unknown[];
}

// ---------------------------------------------------------------------------
// MMP
// ---------------------------------------------------------------------------

export interface MMPPair {
  mol_a_index: number;
  mol_b_index: number;
  core_smiles: string;
  rgroup_a: string;
  rgroup_b: string;
  tanimoto: number;
}

export interface ActivityCliff {
  mol_a_index: number;
  mol_b_index: number;
  sali: number;
  tanimoto: number;
  activity_diff: number;
}

export interface MMPResult {
  pairs: MMPPair[];
  activity_cliffs?: ActivityCliff[];
  lle_values?: Record<string, unknown>[];
}

// ---------------------------------------------------------------------------
// Property statistics
// ---------------------------------------------------------------------------

export interface PropertyStats {
  property_name: string;
  mean: number;
  median: number;
  std: number;
  q1: number;
  q3: number;
  iqr: number;
  min: number;
  max: number;
  count: number;
}

export interface OutlierInfo {
  molecule_index: number;
  property_name: string;
  value: number;
  lower_fence: number;
  upper_fence: number;
}

export interface PropertyCorrelation {
  property_a: string;
  property_b: string;
  pearson_r: number;
}

export interface QualityScore {
  score: number;
  validity_pct: number;
  diversity_pct: number;
  druglikeness_pct: number;
}

export interface StatisticsResult {
  property_stats: PropertyStats[];
  correlations: PropertyCorrelation[];
  outliers: OutlierInfo[];
  quality_score: QualityScore;
}

// ---------------------------------------------------------------------------
// Clustering (Butina)
// ---------------------------------------------------------------------------

export interface ClusterInfo {
  cluster_id: number;
  size: number;
  representative_index: number;
  member_indices: number[];
}

export interface ClusteringResult {
  clusters: ClusterInfo[];
  cluster_count: number;
  singleton_count: number;
  largest_cluster_size: number;
  distance_cutoff: number;
  /** Maps batch result index (as string key) to SMILES for molecule rendering */
  smiles_map?: Record<string, string>;
}

// ---------------------------------------------------------------------------
// Chemical Taxonomy
// ---------------------------------------------------------------------------

export interface TaxonomyCategory {
  name: string;
  category: string;
  description: string;
}

export interface TaxonomyMoleculeResult {
  index: number;
  smiles: string;
  categories: TaxonomyCategory[];
}

export interface TaxonomyResult {
  per_molecule: TaxonomyMoleculeResult[];
  category_counts: Record<string, number>;
  total_molecules: number;
  classified_molecules: number;
  unclassified_molecules: number;
}

// ---------------------------------------------------------------------------
// Registration Hash
// ---------------------------------------------------------------------------

export interface RegistrationHashMolecule {
  index: number;
  smiles: string;
  hash: string;
  canonical_smiles: string;
  formula: string;
}

export interface RegistrationHashCollisionGroup {
  hash: string;
  molecule_indices: number[];
  count: number;
}

export interface RegistrationHashResult {
  per_molecule: RegistrationHashMolecule[];
  unique_count: number;
  total_count: number;
  collision_groups: RegistrationHashCollisionGroup[];
  rdkit_version: string;
  tautomer_hash_v2: boolean;
}

// ---------------------------------------------------------------------------
// MCS Comparison
// ---------------------------------------------------------------------------

export interface MCSPropertyDelta {
  property: string;
  mol_a: number;
  mol_b: number;
  delta: number;
}

export interface MCSComparisonResult {
  mcs_smarts: string;
  num_atoms: number;
  num_bonds: number;
  timed_out: boolean;
  tanimoto: number;
  property_deltas: MCSPropertyDelta[];
  smiles_a: string;
  smiles_b: string;
}

// ---------------------------------------------------------------------------
// Top-level response schemas
// ---------------------------------------------------------------------------

export interface BatchAnalyticsResponse {
  job_id: string;
  status: Record<string, AnalysisStatus>;
  deduplication?: DeduplicationResult;
  scaffold?: ScaffoldResult;
  chemical_space?: ChemSpaceCoordinates;
  similarity_matrix?: SimilarityMatrixResult;
  mmp?: MMPResult;
  statistics?: StatisticsResult;
  clustering?: ClusteringResult;
  taxonomy?: TaxonomyResult;
  registration?: RegistrationHashResult;
}

export interface AnalyticsTriggerResponse {
  job_id: string;
  analysis_type: string;
  status: string;
}
