/**
 * TypeScript interfaces for the Dataset Intelligence API (Phase 12).
 * Covers all endpoints under /api/v1/dataset/.
 *
 * All shapes match backend Pydantic schemas in backend/app/schemas/dataset_intelligence.py.
 */

// ====================== Health Audit ======================

/** Sub-score detail from health audit. */
export interface SubScoreDetail {
  name: string;
  score: number; // 0-1
  weight: number;
  count: number; // affected count
  total: number;
}

/** Single issue entry from health audit. */
export interface IssueEntry {
  row_index: number;
  smiles: string;
  issue_type: string;
  severity: string;
  description: string;
}

/** Histogram data for property distributions. */
export interface HistogramData {
  bins: number[];
  counts: number[];
}

/** Property distribution histograms (MW, LogP, TPSA). */
export interface PropertyDistributions {
  mw: HistogramData;
  logp: HistogramData;
  tpsa: HistogramData;
}

/** Duplicate group identified by InChIKey. */
export interface DedupGroup {
  inchikey: string;
  rows: number[];
}

/** Full health audit result. */
export interface HealthAuditResult {
  overall_score: number; // 0-100
  sub_scores: SubScoreDetail[];
  weights: Record<string, number>;
  molecule_count: number;
  issues: IssueEntry[];
  property_distributions: PropertyDistributions;
  std_pipeline_comparison: Record<string, unknown>;
  std_sample_size: number;
  dedup_groups: DedupGroup[];
}

// ====================== Contradictory Labels ======================

/** Single contradictory label entry. */
export interface ContradictionEntry {
  row_index: number;
  smiles: string;
  activity: number;
}

/** Contradictory label result for a single InChIKey. */
export interface ContradictoryLabelResult {
  inchikey: string;
  entries: ContradictionEntry[];
  fold_difference: number;
  entry_count: number;
  smiles: string;
}

// ====================== Numeric Columns ======================

/** Numeric column detected for activity selection. */
export interface NumericColumnInfo {
  name: string;
  priority: number;
}

// ====================== Full Audit Results ======================

/** Full dataset audit results (shared across all tabs). */
export interface DatasetAuditResults {
  job_id: string;
  status: string;
  health_audit: HealthAuditResult | null;
  contradictions: ContradictoryLabelResult[];
  numeric_columns: NumericColumnInfo[];
  curation_report: Record<string, unknown> | null;
  curated_csv_available: boolean;
}

// ====================== Dataset Diff ======================

/** Property change between two dataset versions. */
export interface PropertyChange {
  column: string;
  old_value: unknown;
  new_value: unknown;
}

/** Single molecule in a dataset diff result. */
export interface DiffMolecule {
  inchikey: string;
  smiles: string;
  row_index: number;
  properties: Record<string, unknown>;
  changes: PropertyChange[];
}

/** Full dataset diff results. */
export interface DatasetDiffResults {
  added: DiffMolecule[];
  removed: DiffMolecule[];
  modified: DiffMolecule[];
  added_count: number;
  removed_count: number;
  modified_count: number;
  unchanged_count: number;
  unique_columns_primary: number;
  unique_columns_comparison: number;
}

// ====================== API Response Types ======================

/** Response from POST /api/v1/dataset/upload. */
export interface DatasetUploadResponse {
  job_id: string;
  filename: string;
  file_type: string;
  status: string;
  message: string;
}

/** Response from GET /api/v1/dataset/{job_id}/status. */
export interface DatasetAuditStatusResponse {
  job_id: string;
  status: string;
  progress: number;
  current_stage: string | null;
  eta_seconds: number | null;
}

// ====================== Weight Profiles ======================

/** Saved weight profile for localStorage persistence. */
export interface DatasetWeightProfile {
  id: string;
  name: string;
  weights: Record<string, number>;
}

// ====================== Page State ======================

/** Tab identifiers for the Dataset Audit page. */
export type DatasetAuditTab = 'health' | 'contradictions' | 'diff' | 'report';

/** Status state machine for dataset audit lifecycle. */
export type DatasetAuditStatus = 'idle' | 'uploading' | 'processing' | 'complete' | 'error';

// ====================== Constants ======================

/**
 * Literature-backed default weights (Fourches 2010).
 * parsability=0.25, stereo=0.15, uniqueness=0.20, alerts=0.20, std_consistency=0.20
 */
export const DEFAULT_WEIGHTS: Record<string, number> = {
  parsability: 0.25,
  stereo: 0.15,
  uniqueness: 0.20,
  alerts: 0.20,
  std_consistency: 0.20,
};

/** Fold-difference threshold options for contradictory label filtering. */
export const FOLD_THRESHOLD_OPTIONS = [3, 5, 10, 50, 100] as const;

/** Default fold-difference threshold. */
export const DEFAULT_FOLD_THRESHOLD = 10;

/** Set of built-in preset weight profile IDs (cannot be deleted). */
export const PRESET_WEIGHT_IDS = new Set(['default_literature']);
