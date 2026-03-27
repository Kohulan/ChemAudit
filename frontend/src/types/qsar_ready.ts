/**
 * TypeScript interfaces for the QSAR-Ready Pipeline API.
 * Covers all endpoints under /api/v1/qsar-ready/.
 *
 * All shapes match backend Pydantic schemas in backend/app/schemas/qsar_ready.py.
 */

// ====================== Config ======================

/**
 * Configuration for the QSAR-Ready curation pipeline.
 * Controls which steps are enabled and molecular property filter thresholds.
 */
export interface QSARReadyConfig {
  enable_metals: boolean;
  enable_desalt: boolean;
  enable_normalize: boolean;
  enable_neutralize: boolean;
  enable_tautomer: boolean;
  enable_stereo_strip: boolean;
  enable_isotope_strip: boolean;
  min_heavy_atoms: number;
  max_heavy_atoms: number;
  max_mw: number;
  remove_inorganics: boolean;
}

// ====================== Single Molecule Result ======================

/**
 * Result for a single pipeline step.
 * Tracks SMILES before/after the step for provenance.
 */
export interface QSARStepResult {
  step_name: string;
  step_index: number;
  enabled: boolean;
  status: 'applied' | 'no_change' | 'skipped' | 'error';
  before_smiles: string | null;
  after_smiles: string | null;
  detail: string | null;
}

/**
 * Full curation result for a single molecule.
 * Carries both original and standardized InChIKeys to track identity changes.
 */
export interface QSARReadyResult {
  original_smiles: string;
  original_inchikey: string | null;
  curated_smiles: string | null;
  standardized_inchikey: string | null;
  inchikey_changed: boolean;
  status: 'ok' | 'rejected' | 'duplicate' | 'error';
  rejection_reason: string | null;
  steps: QSARStepResult[];
}

// ====================== Batch ======================

/**
 * Summary statistics for a completed batch QSAR-ready job.
 */
export interface QSARBatchSummary {
  total: number;
  ok: number;
  rejected: number;
  duplicate: number;
  error: number;
  steps_applied_counts: Record<string, number>;
}

/**
 * Response from the batch upload endpoint.
 * Returns job_id for polling status and fetching results.
 */
export interface QSARBatchUploadResponse {
  job_id: string;
  total_molecules: number;
  status: string;
  message: string;
}

/**
 * Real-time status response for a batch job (WebSocket + polling).
 */
export interface QSARBatchStatusResponse {
  job_id: string;
  status: 'pending' | 'processing' | 'complete' | 'failed';
  progress: number;
  processed: number;
  total: number;
}

/**
 * Paginated results response for a completed batch job.
 */
export interface QSARBatchResultsResponse {
  job_id: string;
  status: string;
  config: QSARReadyConfig;
  summary: QSARBatchSummary;
  results: QSARReadyResult[];
  page: number;
  per_page: number;
  total_pages: number;
}

// ====================== Pipeline Step Metadata ======================

/**
 * Metadata for a single configurable pipeline step.
 * Used to render the step toggle grid in the UI.
 */
export interface QSARPipelineStep {
  key: keyof Omit<QSARReadyConfig, 'min_heavy_atoms' | 'max_heavy_atoms' | 'max_mw' | 'remove_inorganics'>;
  index: number;
  name: string;
  description: string;
}

/**
 * Saved pipeline config profile for localStorage persistence.
 */
export interface QSARPipelineProfile {
  id: string;
  name: string;
  config: QSARReadyConfig;
  createdAt: number;
}

// ====================== Constants ======================

/** Canonical step order for UI rendering. Mirrors backend 10-step pipeline. */
export const PIPELINE_STEPS: QSARPipelineStep[] = [
  { key: 'enable_metals', index: 2, name: 'Metal Disconnect', description: 'Remove metal bonds and ions' },
  { key: 'enable_desalt', index: 3, name: 'Desalt', description: 'Keep largest organic fragment' },
  { key: 'enable_normalize', index: 4, name: 'Normalize', description: 'Fix functional group representations' },
  { key: 'enable_neutralize', index: 5, name: 'Neutralize', description: 'Remove charges where possible' },
  { key: 'enable_tautomer', index: 6, name: 'Tautomer', description: 'Canonicalize tautomeric form' },
  { key: 'enable_stereo_strip', index: 7, name: 'Strip Stereo', description: 'Remove stereochemistry information' },
  { key: 'enable_isotope_strip', index: 8, name: 'Strip Isotopes', description: 'Remove isotope labels' },
];

/** Default configs for the three built-in presets. */
export const PRESET_CONFIGS: Record<string, QSARReadyConfig> = {
  'qsar-2d': {
    enable_metals: true,
    enable_desalt: true,
    enable_normalize: true,
    enable_neutralize: true,
    enable_tautomer: true,
    enable_stereo_strip: true,
    enable_isotope_strip: true,
    min_heavy_atoms: 3,
    max_heavy_atoms: 100,
    max_mw: 1500.0,
    remove_inorganics: true,
  },
  'qsar-3d': {
    enable_metals: true,
    enable_desalt: true,
    enable_normalize: true,
    enable_neutralize: true,
    enable_tautomer: true,
    enable_stereo_strip: false,
    enable_isotope_strip: true,
    min_heavy_atoms: 3,
    max_heavy_atoms: 100,
    max_mw: 1500.0,
    remove_inorganics: true,
  },
  'minimal': {
    enable_metals: false,
    enable_desalt: true,
    enable_normalize: true,
    enable_neutralize: false,
    enable_tautomer: false,
    enable_stereo_strip: false,
    enable_isotope_strip: false,
    min_heavy_atoms: 3,
    max_heavy_atoms: 100,
    max_mw: 1500.0,
    remove_inorganics: false,
  },
};
