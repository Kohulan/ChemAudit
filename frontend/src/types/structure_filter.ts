/**
 * TypeScript interfaces for the GenChem Filter API (Phase 11).
 * Covers all endpoints under /api/v1/genchem/.
 *
 * All shapes match backend Pydantic schemas in backend/app/schemas/genchem.py.
 */

// ====================== Config ======================

/**
 * Configuration for the generative chemistry filter funnel.
 * Controls property thresholds, alert catalogs, and scoring weights.
 */
export interface FilterConfig {
  min_mw: number;
  max_mw: number;
  min_logp: number;
  max_logp: number;
  max_tpsa: number;
  max_rot_bonds: number;
  max_rings: number | null; // null = no ring limit
  max_sa_score: number;
  use_pains: boolean;
  use_brenk: boolean;
  use_kazius: boolean;
  use_nibr: boolean;
  enable_novelty: boolean;
  novelty_threshold: number;
  weight_validity: number;
  weight_qed: number;
  weight_alert_free: number;
  weight_sa: number;
}

// ====================== Filter Results ======================

/**
 * Per-stage funnel statistics.
 * Tracks how many molecules enter and pass each stage.
 */
export interface StageResult {
  stage_name: string;
  stage_index: number;
  input_count: number;
  passed_count: number;
  rejected_count: number;
  enabled: boolean;
}

/**
 * Single molecule outcome from the filter funnel.
 */
export interface MoleculeResult {
  smiles: string;
  status: 'passed' | 'rejected' | 'duplicate' | 'error';
  failed_at: string | null;
  rejection_reason: string | null;
}

/**
 * Full funnel result: stage-level statistics + per-molecule outcomes.
 */
export interface FilterResult {
  input_count: number;
  output_count: number;
  stages: StageResult[];
  molecules: MoleculeResult[];
}

// ====================== Scoring ======================

/**
 * Response from the /genchem/score endpoint.
 * Returns a score in [0, 1] per molecule (null on error).
 */
export interface ScoreResponse {
  scores: (number | null)[];
}

// ====================== REINVENT Scoring ======================

/**
 * Single SMILES input item for the REINVENT-compatible scoring endpoint.
 */
export interface REINVENTInput {
  input_string: string;
  query_id: string;
}

/**
 * Response from the /genchem/reinvent-score endpoint.
 * Follows the REINVENT Component API contract (output.successes_list).
 */
export interface REINVENTResponse {
  output: {
    successes_list: Array<{
      query_id: string;
      output_value: number;
    }>;
  };
}

// ====================== Batch ======================

/**
 * Response from the batch upload endpoint.
 * Returns job_id for tracking progress via WebSocket or polling.
 */
export interface GenChemBatchUploadResponse {
  job_id: string;
  total_molecules: number;
  status: string;
}

/**
 * Real-time status response for a batch filter job.
 */
export interface GenChemBatchStatusResponse {
  job_id: string;
  status: string;
  progress: number | null;
  current_stage: string | null;
}

/**
 * Results response for a completed batch filter job.
 */
export interface GenChemBatchResultsResponse {
  job_id: string;
  status: string;
  result: FilterResult | null;
}

// ====================== Profile ======================

/**
 * Saved filter config profile for localStorage persistence.
 */
export interface GenChemProfile {
  id: string;
  name: string;
  config: FilterConfig;
}

// ====================== Constants ======================

/**
 * Stage color map per D-06 locked decision.
 * Maps stage name to hex color for funnel chart rendering.
 */
export const STAGE_COLORS: Record<string, string> = {
  parse: '#22c55e',
  valence: '#4ade80',
  alerts: '#f59e0b',
  property: '#f97316',
  sa: '#ef4444',
  sa_score: '#ef4444',
  dedup: '#dc2626',
  novelty: '#7c3aed',
};

/**
 * Built-in preset configurations (per D-12 locked, D-15 differentiated weight vectors).
 *
 * Weight vectors are intentionally differentiated per D-15:
 * - drug_like: balanced (0.3/0.3/0.2/0.2)
 * - lead_like: emphasize QED for lead optimization (0.2/0.4/0.2/0.2)
 * - fragment_like: emphasize alert-free and SA for clean building blocks (0.2/0.2/0.3/0.3)
 * - permissive: emphasize validity, lower alert_free since alerts off (0.4/0.3/0.1/0.2)
 */
export const PRESET_CONFIGS: Record<string, FilterConfig> = {
  drug_like: {
    min_mw: 200,
    max_mw: 500,
    min_logp: -1,
    max_logp: 5,
    max_tpsa: 140,
    max_rot_bonds: 10,
    max_rings: null,
    max_sa_score: 5.0,
    use_pains: true,
    use_brenk: true,
    use_kazius: true,
    use_nibr: false,
    enable_novelty: false,
    novelty_threshold: 0.85,
    // D-15: balanced default
    weight_validity: 0.3,
    weight_qed: 0.3,
    weight_alert_free: 0.2,
    weight_sa: 0.2,
  },
  lead_like: {
    min_mw: 200,
    max_mw: 350,
    min_logp: -1,
    max_logp: 3.5,
    max_tpsa: 140,
    max_rot_bonds: 7,
    max_rings: null,
    max_sa_score: 4.0,
    use_pains: true,
    use_brenk: true,
    use_kazius: true,
    use_nibr: true,
    enable_novelty: false,
    novelty_threshold: 0.85,
    // D-15: emphasize QED for lead optimization
    weight_validity: 0.2,
    weight_qed: 0.4,
    weight_alert_free: 0.2,
    weight_sa: 0.2,
  },
  fragment_like: {
    min_mw: 100,
    max_mw: 300,
    min_logp: -1,
    max_logp: 3,
    max_tpsa: 100,
    max_rot_bonds: 3,
    max_rings: 3,
    max_sa_score: 3.0,
    use_pains: true,
    use_brenk: false,
    use_kazius: false,
    use_nibr: false,
    enable_novelty: false,
    novelty_threshold: 0.85,
    // D-15: emphasize alert-free and SA for clean building blocks
    weight_validity: 0.2,
    weight_qed: 0.2,
    weight_alert_free: 0.3,
    weight_sa: 0.3,
  },
  permissive: {
    min_mw: 100,
    max_mw: 800,
    min_logp: -5,
    max_logp: 8,
    max_tpsa: 200,
    max_rot_bonds: 15,
    max_rings: null,
    max_sa_score: 7.0,
    use_pains: false,
    use_brenk: false,
    use_kazius: false,
    use_nibr: false,
    enable_novelty: false,
    novelty_threshold: 0.85,
    // D-15: emphasize validity (main filter), lower alert_free (alerts off)
    weight_validity: 0.4,
    weight_qed: 0.3,
    weight_alert_free: 0.1,
    weight_sa: 0.2,
  },
};

/** Set of built-in preset IDs (cannot be deleted). */
export const PRESET_IDS = new Set(['drug_like', 'lead_like', 'fragment_like', 'permissive']);
