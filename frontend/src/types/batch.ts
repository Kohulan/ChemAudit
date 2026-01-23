/**
 * Batch Processing Types
 *
 * TypeScript interfaces for batch upload, progress tracking, and results display.
 */

/**
 * Progress update from WebSocket
 */
export interface BatchProgress {
  job_id: string;
  status: 'pending' | 'processing' | 'complete' | 'failed' | 'cancelled';
  progress: number; // 0-100
  processed: number;
  total: number;
  eta_seconds: number | null;
  error_message?: string | null;
}

/**
 * Single molecule result in batch
 */
export interface BatchResult {
  smiles: string;
  name: string | null;
  index: number;
  status: 'success' | 'error';
  error: string | null;
  validation: {
    overall_score: number;
    issues?: Array<{
      check_name: string;
      passed: boolean;
      severity: string;
      message: string;
    }>;
    error?: string;
  } | null;
  alerts: {
    has_alerts: boolean;
    alert_count: number;
    alerts: Array<{
      catalog: string;
      rule_name: string;
      severity: string;
    }>;
    error?: string;
  } | null;
  scoring: {
    ml_readiness: {
      score: number;
      interpretation: string;
    };
    error?: string;
  } | null;
}

/**
 * Aggregate statistics for batch job
 */
export interface BatchStatistics {
  total: number;
  successful: number;
  errors: number;
  avg_validation_score: number | null;
  avg_ml_readiness_score: number | null;
  score_distribution: {
    excellent: number;
    good: number;
    moderate: number;
    poor: number;
  };
  alert_summary: Record<string, number>;
  processing_time_seconds: number | null;
}

/**
 * Response from batch upload endpoint
 */
export interface BatchUploadResponse {
  job_id: string;
  status: string;
  total_molecules: number;
  message: string;
}

/**
 * Response from batch results endpoint
 */
export interface BatchResultsResponse {
  job_id: string;
  status: string;
  statistics: BatchStatistics | null;
  results: BatchResult[];
  page: number;
  page_size: number;
  total_results: number;
  total_pages: number;
}

/**
 * Response from column detection endpoint
 */
export interface CSVColumnsResponse {
  columns: string[];
  suggested_smiles: string | null;
  row_count_estimate: number;
}

/**
 * Filters for results query
 */
export interface BatchResultsFilters {
  status_filter?: 'success' | 'error';
  min_score?: number;
  max_score?: number;
}

/**
 * Batch validation page state
 */
export type BatchPageState = 'upload' | 'processing' | 'results';

/**
 * Export format options (matches backend ExportFormat enum)
 */
export type ExportFormat = 'csv' | 'excel' | 'sdf' | 'json' | 'pdf';
