import axios, { AxiosError } from 'axios';

/**
 * Parse a single CSV line, handling quoted fields correctly.
 */
function parseCSVLine(line: string): string[] {
  const result: string[] = [];
  let current = '';
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const char = line[i];
    const nextChar = line[i + 1];

    if (inQuotes) {
      if (char === '"' && nextChar === '"') {
        // Escaped quote
        current += '"';
        i++;
      } else if (char === '"') {
        inQuotes = false;
      } else {
        current += char;
      }
    } else {
      if (char === '"') {
        inQuotes = true;
      } else if (char === ',') {
        result.push(current.trim());
        current = '';
      } else {
        current += char;
      }
    }
  }

  result.push(current.trim());
  return result;
}

/**
 * Escape a value for CSV output.
 */
function escapeCSVValue(value: string): string {
  if (value.includes(',') || value.includes('"') || value.includes('\n')) {
    return `"${value.replace(/"/g, '""')}"`;
  }
  return value;
}
import type {
  ValidationRequest,
  ValidationResponse,
  ValidationError,
  ChecksResponse
} from '../types/validation';
import type {
  AlertScreenRequest,
  AlertScreenResponse,
  AlertError,
  CatalogListResponse
} from '../types/alerts';
import type {
  ScoringRequest,
  ScoringResponse,
  ScoringError
} from '../types/scoring';
import type {
  StandardizeRequest,
  StandardizeResponse,
  StandardizeError,
  StandardizeOptionsResponse
} from '../types/standardization';
import type {
  BatchUploadResponse,
  BatchResultsResponse,
  BatchStatistics,
  BatchResultsFilters,
  CSVColumnsResponse
} from '../types/batch';
import type {
  IntegrationRequest,
  PubChemResult,
  ChEMBLResult,
  COCONUTResult,
  IntegrationError
} from '../types/integrations';

// API Configuration
const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8001/api/v1';
export const API_DOCS_URL = import.meta.env.VITE_API_DOCS_URL || 'http://localhost:8001/docs';
export const DEBUG_MODE = import.meta.env.VITE_DEBUG === 'true';

// Debug logging in development
if (DEBUG_MODE) {
  console.log('[ChemVault API] Configuration:', {
    apiUrl: API_BASE_URL,
    docsUrl: API_DOCS_URL,
    environment: import.meta.env.MODE,
  });
}

export const api = axios.create({
  baseURL: API_BASE_URL,
  timeout: 30000,
  headers: {
    'Content-Type': 'application/json',
  },
});

export const validationApi = {
  validate: async (request: ValidationRequest): Promise<ValidationResponse> => {
    try {
      const response = await api.post<ValidationResponse>('/validate', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<ValidationError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  getChecks: async (): Promise<ChecksResponse> => {
    const response = await api.get<ChecksResponse>('/checks');
    return response.data;
  },

  healthCheck: async (): Promise<{ status: string; rdkit_version: string }> => {
    const response = await api.get('/health');
    return response.data;
  }
};

export const alertsApi = {
  /**
   * Screen a molecule for structural alerts (PAINS, BRENK, etc.)
   */
  screenAlerts: async (request: AlertScreenRequest): Promise<AlertScreenResponse> => {
    try {
      const response = await api.post<AlertScreenResponse>('/alerts', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<AlertError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  /**
   * Quick check if molecule has any alerts (faster, no details)
   */
  quickCheck: async (request: AlertScreenRequest): Promise<{ has_alerts: boolean; checked_catalogs: string[] }> => {
    try {
      const response = await api.post('/alerts/quick-check', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<AlertError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  /**
   * List available alert catalogs
   */
  getCatalogs: async (): Promise<CatalogListResponse> => {
    const response = await api.get<CatalogListResponse>('/alerts/catalogs');
    return response.data;
  }
};

export const scoringApi = {
  /**
   * Calculate comprehensive scores for a molecule including ML-readiness,
   * NP-likeness, scaffold, drug-likeness, safety filters, ADMET, and aggregator.
   */
  getScoring: async (
    molecule: string,
    format: string = 'auto',
    include?: ('ml_readiness' | 'np_likeness' | 'scaffold' | 'druglikeness' | 'safety_filters' | 'admet' | 'aggregator')[]
  ): Promise<ScoringResponse> => {
    try {
      const request: ScoringRequest = {
        molecule,
        format: format as ScoringRequest['format'],
        include: include || [
          'ml_readiness',
          'np_likeness',
          'scaffold',
          'druglikeness',
          'safety_filters',
          'admet',
          'aggregator'
        ]
      };
      const response = await api.post<ScoringResponse>('/score', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<ScoringError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  }
};

export const standardizationApi = {
  /**
   * Standardize a molecule using ChEMBL-compatible pipeline.
   */
  standardize: async (request: StandardizeRequest): Promise<StandardizeResponse> => {
    try {
      const response = await api.post<StandardizeResponse>('/standardize', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<StandardizeError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  /**
   * Get available standardization options and pipeline steps.
   */
  getOptions: async (): Promise<StandardizeOptionsResponse> => {
    const response = await api.get<StandardizeOptionsResponse>('/standardize/options');
    return response.data;
  }
};

export const batchApi = {
  /**
   * Upload a file for batch processing.
   *
   * @param file - SDF or CSV file to upload
   * @param smilesColumn - Column name for SMILES (CSV only)
   * @param nameColumn - Column name for molecule Name/ID (CSV only)
   * @param onUploadProgress - Optional callback for upload progress
   * @returns Job ID and initial status
   */
  uploadBatch: async (
    file: File,
    smilesColumn?: string,
    nameColumn?: string,
    onUploadProgress?: (progress: number) => void
  ): Promise<BatchUploadResponse> => {
    const formData = new FormData();
    formData.append('file', file);
    if (smilesColumn) {
      formData.append('smiles_column', smilesColumn);
    }
    if (nameColumn) {
      formData.append('name_column', nameColumn);
    }

    const response = await api.post<BatchUploadResponse>('/batch/upload', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      timeout: 600000, // 10 minute timeout for large file upload
      onUploadProgress: (progressEvent) => {
        if (onUploadProgress && progressEvent.total) {
          const percent = Math.round((progressEvent.loaded * 100) / progressEvent.total);
          onUploadProgress(percent);
        }
      },
    });
    return response.data;
  },

  /**
   * Get batch job results with optional filtering and pagination.
   */
  getBatchResults: async (
    jobId: string,
    page: number = 1,
    pageSize: number = 50,
    filters?: BatchResultsFilters
  ): Promise<BatchResultsResponse> => {
    const params = new URLSearchParams({
      page: String(page),
      page_size: String(pageSize),
    });

    if (filters?.status_filter) {
      params.append('status_filter', filters.status_filter);
    }
    if (filters?.min_score !== undefined) {
      params.append('min_score', String(filters.min_score));
    }
    if (filters?.max_score !== undefined) {
      params.append('max_score', String(filters.max_score));
    }

    const response = await api.get<BatchResultsResponse>(
      `/batch/${jobId}?${params.toString()}`
    );
    return response.data;
  },

  /**
   * Get statistics only for a batch job.
   */
  getBatchStats: async (jobId: string): Promise<BatchStatistics> => {
    const response = await api.get<BatchStatistics>(`/batch/${jobId}/stats`);
    return response.data;
  },

  /**
   * Cancel a running batch job.
   */
  cancelBatch: async (jobId: string): Promise<void> => {
    await api.delete(`/batch/${jobId}`);
  },

  /**
   * Detect columns in a CSV file locally (no server call).
   * Reads only the first few KB to get header and samples.
   */
  detectColumnsLocal: async (file: File): Promise<CSVColumnsResponse> => {
    return new Promise((resolve, reject) => {
      // Read only first 50KB for column detection
      const SAMPLE_SIZE = 50 * 1024;
      const slice = file.slice(0, SAMPLE_SIZE);
      const reader = new FileReader();

      reader.onload = (e) => {
        try {
          const text = e.target?.result as string;
          const lines = text.split('\n').filter(line => line.trim());

          if (lines.length < 1) {
            reject(new Error('CSV file appears to be empty'));
            return;
          }

          // Parse header (handle quoted columns)
          const header = parseCSVLine(lines[0]);
          const columns = header.filter(col => col && col.length <= 256);

          if (columns.length === 0) {
            reject(new Error('No valid columns found'));
            return;
          }

          // Get sample values from first few data rows
          const columnSamples: Record<string, string> = {};
          const sampleRows = lines.slice(1, 6); // Up to 5 sample rows

          for (const col of columns.slice(0, 20)) {
            const colIndex = header.indexOf(col);
            for (const row of sampleRows) {
              const values = parseCSVLine(row);
              if (values[colIndex] && values[colIndex].trim()) {
                columnSamples[col] = values[colIndex].trim().slice(0, 100);
                break;
              }
            }
          }

          // Detect SMILES column
          const smilesKeywords = ['smiles', 'smi', 'canonical_smiles', 'isomeric_smiles', 'structure'];
          let suggestedSmiles: string | null = null;
          for (const col of columns) {
            const colLower = col.toLowerCase();
            if (smilesKeywords.some(kw => colLower.includes(kw))) {
              suggestedSmiles = col;
              break;
            }
          }

          // Detect Name/ID column
          const nameKeywords = ['name', 'id', 'compound', 'molecule', 'title', 'identifier'];
          let suggestedName: string | null = null;
          for (const col of columns) {
            if (col === suggestedSmiles) continue;
            const colLower = col.toLowerCase();
            if (nameKeywords.some(kw => colLower.includes(kw))) {
              suggestedName = col;
              break;
            }
          }

          // Estimate row count from file size
          const fileSizeMb = file.size / (1024 * 1024);
          let rowCountEstimate: number;
          if (lines.length > 2 && text.length < SAMPLE_SIZE) {
            // Small file - we have all rows
            rowCountEstimate = lines.length - 1;
          } else {
            // Estimate based on average line length from sample
            const avgLineLength = text.length / lines.length;
            rowCountEstimate = Math.max(1, Math.floor(file.size / avgLineLength) - 1);
          }

          resolve({
            columns,
            suggested_smiles: suggestedSmiles,
            suggested_name: suggestedName,
            column_samples: columnSamples,
            row_count_estimate: rowCountEstimate,
            file_size_mb: Math.round(fileSizeMb * 100) / 100,
          });
        } catch (err) {
          reject(new Error('Failed to parse CSV file'));
        }
      };

      reader.onerror = () => reject(new Error('Failed to read file'));
      reader.readAsText(slice);
    });
  },

  /**
   * Extract only selected columns from CSV and create a new file.
   * Uses streaming to handle large files efficiently.
   */
  extractColumns: async (
    file: File,
    smilesColumn: string,
    nameColumn?: string
  ): Promise<File> => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();

      reader.onload = (e) => {
        try {
          const text = e.target?.result as string;
          const lines = text.split('\n');

          if (lines.length < 1) {
            reject(new Error('CSV file is empty'));
            return;
          }

          // Parse header to find column indices
          const header = parseCSVLine(lines[0]);
          const smilesIndex = header.findIndex(
            col => col.toLowerCase() === smilesColumn.toLowerCase()
          );

          if (smilesIndex === -1) {
            reject(new Error(`SMILES column '${smilesColumn}' not found`));
            return;
          }

          const nameIndex = nameColumn
            ? header.findIndex(col => col.toLowerCase() === nameColumn.toLowerCase())
            : -1;

          // Build new CSV with only selected columns
          const outputLines: string[] = [];

          // Header
          const newHeader = nameIndex !== -1
            ? `"SMILES","Name"`
            : `"SMILES"`;
          outputLines.push(newHeader);

          // Data rows
          for (let i = 1; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line) continue;

            const values = parseCSVLine(line);
            const smiles = values[smilesIndex] || '';

            if (!smiles.trim()) continue; // Skip empty SMILES

            const name = nameIndex !== -1 ? (values[nameIndex] || '') : '';

            // Escape values for CSV
            const escapedSmiles = escapeCSVValue(smiles);
            const newLine = nameIndex !== -1
              ? `${escapedSmiles},${escapeCSVValue(name)}`
              : escapedSmiles;

            outputLines.push(newLine);
          }

          const csvContent = outputLines.join('\n');
          const blob = new Blob([csvContent], { type: 'text/csv' });
          const newFile = new File([blob], file.name, { type: 'text/csv' });

          resolve(newFile);
        } catch (err) {
          reject(new Error('Failed to process CSV file'));
        }
      };

      reader.onerror = () => reject(new Error('Failed to read file'));
      reader.readAsText(file);
    });
  },

  /**
   * Detect columns in a CSV file for SMILES selection (server-side fallback).
   */
  detectColumns: async (file: File): Promise<CSVColumnsResponse> => {
    const formData = new FormData();
    formData.append('file', file);

    const response = await api.post<CSVColumnsResponse>(
      '/batch/detect-columns',
      formData,
      {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
        timeout: 120000, // 2 minute timeout for large file analysis
      }
    );
    return response.data;
  },
};

export const integrationsApi = {
  /**
   * Look up molecule in PubChem database.
   */
  lookupPubChem: async (request: IntegrationRequest): Promise<PubChemResult> => {
    try {
      const response = await api.post<PubChemResult>('/integrations/pubchem/lookup', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<IntegrationError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  /**
   * Look up molecule in ChEMBL for bioactivity data.
   */
  lookupChEMBL: async (request: IntegrationRequest): Promise<ChEMBLResult> => {
    try {
      const response = await api.post<ChEMBLResult>('/integrations/chembl/bioactivity', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<IntegrationError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  /**
   * Look up molecule in COCONUT natural products database.
   */
  lookupCOCONUT: async (request: IntegrationRequest): Promise<COCONUTResult> => {
    try {
      const response = await api.post<COCONUTResult>('/integrations/coconut/lookup', request);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        const axiosError = error as AxiosError<IntegrationError>;
        throw axiosError.response?.data || { error: 'Network error' };
      }
      throw error;
    }
  },

  /**
   * Look up molecule in all databases in parallel.
   */
  lookupAll: async (request: IntegrationRequest): Promise<{
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
  }> => {
    const [pubchem, chembl, coconut] = await Promise.allSettled([
      integrationsApi.lookupPubChem(request),
      integrationsApi.lookupChEMBL(request),
      integrationsApi.lookupCOCONUT(request),
    ]);

    return {
      pubchem: pubchem.status === 'fulfilled' ? pubchem.value : null,
      chembl: chembl.status === 'fulfilled' ? chembl.value : null,
      coconut: coconut.status === 'fulfilled' ? coconut.value : null,
    };
  },
};

export type { ValidationRequest, ValidationResponse, ValidationError, ChecksResponse };
export type { AlertScreenRequest, AlertScreenResponse, AlertError, CatalogListResponse };
export type { ScoringRequest, ScoringResponse, ScoringError };
export type { StandardizeRequest, StandardizeResponse, StandardizeError, StandardizeOptionsResponse };
export type { BatchUploadResponse, BatchResultsResponse, BatchStatistics, BatchResultsFilters, CSVColumnsResponse };
export type { IntegrationRequest, PubChemResult, ChEMBLResult, COCONUTResult, IntegrationError };
