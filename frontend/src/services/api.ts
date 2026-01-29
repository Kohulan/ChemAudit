import axios, { AxiosError } from 'axios';
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
const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api/v1';
export const API_DOCS_URL = import.meta.env.VITE_API_DOCS_URL || 'http://localhost:8000/docs';
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
   * Calculate ML-readiness, NP-likeness, and scaffold scores for a molecule.
   */
  getScoring: async (
    molecule: string,
    format: string = 'auto',
    include?: ('ml_readiness' | 'np_likeness' | 'scaffold')[]
  ): Promise<ScoringResponse> => {
    try {
      const request: ScoringRequest = {
        molecule,
        format: format as ScoringRequest['format'],
        include: include || ['ml_readiness', 'np_likeness', 'scaffold']
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
   * @returns Job ID and initial status
   */
  uploadBatch: async (file: File, smilesColumn?: string): Promise<BatchUploadResponse> => {
    const formData = new FormData();
    formData.append('file', file);
    if (smilesColumn) {
      formData.append('smiles_column', smilesColumn);
    }

    const response = await api.post<BatchUploadResponse>('/batch/upload', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      timeout: 600000, // 10 minute timeout for large file upload
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
   * Detect columns in a CSV file for SMILES selection.
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
