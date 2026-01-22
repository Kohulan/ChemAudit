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

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api/v1';

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

export type { ValidationRequest, ValidationResponse, ValidationError, ChecksResponse };
export type { AlertScreenRequest, AlertScreenResponse, AlertError, CatalogListResponse };
export type { ScoringRequest, ScoringResponse, ScoringError };
