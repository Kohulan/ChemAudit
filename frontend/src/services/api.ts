import axios, { AxiosError } from 'axios';
import type {
  ValidationRequest,
  ValidationResponse,
  ValidationError,
  ChecksResponse
} from '../types/validation';

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

export type { ValidationRequest, ValidationResponse, ValidationError, ChecksResponse };
