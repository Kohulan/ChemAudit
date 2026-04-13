import { useState, useCallback } from 'react';
import axios from 'axios';
import { api } from '../services/api';
import type {
  ProfileResponse,
  Shape3DResult,
  LEResult,
  CustomMPOResult,
  MPOProperty,
} from '../types/profiler';

/**
 * Hook for compound profiling API calls.
 *
 * Provides:
 * - profileCompound: full profile via /profiler/full
 * - compute3DShape: on-demand 3D shape descriptors
 * - computeEfficiency: ligand efficiency from activity value
 * - computeCustomMPO: custom MPO score from user-defined property profile
 */
export function useProfiler() {
  const [profile, setProfile] = useState<ProfileResponse | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const profileCompound = useCallback(async (smiles: string) => {
    setIsLoading(true);
    setError(null);
    try {
      const { data } = await api.post<ProfileResponse>(
        '/profiler/full',
        { smiles }
      );
      setProfile(data);
    } catch (err) {
      setError(
        axios.isAxiosError(err)
          ? (err.response?.data?.detail?.error ?? err.response?.data?.detail ?? 'Profile failed')
          : 'Profile failed'
      );
    } finally {
      setIsLoading(false);
    }
  }, []);

  const compute3DShape = useCallback(async (smiles: string): Promise<Shape3DResult | null> => {
    try {
      const { data } = await api.post<Shape3DResult>(
        '/profiler/shape-3d',
        { smiles }
      );
      return data;
    } catch {
      return null;
    }
  }, []);

  const computeEfficiency = useCallback(async (
    smiles: string,
    activityValue: number,
    activityType: string
  ): Promise<LEResult | null> => {
    try {
      const { data } = await api.post<LEResult>(
        '/profiler/efficiency',
        { smiles, activity_value: activityValue, activity_type: activityType }
      );
      return data;
    } catch {
      return null;
    }
  }, []);

  const computeCustomMPO = useCallback(async (
    smiles: string,
    mpoProfile: MPOProperty[]
  ): Promise<CustomMPOResult | null> => {
    try {
      const { data } = await api.post<CustomMPOResult>(
        '/profiler/mpo',
        { smiles, profile: mpoProfile }
      );
      return data;
    } catch {
      return null;
    }
  }, []);

  return {
    profile,
    isLoading,
    error,
    profileCompound,
    compute3DShape,
    computeEfficiency,
    computeCustomMPO,
  };
}
