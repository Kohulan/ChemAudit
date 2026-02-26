/**
 * Hook for managing scoring profiles with caching.
 *
 * Fetches profiles on mount, separates presets from user profiles,
 * and provides CRUD + duplicate actions with automatic re-fetch.
 */

import { useState, useEffect, useCallback, useRef } from 'react';
import { profilesApi } from '../services/api';
import type {
  ScoringProfile,
  ScoringProfileCreate,
  ScoringProfileUpdate,
  ScoringProfileExport,
} from '../types/workflow';

const CACHE_KEY = 'chemaudit_profiles_cache';
const CACHE_TTL_MS = 5 * 60 * 1000; // 5 minutes

interface ProfileCache {
  data: ScoringProfile[];
  timestamp: number;
}

function getCachedProfiles(): ScoringProfile[] | null {
  try {
    const raw = localStorage.getItem(CACHE_KEY);
    if (!raw) return null;
    const cache: ProfileCache = JSON.parse(raw);
    if (Date.now() - cache.timestamp > CACHE_TTL_MS) {
      localStorage.removeItem(CACHE_KEY);
      return null;
    }
    return cache.data;
  } catch {
    return null;
  }
}

function setCacheProfiles(data: ScoringProfile[]) {
  try {
    const cache: ProfileCache = { data, timestamp: Date.now() };
    localStorage.setItem(CACHE_KEY, JSON.stringify(cache));
  } catch {
    // localStorage quota exceeded — ignore
  }
}

export interface UseProfilesReturn {
  profiles: ScoringProfile[];
  presets: ScoringProfile[];
  userProfiles: ScoringProfile[];
  activeProfile: ScoringProfile | null;
  isLoading: boolean;
  error: string | null;
  refetch: () => Promise<void>;
  createProfile: (data: ScoringProfileCreate) => Promise<ScoringProfile>;
  updateProfile: (id: number, data: ScoringProfileUpdate) => Promise<ScoringProfile>;
  deleteProfile: (id: number) => Promise<void>;
  duplicateProfile: (id: number, newName: string) => Promise<ScoringProfile>;
  exportProfile: (id: number) => Promise<ScoringProfileExport>;
  importProfile: (data: ScoringProfileExport) => Promise<ScoringProfile>;
}

export function useProfiles(): UseProfilesReturn {
  const [profiles, setProfiles] = useState<ScoringProfile[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const mountedRef = useRef(true);

  const fetchProfiles = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    try {
      const data = await profilesApi.getProfiles();
      if (mountedRef.current) {
        setProfiles(data);
        setCacheProfiles(data);
      }
    } catch (err) {
      if (mountedRef.current) {
        setError(err instanceof Error ? err.message : 'Failed to load profiles');
      }
    } finally {
      if (mountedRef.current) {
        setIsLoading(false);
      }
    }
  }, []);

  // Load from cache first, then fetch fresh data
  useEffect(() => {
    mountedRef.current = true;
    const cached = getCachedProfiles();
    if (cached) {
      setProfiles(cached);
      setIsLoading(false);
    }
    fetchProfiles();
    return () => {
      mountedRef.current = false;
    };
  }, [fetchProfiles]);

  const presets = profiles.filter((p) => p.is_preset);
  const userProfiles = profiles.filter((p) => !p.is_preset);
  const activeProfile = profiles.find((p) => p.is_active) ?? null;

  const createProfile = useCallback(
    async (data: ScoringProfileCreate) => {
      const created = await profilesApi.createProfile(data);
      await fetchProfiles();
      return created;
    },
    [fetchProfiles]
  );

  const updateProfile = useCallback(
    async (id: number, data: ScoringProfileUpdate) => {
      const updated = await profilesApi.updateProfile(id, data);
      await fetchProfiles();
      return updated;
    },
    [fetchProfiles]
  );

  const deleteProfile = useCallback(
    async (id: number) => {
      await profilesApi.deleteProfile(id);
      await fetchProfiles();
    },
    [fetchProfiles]
  );

  const duplicateProfile = useCallback(
    async (id: number, newName: string) => {
      const dup = await profilesApi.duplicateProfile(id, newName);
      await fetchProfiles();
      return dup;
    },
    [fetchProfiles]
  );

  const exportProfile = useCallback(async (id: number) => {
    return profilesApi.exportProfile(id);
  }, []);

  const importProfile = useCallback(
    async (data: ScoringProfileExport) => {
      const imported = await profilesApi.importProfile(data);
      await fetchProfiles();
      return imported;
    },
    [fetchProfiles]
  );

  return {
    profiles,
    presets,
    userProfiles,
    activeProfile,
    isLoading,
    error,
    refetch: fetchProfiles,
    createProfile,
    updateProfile,
    deleteProfile,
    duplicateProfile,
    exportProfile,
    importProfile,
  };
}
