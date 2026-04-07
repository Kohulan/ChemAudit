import { useState, useCallback } from 'react';
import type { MPOProperty } from '../types/profiler';

/** An MPO profile stored in localStorage. */
export interface MPOProfile {
  /** Unique identifier. Preset IDs are fixed strings; user profiles use Date.now() string. */
  id: string;
  name: string;
  properties: MPOProperty[];
  createdAt: number;
}

// =============================================================================
// Built-in presets (cannot be deleted)
// =============================================================================

/** CNS MPO preset (Wager 2010) — 4-component piecewise linear, max_score=4.0 */
const CNS_MPO_PRESET: MPOProfile = {
  id: 'cns-mpo',
  name: 'CNS MPO (Wager 2010)',
  properties: [
    { property: 'LogP', low: 1, high: 5, weight: 1, shape: 'ramp' as const },
    { property: 'TPSA', low: 20, high: 120, weight: 1, shape: 'ramp' as const },
    { property: 'MW', low: 360, high: 500, weight: 1, shape: 'ramp' as const },
    { property: 'HBD', low: 0, high: 4, weight: 1, shape: 'ramp' as const },
  ],
  createdAt: 0,
};

/** Oral Drug MPO preset — 5-component piecewise linear. */
const ORAL_DRUG_MPO_PRESET: MPOProfile = {
  id: 'oral-drug',
  name: 'Oral Drug MPO',
  properties: [
    { property: 'MW', low: 200, high: 500, weight: 1, shape: 'ramp' as const },
    { property: 'LogP', low: -0.4, high: 5.6, weight: 1, shape: 'ramp' as const },
    { property: 'TPSA', low: 20, high: 130, weight: 1, shape: 'ramp' as const },
    { property: 'HBD', low: 0, high: 5, weight: 1, shape: 'ramp' as const },
    { property: 'RotBonds', low: 0, high: 10, weight: 1, shape: 'ramp' as const },
  ],
  createdAt: 0,
};

const PRESET_IDS = new Set([CNS_MPO_PRESET.id, ORAL_DRUG_MPO_PRESET.id]);
const STORAGE_KEY = 'chemaudit_mpo_profiles';

// =============================================================================
// localStorage helpers
// =============================================================================

function loadCustomProfiles(): MPOProfile[] {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed as MPOProfile[];
  } catch {
    return [];
  }
}

function saveCustomProfiles(profiles: MPOProfile[]): void {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(profiles));
  } catch {
    // localStorage quota exceeded — ignore
  }
}

// =============================================================================
// Hook
// =============================================================================

export interface UseMPOProfilesReturn {
  /** All profiles: presets first, then user-created profiles ordered by createdAt. */
  profiles: MPOProfile[];
  /** Save or update a profile. Preset profiles cannot be overwritten. */
  saveProfile: (profile: MPOProfile) => void;
  /** Delete a profile by id. Preset profiles cannot be deleted. */
  deleteProfile: (id: string) => void;
  /** Whether a given id is a preset (cannot be modified or deleted). */
  isPreset: (id: string) => boolean;
}

/**
 * Hook for managing custom MPO profiles with localStorage persistence.
 *
 * Per D-13: uses localStorage key `chemaudit_mpo_profiles` for custom profiles.
 * Built-in presets (CNS MPO, Oral Drug) are always prepended and cannot be deleted.
 */
export function useMPOProfiles(): UseMPOProfilesReturn {
  const [customProfiles, setCustomProfiles] = useState<MPOProfile[]>(loadCustomProfiles);

  const profiles: MPOProfile[] = [CNS_MPO_PRESET, ORAL_DRUG_MPO_PRESET, ...customProfiles];

  const saveProfile = useCallback((profile: MPOProfile) => {
    if (PRESET_IDS.has(profile.id)) return; // Cannot overwrite presets

    setCustomProfiles((prev) => {
      const existingIndex = prev.findIndex((p) => p.id === profile.id);
      let next: MPOProfile[];
      if (existingIndex >= 0) {
        next = [...prev];
        next[existingIndex] = profile;
      } else {
        next = [...prev, profile];
      }
      saveCustomProfiles(next);
      return next;
    });
  }, []);

  const deleteProfile = useCallback((id: string) => {
    if (PRESET_IDS.has(id)) return; // Cannot delete presets

    setCustomProfiles((prev) => {
      const next = prev.filter((p) => p.id !== id);
      saveCustomProfiles(next);
      return next;
    });
  }, []);

  const isPreset = useCallback((id: string) => PRESET_IDS.has(id), []);

  return { profiles, saveProfile, deleteProfile, isPreset };
}
