import { useState, useCallback } from 'react';
import type { SubScoreDetail, DatasetWeightProfile } from '../types/dataset_intelligence';
import { DEFAULT_WEIGHTS, PRESET_WEIGHT_IDS } from '../types/dataset_intelligence';

// =============================================================================
// Constants
// =============================================================================

const STORAGE_KEY = 'chemaudit_dataset_weight_profiles';

const DEFAULT_PROFILE: DatasetWeightProfile = {
  id: 'default_literature',
  name: 'Literature Default (Fourches 2010)',
  weights: { ...DEFAULT_WEIGHTS },
};

// =============================================================================
// localStorage helpers
// =============================================================================

function loadCustomProfiles(): DatasetWeightProfile[] {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed as DatasetWeightProfile[];
  } catch {
    return [];
  }
}

function saveCustomProfiles(profiles: DatasetWeightProfile[]): void {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(profiles));
  } catch {
    // localStorage quota exceeded — ignore
  }
}

// =============================================================================
// Types
// =============================================================================

export interface UseDatasetWeightsReturn {
  /** Current weight values (always sum to 1.0). */
  weights: Record<string, number>;
  /** ID of the active weight profile, or null for custom. */
  activeProfileId: string | null;
  /** All saved profiles (preset first, then custom). */
  savedProfiles: DatasetWeightProfile[];
  /**
   * Update a single weight and normalize others proportionally.
   * Per Pitfall 8: synchronous in onChange, not useEffect.
   */
  updateWeight: (key: string, value: number) => void;
  /** Compute overall score from sub-scores using current weights. */
  computeOverallScore: (subScores: SubScoreDetail[]) => number;
  /** Save current weights as a named profile. */
  saveProfile: (name: string) => void;
  /** Delete a custom profile by id. No-op for preset IDs. */
  deleteProfile: (id: string) => void;
  /** Load weights from a saved profile. */
  loadProfile: (id: string) => void;
  /** Reset weights to literature defaults. */
  resetToDefaults: () => void;
}

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook for managing dataset health audit weight sliders.
 *
 * Features:
 * - Weight normalization: changing one slider redistributes remaining 1.0
 *   proportionally across other sliders (Pitfall 8: synchronous in onChange)
 * - localStorage persistence for custom weight profiles
 * - Preset guard via PRESET_WEIGHT_IDS Set (same pattern as useMPOProfiles)
 * - Client-side score recomputation without server call
 */
export function useDatasetWeights(): UseDatasetWeightsReturn {
  const [weights, setWeights] = useState<Record<string, number>>({ ...DEFAULT_WEIGHTS });
  const [customProfiles, setCustomProfiles] = useState<DatasetWeightProfile[]>(loadCustomProfiles);
  const [activeProfileId, setActiveProfileId] = useState<string | null>('default_literature');

  // All profiles: preset first, then custom
  const savedProfiles: DatasetWeightProfile[] = [DEFAULT_PROFILE, ...customProfiles];

  // ==========================================================================
  // updateWeight — with proportional normalization
  // ==========================================================================

  const updateWeight = useCallback((key: string, value: number) => {
    setWeights((prev) => {
      const clamped = Math.max(0, Math.min(1, value));
      const otherKeys = Object.keys(prev).filter((k) => k !== key);
      const othersSum = otherKeys.reduce((sum, k) => sum + prev[k], 0);
      const remaining = 1.0 - clamped;

      const next: Record<string, number> = { ...prev, [key]: clamped };

      if (othersSum === 0) {
        // All others are zero — distribute evenly
        const even = remaining / otherKeys.length;
        for (const k of otherKeys) {
          next[k] = even;
        }
      } else {
        // Redistribute proportionally
        for (const k of otherKeys) {
          next[k] = (prev[k] / othersSum) * remaining;
        }
      }

      return next;
    });
    // Custom weights — no longer tracking a preset
    setActiveProfileId(null);
  }, []);

  // ==========================================================================
  // computeOverallScore
  // ==========================================================================

  const computeOverallScore = useCallback(
    (subScores: SubScoreDetail[]): number => {
      let total = 0;
      for (const sub of subScores) {
        const w = weights[sub.name] ?? 0;
        total += sub.score * w;
      }
      // sub.score is 0-1, multiply by 100 for 0-100 scale
      return Math.max(0, Math.min(100, total * 100));
    },
    [weights],
  );

  // ==========================================================================
  // Profile management
  // ==========================================================================

  const saveProfile = useCallback(
    (name: string) => {
      const id = `custom_${Date.now()}`;
      const newProfile: DatasetWeightProfile = {
        id,
        name,
        weights: { ...weights },
      };
      setCustomProfiles((prev) => {
        const next = [...prev, newProfile];
        saveCustomProfiles(next);
        return next;
      });
      setActiveProfileId(id);
    },
    [weights],
  );

  const deleteProfile = useCallback((id: string) => {
    if (PRESET_WEIGHT_IDS.has(id)) return; // Cannot delete presets

    setCustomProfiles((prev) => {
      const next = prev.filter((p) => p.id !== id);
      saveCustomProfiles(next);
      return next;
    });
    setActiveProfileId((prev) => (prev === id ? null : prev));
  }, []);

  const loadProfile = useCallback(
    (id: string) => {
      const profile = savedProfiles.find((p) => p.id === id);
      if (!profile) return;
      setWeights({ ...profile.weights });
      setActiveProfileId(id);
    },
    [savedProfiles],
  );

  const resetToDefaults = useCallback(() => {
    setWeights({ ...DEFAULT_WEIGHTS });
    setActiveProfileId('default_literature');
  }, []);

  return {
    weights,
    activeProfileId,
    savedProfiles,
    updateWeight,
    computeOverallScore,
    saveProfile,
    deleteProfile,
    loadProfile,
    resetToDefaults,
  };
}
