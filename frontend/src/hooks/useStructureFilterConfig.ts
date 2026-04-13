import { useState, useCallback } from 'react';
import type { FilterConfig, StructureFilterProfile } from '../types/structure_filter';
import { PRESET_CONFIGS, PRESET_IDS } from '../types/structure_filter';

// =============================================================================
// Constants
// =============================================================================

const STORAGE_KEY = 'chemaudit_structure_filter_configs';

// =============================================================================
// localStorage helpers
// =============================================================================

function loadCustomProfiles(): StructureFilterProfile[] {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed as StructureFilterProfile[];
  } catch {
    return [];
  }
}

function saveCustomProfiles(profiles: StructureFilterProfile[]): void {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(profiles));
  } catch {
    // localStorage quota exceeded — ignore
  }
}

// =============================================================================
// Return type
// =============================================================================

export interface UseStructureFilterConfigReturn {
  /** Current active filter configuration. */
  config: FilterConfig;
  /** ID of the currently active preset, or null if a custom config is active. */
  activePreset: string | null;
  /** Name of the preset this config was modified from (if any). */
  modifiedFrom: string | null;
  /** All custom profiles saved to localStorage. */
  profiles: StructureFilterProfile[];
  /** Activate a preset by ID (drug_like, lead_like, fragment_like, permissive, or custom). */
  selectPreset: (id: string) => void;
  /** Apply a partial config update. Clears activePreset and tracks modifiedFrom. */
  updateConfig: (partial: Partial<FilterConfig>) => void;
  /** Save current config as a named custom profile. */
  saveProfile: (name: string) => void;
  /** Delete a custom profile by ID. No-op for built-in preset IDs. */
  deleteProfile: (id: string) => void;
}

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook for managing Structure Filter configuration with localStorage
 * persistence for custom profiles.
 *
 * Four built-in presets (drug_like, lead_like, fragment_like, permissive)
 * are always available and cannot be deleted. Custom profiles are stored
 * under `chemaudit_structure_filter_configs` in localStorage.
 *
 * When a field changes while a preset is active, the preset is cleared and
 * `modifiedFrom` tracks which preset was the base — matching the
 * useQSARPipelineConfig pattern from Phase 10.
 */
export function useStructureFilterConfig(): UseStructureFilterConfigReturn {
  // Initialize with drug_like preset (D-12 default)
  const [config, setConfig] = useState<FilterConfig>(() => ({ ...PRESET_CONFIGS.drug_like }));
  const [activePreset, setActivePreset] = useState<string | null>('drug_like');
  const [modifiedFrom, setModifiedFrom] = useState<string | null>(null);
  const [profiles, setProfiles] = useState<StructureFilterProfile[]>(loadCustomProfiles);

  const selectPreset = useCallback(
    (id: string) => {
      if (PRESET_IDS.has(id)) {
        // Built-in preset
        setConfig({ ...PRESET_CONFIGS[id] });
        setActivePreset(id);
        setModifiedFrom(null);
      } else {
        // Custom profile — load from state
        const found = profiles.find((p) => p.id === id);
        if (found) {
          setConfig({ ...found.config });
          setActivePreset(id);
          setModifiedFrom(null);
        }
      }
    },
    [profiles],
  );

  const updateConfig = useCallback(
    (partial: Partial<FilterConfig>) => {
      setConfig((prev) => ({ ...prev, ...partial }));
      // When editing while on a preset, track what it was modified from
      setActivePreset((prevPreset) => {
        if (prevPreset !== null) {
          setModifiedFrom(prevPreset);
        }
        return null;
      });
    },
    [],
  );

  const saveProfile = useCallback(
    (name: string) => {
      const id = name.toLowerCase().replace(/\s+/g, '_') + '_' + Date.now();
      const newProfile: StructureFilterProfile = { id, name, config: { ...config } };
      setProfiles((prev) => {
        const next = [...prev, newProfile];
        saveCustomProfiles(next);
        return next;
      });
      setActivePreset(id);
      setModifiedFrom(null);
    },
    [config],
  );

  const deleteProfile = useCallback(
    (id: string) => {
      if (PRESET_IDS.has(id)) return; // No-op for built-in presets (Phase 10 pattern)

      setProfiles((prev) => {
        const next = prev.filter((p) => p.id !== id);
        saveCustomProfiles(next);
        return next;
      });

      // If the deleted profile was active, reset to default preset
      setActivePreset((prev) => {
        if (prev === id) {
          setConfig({ ...PRESET_CONFIGS.drug_like });
          setModifiedFrom(null);
          return 'drug_like';
        }
        return prev;
      });
    },
    [],
  );

  return {
    config,
    activePreset,
    modifiedFrom,
    profiles,
    selectPreset,
    updateConfig,
    saveProfile,
    deleteProfile,
  };
}
