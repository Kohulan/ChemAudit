import { useState, useCallback } from 'react';
import type { QSARReadyConfig, QSARPipelineProfile } from '../types/qsar_ready';
import { PRESET_CONFIGS } from '../types/qsar_ready';

// =============================================================================
// Constants
// =============================================================================

const STORAGE_KEY = 'chemaudit_qsar_configs';
const PRESET_IDS = new Set(['qsar-2d', 'qsar-3d', 'minimal']);

// =============================================================================
// localStorage helpers
// =============================================================================

function loadCustomProfiles(): QSARPipelineProfile[] {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed as QSARPipelineProfile[];
  } catch {
    return [];
  }
}

function saveCustomProfiles(profiles: QSARPipelineProfile[]): void {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(profiles));
  } catch {
    // localStorage quota exceeded — ignore
  }
}

// =============================================================================
// Return type
// =============================================================================

export interface UseQSARPipelineConfigReturn {
  /** Current active pipeline configuration. */
  config: QSARReadyConfig;
  /** ID of the currently active preset, or null if a custom config is active. */
  activePreset: string | null;
  /** Name of the preset this config was modified from (if any). */
  modifiedFrom: string | null;
  /** All custom profiles saved to localStorage. */
  profiles: QSARPipelineProfile[];
  /** Activate a preset by ID (qsar-2d, qsar-3d, minimal, or custom profile). */
  selectPreset: (id: string) => void;
  /** Apply a partial config update. Clears activePreset and tracks modifiedFrom. */
  updateConfig: (partial: Partial<QSARReadyConfig>) => void;
  /** Save current config as a named custom profile. */
  saveProfile: (name: string) => void;
  /** Delete a custom profile by ID. No-op for preset IDs. */
  deleteProfile: (id: string) => void;
  /** Load a custom profile config by ID. */
  loadProfile: (id: string) => void;
  /** Whether a given ID is a built-in preset (cannot be deleted). */
  isPreset: (id: string) => boolean;
}

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook for managing QSAR-ready pipeline configuration with localStorage
 * persistence for custom profiles.
 *
 * Three built-in presets (qsar-2d, qsar-3d, minimal) are always available
 * and cannot be deleted. Custom profiles are stored under
 * `chemaudit_qsar_configs` in localStorage.
 *
 * When a step toggle changes while a preset is active, the preset is cleared
 * and `modifiedFrom` tracks which preset was the base — matching the
 * useMPOProfiles pattern from Phase 07.
 */
export function useQSARPipelineConfig(): UseQSARPipelineConfigReturn {
  const [config, setConfig] = useState<QSARReadyConfig>(() => PRESET_CONFIGS['qsar-2d']);
  const [activePreset, setActivePreset] = useState<string | null>('qsar-2d');
  const [modifiedFrom, setModifiedFrom] = useState<string | null>(null);
  const [customProfiles, setCustomProfiles] = useState<QSARPipelineProfile[]>(loadCustomProfiles);

  const selectPreset = useCallback((id: string) => {
    if (PRESET_IDS.has(id)) {
      // Built-in preset
      setConfig(PRESET_CONFIGS[id]);
      setActivePreset(id);
      setModifiedFrom(null);
    } else {
      // Custom profile — load from state
      const found = customProfiles.find((p) => p.id === id);
      if (found) {
        setConfig(found.config);
        setActivePreset(id);
        setModifiedFrom(null);
      }
    }
  }, [customProfiles]);

  const updateConfig = useCallback((partial: Partial<QSARReadyConfig>) => {
    setConfig((prev) => ({ ...prev, ...partial }));
    // When editing while on a preset, track what it was modified from
    setActivePreset((prevPreset) => {
      if (prevPreset !== null) {
        setModifiedFrom(prevPreset);
      }
      return null;
    });
  }, []);

  const saveProfile = useCallback((name: string) => {
    const newProfile: QSARPipelineProfile = {
      id: String(Date.now()),
      name,
      config,
      createdAt: Date.now(),
    };
    setCustomProfiles((prev) => {
      const next = [...prev, newProfile];
      saveCustomProfiles(next);
      return next;
    });
    // Switch to the newly created profile
    setActivePreset(newProfile.id);
    setModifiedFrom(null);
  }, [config]);

  const deleteProfile = useCallback((id: string) => {
    if (PRESET_IDS.has(id)) return; // Cannot delete built-in presets

    setCustomProfiles((prev) => {
      const next = prev.filter((p) => p.id !== id);
      saveCustomProfiles(next);
      return next;
    });

    // If the deleted profile was active, reset to default preset
    setActivePreset((prev) => {
      if (prev === id) {
        setConfig(PRESET_CONFIGS['qsar-2d']);
        setModifiedFrom(null);
        return 'qsar-2d';
      }
      return prev;
    });
  }, []);

  const loadProfile = useCallback((id: string) => {
    selectPreset(id);
  }, [selectPreset]);

  const isPreset = useCallback((id: string) => PRESET_IDS.has(id), []);

  return {
    config,
    activePreset,
    modifiedFrom,
    profiles: customProfiles,
    selectPreset,
    updateConfig,
    saveProfile,
    deleteProfile,
    loadProfile,
    isPreset,
  };
}
