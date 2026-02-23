import { useLocalStorage } from './useLocalStorage';
import type { DeepValidationConfig, SeverityOverride } from '../types/validation';

const DEFAULT_CONFIG: DeepValidationConfig = {
  severityOverrides: {},
};

/**
 * Hook for managing per-check deep validation severity overrides.
 * Persists configuration to localStorage under 'chemaudit:deep-validation-config'.
 */
export function useDeepValidationConfig() {
  const [config, setConfig, clearConfig] = useLocalStorage<DeepValidationConfig>(
    'chemaudit:deep-validation-config',
    DEFAULT_CONFIG
  );

  const setSeverityOverride = (checkName: string, severity: SeverityOverride) => {
    setConfig((prev: DeepValidationConfig) => ({
      ...prev,
      severityOverrides: { ...prev.severityOverrides, [checkName]: severity },
    }));
  };

  const removeSeverityOverride = (checkName: string) => {
    setConfig((prev: DeepValidationConfig) => {
      const { [checkName]: _, ...rest } = prev.severityOverrides;
      return { ...prev, severityOverrides: rest };
    });
  };

  const resetAllOverrides = () => {
    setConfig(DEFAULT_CONFIG);
  };

  const getEffectiveSeverity = (checkName: string, defaultSeverity: string): string => {
    return config.severityOverrides[checkName] || defaultSeverity;
  };

  return {
    config,
    setSeverityOverride,
    removeSeverityOverride,
    resetAllOverrides,
    getEffectiveSeverity,
    clearConfig,
  };
}
