import { useState, useCallback } from 'react';
import { safetyApi } from '../services/api';
import type { AlertScreenResponse, SafetyAssessResponse } from '../types/safety';

/**
 * Hook for safety screening API calls.
 *
 * Fires both alert screening and safety assessment in parallel via
 * Promise.allSettled so errors are independent — one can fail without
 * blocking the other (per UI-SPEC).
 *
 * Uses safetyApi from services/api.ts (which uses the project's axios
 * instance with CSRF tokens, auth headers, and error interceptors).
 */
export function useSafety() {
  const [alertResult, setAlertResult] = useState<AlertScreenResponse | null>(null);
  const [safetyResult, setSafetyResult] = useState<SafetyAssessResponse | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [alertError, setAlertError] = useState<string | null>(null);
  const [safetyError, setSafetyError] = useState<string | null>(null);

  const screenMolecule = useCallback(async (smiles: string) => {
    setIsLoading(true);
    setAlertError(null);
    setSafetyError(null);
    setAlertResult(null);
    setSafetyResult(null);

    // Fire both requests in parallel — errors are independent (per UI-SPEC)
    const [alertRes, safetyRes] = await Promise.allSettled([
      safetyApi.screen(smiles),
      safetyApi.assess(smiles),
    ]);

    if (alertRes.status === 'fulfilled') {
      setAlertResult(alertRes.value);
    } else {
      const reason = alertRes.reason as { error?: string; detail?: string };
      setAlertError(reason?.error ?? reason?.detail ?? 'Alert screening failed');
    }

    if (safetyRes.status === 'fulfilled') {
      setSafetyResult(safetyRes.value);
    } else {
      const reason = safetyRes.reason as { error?: string; detail?: string };
      setSafetyError(reason?.error ?? reason?.detail ?? 'Safety assessment failed');
    }

    setIsLoading(false);
  }, []);

  return {
    alertResult,
    safetyResult,
    isLoading,
    alertError,
    safetyError,
    screenMolecule,
  };
}
