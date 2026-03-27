import { useState, useCallback, useRef, useEffect } from 'react';
import type { FilterConfig, FilterResult } from '../types/genchem';
import { genchemApi } from '../services/api';

// =============================================================================
// Types
// =============================================================================

/**
 * State machine for the filter execution lifecycle.
 * - idle: no filter running or completed
 * - loading-sync: synchronous filter in progress (<=1000 SMILES)
 * - loading-async: async batch upload in progress (>1000 SMILES)
 * - processing: batch job running on server (WebSocket/polling active)
 * - success: filter completed with results available
 * - error: filter failed with error message
 */
export type FilterState = 'idle' | 'loading-sync' | 'loading-async' | 'processing' | 'success' | 'error';

export interface UseGenChemFilterReturn {
  state: FilterState;
  result: FilterResult | null;
  error: string | null;
  jobId: string | null;
  progress: number | null;
  currentStage: string | null;
  selectedStage: string | null;
  setSelectedStage: (stage: string | null) => void;
  runFilter: (smilesList: string[], config: FilterConfig, preset?: string) => Promise<void>;
  reset: () => void;
}

// =============================================================================
// Constants
// =============================================================================

/** Threshold for sync vs async execution per D-22. */
const ASYNC_THRESHOLD = 1000;

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook for GenChem Filter execution state management.
 *
 * Manages:
 * - Synchronous filter for <=1000 SMILES (runFilter → genchemApi.filter)
 * - Async batch for >1000 SMILES (runFilter → genchemApi.batchUpload → WebSocket/poll)
 * - WebSocket progress tracking with polling fallback on ws.onerror
 * - Stage selection for funnel chart drill-down
 *
 * The sync/async split at ASYNC_THRESHOLD=1000 follows the D-22 locked decision.
 * WebSocket + polling fallback follows the useQSARReady pattern from Phase 10.
 */
export function useGenChemFilter(): UseGenChemFilterReturn {
  const [state, setState] = useState<FilterState>('idle');
  const [result, setResult] = useState<FilterResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);
  const [progress, setProgress] = useState<number | null>(null);
  const [currentStage, setCurrentStage] = useState<string | null>(null);
  const [selectedStage, setSelectedStage] = useState<string | null>(null);
  const wsRef = useRef<WebSocket | null>(null);

  // WebSocket cleanup on unmount
  useEffect(() => {
    return () => {
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }
    };
  }, []);

  // WebSocket connection for async batch progress
  useEffect(() => {
    if (!jobId || state !== 'loading-async') return;

    // Close any existing connection before opening a new one
    if (wsRef.current) {
      wsRef.current.close();
      wsRef.current = null;
    }

    const wsBase = window.location.origin.replace(/^http/, 'ws');
    const wsUrl = `${wsBase}/ws/genchem/${jobId}`;
    const ws = new WebSocket(wsUrl);
    wsRef.current = ws;

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data) as {
          progress?: number;
          current_stage?: string;
          status?: string;
        };
        if (data.progress !== undefined) setProgress(data.progress);
        if (data.current_stage) setCurrentStage(data.current_stage);
        if (data.status === 'complete') {
          setState('processing');
          ws.close();
          wsRef.current = null;
          // Fetch results
          genchemApi
            .batchResults(jobId)
            .then((resp) => {
              if (resp.result) {
                setResult(resp.result);
                setState('success');
              } else {
                setError('Batch completed but returned no results');
                setState('error');
              }
            })
            .catch((err: unknown) => {
              const e = err as { error?: string; detail?: string };
              setError(e?.error ?? e?.detail ?? 'Failed to fetch batch results');
              setState('error');
            });
        } else if (data.status === 'failed') {
          ws.close();
          wsRef.current = null;
          setError('Batch processing failed');
          setState('error');
        }
      } catch {
        // Ignore malformed WebSocket messages
      }
    };

    ws.onerror = () => {
      // WebSocket failed — fall back to polling (Phase 10 pattern)
      ws.close();
      wsRef.current = null;

      const poll = async () => {
        try {
          const status = await genchemApi.batchStatus(jobId);
          if (status.progress !== null) setProgress(status.progress);
          if (status.current_stage) setCurrentStage(status.current_stage);

          if (status.status === 'complete') {
            const resp = await genchemApi.batchResults(jobId);
            if (resp.result) {
              setResult(resp.result);
              setState('success');
            } else {
              setError('Batch completed but returned no results');
              setState('error');
            }
          } else if (status.status === 'failed') {
            setError('Batch processing failed');
            setState('error');
          } else {
            // Still processing — poll again in 2s
            setTimeout(poll, 2000);
          }
        } catch (err: unknown) {
          const e = err as { error?: string; detail?: string };
          setError(e?.error ?? e?.detail ?? 'Status polling failed');
          setState('error');
        }
      };

      setTimeout(poll, 2000);
    };

    ws.onclose = () => {
      if (wsRef.current === ws) {
        wsRef.current = null;
      }
    };
  }, [jobId, state]);

  // ==========================================================================
  // runFilter
  // ==========================================================================

  const runFilter = useCallback(async (smilesList: string[], config: FilterConfig, preset?: string) => {
    setError(null);
    setResult(null);
    setSelectedStage(null);
    setProgress(null);
    setCurrentStage(null);

    if (smilesList.length <= ASYNC_THRESHOLD) {
      // Synchronous path: direct API call
      setState('loading-sync');
      try {
        const res = await genchemApi.filter(smilesList, preset, config);
        setResult(res);
        setState('success');
      } catch (err: unknown) {
        const e = err as { error?: string; detail?: string };
        setError(e?.error ?? e?.detail ?? 'Filtering failed');
        setState('error');
      }
    } else {
      // Async path: upload file, track via WebSocket/polling
      setState('loading-async');
      // Convert SMILES list to a text blob for file upload
      const blob = new Blob([smilesList.join('\n')], { type: 'text/plain' });
      const file = new File([blob], 'genchem_input.txt', { type: 'text/plain' });
      try {
        const resp = await genchemApi.batchUpload(file, preset, config);
        setJobId(resp.job_id);
        // WebSocket useEffect will trigger on jobId + state='loading-async'
      } catch (err: unknown) {
        const e = err as { error?: string; detail?: string };
        setError(e?.error ?? e?.detail ?? 'Batch upload failed');
        setState('error');
      }
    }
  }, []);

  // ==========================================================================
  // reset
  // ==========================================================================

  const reset = useCallback(() => {
    if (wsRef.current) {
      wsRef.current.close();
      wsRef.current = null;
    }
    setState('idle');
    setResult(null);
    setError(null);
    setJobId(null);
    setProgress(null);
    setCurrentStage(null);
    setSelectedStage(null);
  }, []);

  return {
    state,
    result,
    error,
    jobId,
    progress,
    currentStage,
    selectedStage,
    setSelectedStage,
    runFilter,
    reset,
  };
}
