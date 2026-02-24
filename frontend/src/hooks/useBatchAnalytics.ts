/**
 * useBatchAnalytics Hook
 *
 * Triggers analytics computation and polls until all requested types complete.
 * Returns analytics data, status, error, and a retrigger function.
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { batchApi } from '../services/api';
import type { BatchAnalyticsResponse, AnalysisStatus } from '../types/analytics';

export type AnalyticsHookStatus = 'idle' | 'computing' | 'complete' | 'error';

interface UseBatchAnalyticsReturn {
  data: BatchAnalyticsResponse | null;
  status: AnalyticsHookStatus;
  error: string | null;
  retrigger: (type: string) => void;
}

const POLL_INTERVAL_MS = 2000;

/**
 * Custom hook for polling batch analytics data.
 *
 * @param jobId - The batch job ID to fetch analytics for
 * @param types - Analytics types to trigger (e.g., ['scaffold', 'chemical_space', 'statistics'])
 */
export function useBatchAnalytics(
  jobId: string | null,
  types: string[]
): UseBatchAnalyticsReturn {
  const [data, setData] = useState<BatchAnalyticsResponse | null>(null);
  const [status, setStatus] = useState<AnalyticsHookStatus>('idle');
  const [error, setError] = useState<string | null>(null);
  const intervalRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const triggeredRef = useRef(false);

  // Trigger analytics on mount
  useEffect(() => {
    if (!jobId || types.length === 0) return;

    // Prevent duplicate triggers on strict-mode re-renders
    if (triggeredRef.current) return;
    triggeredRef.current = true;

    setStatus('computing');
    setError(null);

    // Fire-and-forget trigger for each type
    for (const type of types) {
      batchApi.triggerAnalytics(jobId, type).catch(() => {
        // Ignore trigger errors — polling will pick up the status
      });
    }

    // Start polling
    const poll = async () => {
      try {
        const response = await batchApi.getAnalytics(jobId);
        setData(response);

        // Check if all requested types are terminal (complete or failed)
        const allTerminal = types.every((type) => {
          const st: AnalysisStatus | undefined = response.status[type];
          return st && (st.status === 'complete' || st.status === 'failed');
        });

        if (allTerminal) {
          // Check if any failed
          const anyFailed = types.some((type) => {
            const st = response.status[type];
            return st?.status === 'failed';
          });

          if (anyFailed) {
            const failedTypes = types.filter(
              (type) => response.status[type]?.status === 'failed'
            );
            setError(`Analytics failed for: ${failedTypes.join(', ')}`);
          }

          setStatus(anyFailed ? 'error' : 'complete');

          if (intervalRef.current) {
            clearInterval(intervalRef.current);
            intervalRef.current = null;
          }
        }
      } catch (e: unknown) {
        const errMsg = e instanceof Error ? e.message : 'Failed to fetch analytics';
        setError(errMsg);
        setStatus('error');
        if (intervalRef.current) {
          clearInterval(intervalRef.current);
          intervalRef.current = null;
        }
      }
    };

    // Initial fetch
    poll();
    intervalRef.current = setInterval(poll, POLL_INTERVAL_MS);

    return () => {
      if (intervalRef.current) {
        clearInterval(intervalRef.current);
        intervalRef.current = null;
      }
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId]);

  // Reset triggered ref when jobId changes
  useEffect(() => {
    triggeredRef.current = false;
  }, [jobId]);

  const retrigger = useCallback(
    (type: string) => {
      if (!jobId) return;

      setStatus('computing');
      setError(null);

      batchApi.triggerAnalytics(jobId, type).catch(() => {
        // Ignore — polling picks up status
      });

      // Restart polling if it stopped
      if (!intervalRef.current) {
        const poll = async () => {
          try {
            const response = await batchApi.getAnalytics(jobId);
            setData(response);

            const allTerminal = types.every((t) => {
              const st = response.status[t];
              return st && (st.status === 'complete' || st.status === 'failed');
            });

            if (allTerminal) {
              setStatus('complete');
              if (intervalRef.current) {
                clearInterval(intervalRef.current);
                intervalRef.current = null;
              }
            }
          } catch {
            // Silently continue polling
          }
        };

        poll();
        intervalRef.current = setInterval(poll, POLL_INTERVAL_MS);
      }
    },
    [jobId, types]
  );

  return { data, status, error, retrigger };
}
