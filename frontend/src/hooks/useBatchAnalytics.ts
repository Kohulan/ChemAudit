/**
 * useBatchAnalytics Hook
 *
 * Single-instance analytics polling with exponential backoff on errors.
 * Triggers expensive analytics, polls until all requested types complete.
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

const BASE_POLL_MS = 3000;
const MAX_POLL_MS = 30000;

export function useBatchAnalytics(
  jobId: string | null,
  types: string[]
): UseBatchAnalyticsReturn {
  const [data, setData] = useState<BatchAnalyticsResponse | null>(null);
  const [status, setStatus] = useState<AnalyticsHookStatus>('idle');
  const [error, setError] = useState<string | null>(null);
  const timeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const triggeredRef = useRef(false);
  const backoffRef = useRef(BASE_POLL_MS);
  const mountedRef = useRef(true);

  useEffect(() => {
    mountedRef.current = true;
    return () => { mountedRef.current = false; };
  }, []);

  useEffect(() => {
    if (!jobId || types.length === 0) return;

    // Prevent duplicate triggers on strict-mode re-renders
    if (triggeredRef.current) return;
    triggeredRef.current = true;

    setStatus('computing');
    setError(null);
    backoffRef.current = BASE_POLL_MS;

    // Fire-and-forget trigger for each expensive type
    for (const type of types) {
      batchApi.triggerAnalytics(jobId, type).catch(() => {});
    }

    const poll = async () => {
      if (!mountedRef.current) return;

      try {
        const response = await batchApi.getAnalytics(jobId);
        if (!mountedRef.current) return;

        setData(response);
        backoffRef.current = BASE_POLL_MS; // Reset backoff on success

        // Check if all requested types are terminal
        const allTerminal = types.every((type) => {
          const st: AnalysisStatus | undefined = response.status[type];
          return st && (st.status === 'complete' || st.status === 'failed' || st.status === 'skipped');
        });

        if (allTerminal) {
          const anyFailed = types.some((type) => response.status[type]?.status === 'failed');
          if (anyFailed) {
            const failedTypes = types.filter((type) => response.status[type]?.status === 'failed');
            setError(`Analytics failed for: ${failedTypes.join(', ')}`);
          }
          setStatus(anyFailed ? 'error' : 'complete');
          return; // Stop polling
        }

        // Schedule next poll
        timeoutRef.current = setTimeout(poll, backoffRef.current);
      } catch (e: unknown) {
        if (!mountedRef.current) return;

        const axiosResponse = (e as { response?: { status?: number } })?.response;
        const httpStatus = axiosResponse?.status;

        // 404: analytics not initialized yet — keep polling
        if (httpStatus === 404) {
          timeoutRef.current = setTimeout(poll, backoffRef.current);
          return;
        }

        // 429: rate limited — exponential backoff, keep polling
        if (httpStatus === 429) {
          backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);
          timeoutRef.current = setTimeout(poll, backoffRef.current);
          return;
        }

        // Other errors — exponential backoff up to 3 retries then give up
        backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);
        if (backoffRef.current < MAX_POLL_MS) {
          timeoutRef.current = setTimeout(poll, backoffRef.current);
        } else {
          const errMsg = e instanceof Error ? e.message : 'Failed to fetch analytics';
          setError(errMsg);
          setStatus('error');
        }
      }
    };

    // Start polling
    poll();

    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
        timeoutRef.current = null;
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
      backoffRef.current = BASE_POLL_MS;

      batchApi.triggerAnalytics(jobId, type).catch(() => {});

      // Restart polling if stopped
      if (!timeoutRef.current) {
        const poll = async () => {
          if (!mountedRef.current) return;
          try {
            const response = await batchApi.getAnalytics(jobId);
            if (!mountedRef.current) return;
            setData(response);

            const allTerminal = types.every((t) => {
              const st = response.status[t];
              return st && (st.status === 'complete' || st.status === 'failed' || st.status === 'skipped');
            });

            if (allTerminal) {
              setStatus('complete');
              return;
            }
            timeoutRef.current = setTimeout(poll, backoffRef.current);
          } catch {
            backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);
            timeoutRef.current = setTimeout(poll, backoffRef.current);
          }
        };
        poll();
      }
    },
    [jobId, types]
  );

  return { data, status, error, retrigger };
}
