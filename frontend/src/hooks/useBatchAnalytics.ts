/**
 * useBatchAnalytics Hook
 *
 * Single-instance analytics polling with exponential backoff on errors.
 * Triggers expensive analytics, polls until all requested types complete.
 * Includes timeout detection and per-type progress tracking.
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { batchApi } from '../services/api';
import type { BatchAnalyticsResponse, AnalysisStatus } from '../types/analytics';

export type AnalyticsHookStatus = 'idle' | 'computing' | 'complete' | 'error';

export interface AnalyticsProgressInfo {
  /** How many of the requested types have reached a terminal state */
  completedCount: number;
  /** Total requested types */
  totalCount: number;
  /** Per-type status from the last successful poll */
  typeStatuses: Record<string, AnalysisStatus>;
  /** Seconds since polling started */
  elapsedSeconds: number;
}

interface UseBatchAnalyticsReturn {
  data: BatchAnalyticsResponse | null;
  status: AnalyticsHookStatus;
  error: string | null;
  progress: AnalyticsProgressInfo;
  retrigger: (type: string) => void;
}

const BASE_POLL_MS = 3000;
const MAX_POLL_MS = 30000;
/** After this many seconds stuck in "computing", treat as timed out */
const COMPUTING_TIMEOUT_S = 300; // 5 minutes

export function useBatchAnalytics(
  jobId: string | null,
  types: string[]
): UseBatchAnalyticsReturn {
  const [data, setData] = useState<BatchAnalyticsResponse | null>(null);
  const [status, setStatus] = useState<AnalyticsHookStatus>('idle');
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState<AnalyticsProgressInfo>({
    completedCount: 0,
    totalCount: types.length,
    typeStatuses: {},
    elapsedSeconds: 0,
  });
  const timeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const triggeredRef = useRef(false);
  const backoffRef = useRef(BASE_POLL_MS);
  const mountedRef = useRef(true);
  const startTimeRef = useRef<number>(0);
  const consecutiveErrorsRef = useRef(0);

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
    startTimeRef.current = Date.now();
    consecutiveErrorsRef.current = 0;
    setProgress({
      completedCount: 0,
      totalCount: types.length,
      typeStatuses: {},
      elapsedSeconds: 0,
    });

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
        consecutiveErrorsRef.current = 0;

        // Update progress info
        const elapsedSeconds = Math.round((Date.now() - startTimeRef.current) / 1000);
        const terminalStatuses = ['complete', 'failed', 'skipped'];
        let completedCount = 0;
        const typeStatuses: Record<string, AnalysisStatus> = {};

        for (const type of types) {
          const st: AnalysisStatus | undefined = response.status[type];
          typeStatuses[type] = st || { status: 'pending', computed_at: null, error: null };
          if (st && terminalStatuses.includes(st.status)) {
            completedCount++;
          }
        }

        setProgress({
          completedCount,
          totalCount: types.length,
          typeStatuses,
          elapsedSeconds,
        });

        // Check if all requested types are terminal
        const allTerminal = completedCount === types.length;

        if (allTerminal) {
          const anyFailed = types.some((type) => response.status[type]?.status === 'failed');
          if (anyFailed) {
            const failedTypes = types.filter((type) => response.status[type]?.status === 'failed');
            setError(`Analytics failed for: ${failedTypes.join(', ')}`);
          }
          setStatus(anyFailed ? 'error' : 'complete');
          return; // Stop polling
        }

        // Check for timeout: if stuck in "computing" too long
        if (elapsedSeconds > COMPUTING_TIMEOUT_S) {
          const stuckTypes = types.filter((type) => {
            const st = response.status[type];
            return st && (st.status === 'computing' || st.status === 'pending');
          });
          if (stuckTypes.length > 0) {
            setError(`Analytics timed out after ${Math.round(COMPUTING_TIMEOUT_S / 60)} minutes for: ${stuckTypes.join(', ')}. Click retry to try again.`);
            setStatus('error');
            return; // Stop polling
          }
        }

        // Schedule next poll
        timeoutRef.current = setTimeout(poll, backoffRef.current);
      } catch (e: unknown) {
        if (!mountedRef.current) return;

        const axiosResponse = (e as { response?: { status?: number } })?.response;
        const httpStatus = axiosResponse?.status;

        // 404: analytics not initialized yet — keep polling
        if (httpStatus === 404) {
          // Check timeout even for 404
          const elapsedSeconds = Math.round((Date.now() - startTimeRef.current) / 1000);
          setProgress(prev => ({ ...prev, elapsedSeconds }));

          if (elapsedSeconds > COMPUTING_TIMEOUT_S) {
            setError('Analytics initialization timed out. The batch may still be processing.');
            setStatus('error');
            return;
          }

          timeoutRef.current = setTimeout(poll, backoffRef.current);
          return;
        }

        // 429: rate limited — exponential backoff, keep polling
        if (httpStatus === 429) {
          backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);
          timeoutRef.current = setTimeout(poll, backoffRef.current);
          return;
        }

        // Other errors — exponential backoff, max 5 consecutive errors then give up
        consecutiveErrorsRef.current++;
        backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);

        if (consecutiveErrorsRef.current >= 5) {
          const errMsg = e instanceof Error ? e.message : 'Failed to fetch analytics';
          setError(errMsg);
          setStatus('error');
        } else {
          timeoutRef.current = setTimeout(poll, backoffRef.current);
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
      startTimeRef.current = Date.now();
      consecutiveErrorsRef.current = 0;

      batchApi.triggerAnalytics(jobId, type).catch(() => {});

      // Restart polling if stopped
      if (!timeoutRef.current) {
        const poll = async () => {
          if (!mountedRef.current) return;
          try {
            const response = await batchApi.getAnalytics(jobId);
            if (!mountedRef.current) return;
            setData(response);

            const elapsedSeconds = Math.round((Date.now() - startTimeRef.current) / 1000);
            const terminalStatuses = ['complete', 'failed', 'skipped'];
            let completedCount = 0;
            const typeStatuses: Record<string, AnalysisStatus> = {};

            for (const t of types) {
              const st = response.status[t];
              typeStatuses[t] = st || { status: 'pending', computed_at: null, error: null };
              if (st && terminalStatuses.includes(st.status)) {
                completedCount++;
              }
            }

            setProgress({
              completedCount,
              totalCount: types.length,
              typeStatuses,
              elapsedSeconds,
            });

            const allTerminal = completedCount === types.length;
            if (allTerminal) {
              const anyFailed = types.some((t) => response.status[t]?.status === 'failed');
              setStatus(anyFailed ? 'error' : 'complete');
              return;
            }

            if (elapsedSeconds > COMPUTING_TIMEOUT_S) {
              setError(`Analytics timed out. Click retry to try again.`);
              setStatus('error');
              return;
            }

            timeoutRef.current = setTimeout(poll, backoffRef.current);
          } catch {
            consecutiveErrorsRef.current++;
            if (consecutiveErrorsRef.current >= 5) {
              setError('Failed to fetch analytics after multiple retries');
              setStatus('error');
            } else {
              backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);
              timeoutRef.current = setTimeout(poll, backoffRef.current);
            }
          }
        };
        poll();
      }
    },
    [jobId, types]
  );

  return { data, status, error, progress, retrigger };
}
