/**
 * useBatchAnalytics Hook
 *
 * Triggers expensive analytics and polls until all active analytics reach a
 * terminal state. "Active" means the type was either explicitly triggered
 * (triggerTypes) or has moved beyond "pending" in the backend response
 * (auto-started by cheap analytics like deduplication and statistics).
 *
 * This decouples "what to trigger" from "what to wait for", ensuring the hook
 * keeps polling until every started analysis completes — not just the
 * user-triggered ones.
 */

import { useState, useEffect, useRef, useCallback } from 'react';
import { batchApi } from '../services/api';
import type { BatchAnalyticsResponse, AnalysisStatus } from '../types/analytics';

export type AnalyticsHookStatus = 'idle' | 'computing' | 'complete' | 'error';

export interface AnalyticsProgressInfo {
  /** How many active types have reached a terminal state */
  completedCount: number;
  /** Total active types (triggered + auto-started) */
  totalCount: number;
  /** Per-type status for all active types from the last successful poll */
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

const TERMINAL_STATUSES = ['complete', 'failed', 'skipped'];
const BASE_POLL_MS = 3000;
const MAX_POLL_MS = 30000;
/** After this many seconds stuck in "computing", treat as timed out */
const COMPUTING_TIMEOUT_S = 300; // 5 minutes

/**
 * Expand the set of active types by adding any types from the response that
 * have moved beyond "pending" (auto-started by the backend). The set is
 * monotonically growing — once a type is tracked, it stays tracked.
 */
function expandActiveTypes(
  responseStatus: Record<string, AnalysisStatus>,
  accumulated: Set<string>,
): void {
  for (const [type, st] of Object.entries(responseStatus)) {
    if (st.status !== 'pending') {
      accumulated.add(type);
    }
  }
}

/**
 * Count how many active types have reached a terminal state.
 */
function countTerminal(
  activeTypes: Set<string>,
  responseStatus: Record<string, AnalysisStatus>,
): { completedCount: number; totalCount: number } {
  let completedCount = 0;
  for (const type of activeTypes) {
    const st = responseStatus[type];
    if (st && TERMINAL_STATUSES.includes(st.status)) {
      completedCount++;
    }
  }
  return { completedCount, totalCount: activeTypes.size };
}

/**
 * Build typeStatuses record for all active types.
 */
function buildTypeStatuses(
  activeTypes: Set<string>,
  responseStatus: Record<string, AnalysisStatus>,
): Record<string, AnalysisStatus> {
  const result: Record<string, AnalysisStatus> = {};
  for (const type of activeTypes) {
    result[type] = responseStatus[type] || { status: 'pending', computed_at: null, error: null };
  }
  return result;
}

/**
 * Find active types that are still stuck in a non-terminal state.
 */
function findStuckTypes(
  activeTypes: Set<string>,
  responseStatus: Record<string, AnalysisStatus>,
): string[] {
  return Array.from(activeTypes).filter((type) => {
    const st = responseStatus[type];
    return st && !TERMINAL_STATUSES.includes(st.status);
  });
}

/**
 * Process a successful poll response: update progress, determine if polling
 * should continue, and return the next hook status (or null to keep polling).
 */
function processPollResponse(
  response: BatchAnalyticsResponse,
  triggerTypes: string[],
  activeTypesRef: React.MutableRefObject<Set<string>>,
  startTime: number,
): {
  completedCount: number;
  totalCount: number;
  typeStatuses: Record<string, AnalysisStatus>;
  elapsedSeconds: number;
  nextStatus: AnalyticsHookStatus | null;
  errorMessage: string | null;
} {
  expandActiveTypes(response.status, activeTypesRef.current);
  // Always include trigger types even if still pending in the response
  for (const t of triggerTypes) activeTypesRef.current.add(t);

  const activeTypes = activeTypesRef.current;
  const { completedCount, totalCount } = countTerminal(activeTypes, response.status);
  const typeStatuses = buildTypeStatuses(activeTypes, response.status);
  const elapsedSeconds = Math.round((Date.now() - startTime) / 1000);

  // All active types terminal?
  if (completedCount === totalCount && totalCount > 0) {
    const failedTypes = Array.from(activeTypes).filter(
      (type) => response.status[type]?.status === 'failed'
    );
    if (failedTypes.length > 0) {
      return {
        completedCount, totalCount, typeStatuses, elapsedSeconds,
        nextStatus: 'error',
        errorMessage: `Analytics failed for: ${failedTypes.join(', ')}`,
      };
    }
    return {
      completedCount, totalCount, typeStatuses, elapsedSeconds,
      nextStatus: 'complete',
      errorMessage: null,
    };
  }

  // Timeout?
  if (elapsedSeconds > COMPUTING_TIMEOUT_S) {
    const stuck = findStuckTypes(activeTypes, response.status);
    if (stuck.length > 0) {
      return {
        completedCount, totalCount, typeStatuses, elapsedSeconds,
        nextStatus: 'error',
        errorMessage: `Analytics timed out after ${Math.round(COMPUTING_TIMEOUT_S / 60)} minutes for: ${stuck.join(', ')}. Click retry to try again.`,
      };
    }
  }

  // Keep polling
  return {
    completedCount, totalCount, typeStatuses, elapsedSeconds,
    nextStatus: null,
    errorMessage: null,
  };
}

export function useBatchAnalytics(
  jobId: string | null,
  triggerTypes: string[]
): UseBatchAnalyticsReturn {
  const [data, setData] = useState<BatchAnalyticsResponse | null>(null);
  const [status, setStatus] = useState<AnalyticsHookStatus>('idle');
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState<AnalyticsProgressInfo>({
    completedCount: 0,
    totalCount: triggerTypes.length,
    typeStatuses: {},
    elapsedSeconds: 0,
  });
  const timeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const triggeredRef = useRef(false);
  const backoffRef = useRef(BASE_POLL_MS);
  const mountedRef = useRef(true);
  const startTimeRef = useRef<number>(0);
  const consecutiveErrorsRef = useRef(0);
  // Accumulates all types we're tracking: trigger types + auto-started types.
  // Monotonically growing — once tracked, stays tracked.
  const activeTypesRef = useRef<Set<string>>(new Set(triggerTypes));

  useEffect(() => {
    mountedRef.current = true;
    return () => { mountedRef.current = false; };
  }, []);

  useEffect(() => {
    if (!jobId || triggerTypes.length === 0) return;

    // Prevent duplicate triggers on strict-mode re-renders
    if (triggeredRef.current) return;
    triggeredRef.current = true;

    setStatus('computing');
    setError(null);
    backoffRef.current = BASE_POLL_MS;
    startTimeRef.current = Date.now();
    consecutiveErrorsRef.current = 0;
    activeTypesRef.current = new Set(triggerTypes);
    setProgress({
      completedCount: 0,
      totalCount: triggerTypes.length,
      typeStatuses: {},
      elapsedSeconds: 0,
    });

    // Fire-and-forget trigger for each expensive type
    for (const type of triggerTypes) {
      batchApi.triggerAnalytics(jobId, type).catch(() => {});
    }

    const poll = async () => {
      if (!mountedRef.current) return;

      try {
        const response = await batchApi.getAnalytics(jobId);
        if (!mountedRef.current) return;

        setData(response);
        backoffRef.current = BASE_POLL_MS;
        consecutiveErrorsRef.current = 0;

        const result = processPollResponse(
          response, triggerTypes, activeTypesRef, startTimeRef.current
        );

        setProgress({
          completedCount: result.completedCount,
          totalCount: result.totalCount,
          typeStatuses: result.typeStatuses,
          elapsedSeconds: result.elapsedSeconds,
        });

        if (result.nextStatus) {
          if (result.errorMessage) setError(result.errorMessage);
          setStatus(result.nextStatus);
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

            const result = processPollResponse(
              response, triggerTypes, activeTypesRef, startTimeRef.current
            );

            setProgress({
              completedCount: result.completedCount,
              totalCount: result.totalCount,
              typeStatuses: result.typeStatuses,
              elapsedSeconds: result.elapsedSeconds,
            });

            if (result.nextStatus) {
              if (result.errorMessage) setError(result.errorMessage);
              setStatus(result.nextStatus);
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
    [jobId, triggerTypes]
  );

  return { data, status, error, progress, retrigger };
}
