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
  retrigger: (type: string, params?: Record<string, string>) => void;
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
 * Merge a fresh analytics response with the previous state, preserving result
 * fields the backend omits for types still "computing". This prevents a
 * retrigger (e.g. t-SNE) from wiping previously loaded data (e.g. PCA).
 */
function mergeAnalyticsData(
  prev: BatchAnalyticsResponse | null,
  next: BatchAnalyticsResponse,
): BatchAnalyticsResponse {
  if (!prev) return next;
  return {
    ...next,
    deduplication: next.deduplication ?? prev.deduplication,
    scaffold: next.scaffold ?? prev.scaffold,
    chemical_space: next.chemical_space ?? prev.chemical_space,
    similarity_matrix: next.similarity_matrix ?? prev.similarity_matrix,
    mmp: next.mmp ?? prev.mmp,
    statistics: next.statistics ?? prev.statistics,
  };
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

    // Trigger each expensive type and wait for backend acknowledgement
    // before polling. Without this, the first GET can arrive before any
    // POST is processed, yielding a 404 (analytics storage not yet init).
    const triggers = triggerTypes.map(type =>
      batchApi.triggerAnalytics(jobId, type).catch(() => {})
    );

    const poll = async () => {
      if (!mountedRef.current) return;

      try {
        const response = await batchApi.getAnalytics(jobId);
        if (!mountedRef.current) return;

        setData(prev => mergeAnalyticsData(prev, response));
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
          timeoutRef.current = null; // Clear so retrigger can restart polling
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
            timeoutRef.current = null;
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
          timeoutRef.current = null;
        } else {
          timeoutRef.current = setTimeout(poll, backoffRef.current);
        }
      }
    };

    // Start polling only after triggers are acknowledged
    Promise.all(triggers).then(() => {
      if (mountedRef.current) poll();
    });

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
    async (type: string, params?: Record<string, string>) => {
      if (!jobId) return;

      setStatus('computing');
      setError(null);
      backoffRef.current = BASE_POLL_MS;
      startTimeRef.current = Date.now();
      consecutiveErrorsRef.current = 0;

      // Await the trigger so the backend status transitions to "computing"
      // before we start polling. Without this, the first poll can see the
      // OLD "complete" status (e.g. PCA) and exit immediately.
      try {
        await batchApi.triggerAnalytics(jobId, type, params);
      } catch {
        // Trigger failed — still start polling in case backend processed it
      }

      // Always start a fresh poll loop for the retriggered type.
      // If a previous poll is running (e.g. initial mount poll tracking
      // multiple types), cancel it — otherwise the retrigger has no poll
      // watching for the new computation to complete.
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
        timeoutRef.current = null;
      }

      const poll = async () => {
        if (!mountedRef.current) return;
        try {
          const response = await batchApi.getAnalytics(jobId);
          if (!mountedRef.current) return;

          setData(prev => mergeAnalyticsData(prev, response));

          const st = response.status[type];
          const elapsedSeconds = Math.round((Date.now() - startTimeRef.current) / 1000);
          const isTerminal = st && TERMINAL_STATUSES.includes(st.status);

          setProgress({
            completedCount: isTerminal ? 1 : 0,
            totalCount: 1,
            typeStatuses: { [type]: st || { status: 'pending', computed_at: null, error: null } },
            elapsedSeconds,
          });

          if (isTerminal) {
            if (st.status === 'failed') {
              setError(`Analytics failed for: ${type}`);
              setStatus('error');
            } else {
              setStatus('complete');
            }
            timeoutRef.current = null;
            return;
          }

          if (elapsedSeconds > COMPUTING_TIMEOUT_S) {
            setError(`Analytics timed out for: ${type}. Click retry to try again.`);
            setStatus('error');
            timeoutRef.current = null;
            return;
          }

          timeoutRef.current = setTimeout(poll, backoffRef.current);
        } catch {
          consecutiveErrorsRef.current++;
          if (consecutiveErrorsRef.current >= 5) {
            setError('Failed to fetch analytics after multiple retries');
            setStatus('error');
            timeoutRef.current = null;
          } else {
            backoffRef.current = Math.min(backoffRef.current * 2, MAX_POLL_MS);
            timeoutRef.current = setTimeout(poll, backoffRef.current);
          }
        }
      };
      poll();
    },
    [jobId]
  );

  return { data, status, error, progress, retrigger };
}
