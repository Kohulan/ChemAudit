import { useState, useCallback, useRef, useEffect } from 'react';
import type {
  DatasetAuditStatus,
  DatasetAuditResults,
  DatasetDiffResults,
} from '../types/dataset_intelligence';
import { datasetApi } from '../services/api';

// =============================================================================
// Types
// =============================================================================

export interface UseDatasetAuditReturn {
  /** Currently selected file (null before first upload). */
  file: File | null;
  /** Job ID returned from the upload endpoint. */
  jobId: string | null;
  /** Current lifecycle status. */
  status: DatasetAuditStatus;
  /** Processing progress 0-100. */
  progress: number;
  /** Current processing stage label. */
  currentStage: string | null;
  /** Full audit results after completion. */
  results: DatasetAuditResults | null;
  /** Diff results from comparison file upload. */
  diffResults: DatasetDiffResults | null;
  /** Error message (upload or processing failure). */
  error: string | null;
  /** Error message for diff upload. */
  diffError: string | null;
  /** Whether a diff upload is in progress. */
  diffLoading: boolean;
  /** Upload a CSV/SDF file and start the audit pipeline. */
  uploadFile: (file: File) => Promise<void>;
  /** Upload a comparison file for dataset diff. */
  uploadDiffFile: (file: File) => Promise<void>;
  /** Reset all state to initial values. */
  resetAll: () => void;
}

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook for Dataset Audit state management.
 *
 * Manages the full lifecycle:
 * 1. Upload file via datasetApi.upload
 * 2. Connect WebSocket at ws/dataset/{job_id} for progress
 * 3. On complete, fetch results via datasetApi.getResults
 * 4. Optionally upload a comparison file for dataset diff
 *
 * WebSocket follows the established pattern from useStructureFilter / useQSARReady
 * (Phase 10/11), with polling fallback on ws.onerror (2s interval).
 */
export function useDatasetAudit(): UseDatasetAuditReturn {
  const [file, setFile] = useState<File | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);
  const [status, setStatus] = useState<DatasetAuditStatus>('idle');
  const [progress, setProgress] = useState<number>(0);
  const [currentStage, setCurrentStage] = useState<string | null>(null);
  const [results, setResults] = useState<DatasetAuditResults | null>(null);
  const [diffResults, setDiffResults] = useState<DatasetDiffResults | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [diffError, setDiffError] = useState<string | null>(null);
  const [diffLoading, setDiffLoading] = useState(false);
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

  // WebSocket connection for processing progress
  useEffect(() => {
    if (!jobId || status !== 'processing') return;

    // Close any existing connection before opening a new one
    if (wsRef.current) {
      wsRef.current.close();
      wsRef.current = null;
    }

    const wsBase = window.location.origin.replace(/^http/, 'ws');
    const wsUrl = `${wsBase}/ws/dataset/${jobId}`;
    const ws = new WebSocket(wsUrl);
    wsRef.current = ws;

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data) as {
          progress?: number;
          current_stage?: string;
          status?: string;
          eta_seconds?: number;
        };
        if (data.progress !== undefined) setProgress(data.progress);
        if (data.current_stage) setCurrentStage(data.current_stage);
        if (data.status === 'complete') {
          ws.close();
          wsRef.current = null;
          // Fetch full results
          datasetApi
            .getResults(jobId)
            .then((res) => {
              setResults(res);
              setStatus('complete');
            })
            .catch((err: unknown) => {
              const e = err as { message?: string; detail?: string };
              setError(e?.message ?? e?.detail ?? 'Failed to fetch audit results');
              setStatus('error');
            });
        } else if (data.status === 'error' || data.status === 'failed') {
          ws.close();
          wsRef.current = null;
          setError('Dataset audit processing failed');
          setStatus('error');
        }
      } catch {
        // Ignore malformed WebSocket messages
      }
    };

    ws.onerror = () => {
      // WebSocket failed — fall back to polling (Phase 10/11 pattern)
      ws.close();
      wsRef.current = null;

      const poll = async () => {
        try {
          const statusResp = await datasetApi.getStatus(jobId);
          if (statusResp.progress !== undefined) setProgress(statusResp.progress);
          if (statusResp.current_stage) setCurrentStage(statusResp.current_stage);

          if (statusResp.status === 'complete') {
            const res = await datasetApi.getResults(jobId);
            setResults(res);
            setStatus('complete');
          } else if (statusResp.status === 'error' || statusResp.status === 'failed') {
            setError('Dataset audit processing failed');
            setStatus('error');
          } else {
            // Still processing — poll again in 2s
            setTimeout(poll, 2000);
          }
        } catch (err: unknown) {
          const e = err as { message?: string; detail?: string };
          setError(e?.message ?? e?.detail ?? 'Status polling failed');
          setStatus('error');
        }
      };

      setTimeout(poll, 2000);
    };

    ws.onclose = () => {
      if (wsRef.current === ws) {
        wsRef.current = null;
      }
    };
  }, [jobId, status]);

  // ==========================================================================
  // uploadFile
  // ==========================================================================

  const uploadFile = useCallback(async (newFile: File) => {
    setFile(newFile);
    setError(null);
    setResults(null);
    setDiffResults(null);
    setDiffError(null);
    setProgress(0);
    setCurrentStage(null);
    setStatus('uploading');

    try {
      const resp = await datasetApi.upload(newFile);
      setJobId(resp.job_id);
      setStatus('processing');
      // WebSocket useEffect will trigger on jobId + status='processing'
    } catch (err: unknown) {
      const e = err as { response?: { data?: { detail?: string } }; message?: string };
      setError(
        e?.response?.data?.detail ?? e?.message ?? 'Upload failed',
      );
      setStatus('error');
    }
  }, []);

  // ==========================================================================
  // uploadDiffFile
  // ==========================================================================

  const uploadDiffFile = useCallback(
    async (diffFile: File) => {
      if (!jobId) return;

      setDiffLoading(true);
      setDiffError(null);
      setDiffResults(null);

      try {
        const res = await datasetApi.uploadDiff(jobId, diffFile);
        setDiffResults(res);
      } catch (err: unknown) {
        const e = err as { response?: { data?: { detail?: string } }; message?: string };
        setDiffError(
          e?.response?.data?.detail ?? e?.message ?? 'Diff upload failed',
        );
      } finally {
        setDiffLoading(false);
      }
    },
    [jobId],
  );

  // ==========================================================================
  // resetAll
  // ==========================================================================

  const resetAll = useCallback(() => {
    if (wsRef.current) {
      wsRef.current.close();
      wsRef.current = null;
    }
    setFile(null);
    setJobId(null);
    setStatus('idle');
    setProgress(0);
    setCurrentStage(null);
    setResults(null);
    setDiffResults(null);
    setError(null);
    setDiffError(null);
    setDiffLoading(false);
  }, []);

  return {
    file,
    jobId,
    status,
    progress,
    currentStage,
    results,
    diffResults,
    error,
    diffError,
    diffLoading,
    uploadFile,
    uploadDiffFile,
    resetAll,
  };
}
