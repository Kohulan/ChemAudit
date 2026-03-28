import { useState, useCallback, useRef, useEffect } from 'react';
import { qsarReadyApi } from '../services/api';
import type {
  QSARReadyConfig,
  QSARReadyResult,
  QSARBatchStatusResponse,
  QSARBatchResultsResponse,
} from '../types/qsar_ready';

// =============================================================================
// Return type
// =============================================================================

export interface UseQSARReadyReturn {
  // Single molecule state
  singleResult: QSARReadyResult | null;
  singleLoading: boolean;
  singleError: string | null;

  // Batch state
  batchJobId: string | null;
  batchStatus: QSARBatchStatusResponse | null;
  batchResults: QSARBatchResultsResponse | null;
  batchLoading: boolean;
  batchError: string | null;
  batchPage: number;

  // Actions
  runSingle: (smiles: string, config: QSARReadyConfig) => Promise<void>;
  runBatch: (
    file: File | null,
    smilesText: string | null,
    config: QSARReadyConfig,
  ) => Promise<void>;
  fetchBatchPage: (page: number) => Promise<void>;
  downloadBatch: (format: 'csv' | 'sdf' | 'json') => Promise<void>;
  clearSingle: () => void;
  clearBatch: () => void;
}

// =============================================================================
// Hook
// =============================================================================

/**
 * Hook for QSAR-Ready Pipeline state management.
 *
 * Manages:
 * - Single molecule curation (runSingle)
 * - Batch processing with WebSocket progress tracking (runBatch)
 * - Paginated results fetching (fetchBatchPage)
 * - Result download in multiple formats (downloadBatch)
 *
 * Error handling follows the per-tool isolation pattern from useDiagnostics.ts:
 * each state triplet (result/loading/error) is fully independent.
 */
export function useQSARReady(): UseQSARReadyReturn {
  // ── Single molecule state ──
  const [singleResult, setSingleResult] = useState<QSARReadyResult | null>(null);
  const [singleLoading, setSingleLoading] = useState(false);
  const [singleError, setSingleError] = useState<string | null>(null);

  // ── Batch state ──
  const [batchJobId, setBatchJobId] = useState<string | null>(null);
  const [batchStatus, setBatchStatus] = useState<QSARBatchStatusResponse | null>(null);
  const [batchResults, setBatchResults] = useState<QSARBatchResultsResponse | null>(null);
  const [batchLoading, setBatchLoading] = useState(false);
  const [batchError, setBatchError] = useState<string | null>(null);
  const [batchPage, setBatchPage] = useState(1);

  // WebSocket ref for cleanup
  const wsRef = useRef<WebSocket | null>(null);

  // ── WebSocket cleanup on unmount ──
  useEffect(() => {
    return () => {
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }
    };
  }, []);

  // ==========================================================================
  // Internal: open WebSocket and track batch progress
  // ==========================================================================

  const connectBatchWebSocket = useCallback(
    (jobId: string, currentConfig: QSARReadyConfig) => {
      // Close any existing WS before opening a new one
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }

      const wsBase = window.location.origin.replace(/^http/, 'ws');
      const wsUrl = `${wsBase}/ws/qsar/${jobId}`;
      const ws = new WebSocket(wsUrl);
      wsRef.current = ws;

      ws.onmessage = (event) => {
        try {
          const msg = JSON.parse(event.data) as Partial<QSARBatchStatusResponse>;
          setBatchStatus((prev) => ({
            job_id: jobId,
            status: msg.status ?? prev?.status ?? 'processing',
            progress: msg.progress ?? prev?.progress ?? 0,
            processed: msg.processed ?? prev?.processed ?? 0,
            total: msg.total ?? prev?.total ?? 0,
            eta_seconds: msg.eta_seconds ?? prev?.eta_seconds ?? null,
          }));

          // On completion, fetch the full results
          if (msg.status === 'complete') {
            ws.close();
            wsRef.current = null;
            setBatchLoading(false);
            // Fetch first page of results
            qsarReadyApi
              .batchResults(jobId, 1, 50)
              .then((results) => {
                setBatchResults(results);
                setBatchPage(1);
              })
              .catch((err: unknown) => {
                const e = err as { error?: string; detail?: string };
                setBatchError(e?.error ?? e?.detail ?? 'Failed to fetch batch results');
              });
          } else if (msg.status === 'failed') {
            ws.close();
            wsRef.current = null;
            setBatchLoading(false);
            setBatchError('Batch processing failed');
          }
        } catch {
          // Ignore malformed WS messages
        }
      };

      ws.onerror = () => {
        // WS failed — fall back to polling via batchStatus endpoint
        ws.close();
        wsRef.current = null;
        // Start polling as fallback
        startPolling(jobId, currentConfig);
      };

      ws.onclose = () => {
        if (wsRef.current === ws) {
          wsRef.current = null;
        }
      };
    },
    [], // eslint-disable-line react-hooks/exhaustive-deps
  );

  // ==========================================================================
  // Internal: polling fallback when WebSocket is unavailable
  // ==========================================================================

  const startPolling = useCallback(
    (jobId: string, _config: QSARReadyConfig) => {
      const poll = async () => {
        try {
          const status = await qsarReadyApi.batchStatus(jobId);
          setBatchStatus(status);

          if (status.status === 'complete') {
            setBatchLoading(false);
            const results = await qsarReadyApi.batchResults(jobId, 1, 50);
            setBatchResults(results);
            setBatchPage(1);
          } else if (status.status === 'failed') {
            setBatchLoading(false);
            setBatchError('Batch processing failed');
          } else {
            // Still processing — poll again in 1.5s
            setTimeout(poll, 1500);
          }
        } catch (err: unknown) {
          const e = err as { error?: string; detail?: string };
          setBatchLoading(false);
          setBatchError(e?.error ?? e?.detail ?? 'Status polling failed');
        }
      };
      setTimeout(poll, 1500);
    },
    [],
  );

  // ==========================================================================
  // runSingle
  // ==========================================================================

  const runSingle = useCallback(async (smiles: string, config: QSARReadyConfig) => {
    setSingleLoading(true);
    setSingleError(null);
    setSingleResult(null);
    try {
      const result = await qsarReadyApi.single(smiles, config);
      setSingleResult(result);
    } catch (err: unknown) {
      const e = err as { error?: string; detail?: string };
      setSingleError(e?.error ?? e?.detail ?? 'QSAR-ready curation failed');
    } finally {
      setSingleLoading(false);
    }
  }, []);

  // ==========================================================================
  // runBatch
  // ==========================================================================

  const runBatch = useCallback(
    async (file: File | null, smilesText: string | null, config: QSARReadyConfig) => {
      setBatchLoading(true);
      setBatchError(null);
      setBatchStatus(null);
      setBatchResults(null);
      setBatchJobId(null);
      setBatchPage(1);

      try {
        const uploadResp = await qsarReadyApi.batchUpload(file, smilesText, config);
        const jobId = uploadResp.job_id;
        setBatchJobId(jobId);

        // Set initial status
        setBatchStatus({
          job_id: jobId,
          status: uploadResp.status as QSARBatchStatusResponse['status'],
          progress: 0,
          processed: 0,
          total: uploadResp.total_molecules,
          eta_seconds: null,
        });

        // Connect WebSocket for progress tracking
        connectBatchWebSocket(jobId, config);
      } catch (err: unknown) {
        const e = err as { error?: string; detail?: string };
        setBatchLoading(false);
        setBatchError(e?.error ?? e?.detail ?? 'Batch upload failed');
      }
    },
    [connectBatchWebSocket],
  );

  // ==========================================================================
  // fetchBatchPage
  // ==========================================================================

  const fetchBatchPage = useCallback(async (page: number) => {
    if (!batchJobId) return;
    setBatchLoading(true);
    setBatchError(null);
    try {
      const results = await qsarReadyApi.batchResults(batchJobId, page, 50);
      setBatchResults(results);
      setBatchPage(page);
    } catch (err: unknown) {
      const e = err as { error?: string; detail?: string };
      setBatchError(e?.error ?? e?.detail ?? 'Failed to fetch page');
    } finally {
      setBatchLoading(false);
    }
  }, [batchJobId]);

  // ==========================================================================
  // downloadBatch
  // ==========================================================================

  const downloadBatch = useCallback(async (format: 'csv' | 'sdf' | 'json') => {
    if (!batchJobId) return;
    try {
      const blob = await qsarReadyApi.batchDownload(batchJobId, format);
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `qsar_ready_${batchJobId}.${format}`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    } catch (err: unknown) {
      const e = err as { error?: string; detail?: string };
      setBatchError(e?.error ?? e?.detail ?? `Download failed (${format})`);
    }
  }, [batchJobId]);

  // ==========================================================================
  // clearSingle / clearBatch
  // ==========================================================================

  const clearSingle = useCallback(() => {
    setSingleResult(null);
    setSingleError(null);
  }, []);

  const clearBatch = useCallback(() => {
    if (wsRef.current) {
      wsRef.current.close();
      wsRef.current = null;
    }
    setBatchJobId(null);
    setBatchStatus(null);
    setBatchResults(null);
    setBatchError(null);
    setBatchPage(1);
  }, []);

  return {
    singleResult,
    singleLoading,
    singleError,
    batchJobId,
    batchStatus,
    batchResults,
    batchLoading,
    batchError,
    batchPage,
    runSingle,
    runBatch,
    fetchBatchPage,
    downloadBatch,
    clearSingle,
    clearBatch,
  };
}
