import { useState, useEffect, useCallback, useRef } from 'react';
import type { BatchProgress } from '../types/batch';

/**
 * WebSocket hook for real-time batch progress updates.
 *
 * Connects to the batch progress WebSocket endpoint and provides
 * reactive progress state updates.
 *
 * @param jobId - Job ID to subscribe to (null to disconnect)
 * @returns Progress state, connection status, and manual close function
 */
export function useBatchProgress(jobId: string | null) {
  const [progress, setProgress] = useState<BatchProgress | null>(null);
  const [isConnected, setIsConnected] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const wsRef = useRef<WebSocket | null>(null);
  const reconnectAttempts = useRef(0);
  const maxReconnectAttempts = 3;

  const connect = useCallback(() => {
    if (!jobId) return;

    // Close existing connection
    if (wsRef.current) {
      wsRef.current.close();
    }

    // Determine WebSocket URL based on current location
    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
    const host = window.location.hostname;
    // Use port 8000 for development (backend), or same port for production
    const port = import.meta.env.DEV ? '8000' : window.location.port;
    const wsUrl = `${protocol}//${host}:${port}/ws/batch/${jobId}`;

    try {
      const ws = new WebSocket(wsUrl);
      wsRef.current = ws;

      ws.onopen = () => {
        setIsConnected(true);
        setError(null);
        reconnectAttempts.current = 0;
      };

      ws.onclose = (event) => {
        setIsConnected(false);

        // If job is complete/failed/cancelled, don't reconnect
        if (
          progress?.status === 'complete' ||
          progress?.status === 'failed' ||
          progress?.status === 'cancelled'
        ) {
          return;
        }

        // Attempt reconnection for unexpected closes
        if (
          !event.wasClean &&
          reconnectAttempts.current < maxReconnectAttempts
        ) {
          reconnectAttempts.current += 1;
          setTimeout(() => connect(), 1000 * reconnectAttempts.current);
        }
      };

      ws.onerror = () => {
        setError('WebSocket connection error');
      };

      ws.onmessage = (event) => {
        try {
          const data: BatchProgress = JSON.parse(event.data);
          setProgress(data);

          // If job completed, we can close the connection
          if (
            data.status === 'complete' ||
            data.status === 'failed' ||
            data.status === 'cancelled'
          ) {
            // Keep connection open briefly to ensure we got the final message
            setTimeout(() => {
              if (wsRef.current) {
                wsRef.current.close();
              }
            }, 500);
          }
        } catch (e) {
          console.error('Failed to parse WebSocket message:', e);
        }
      };
    } catch (e) {
      setError('Failed to create WebSocket connection');
    }
  }, [jobId, progress?.status]);

  // Connect/disconnect when jobId changes
  useEffect(() => {
    if (jobId) {
      connect();
    } else {
      // Clear state when jobId is null
      setProgress(null);
      setIsConnected(false);
      setError(null);
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }
    }

    // Cleanup on unmount
    return () => {
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
      }
    };
  }, [jobId, connect]);

  // Manual close function
  const close = useCallback(() => {
    if (wsRef.current) {
      wsRef.current.close();
      wsRef.current = null;
    }
    setIsConnected(false);
  }, []);

  return {
    progress,
    isConnected,
    error,
    close,
  };
}
