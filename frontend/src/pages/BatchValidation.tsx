import { useState, useCallback } from 'react';
import { BatchUpload } from '../components/batch/BatchUpload';
import { BatchProgress } from '../components/batch/BatchProgress';
import { BatchSummary } from '../components/batch/BatchSummary';
import { BatchResultsTable } from '../components/batch/BatchResultsTable';
import { useBatchProgress } from '../hooks/useBatchProgress';
import { batchApi } from '../services/api';
import type {
  BatchPageState,
  BatchResultsResponse,
  BatchResultsFilters,
} from '../types/batch';

/**
 * Batch validation page with state machine:
 * upload -> processing -> results
 */
export function BatchValidationPage() {
  // Page state machine
  const [pageState, setPageState] = useState<BatchPageState>('upload');
  const [jobId, setJobId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  // Results state
  const [resultsData, setResultsData] = useState<BatchResultsResponse | null>(null);
  const [resultsLoading, setResultsLoading] = useState(false);
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(50);
  const [filters, setFilters] = useState<BatchResultsFilters>({});

  // WebSocket progress
  const { progress, isConnected } = useBatchProgress(
    pageState === 'processing' ? jobId : null
  );

  // Handle successful upload
  const handleUploadSuccess = useCallback((newJobId: string, _molecules: number) => {
    setJobId(newJobId);
    setError(null);
    setPageState('processing');
  }, []);

  // Handle upload error
  const handleUploadError = useCallback((errorMsg: string) => {
    setError(errorMsg);
  }, []);

  // Handle processing complete
  const handleProcessingComplete = useCallback(async () => {
    if (!jobId) return;

    setResultsLoading(true);
    try {
      const data = await batchApi.getBatchResults(jobId, 1, pageSize, {});
      setResultsData(data);
      setPage(1);
      setFilters({});
      setPageState('results');
    } catch (e: any) {
      setError(e.message || 'Failed to load results');
    } finally {
      setResultsLoading(false);
    }
  }, [jobId, pageSize]);

  // Handle cancel
  const handleCancel = useCallback(async () => {
    if (!jobId) return;

    try {
      await batchApi.cancelBatch(jobId);
    } catch (e) {
      // Ignore cancel errors
    }
  }, [jobId]);

  // Handle page change
  const handlePageChange = useCallback(
    async (newPage: number) => {
      if (!jobId) return;

      setResultsLoading(true);
      try {
        const data = await batchApi.getBatchResults(jobId, newPage, pageSize, filters);
        setResultsData(data);
        setPage(newPage);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, pageSize, filters]
  );

  // Handle page size change
  const handlePageSizeChange = useCallback(
    async (newSize: number) => {
      if (!jobId) return;

      setResultsLoading(true);
      try {
        const data = await batchApi.getBatchResults(jobId, 1, newSize, filters);
        setResultsData(data);
        setPage(1);
        setPageSize(newSize);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, filters]
  );

  // Handle filter change
  const handleFiltersChange = useCallback(
    async (newFilters: BatchResultsFilters) => {
      if (!jobId) return;

      setResultsLoading(true);
      try {
        const data = await batchApi.getBatchResults(jobId, 1, pageSize, newFilters);
        setResultsData(data);
        setPage(1);
        setFilters(newFilters);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, pageSize]
  );

  // Reset to upload state
  const handleStartNew = useCallback(() => {
    setPageState('upload');
    setJobId(null);
    setError(null);
    setResultsData(null);
    setPage(1);
    setFilters({});
  }, []);

  return (
    <div className="max-w-7xl mx-auto space-y-6">
      {/* Header */}
      <div className="text-center">
        <h2 className="text-3xl font-bold text-gray-900">Batch Validation</h2>
        <p className="text-gray-500 mt-2">
          Validate up to 10,000 molecules at once from SDF or CSV files
        </p>
      </div>

      {/* Error display */}
      {error && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-4">
          <div className="flex items-start">
            <svg
              className="h-5 w-5 text-red-400 mt-0.5"
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M12 8v4m0 4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"
              />
            </svg>
            <div className="ml-3">
              <p className="text-sm text-red-700">{error}</p>
            </div>
            <button
              onClick={() => setError(null)}
              className="ml-auto text-red-400 hover:text-red-600"
            >
              <svg className="h-5 w-5" fill="currentColor" viewBox="0 0 20 20">
                <path
                  fillRule="evenodd"
                  d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z"
                  clipRule="evenodd"
                />
              </svg>
            </button>
          </div>
        </div>
      )}

      {/* Upload state */}
      {pageState === 'upload' && (
        <div className="bg-white rounded-lg shadow-md p-6">
          <h3 className="font-semibold text-gray-900 mb-4">Upload File</h3>
          <BatchUpload
            onUploadSuccess={handleUploadSuccess}
            onUploadError={handleUploadError}
          />
        </div>
      )}

      {/* Processing state */}
      {pageState === 'processing' && (
        <BatchProgress
          progress={progress}
          isConnected={isConnected}
          onCancel={handleCancel}
          onComplete={handleProcessingComplete}
        />
      )}

      {/* Results state */}
      {pageState === 'results' && resultsData && (
        <div className="space-y-6">
          {/* Start new batch button */}
          <div className="flex justify-end">
            <button
              onClick={handleStartNew}
              className="px-4 py-2 text-sm bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors"
            >
              Start New Batch
            </button>
          </div>

          {/* Summary */}
          {resultsData.statistics && jobId && (
            <div className="bg-white rounded-lg shadow-md p-6">
              <h3 className="text-lg font-semibold text-gray-900 mb-4">Summary</h3>
              <BatchSummary jobId={jobId} statistics={resultsData.statistics} />
            </div>
          )}

          {/* Results table */}
          <div className="bg-white rounded-lg shadow-md p-6">
            <h3 className="text-lg font-semibold text-gray-900 mb-4">Results</h3>
            <BatchResultsTable
              results={resultsData.results}
              page={page}
              pageSize={pageSize}
              totalResults={resultsData.total_results}
              totalPages={resultsData.total_pages}
              filters={filters}
              onPageChange={handlePageChange}
              onPageSizeChange={handlePageSizeChange}
              onFiltersChange={handleFiltersChange}
              isLoading={resultsLoading}
            />
          </div>
        </div>
      )}
    </div>
  );
}
