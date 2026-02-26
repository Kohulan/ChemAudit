import { useState, useCallback, useEffect, useRef } from 'react';
import { useLocation } from 'react-router-dom';
import { useBatchCache } from '../contexts/BatchCacheContext';
import { motion, AnimatePresence } from 'framer-motion';
import { Upload, AlertTriangle, X, ArrowRight, RotateCcw, FileSpreadsheet, Sparkles, Clock, BarChart3, Share2, Check } from 'lucide-react';
import { BatchUpload } from '../components/batch/BatchUpload';
import { BatchProgress } from '../components/batch/BatchProgress';
import { BatchSummary } from '../components/batch/BatchSummary';
import { BatchResultsTable } from '../components/batch/BatchResultsTable';
import { BatchAnalyticsPanel } from '../components/batch/BatchAnalyticsPanel';
import { MoleculeComparisonPanel } from '../components/batch/MoleculeComparisonPanel';
import { SubsetActionPanel } from '../components/batch/SubsetActionPanel';
import { BatchTimeline } from '../components/batch/BatchTimeline';
import { ProfileSidebar } from '../components/batch/ProfileSidebar';
import { ClayButton } from '../components/ui/ClayButton';
import { useBatchProgress } from '../hooks/useBatchProgress';
import { useBatchAnalytics } from '../hooks/useBatchAnalytics';
import { useBrushSelection, setSelection, toggleIndex, clearSelection } from '../hooks/useBrushSelection';
import { useLimits } from '../context/ConfigContext';
import { batchApi, permalinksApi } from '../services/api';
import { cn } from '../lib/utils';
import type {
  BatchPageState,
  BatchResultsResponse,
  BatchResultsFilters,
  SortField,
} from '../types/batch';

/**
 * Premium batch validation page with state machine:
 * upload -> processing -> results
 */
export function BatchValidationPage() {
  const limits = useLimits();
  const location = useLocation();
  const { getCache, setCache, clearCache } = useBatchCache();

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
  const [sortBy, setSortBy] = useState<SortField>('index');
  const [sortDir, setSortDir] = useState<'asc' | 'desc'>('asc');
  const [selectedIndices, selectionDispatch] = useBrushSelection();
  const [compareMode, setCompareMode] = useState(false);
  const [comparisonResults, setComparisonResults] = useState<import('../types/batch').BatchResult[]>([]);
  const [compareLoading, setCompareLoading] = useState(false);
  const [includeAnalytics, setIncludeAnalytics] = useState(true);
  const [selectedProfileId, setSelectedProfileId] = useState<number | null>(null);
  const [subsetPanelOpen, setSubsetPanelOpen] = useState(false);
  const [shareCopied, setShareCopied] = useState(false);

  // Analytics data for timeline and comparison radar
  const { data: analyticsData, status: analyticsStatus, error: analyticsError, progress: analyticsProgress, retrigger: analyticsRetrigger } = useBatchAnalytics(
    pageState === 'results' && includeAnalytics ? jobId : null,
    ['scaffold', 'chemical_space']
  );

  // WebSocket progress
  const { progress, isConnected } = useBatchProgress(
    pageState === 'processing' ? jobId : null
  );

  // Handle successful upload
  const handleUploadSuccess = useCallback((newJobId: string, _molecules: number, options?: { includeAnalytics: boolean }) => {
    setJobId(newJobId);
    setError(null);
    setIncludeAnalytics(options?.includeAnalytics ?? true);
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
        const filtersWithSort = { ...filters, sort_by: sortBy, sort_dir: sortDir };
        const data = await batchApi.getBatchResults(jobId, newPage, pageSize, filtersWithSort);
        setResultsData(data);
        setPage(newPage);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, pageSize, filters, sortBy, sortDir]
  );

  // Handle page size change
  const handlePageSizeChange = useCallback(
    async (newSize: number) => {
      if (!jobId) return;

      setResultsLoading(true);
      try {
        const filtersWithSort = { ...filters, sort_by: sortBy, sort_dir: sortDir };
        const data = await batchApi.getBatchResults(jobId, 1, newSize, filtersWithSort);
        setResultsData(data);
        setPage(1);
        setPageSize(newSize);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, filters, sortBy, sortDir]
  );

  // Handle filter change
  const handleFiltersChange = useCallback(
    async (newFilters: BatchResultsFilters) => {
      if (!jobId) return;

      setResultsLoading(true);
      try {
        const filtersWithSort = { ...newFilters, sort_by: sortBy, sort_dir: sortDir };
        const data = await batchApi.getBatchResults(jobId, 1, pageSize, filtersWithSort);
        setResultsData(data);
        setPage(1);
        setFilters(filtersWithSort);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, pageSize, sortBy, sortDir]
  );

  // Handle sort change
  const handleSortChange = useCallback(
    async (newSortBy: SortField, newSortDir: 'asc' | 'desc') => {
      if (!jobId) return;

      setResultsLoading(true);
      try {
        const newFilters = { ...filters, sort_by: newSortBy, sort_dir: newSortDir };
        const data = await batchApi.getBatchResults(jobId, 1, pageSize, newFilters);
        setResultsData(data);
        setPage(1);
        setSortBy(newSortBy);
        setSortDir(newSortDir);
        setFilters(newFilters);
      } catch (e: any) {
        setError(e.message || 'Failed to load results');
      } finally {
        setResultsLoading(false);
      }
    },
    [jobId, pageSize, filters]
  );

  // Handle selection change
  const handleSelectionChange = useCallback((indices: Set<number>) => {
    selectionDispatch(setSelection(indices));
  }, [selectionDispatch]);

  // Handle clear selection
  const handleClearSelection = useCallback(() => {
    selectionDispatch(clearSelection());
  }, [selectionDispatch]);

  // Handle compare button click — fetch molecules from API (not paginated local data)
  const handleCompare = useCallback(async () => {
    if (!jobId || selectedIndices.size === 0) return;
    const indices = Array.from(selectedIndices).slice(0, 2);

    setCompareLoading(true);
    try {
      const results: import('../types/batch').BatchResult[] = [];
      for (const idx of indices) {
        // Compute the page that contains this molecule index (sort_by=index asc)
        const pageNum = Math.floor(idx / pageSize) + 1;
        const data = await batchApi.getBatchResults(jobId, pageNum, pageSize, {
          sort_by: 'index',
          sort_dir: 'asc',
        });
        const mol = data.results.find((r) => r.index === idx);
        if (mol) results.push(mol);
      }
      setComparisonResults(results);
      setCompareMode(true);
    } catch (e: any) {
      setError(e.message || 'Failed to fetch molecules for comparison');
    } finally {
      setCompareLoading(false);
    }
  }, [jobId, selectedIndices, pageSize]);

  // Handle close comparison panel
  const handleCloseCompare = useCallback(() => {
    setCompareMode(false);
    setComparisonResults([]);
  }, []);

  // Handle removing a molecule from comparison
  const handleRemoveFromCompare = useCallback((index: number) => {
    // index is the position in comparisonResults array (0 or 1)
    const molIndex = comparisonResults[index]?.index;
    if (molIndex !== undefined) {
      selectionDispatch(toggleIndex(molIndex));
      setComparisonResults((prev) => prev.filter((_, i) => i !== index));
    }
    // If removing leaves 0, close the panel
    if (comparisonResults.length <= 1) {
      setCompareMode(false);
    }
  }, [comparisonResults, selectionDispatch]);

  // Handle share permalink
  const handleShare = useCallback(async () => {
    if (!jobId) return;
    try {
      const result = await permalinksApi.createPermalink(jobId);
      const url = `${window.location.origin}/report/${result.short_id}`;
      await navigator.clipboard.writeText(url);
      setShareCopied(true);
      setTimeout(() => setShareCopied(false), 2000);
    } catch (err) {
      console.error('Failed to create permalink:', err);
    }
  }, [jobId]);

  // Reset to upload state
  const handleStartNew = useCallback(() => {
    setPageState('upload');
    setJobId(null);
    setError(null);
    setResultsData(null);
    setPage(1);
    setFilters({});
    setSortBy('index');
    setSortDir('asc');
    setCompareMode(false);
    setComparisonResults([]);
    selectionDispatch(clearSelection());
    clearCache();
  }, [selectionDispatch, clearCache]);

  // Restore cached batch state on mount (when navigating back to this page)
  const didRestore = useRef(false);
  useEffect(() => {
    if (didRestore.current) return;
    didRestore.current = true;
    const cached = getCache();
    if (cached && cached.jobId) {
      setPageState(cached.pageState);
      setJobId(cached.jobId);
      setResultsData(cached.resultsData);
      setPage(cached.page);
      setPageSize(cached.pageSize);
      setFilters(cached.filters);
      setSortBy(cached.sortBy);
      setSortDir(cached.sortDir);
      setIncludeAnalytics(cached.includeAnalytics);
      setSelectedProfileId(cached.selectedProfileId ?? null);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Restore from permalink navigation (location.state.permalinkJobId)
  const didRestorePermalink = useRef(false);
  useEffect(() => {
    if (didRestorePermalink.current) return;
    const state = location.state as { permalinkJobId?: string; permalinkSnapshot?: Record<string, unknown> } | null;
    if (!state?.permalinkJobId) return;
    didRestorePermalink.current = true;
    // Clear the location state so refresh doesn't re-trigger
    window.history.replaceState({}, '', '/batch');

    const permalinkJobId = state.permalinkJobId;
    setJobId(permalinkJobId);
    setPageState('results');

    // Fetch results for the permalink job
    setResultsLoading(true);
    batchApi.getBatchResults(permalinkJobId, 1, pageSize, {}).then((data) => {
      setResultsData(data);
      setPage(1);
      setFilters({});
    }).catch((e: any) => {
      setError(e.message || 'Failed to load permalink results');
      setPageState('upload');
    }).finally(() => {
      setResultsLoading(false);
    });
  }, [location.state, pageSize]);

  // Persist batch state to cache whenever it changes
  useEffect(() => {
    if (!jobId) return;
    setCache({
      pageState,
      jobId,
      resultsData,
      page,
      pageSize,
      filters,
      sortBy,
      sortDir,
      includeAnalytics,
      selectedProfileId,
    });
  }, [pageState, jobId, resultsData, page, pageSize, filters, sortBy, sortDir, includeAnalytics, selectedProfileId, setCache]);

  return (
    <div className={cn(
      "mx-auto space-y-8 px-4 sm:px-6",
      pageState === 'results' ? 'max-w-[1600px]' : 'max-w-7xl'
    )}>
      {/* Hero Header */}
      <motion.div
        className="text-center pt-4"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: [0.25, 0.46, 0.45, 0.94] }}
      >
        <div className="inline-flex items-center gap-2 px-4 py-1.5 rounded-full bg-[var(--color-primary)]/10 border border-[var(--color-primary)]/20 mb-4">
          <FileSpreadsheet className="w-4 h-4 text-[var(--color-primary)]" />
          <span className="text-sm font-medium text-[var(--color-primary)]">High-throughput Processing</span>
        </div>
        <h1 className="text-3xl sm:text-4xl font-bold text-gradient tracking-tight font-display">
          Batch Validation
        </h1>
        <p className="text-[var(--color-text-secondary)] mt-3 text-base sm:text-lg max-w-2xl mx-auto">
          Validate up to <span className="font-semibold text-[var(--color-text-primary)]">{limits.max_batch_size.toLocaleString()} molecules</span> at once from SDF or CSV files
        </p>

        {/* Large batch warning for xl/coconut profiles */}
        {limits.max_batch_size >= 50000 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.3 }}
            className="inline-flex items-center gap-2 mt-4 px-4 py-2 rounded-lg bg-amber-500/10 border border-amber-500/20 text-amber-600 dark:text-amber-400"
          >
            <Clock className="w-4 h-4" />
            <span className="text-sm">Large batches may take several minutes to process</span>
          </motion.div>
        )}
      </motion.div>

      {/* Error display */}
      <AnimatePresence>
        {error && (
          <motion.div
            initial={{ opacity: 0, y: -10, scale: 0.98 }}
            animate={{ opacity: 1, y: 0, scale: 1 }}
            exit={{ opacity: 0, y: -10, scale: 0.98 }}
            className="card p-5 border-red-500/30 bg-red-500/5"
          >
            <div className="flex items-start gap-4">
              <div className="w-10 h-10 rounded-xl bg-red-500/10 flex items-center justify-center flex-shrink-0">
                <AlertTriangle className="w-5 h-5 text-red-500" />
              </div>
              <div className="flex-1">
                <h3 className="font-semibold text-red-600 dark:text-red-400 mb-1 font-display">
                  Error Occurred
                </h3>
                <p className="text-sm text-[var(--color-text-secondary)]">{error}</p>
              </div>
              <button
                onClick={() => setError(null)}
                className="p-2 rounded-lg hover:bg-red-500/10 transition-colors"
              >
                <X className="w-5 h-5 text-red-500" />
              </button>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Upload state */}
      <AnimatePresence mode="wait">
        {pageState === 'upload' && (
          <motion.div
            key="upload"
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
            transition={{ duration: 0.4 }}
          >
            <div className="card-glass p-8">
              <div className="flex items-center gap-4 mb-6">
                <div className="w-12 h-12 rounded-2xl bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                  <Upload className="w-6 h-6" />
                </div>
                <div>
                  <h3 className="font-semibold text-[var(--color-text-primary)] text-lg font-display">
                    Upload Your File
                  </h3>
                  <p className="text-sm text-[var(--color-text-muted)]">
                    Supported formats: SDF, CSV (with SMILES column)
                  </p>
                </div>
              </div>
              <ProfileSidebar
                selectedProfileId={selectedProfileId}
                onProfileChange={setSelectedProfileId}
                disabled={pageState !== 'upload'}
              />
              <BatchUpload
                onUploadSuccess={handleUploadSuccess}
                onUploadError={handleUploadError}
                profileId={selectedProfileId}
              />
            </div>

            {/* Features Grid */}
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mt-6">
              <FeatureCard
                icon={<Sparkles className="w-5 h-5" />}
                title="ML-Readiness Scoring"
                description="Automated quality assessment for ML datasets"
              />
              <FeatureCard
                icon={<FileSpreadsheet className="w-5 h-5" />}
                title="Multiple Export Formats"
                description="CSV, JSON, Excel, SDF, and PDF reports"
              />
              <FeatureCard
                icon={<ArrowRight className="w-5 h-5" />}
                title="Real-time Progress"
                description="WebSocket-powered live updates"
              />
            </div>
          </motion.div>
        )}

        {/* Processing state */}
        {pageState === 'processing' && (
          <motion.div
            key="processing"
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
            transition={{ duration: 0.4 }}
          >
            <BatchProgress
              progress={progress}
              isConnected={isConnected}
              onCancel={handleCancel}
              onComplete={handleProcessingComplete}
            />
          </motion.div>
        )}

        {/* Results state */}
        {pageState === 'results' && resultsData && (
          <motion.div
            key="results"
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
            transition={{ duration: 0.4 }}
            className="space-y-6"
          >
            {/* Top bar: quick-nav + Start New Batch */}
            <div className="sticky top-0 z-20 -mx-4 sm:-mx-6 px-4 sm:px-6 py-3 bg-[var(--color-bg)]/80 backdrop-blur-md border-b border-[var(--color-border)]/50">
              <div className="flex items-center justify-between gap-4">
                {/* Quick-nav clay pills */}
                <nav className="flex gap-2.5" aria-label="Jump to section">
                  {includeAnalytics && (
                    <ClayButton
                      variant="accent"
                      size="sm"
                      onClick={() => document.getElementById('section-analytics')?.scrollIntoView({ behavior: 'smooth', block: 'start' })}
                      leftIcon={<BarChart3 className="w-3.5 h-3.5" />}
                    >
                      Analytics
                    </ClayButton>
                  )}
                  <ClayButton
                    variant="stone"
                    size="sm"
                    onClick={() => document.getElementById('section-results')?.scrollIntoView({ behavior: 'smooth', block: 'start' })}
                    leftIcon={<FileSpreadsheet className="w-3.5 h-3.5" />}
                  >
                    Detailed Results
                  </ClayButton>
                </nav>

                <div className="flex items-center gap-2">
                  <ClayButton
                    size="sm"
                    onClick={handleShare}
                    leftIcon={shareCopied ? <Check className="w-3.5 h-3.5 text-green-500" /> : <Share2 className="w-3.5 h-3.5" />}
                  >
                    {shareCopied ? 'Link copied!' : 'Share'}
                  </ClayButton>
                  <ClayButton
                    variant="primary"
                    onClick={handleStartNew}
                    leftIcon={<RotateCcw className="w-4 h-4" />}
                  >
                    Start New Batch
                  </ClayButton>
                </div>
              </div>
            </div>

            {/* Summary */}
            {resultsData.statistics && jobId && (
              <div className="card p-6">
                <div className="flex items-center gap-3 mb-5">
                  <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                    <Sparkles className="w-5 h-5" />
                  </div>
                  <div>
                    <h3 className="text-lg font-semibold text-[var(--color-text-primary)] font-display">
                      Validation Summary
                    </h3>
                    <p className="text-xs text-[var(--color-text-muted)]">
                      Overview of batch processing results
                    </p>
                  </div>
                </div>
                <BatchSummary
                  jobId={jobId}
                  statistics={resultsData.statistics}
                  selectedIndices={selectedIndices}
                  onClearSelection={handleClearSelection}
                />
              </div>
            )}

            {/* Analytics & Visualizations (overview first) */}
            {jobId && includeAnalytics && (
              <div id="section-analytics" className="card p-6 scroll-mt-20">
                <div className="flex items-center gap-3 mb-5">
                  <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
                    <BarChart3 className="w-5 h-5" />
                  </div>
                  <div>
                    <h3 className="text-lg font-semibold text-[var(--color-text-primary)] font-display">
                      Analytics &amp; Visualizations
                    </h3>
                    <p className="text-xs text-[var(--color-text-muted)]">
                      Interactive charts for batch analysis
                    </p>
                  </div>
                </div>

                {/* Batch Timeline (VIZ-09) */}
                {resultsData.statistics && (
                  <div className="mb-6">
                    <BatchTimeline
                      statistics={resultsData.statistics}
                      analyticsStatus={analyticsData?.status ?? null}
                    />
                  </div>
                )}

                <BatchAnalyticsPanel
                  jobId={jobId}
                  statistics={resultsData.statistics}
                  results={resultsData.results}
                  selectedIndices={selectedIndices}
                  onSelectionChange={handleSelectionChange}
                  analyticsData={analyticsData}
                  analyticsStatus={analyticsStatus}
                  analyticsError={analyticsError}
                  analyticsProgress={analyticsProgress}
                  onRetrigger={analyticsRetrigger}
                  onCompare={handleCompare}
                  onOpenActions={() => setSubsetPanelOpen(true)}
                />
              </div>
            )}

            {/* Detailed Results table (drill-down) */}
            <div id="section-results" className="card p-6 scroll-mt-20">
              <div className="flex items-center gap-3 mb-5">
                <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-accent)]/10 to-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-accent)]">
                  <FileSpreadsheet className="w-5 h-5" />
                </div>
                <div>
                  <h3 className="text-lg font-semibold text-[var(--color-text-primary)] font-display">
                    Detailed Results
                  </h3>
                  <p className="text-xs text-[var(--color-text-muted)]">
                    {resultsData.total_results.toLocaleString()} molecules processed
                  </p>
                </div>
              </div>
              <BatchResultsTable
                results={resultsData.results}
                page={page}
                pageSize={pageSize}
                totalResults={resultsData.total_results}
                totalPages={resultsData.total_pages}
                filters={filters}
                sortBy={sortBy}
                sortDir={sortDir}
                onPageChange={handlePageChange}
                onPageSizeChange={handlePageSizeChange}
                onFiltersChange={handleFiltersChange}
                onSortChange={handleSortChange}
                isLoading={resultsLoading}
                selectedIndices={selectedIndices}
                onSelectionChange={handleSelectionChange}
              />
            </div>

            {/* Molecule Comparison Panel (VIZ-07) */}
            {compareLoading && (
              <div className="card p-8 text-center">
                <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-[var(--color-primary)] mx-auto mb-3" />
                <p className="text-sm text-[var(--color-text-muted)]">Fetching molecules for comparison...</p>
              </div>
            )}
            {compareMode && comparisonResults.length > 0 && (
              <MoleculeComparisonPanel
                molecules={comparisonResults}
                datasetStats={analyticsData?.statistics?.property_stats ?? null}
                onClose={handleCloseCompare}
                onRemoveMolecule={handleRemoveFromCompare}
              />
            )}

            {/* Subset Action Panel */}
            {jobId && (
              <SubsetActionPanel
                jobId={jobId}
                selectedIndices={selectedIndices}
                isOpen={subsetPanelOpen}
                onClose={() => setSubsetPanelOpen(false)}
              />
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

interface FeatureCardProps {
  icon: React.ReactNode;
  title: string;
  description: string;
}

function FeatureCard({ icon, title, description }: FeatureCardProps) {
  return (
    <motion.div
      className={cn(
        'p-5 rounded-2xl',
        'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
        'hover:border-[var(--color-primary)]/30 hover:shadow-[0_0_20px_var(--glow-soft)]',
        'transition-all duration-300'
      )}
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      whileHover={{ y: -2 }}
    >
      <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)] mb-3">
        {icon}
      </div>
      <h4 className="font-semibold text-[var(--color-text-primary)] text-sm font-display mb-1">
        {title}
      </h4>
      <p className="text-xs text-[var(--color-text-muted)]">
        {description}
      </p>
    </motion.div>
  );
}
