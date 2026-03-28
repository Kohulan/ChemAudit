import { useState, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Database, AlertTriangle, GitCompare, FileText } from 'lucide-react';
import type { DatasetAuditTab } from '../types/dataset_intelligence';
import { useDatasetAudit } from '../hooks/useDatasetAudit';
import { useDatasetWeights } from '../hooks/useDatasetWeights';
import { DatasetUploadZone } from '../components/dataset-audit/DatasetUploadZone';
import { DatasetProgressBar } from '../components/dataset-audit/DatasetProgressBar';
import { cn } from '../lib/utils';

// =============================================================================
// Tab configuration
// =============================================================================

const TABS: Array<{ id: DatasetAuditTab; label: string; icon: typeof Database }> = [
  { id: 'health', label: 'Health Audit', icon: Database },
  { id: 'contradictions', label: 'Contradictory Labels', icon: AlertTriangle },
  { id: 'diff', label: 'Dataset Diff', icon: GitCompare },
  { id: 'report', label: 'Curation Report', icon: FileText },
];

// =============================================================================
// Page component
// =============================================================================

/**
 * Dataset Audit page -- /dataset-audit
 *
 * Phase 12 Plan 03:
 * 1. Page heading + subtitle
 * 2. DatasetUploadZone (sticky above tabs)
 * 3. DatasetProgressBar (shown during processing)
 * 4. 4-tab layout: Health Audit, Contradictory Labels, Dataset Diff, Curation Report
 * 5. Empty states before upload; placeholder divs after upload for Plan 04/05
 *
 * Both hooks wired at top level so Plan 04/05 can add tab content
 * components without prop drilling restructuring.
 */
export default function DatasetAudit() {
  const [activeTab, setActiveTab] = useState<DatasetAuditTab>('health');

  // Wire both hooks at page top level
  const auditState = useDatasetAudit();
  const weightState = useDatasetWeights();

  // Suppress unused variable warnings for Plan 04/05
  void weightState;

  // Handle file selection: reset + upload
  const handleFileSelect = useCallback(
    (file: File) => {
      if (auditState.status === 'complete' || auditState.status === 'error') {
        auditState.resetAll();
      }
      void auditState.uploadFile(file);
    },
    [auditState],
  );

  const isComplete = auditState.status === 'complete';

  return (
    <div className="max-w-[1200px] mx-auto px-4 pt-16 pb-16">
      {/* Page heading */}
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.3 }}
      >
        <h1 className="text-2xl font-semibold font-display mb-2">
          Dataset Audit
        </h1>
        <p className="text-sm text-[var(--color-text-secondary)] mb-8">
          Upload a compound dataset to assess quality, detect contradictions,
          compare versions, and generate reproducible curation reports.
        </p>
      </motion.div>

      {/* Upload zone -- sticky above tabs */}
      <div className="mb-6">
        <DatasetUploadZone
          onFileSelect={handleFileSelect}
          disabled={auditState.status === 'uploading' || auditState.status === 'processing'}
          file={auditState.file}
          status={auditState.status}
          progress={auditState.progress}
        />
      </div>

      {/* Error display */}
      {auditState.error && (
        <div className="mb-6 p-4 rounded-xl bg-red-500/10 border border-red-500/20">
          <p className="text-sm text-red-600 dark:text-red-400">{auditState.error}</p>
        </div>
      )}

      {/* Progress bar -- shown during processing */}
      {(auditState.status === 'processing' || auditState.status === 'uploading') && (
        <div className="mb-6">
          <DatasetProgressBar
            progress={auditState.progress}
            currentStage={auditState.currentStage}
            status={auditState.status}
          />
        </div>
      )}

      {/* Tab bar */}
      <div
        className="border-b border-[var(--color-border)] mb-6"
        role="tablist"
        aria-label="Dataset audit tabs"
      >
        <div className="flex gap-0">
          {TABS.map((tab) => {
            const isActive = activeTab === tab.id;
            return (
              <button
                key={tab.id}
                role="tab"
                id={`tab-${tab.id}`}
                aria-selected={isActive}
                aria-controls={`tabpanel-${tab.id}`}
                onClick={() => setActiveTab(tab.id)}
                className={cn(
                  'relative px-4 py-3 text-sm font-medium transition-colors duration-200',
                  'flex items-center gap-2',
                  'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-inset',
                  isActive
                    ? 'text-[var(--color-primary)] font-semibold'
                    : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)]',
                )}
              >
                <tab.icon className="w-4 h-4" />
                <span className="font-display">{tab.label}</span>
                {/* Active underline */}
                {isActive && (
                  <motion.div
                    layoutId="datasetTabUnderline"
                    className="absolute bottom-0 left-0 right-0 h-0.5 bg-[var(--color-primary)]"
                    transition={{ type: 'spring', stiffness: 380, damping: 28, mass: 0.8 }}
                  />
                )}
              </button>
            );
          })}
        </div>
      </div>

      {/* Tab content */}
      <AnimatePresence mode="wait">
        <motion.div
          key={activeTab}
          initial={{ opacity: 0, y: 8 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -8 }}
          transition={{ duration: 0.2, ease: 'easeOut' }}
          role="tabpanel"
          id={`tabpanel-${activeTab}`}
          aria-labelledby={`tab-${activeTab}`}
        >
          {/* Before upload: empty state */}
          {!isComplete && auditState.status !== 'processing' && auditState.status !== 'uploading' && (
            <EmptyState />
          )}

          {/* After upload complete: tab-specific content */}
          {isComplete && activeTab === 'health' && (
            <div className="min-h-[200px]">
              {/* Plan 04: HealthAuditTab content */}
              <PlaceholderTab label="Health Audit results will appear here" />
            </div>
          )}
          {isComplete && activeTab === 'contradictions' && (
            <div className="min-h-[200px]">
              {/* Plan 05: ContradictoryLabelsTab */}
              <PlaceholderTab label="Contradictory Labels results will appear here" />
            </div>
          )}
          {isComplete && activeTab === 'diff' && (
            <div className="min-h-[200px]">
              {/* Plan 05: DatasetDiffTab */}
              <PlaceholderTab label="Dataset Diff tool will appear here" />
            </div>
          )}
          {isComplete && activeTab === 'report' && (
            <div className="min-h-[200px]">
              {/* Plan 05: CurationReportTab */}
              <PlaceholderTab label="Curation Report will appear here" />
            </div>
          )}
        </motion.div>
      </AnimatePresence>
    </div>
  );
}

// =============================================================================
// Empty state
// =============================================================================

/** Empty state shown in all tabs before a dataset is uploaded. */
function EmptyState() {
  return (
    <div className="flex flex-col items-center justify-center py-16 text-center">
      <div className="w-16 h-16 rounded-2xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center mb-4">
        <Database className="w-8 h-8 text-[var(--color-text-muted)]" />
      </div>
      <h3 className="text-lg font-semibold font-display text-[var(--color-text-primary)] mb-2">
        No dataset loaded
      </h3>
      <p className="text-sm text-[var(--color-text-secondary)] max-w-md">
        Drop a CSV or SDF file above to assess dataset health, detect
        contradictory labels, and generate a curation report.
      </p>
    </div>
  );
}

// =============================================================================
// Placeholder tab content (Plan 04/05 will replace these)
// =============================================================================

function PlaceholderTab({ label }: { label: string }) {
  return (
    <div className="flex items-center justify-center py-12 text-[var(--color-text-muted)] text-sm">
      {label}
    </div>
  );
}
