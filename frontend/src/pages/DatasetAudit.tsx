import { useState, useCallback, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Database, AlertTriangle, GitCompare, FileText } from 'lucide-react';
import type { DatasetAuditTab } from '../types/dataset_intelligence';
import { useDatasetAudit } from '../hooks/useDatasetAudit';
import { useDatasetWeights } from '../hooks/useDatasetWeights';
import { DatasetUploadZone } from '../components/dataset-audit/DatasetUploadZone';
import { DatasetProgressBar } from '../components/dataset-audit/DatasetProgressBar';
import { HealthScoreGauge } from '../components/dataset-audit/HealthScoreGauge';
import { SubScoreCards } from '../components/dataset-audit/SubScoreCards';
import { WeightSliders } from '../components/dataset-audit/WeightSliders';
import { IssueTreemap } from '../components/dataset-audit/IssueTreemap';
import { TreemapDrillDown } from '../components/dataset-audit/TreemapDrillDown';
import { PropertyDistOverlay } from '../components/dataset-audit/PropertyDistOverlay';
import { StdConsistencyPanel } from '../components/dataset-audit/StdConsistencyPanel';
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
 * Phase 12 Plan 03 (skeleton) + Plan 04 (Health Audit tab):
 * 1. Page heading + subtitle
 * 2. DatasetUploadZone (sticky above tabs)
 * 3. DatasetProgressBar (shown during processing)
 * 4. 4-tab layout: Health Audit, Contradictory Labels, Dataset Diff, Curation Report
 * 5. Health Audit tab: 7 components (gauge, sub-scores, weights, treemap, drill-down,
 *    property distributions, std consistency panel)
 *
 * Both hooks wired at top level so Plan 05 can add tab content
 * components without prop drilling restructuring.
 */
export default function DatasetAudit() {
  const [activeTab, setActiveTab] = useState<DatasetAuditTab>('health');
  const [selectedIssueType, setSelectedIssueType] = useState<string | null>(null);

  // Wire both hooks at page top level
  const auditState = useDatasetAudit();
  const weightState = useDatasetWeights();

  // Ref for scrolling to treemap on sub-score card click
  const treemapRef = useRef<HTMLDivElement>(null);

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

  // Handle sub-score card click: scroll to treemap and highlight matching category
  const handleSubScoreCardClick = useCallback((name: string) => {
    // Map sub-score name to closest issue type category for highlighting
    const nameToCategory: Record<string, string> = {
      parsability: 'Parse failures',
      stereo: 'Undefined stereo',
      uniqueness: 'Duplicates',
      alerts: 'Alert hits',
      std_consistency: 'Std. inconsistency',
    };
    const category = nameToCategory[name] ?? null;
    setSelectedIssueType(category);
    // Scroll treemap into view
    treemapRef.current?.scrollIntoView({ behavior: 'smooth', block: 'center' });
  }, []);

  // Handle treemap cell click: toggle drill-down
  const handleTreemapCellClick = useCallback((issueType: string) => {
    setSelectedIssueType((prev) => {
      // Empty string from toggle = clear selection
      if (!issueType) return null;
      return prev === issueType ? null : issueType;
    });
  }, []);

  const isComplete = auditState.status === 'complete';
  const healthAudit = auditState.results?.health_audit ?? null;

  // Compute overall score using current weights
  const overallScore = healthAudit
    ? weightState.computeOverallScore(healthAudit.sub_scores)
    : 0;

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

          {/* ================================================================
              Health Audit Tab (Plan 04)
              ================================================================ */}
          {isComplete && activeTab === 'health' && healthAudit && (
            <div className="min-h-[200px] space-y-8">
              {/* Row 1: Health Score Gauge (centered) */}
              <div className="flex justify-center">
                <HealthScoreGauge score={overallScore} />
              </div>

              {/* Row 2: Sub-Score Cards (5 cards) */}
              <SubScoreCards
                subScores={healthAudit.sub_scores}
                onCardClick={handleSubScoreCardClick}
              />

              {/* Row 3: Weight Sliders (30%) + Issue Treemap (70%) */}
              <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                <div className="lg:col-span-1">
                  <WeightSliders
                    weights={weightState.weights}
                    activeProfileId={weightState.activeProfileId}
                    savedProfiles={weightState.savedProfiles}
                    onUpdateWeight={weightState.updateWeight}
                    onSaveProfile={weightState.saveProfile}
                    onDeleteProfile={weightState.deleteProfile}
                    onLoadProfile={weightState.loadProfile}
                    onResetToDefaults={weightState.resetToDefaults}
                  />
                </div>
                <div className="lg:col-span-2" ref={treemapRef}>
                  <IssueTreemap
                    issues={healthAudit.issues}
                    onCellClick={handleTreemapCellClick}
                    selectedIssueType={selectedIssueType}
                  />
                </div>
              </div>

              {/* Row 4: Treemap Drill-Down (full width, expandable) */}
              <TreemapDrillDown
                issueType={selectedIssueType}
                issues={healthAudit.issues}
                isOpen={!!selectedIssueType}
              />

              {/* Row 5: Property Distribution Overlay (3 histograms) */}
              <PropertyDistOverlay distributions={healthAudit.property_distributions} />

              {/* Row 6: Standardization Consistency Panel (full width) */}
              <StdConsistencyPanel
                comparison={healthAudit.std_pipeline_comparison}
                sampleSize={healthAudit.std_sample_size}
              />
            </div>
          )}

          {/* Health tab with no health_audit data */}
          {isComplete && activeTab === 'health' && !healthAudit && (
            <div className="min-h-[200px]">
              <PlaceholderTab label="Health audit data is not available for this dataset." />
            </div>
          )}

          {/* Other tabs (Plan 05 will replace placeholders) */}
          {isComplete && activeTab === 'contradictions' && (
            <div className="min-h-[200px]">
              <PlaceholderTab label="Contradictory Labels results will appear here" />
            </div>
          )}
          {isComplete && activeTab === 'diff' && (
            <div className="min-h-[200px]">
              <PlaceholderTab label="Dataset Diff tool will appear here" />
            </div>
          )}
          {isComplete && activeTab === 'report' && (
            <div className="min-h-[200px]">
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
// Placeholder tab content (Plan 05 will replace these)
// =============================================================================

function PlaceholderTab({ label }: { label: string }) {
  return (
    <div className="flex items-center justify-center py-12 text-[var(--color-text-muted)] text-sm">
      {label}
    </div>
  );
}
