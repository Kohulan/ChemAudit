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
import { ContradictoryLabelsTab } from '../components/dataset-audit/ContradictoryLabelsTab';
import { DatasetDiffTab } from '../components/dataset-audit/DatasetDiffTab';
import { CurationReportTab } from '../components/dataset-audit/CurationReportTab';
import { InfoTooltip } from '../components/ui/Tooltip';
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
 * Phase 12 Plans 03-05 (complete):
 * 1. Page heading + subtitle
 * 2. DatasetUploadZone (sticky above tabs)
 * 3. DatasetProgressBar (shown during processing)
 * 4. 4-tab layout: Health Audit, Contradictory Labels, Dataset Diff, Curation Report
 * 5. Health Audit tab (Plan 04): 7 components (gauge, sub-scores, weights, treemap,
 *    drill-down, property distributions, std consistency panel)
 * 6. Contradictory Labels tab (Plan 05): activity column selector, fold-threshold filter,
 *    paginated contradiction cards with 2D structures
 * 7. Dataset Diff tab (Plan 05): comparison file upload, summary badges, filterable
 *    molecule table with expandable modified rows
 * 8. Curation Report tab (Plan 05): collapsible preview sections, JSON + CSV downloads
 *
 * Both hooks wired at top level for all tab content components.
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
    [auditState.status, auditState.resetAll, auditState.uploadFile],
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

      {/* Claymorphic tab bar */}
      <div
        className="mb-8 p-1.5 rounded-2xl"
        style={{
          backgroundColor: 'var(--color-surface-sunken)',
          boxShadow: `
            2px 4px 8px 0 rgba(var(--color-primary-rgb), 0.08),
            inset 2px 2px 6px rgba(255, 255, 255, 0.5),
            inset -2px -2px 6px rgba(0, 0, 0, 0.06)
          `,
        }}
        role="tablist"
        aria-label="Dataset audit tabs"
      >
        <div className="flex gap-1">
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
                  'relative flex-1 px-3 py-2.5 text-sm font-medium rounded-xl',
                  'flex items-center justify-center gap-2',
                  'transition-all duration-400 ease-out',
                  'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-offset-2',
                  isActive
                    ? 'text-[var(--color-primary-dark)] font-semibold'
                    : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)]',
                )}
                style={isActive ? {
                  backgroundColor: 'var(--color-surface-elevated)',
                  boxShadow: `
                    1px 2px 4px 0 rgba(var(--color-primary-rgb), 0.12),
                    2px 4px 10px 0 rgba(var(--color-primary-rgb), 0.08),
                    inset 1px 1px 4px rgba(255, 255, 255, 0.6),
                    inset -1px -1px 4px rgba(var(--color-primary-rgb), 0.06)
                  `,
                } : {
                  backgroundColor: 'transparent',
                  boxShadow: 'none',
                }}
              >
                <tab.icon className={cn(
                  'w-4 h-4 transition-transform duration-300',
                  isActive && 'scale-110',
                )} />
                <span className="font-display">{tab.label}</span>
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
              <div className="flex flex-col items-center">
                <div className="flex items-center gap-1.5 mb-2">
                  <InfoTooltip
                    title="Dataset Health Score"
                    content={
                      <div className="text-xs space-y-1">
                        <p>Weighted composite score (0-100) assessing five quality dimensions of your chemical dataset.</p>
                        <ul className="mt-1 text-white/70 space-y-0.5">
                          <li>&ge;80: Excellent &mdash; ready for ML workflows</li>
                          <li>50-79: Fair &mdash; review sub-scores for targeted fixes</li>
                          <li>&lt;50: Poor &mdash; curation strongly recommended</li>
                        </ul>
                        <p className="mt-1 text-white/60">Adjust the weight sliders below to emphasize the dimensions most important to your use case.</p>
                      </div>
                    }
                    size="small"
                  />
                </div>
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

          {/* ================================================================
              Contradictory Labels Tab (Plan 05)
              ================================================================ */}
          {isComplete && activeTab === 'contradictions' && (
            <ContradictoryLabelsTab
              contradictions={auditState.results?.contradictions ?? []}
              numericColumns={auditState.results?.numeric_columns ?? []}
            />
          )}

          {/* ================================================================
              Dataset Diff Tab (Plan 05)
              ================================================================ */}
          {isComplete && activeTab === 'diff' && (
            <DatasetDiffTab
              diffResults={auditState.diffResults}
              diffLoading={auditState.diffLoading}
              diffError={auditState.diffError}
              onUploadDiffFile={auditState.uploadDiffFile}
              primaryJobId={auditState.jobId}
            />
          )}

          {/* ================================================================
              Curation Report Tab (Plan 05)
              ================================================================ */}
          {isComplete && activeTab === 'report' && (
            <CurationReportTab
              report={auditState.results?.curation_report ?? null}
              jobId={auditState.jobId}
              curatedCsvAvailable={auditState.results?.curated_csv_available ?? false}
            />
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
// Placeholder (used for health tab no-data fallback only)
// =============================================================================

function PlaceholderTab({ label }: { label: string }) {
  return (
    <div className="flex items-center justify-center py-12 text-[var(--color-text-muted)] text-sm">
      {label}
    </div>
  );
}
