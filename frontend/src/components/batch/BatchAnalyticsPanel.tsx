/**
 * BatchAnalyticsPanel
 *
 * Tabbed container for batch analytics charts.
 * Two tabs: "Distributions" and "Chemical Space" with skeleton loading,
 * error states, and summary badges.
 */

import React, { useState, useCallback, useRef, useMemo } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Download, AlertTriangle, RotateCcw, Loader2, GitCompare, X, Beaker } from 'lucide-react';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  Tooltip as RechartsTooltip,
  ResponsiveContainer,
  CartesianGrid,
  Cell,
} from 'recharts';
import { ScoreHistogram } from './charts/ScoreHistogram';
import { PropertyScatterPlot } from './charts/PropertyScatterPlot';
import { AlertFrequencyChart } from './charts/AlertFrequencyChart';
import { ValidationTreemap } from './charts/ValidationTreemap';
import { ScaffoldTreemap } from './charts/ScaffoldTreemap';
import { ChemicalSpaceScatter } from './charts/ChemicalSpaceScatter';
import { ClayButton } from '../ui/ClayButton';
import type { BatchStatistics, BatchResult } from '../../types/batch';
import type { AnalyticsHookStatus, AnalyticsProgressInfo } from '../../hooks/useBatchAnalytics';
import { cn } from '../../lib/utils';

interface BatchAnalyticsPanelProps {
  jobId: string;
  statistics: BatchStatistics | null;
  results: BatchResult[];
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
  analyticsData: import('../../types/analytics').BatchAnalyticsResponse | null;
  analyticsStatus?: AnalyticsHookStatus;
  analyticsError?: string | null;
  analyticsProgress?: AnalyticsProgressInfo;
  onRetrigger: (type: string) => void;
  onCompare?: () => void;
}

const TABS = ['Distributions', 'Chemical Space'] as const;
type TabName = (typeof TABS)[number];

function ChartSkeleton() {
  return (
    <div className="animate-pulse space-y-3">
      <div className="h-4 w-32 bg-[var(--color-border)] rounded" />
      <div className="h-[240px] bg-[var(--color-border)]/50 rounded-xl" />
    </div>
  );
}

interface ChartCardProps {
  title: string;
  icon: React.ReactNode;
  onDownload?: () => void;
  ariaLabel: string;
  caption: string;
  children: React.ReactNode;
}

function ChartCard({ title, icon, onDownload, ariaLabel, caption, children }: ChartCardProps) {
  return (
    <figure
      aria-label={ariaLabel}
      className="rounded-2xl p-5 bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)] border border-[var(--color-border)]"
    >
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center gap-2">
          <span className="text-[var(--color-primary)]">{icon}</span>
          <h4 className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
            {title}
          </h4>
        </div>
        {onDownload && (
          <button
            onClick={onDownload}
            className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]"
            title="Download as PNG"
          >
            <Download className="w-4 h-4" />
          </button>
        )}
      </div>
      {children}
      <figcaption className="sr-only">{caption}</figcaption>
    </figure>
  );
}

interface ErrorCardProps {
  message: string;
  onRetry: () => void;
}

function ErrorCard({ message, onRetry }: ErrorCardProps) {
  return (
    <div className="rounded-2xl p-5 bg-amber-500/5 border border-amber-500/20">
      <div className="flex items-center gap-3">
        <AlertTriangle className="w-5 h-5 text-amber-500 flex-shrink-0" />
        <div className="flex-1">
          <p className="text-sm text-amber-700 dark:text-amber-300">{message}</p>
        </div>
        <button
          onClick={onRetry}
          className="p-1.5 rounded-lg bg-amber-500/10 hover:bg-amber-500/20 transition-colors text-amber-600 dark:text-amber-400"
          title="Retry"
        >
          <RotateCcw className="w-4 h-4" />
        </button>
      </div>
    </div>
  );
}

/**
 * Download SVG chart as PNG by serializing the SVG element inside a container.
 */
function downloadSvgAsPng(containerRef: React.RefObject<HTMLDivElement | null>, filename: string) {
  const container = containerRef.current;
  if (!container) return;

  const svgEl = container.querySelector('svg');
  if (!svgEl) return;

  const serializer = new XMLSerializer();
  const svgStr = serializer.serializeToString(svgEl);
  const svgBlob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(svgBlob);

  const img = new Image();
  img.onload = () => {
    const canvas = document.createElement('canvas');
    const dpr = window.devicePixelRatio || 1;
    canvas.width = img.width * dpr;
    canvas.height = img.height * dpr;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    ctx.scale(dpr, dpr);
    ctx.drawImage(img, 0, 0);

    canvas.toBlob((blob) => {
      if (!blob) return;
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = filename;
      a.click();
      URL.revokeObjectURL(a.href);
    }, 'image/png');

    URL.revokeObjectURL(url);
  };
  img.src = url;
}

function ProfileScoreHistogram({ results }: { results: BatchResult[] }) {
  const buckets = { Excellent: 0, Good: 0, Moderate: 0, Poor: 0 };
  for (const r of results) {
    const score = r.scoring?.profile?.score;
    if (score == null) continue;
    if (score >= 80) buckets.Excellent++;
    else if (score >= 50) buckets.Good++;
    else if (score >= 20) buckets.Moderate++;
    else buckets.Poor++;
  }
  const data = [
    { name: 'Excellent', count: buckets.Excellent, fill: '#22c55e' },
    { name: 'Good', count: buckets.Good, fill: '#eab308' },
    { name: 'Moderate', count: buckets.Moderate, fill: '#f97316' },
    { name: 'Poor', count: buckets.Poor, fill: '#ef4444' },
  ];
  return (
    <ResponsiveContainer width="100%" height={250}>
      <BarChart data={data} margin={{ top: 5, right: 20, bottom: 5, left: 0 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="var(--color-border)" />
        <XAxis dataKey="name" tick={{ fill: 'var(--color-text-muted)', fontSize: 12 }} />
        <YAxis allowDecimals={false} tick={{ fill: 'var(--color-text-muted)', fontSize: 12 }} />
        <RechartsTooltip
          contentStyle={{
            backgroundColor: 'var(--color-surface-elevated)',
            border: '1px solid var(--color-border)',
            borderRadius: '8px',
            color: 'var(--color-text-primary)',
          }}
        />
        <Bar dataKey="count" radius={[4, 4, 0, 0]}>
          {data.map((entry, index) => (
            <Cell key={index} fill={entry.fill} />
          ))}
        </Bar>
      </BarChart>
    </ResponsiveContainer>
  );
}

function AnalyticsProgressBar({ progress, status, error, onRetrigger }: {
  progress?: AnalyticsProgressInfo;
  status?: AnalyticsHookStatus;
  error?: string | null;
  onRetrigger: (type: string) => void;
}) {
  if (!progress || status === 'idle' || status === 'complete') return null;

  const pct = progress.totalCount > 0
    ? Math.round((progress.completedCount / progress.totalCount) * 100)
    : 0;

  const formatTime = (s: number) => {
    if (s < 60) return `${s}s`;
    return `${Math.floor(s / 60)}m ${s % 60}s`;
  };

  const typeLabels: Record<string, string> = {
    deduplication: 'Deduplication',
    statistics: 'Statistics',
    scaffold: 'Scaffold Diversity',
    chemical_space: 'Chemical Space (PCA)',
  };

  return (
    <div className="mb-4 rounded-xl border border-[var(--color-border)] bg-[var(--color-surface-sunken)]/50 p-4 space-y-3">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          {status === 'computing' && <Loader2 className="w-4 h-4 text-[var(--color-primary)] animate-spin" />}
          {status === 'error' && <AlertTriangle className="w-4 h-4 text-amber-500" />}
          <span className="text-sm font-medium text-[var(--color-text-primary)]">
            {status === 'computing' ? 'Computing Analytics...' : 'Analytics Error'}
          </span>
        </div>
        <span className="text-xs text-[var(--color-text-muted)]">
          {formatTime(progress.elapsedSeconds)} elapsed
        </span>
      </div>

      {/* Progress bar */}
      <div className="h-2 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden">
        <div
          className="h-full rounded-full transition-all duration-500 ease-out"
          style={{
            width: `${Math.max(status === 'computing' ? 5 : 0, pct)}%`,
            background: status === 'error'
              ? 'linear-gradient(90deg, #f59e0b, #ef4444)'
              : 'linear-gradient(90deg, #c41e3a, #e11d48, #d97706)',
          }}
        />
      </div>

      {/* Per-type status */}
      <div className="flex flex-wrap gap-2">
        {Object.entries(progress.typeStatuses).map(([type, st]) => (
          <span
            key={type}
            className={cn(
              'inline-flex items-center gap-1.5 px-2.5 py-1 rounded-md text-xs font-medium',
              st.status === 'complete' && 'bg-green-500/10 text-green-600 dark:text-green-400',
              st.status === 'computing' && 'bg-[var(--color-primary)]/10 text-[var(--color-primary)]',
              st.status === 'pending' && 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]',
              st.status === 'failed' && 'bg-red-500/10 text-red-500',
            )}
          >
            {st.status === 'complete' && <span className="w-1.5 h-1.5 rounded-full bg-green-500" />}
            {st.status === 'computing' && <Loader2 className="w-3 h-3 animate-spin" />}
            {st.status === 'failed' && <AlertTriangle className="w-3 h-3" />}
            {typeLabels[type] || type}
          </span>
        ))}
      </div>

      {/* Error message with retry */}
      {error && (
        <div className="flex items-center justify-between gap-2 p-2 rounded-lg bg-amber-500/5 border border-amber-500/20">
          <p className="text-xs text-amber-600 dark:text-amber-400">{error}</p>
          <button
            onClick={() => {
              for (const type of Object.keys(progress.typeStatuses)) {
                const st = progress.typeStatuses[type];
                if (st.status !== 'complete') {
                  onRetrigger(type);
                }
              }
            }}
            className="flex-shrink-0 p-1.5 rounded-lg bg-amber-500/10 hover:bg-amber-500/20 transition-colors text-amber-600 dark:text-amber-400"
            title="Retry failed analytics"
          >
            <RotateCcw className="w-3.5 h-3.5" />
          </button>
        </div>
      )}
    </div>
  );
}

export const BatchAnalyticsPanel = React.memo(function BatchAnalyticsPanel({
  jobId,
  statistics,
  results,
  selectedIndices,
  onSelectionChange,
  analyticsData,
  analyticsStatus,
  analyticsError,
  analyticsProgress,
  onRetrigger,
  onCompare,
}: BatchAnalyticsPanelProps) {
  const [activeTab, setActiveTab] = useState<TabName>('Distributions');
  const [scatterXProp, setScatterXProp] = useState('MW');
  const [scatterYProp, setScatterYProp] = useState('LogP');
  const [scatterColorProp, setScatterColorProp] = useState('overall_score');
  const [chemSpaceColorProp, setChemSpaceColorProp] = useState('overall_score');

  // Chart container refs for SVG download
  const scoreHistRef = useRef<HTMLDivElement>(null);
  const alertFreqRef = useRef<HTMLDivElement>(null);
  const validTreemapRef = useRef<HTMLDivElement>(null);
  const propScatterRef = useRef<HTMLDivElement>(null);
  const scaffoldTreemapRef = useRef<HTMLDivElement>(null);
  const profileHistRef = useRef<HTMLDivElement>(null);

  // Badge counts
  const outlierCount = analyticsData?.statistics?.outliers?.length ?? 0;
  const scaffoldCount = analyticsData?.scaffold?.unique_scaffold_count ?? 0;
  const [outlierPopoverOpen, setOutlierPopoverOpen] = useState(false);
  const outlierBadgeRef = useRef<HTMLSpanElement>(null);

  // Group outliers by property for the popover
  const outlierBreakdown = useMemo(() => {
    const outliers = analyticsData?.statistics?.outliers;
    if (!outliers || outliers.length === 0) return null;
    const byProp: Record<string, number> = {};
    const uniqueMols = new Set<number>();
    for (const o of outliers) {
      byProp[o.property_name] = (byProp[o.property_name] || 0) + 1;
      uniqueMols.add(o.molecule_index);
    }
    return {
      byProperty: Object.entries(byProp).sort((a, b) => b[1] - a[1]),
      uniqueMoleculeCount: uniqueMols.size,
    };
  }, [analyticsData?.statistics?.outliers]);

  const handlePropertyChange = useCallback(
    (axis: 'x' | 'y' | 'color', property: string) => {
      if (axis === 'x') setScatterXProp(property);
      else if (axis === 'y') setScatterYProp(property);
      else setScatterColorProp(property);
    },
    []
  );

  const isChemSpaceStatus = (type: string) =>
    analyticsData?.status?.[type];

  return (
    <div className="space-y-4">
      {/* Tab bar */}
      <div className="flex gap-1 p-1 rounded-xl bg-[var(--color-surface-sunken)] w-fit">
        {TABS.map((tab) => (
          <button
            key={tab}
            onClick={() => setActiveTab(tab)}
            className={cn(
              'px-4 py-2 rounded-lg text-sm font-medium transition-all duration-200 flex items-center gap-2',
              activeTab === tab
                ? 'bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] shadow-sm'
                : 'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]'
            )}
          >
            {tab}
            {/* Outlier badge with hover popover */}
            {tab === 'Distributions' && outlierCount > 0 && (
              <span
                ref={outlierBadgeRef}
                className="relative text-xs px-1.5 py-0.5 rounded-full bg-amber-500/10 text-amber-600 dark:text-amber-400 cursor-help"
                onMouseEnter={() => setOutlierPopoverOpen(true)}
                onMouseLeave={() => setOutlierPopoverOpen(false)}
              >
                {outlierCount} outlier{outlierCount !== 1 ? 's' : ''}
                <AnimatePresence>
                  {outlierPopoverOpen && outlierBreakdown && (
                    <motion.div
                      initial={{ opacity: 0, y: 4, scale: 0.96 }}
                      animate={{ opacity: 1, y: 0, scale: 1 }}
                      exit={{ opacity: 0, y: 4, scale: 0.96 }}
                      transition={{ duration: 0.15 }}
                      className="absolute top-full left-1/2 -translate-x-1/2 mt-2 z-50 w-64 p-3 rounded-xl bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-xl"
                      onClick={(e) => e.stopPropagation()}
                    >
                      <p className="text-xs font-semibold text-[var(--color-text-primary)] mb-1.5">
                        {outlierBreakdown.uniqueMoleculeCount} molecule{outlierBreakdown.uniqueMoleculeCount !== 1 ? 's' : ''} with values outside statistical bounds
                      </p>
                      <div className="space-y-1 mb-2">
                        {outlierBreakdown.byProperty.map(([prop, count]) => (
                          <div key={prop} className="flex items-center justify-between text-xs">
                            <span className="text-[var(--color-text-secondary)] font-medium">{prop}</span>
                            <span className="text-amber-600 dark:text-amber-400 font-mono tabular-nums">
                              {count} outlier{count !== 1 ? 's' : ''}
                            </span>
                          </div>
                        ))}
                      </div>
                      <p className="text-[10px] text-[var(--color-text-muted)] leading-snug border-t border-[var(--color-border)] pt-2">
                        Use the Property Scatter Plot to identify and select these molecules
                      </p>
                    </motion.div>
                  )}
                </AnimatePresence>
              </span>
            )}
            {tab === 'Chemical Space' && scaffoldCount > 0 && (
              <span className="text-xs px-1.5 py-0.5 rounded-full bg-[var(--color-primary)]/10 text-[var(--color-primary)]">
                {scaffoldCount} scaffold{scaffoldCount !== 1 ? 's' : ''}
              </span>
            )}
          </button>
        ))}
      </div>

      {/* Selection toolbar — contextual compare action */}
      <AnimatePresence>
        {selectedIndices.size > 0 && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: 'auto' }}
            exit={{ opacity: 0, height: 0 }}
            transition={{ duration: 0.2 }}
            className="overflow-hidden"
          >
            <div className="flex items-center gap-3 px-4 py-2.5 rounded-xl bg-[var(--color-primary)]/5 border border-[var(--color-primary)]/20">
              <div className="flex items-center gap-2 flex-1 min-w-0">
                <GitCompare className="w-4 h-4 text-[var(--color-primary)] flex-shrink-0" />
                <span className="text-sm font-medium text-[var(--color-text-primary)]">
                  {selectedIndices.size} molecule{selectedIndices.size !== 1 ? 's' : ''} selected
                </span>
                {selectedIndices.size > 2 && (
                  <span className="text-xs text-[var(--color-text-muted)]">
                    (select up to 2 to compare)
                  </span>
                )}
              </div>
              <div className="flex items-center gap-2 flex-shrink-0">
                {onCompare && selectedIndices.size >= 1 && selectedIndices.size <= 2 && (
                  <ClayButton
                    variant="primary"
                    size="sm"
                    onClick={onCompare}
                    leftIcon={<GitCompare className="w-3.5 h-3.5" />}
                  >
                    Compare
                  </ClayButton>
                )}
                <ClayButton
                  variant="ghost"
                  size="sm"
                  onClick={() => onSelectionChange(new Set())}
                  leftIcon={<X className="w-3.5 h-3.5" />}
                >
                  Clear
                </ClayButton>
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Analytics progress bar */}
      <AnalyticsProgressBar
        progress={analyticsProgress}
        status={analyticsStatus}
        error={analyticsError}
        onRetrigger={onRetrigger}
      />

      {/* Tab content */}
      <AnimatePresence mode="wait">
        {activeTab === 'Distributions' && (
          <motion.div
            key="distributions"
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -10 }}
            transition={{ duration: 0.2 }}
            className="grid grid-cols-1 lg:grid-cols-2 gap-4"
          >
            {/* VIZ-01: Score Histogram */}
            <ChartCard
              title="Score Distribution"
              icon={<svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24"><rect x="3" y="12" width="4" height="9" rx="1" /><rect x="10" y="6" width="4" height="15" rx="1" /><rect x="17" y="3" width="4" height="18" rx="1" /></svg>}
              onDownload={() => downloadSvgAsPng(scoreHistRef, 'score-distribution.png')}
              ariaLabel="Score distribution histogram showing excellent, good, moderate, and poor categories"
              caption="Distribution of validation scores across the batch"
            >
              <div ref={scoreHistRef}>
                {statistics ? (
                  <ScoreHistogram
                    data={statistics.score_distribution}
                    results={results}
                    selectedIndices={selectedIndices}
                    onSelectionChange={onSelectionChange}
                  />
                ) : (
                  <ChartSkeleton />
                )}
              </div>
            </ChartCard>

            {/* VIZ-03: Alert Frequency */}
            <ChartCard
              title="Alert Frequency"
              icon={<AlertTriangle className="w-4 h-4" />}
              onDownload={() => downloadSvgAsPng(alertFreqRef, 'alert-frequency.png')}
              ariaLabel="Alert frequency bar chart showing alert type counts"
              caption="Frequency of structural alerts detected in the batch"
            >
              <div ref={alertFreqRef}>
                {statistics ? (
                  <AlertFrequencyChart
                    data={statistics.alert_summary}
                    selectedIndices={selectedIndices}
                    onSelectionChange={onSelectionChange}
                  />
                ) : (
                  <ChartSkeleton />
                )}
              </div>
            </ChartCard>

            {/* VIZ: Profile Score Distribution (conditional) */}
            {results.some(r => r.scoring?.profile?.score != null) && (
              <ChartCard
                title="Profile Score Distribution"
                icon={<Beaker className="w-4 h-4" />}
                onDownload={() => downloadSvgAsPng(profileHistRef, 'profile-score-distribution.png')}
                ariaLabel="Profile score distribution histogram"
                caption="Distribution of profile desirability scores across the batch"
              >
                <div ref={profileHistRef}>
                  <ProfileScoreHistogram results={results} />
                </div>
              </ChartCard>
            )}

            {/* VIZ-04: Validation Treemap */}
            <ChartCard
              title="Validation Issues"
              icon={<svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24"><rect x="2" y="2" width="9" height="9" rx="1" /><rect x="13" y="2" width="9" height="5" rx="1" /><rect x="13" y="9" width="9" height="6" rx="1" /><rect x="2" y="13" width="9" height="9" rx="1" /><rect x="13" y="17" width="9" height="5" rx="1" /></svg>}
              onDownload={() => downloadSvgAsPng(validTreemapRef, 'validation-issues.png')}
              ariaLabel="Validation issue treemap showing category-to-issue counts"
              caption="Hierarchical view of validation issues by category"
            >
              <div ref={validTreemapRef}>
                {statistics ? (
                  <ValidationTreemap
                    data={statistics.issue_summary}
                    selectedIndices={selectedIndices}
                    onSelectionChange={onSelectionChange}
                  />
                ) : (
                  <ChartSkeleton />
                )}
              </div>
            </ChartCard>

            {/* VIZ-02: Property Scatter (in Distributions tab) */}
            <ChartCard
              title="Property Scatter Plot"
              icon={<svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24"><circle cx="6" cy="18" r="2" /><circle cx="12" cy="8" r="2" /><circle cx="18" cy="14" r="2" /><circle cx="8" cy="12" r="2" /><circle cx="16" cy="6" r="2" /></svg>}
              onDownload={() => downloadSvgAsPng(propScatterRef, 'property-scatter.png')}
              ariaLabel="Property scatter plot showing relationship between two selected molecular properties"
              caption="Scatter plot of selected molecular properties"
            >
              <div ref={propScatterRef}>
                {results.length > 0 ? (
                  <PropertyScatterPlot
                    results={results}
                    xProperty={scatterXProp}
                    yProperty={scatterYProp}
                    colorByProperty={scatterColorProp}
                    selectedIndices={selectedIndices}
                    onSelectionChange={onSelectionChange}
                    onPropertyChange={handlePropertyChange}
                  />
                ) : (
                  <ChartSkeleton />
                )}
              </div>
            </ChartCard>
          </motion.div>
        )}

        {activeTab === 'Chemical Space' && (
          <motion.div
            key="chemical-space"
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -10 }}
            transition={{ duration: 0.2 }}
            className="grid grid-cols-1 lg:grid-cols-2 gap-4"
          >
            {/* VIZ-05: Scaffold Treemap */}
            <ChartCard
              title="Scaffold Diversity"
              icon={<svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24"><rect x="2" y="2" width="9" height="9" rx="1" /><rect x="13" y="2" width="9" height="5" rx="1" /><rect x="13" y="9" width="9" height="6" rx="1" /><rect x="2" y="13" width="9" height="9" rx="1" /><rect x="13" y="17" width="9" height="5" rx="1" /></svg>}
              onDownload={() => downloadSvgAsPng(scaffoldTreemapRef, 'scaffold-treemap.png')}
              ariaLabel="Scaffold treemap showing scaffold frequency distribution"
              caption="Top scaffolds in the batch by frequency"
            >
              <div ref={scaffoldTreemapRef}>
                {isChemSpaceStatus('scaffold')?.status === 'failed' ? (
                  <ErrorCard
                    message="Scaffold analysis failed"
                    onRetry={() => onRetrigger('scaffold')}
                  />
                ) : analyticsData?.scaffold ? (
                  <ScaffoldTreemap
                    data={analyticsData.scaffold}
                    selectedIndices={selectedIndices}
                    onSelectionChange={onSelectionChange}
                  />
                ) : (
                  <ChartSkeleton />
                )}
              </div>
            </ChartCard>

            {/* VIZ-06: Chemical Space Scatter */}
            <ChartCard
              title="Chemical Space"
              icon={<svg className="w-4 h-4" fill="currentColor" viewBox="0 0 24 24"><circle cx="6" cy="18" r="2" /><circle cx="12" cy="8" r="2" /><circle cx="18" cy="14" r="2" /><circle cx="8" cy="12" r="2" /><circle cx="16" cy="6" r="2" /></svg>}
              ariaLabel="Chemical space scatter plot showing PCA or t-SNE projection"
              caption="2D chemical space embedding of batch molecules"
            >
              {isChemSpaceStatus('chemical_space')?.status === 'failed' ? (
                <ErrorCard
                  message="Chemical space analysis failed"
                  onRetry={() => onRetrigger('chemical_space')}
                />
              ) : (
                <ChemicalSpaceScatter
                  data={analyticsData?.chemical_space ?? null}
                  results={results}
                  colorByProperty={chemSpaceColorProp}
                  onColorByChange={setChemSpaceColorProp}
                  selectedIndices={selectedIndices}
                  onSelectionChange={onSelectionChange}
                  jobId={jobId}
                />
              )}
            </ChartCard>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
});
