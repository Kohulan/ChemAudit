/**
 * BatchAnalyticsPanel
 *
 * Tabbed container for batch analytics charts.
 * Two tabs: "Distributions" and "Chemical Space" with skeleton loading,
 * error states, and summary badges.
 */

import React, { useState, useCallback, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Download, AlertTriangle, RotateCcw } from 'lucide-react';
import { ScoreHistogram } from './charts/ScoreHistogram';
import { PropertyScatterPlot } from './charts/PropertyScatterPlot';
import { AlertFrequencyChart } from './charts/AlertFrequencyChart';
import { ValidationTreemap } from './charts/ValidationTreemap';
import { ScaffoldTreemap } from './charts/ScaffoldTreemap';
import { ChemicalSpaceScatter } from './charts/ChemicalSpaceScatter';
import type { BatchStatistics, BatchResult } from '../../types/batch';
import { cn } from '../../lib/utils';

interface BatchAnalyticsPanelProps {
  jobId: string;
  statistics: BatchStatistics | null;
  results: BatchResult[];
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
  analyticsData: import('../../types/analytics').BatchAnalyticsResponse | null;
  onRetrigger: (type: string) => void;
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

export const BatchAnalyticsPanel = React.memo(function BatchAnalyticsPanel({
  jobId,
  statistics,
  results,
  selectedIndices,
  onSelectionChange,
  analyticsData,
  onRetrigger,
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

  // Badge counts
  const outlierCount = analyticsData?.statistics?.outliers?.length ?? 0;
  const scaffoldCount = analyticsData?.scaffold?.unique_scaffold_count ?? 0;

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
            {/* Summary badge */}
            {tab === 'Distributions' && outlierCount > 0 && (
              <span className="text-xs px-1.5 py-0.5 rounded-full bg-amber-500/10 text-amber-600 dark:text-amber-400">
                {outlierCount} outlier{outlierCount !== 1 ? 's' : ''}
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
