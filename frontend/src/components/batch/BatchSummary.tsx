import { useState } from 'react';
import type { BatchStatistics } from '../../types/batch';
import { ExportDialog } from './ExportDialog';

interface BatchSummaryProps {
  jobId: string;
  statistics: BatchStatistics;
}

/**
 * Summary statistics display with cards and charts.
 */
export function BatchSummary({ jobId, statistics }: BatchSummaryProps) {
  const [isExportDialogOpen, setIsExportDialogOpen] = useState(false);
  const formatTime = (seconds: number | null): string => {
    if (seconds === null) return '-';
    if (seconds < 60) return `${seconds.toFixed(1)}s`;
    const mins = Math.floor(seconds / 60);
    const secs = Math.round(seconds % 60);
    return `${mins}m ${secs}s`;
  };

  const getScoreColor = (score: number | null): string => {
    if (score === null) return 'text-[var(--color-text-muted)]';
    if (score >= 80) return 'text-amber-600 dark:text-yellow-400';
    if (score >= 50) return 'text-orange-600 dark:text-orange-400';
    return 'text-red-600 dark:text-red-400';
  };

  // Calculate percentages for distribution chart
  const total = statistics.successful;
  const distributionPct = {
    excellent: total > 0 ? (statistics.score_distribution.excellent / total) * 100 : 0,
    good: total > 0 ? (statistics.score_distribution.good / total) * 100 : 0,
    moderate: total > 0 ? (statistics.score_distribution.moderate / total) * 100 : 0,
    poor: total > 0 ? (statistics.score_distribution.poor / total) * 100 : 0,
  };

  return (
    <div className="space-y-6">
      {/* Main stats cards */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        {/* Total */}
        <div className="bg-blue-500/10 dark:bg-blue-500/20 rounded-lg p-4">
          <p className="text-3xl font-bold text-blue-700 dark:text-blue-400">{statistics.total}</p>
          <p className="text-sm text-blue-600 dark:text-blue-500">Total Molecules</p>
        </div>

        {/* Successful */}
        <div className="bg-yellow-500/10 dark:bg-yellow-500/20 rounded-lg p-4">
          <p className="text-3xl font-bold text-amber-700 dark:text-yellow-400">
            {statistics.successful}
          </p>
          <p className="text-sm text-amber-600 dark:text-yellow-500">Successful</p>
        </div>

        {/* Errors */}
        <div className="bg-red-500/10 dark:bg-red-500/20 rounded-lg p-4">
          <p className="text-3xl font-bold text-red-700 dark:text-red-400">{statistics.errors}</p>
          <p className="text-sm text-red-600 dark:text-red-500">Errors</p>
        </div>

        {/* Processing time */}
        <div className="bg-purple-500/10 dark:bg-purple-500/20 rounded-lg p-4">
          <p className="text-3xl font-bold text-purple-700 dark:text-purple-400">
            {formatTime(statistics.processing_time_seconds)}
          </p>
          <p className="text-sm text-purple-600 dark:text-purple-500">Processing Time</p>
        </div>
      </div>

      {/* Score averages */}
      <div className="grid grid-cols-2 gap-4">
        <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
          <p className="text-sm text-[var(--color-text-muted)] mb-1">Avg Validation Score</p>
          <p
            className={`text-4xl font-bold ${getScoreColor(
              statistics.avg_validation_score
            )}`}
          >
            {statistics.avg_validation_score?.toFixed(1) ?? '-'}
          </p>
        </div>
        <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
          <p className="text-sm text-[var(--color-text-muted)] mb-1">Avg ML-Readiness Score</p>
          <p
            className={`text-4xl font-bold ${getScoreColor(
              statistics.avg_ml_readiness_score
            )}`}
          >
            {statistics.avg_ml_readiness_score?.toFixed(1) ?? '-'}
          </p>
        </div>
      </div>

      {/* Score distribution chart */}
      <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
        <h4 className="font-medium text-[var(--color-text-primary)] mb-4">Score Distribution</h4>
        <div className="space-y-3">
          {/* Excellent */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-[var(--color-text-secondary)]">Excellent</span>
            <div className="flex-1 h-6 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-yellow-500 transition-all duration-500"
                style={{ width: `${distributionPct.excellent}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-[var(--color-text-secondary)]">
              {statistics.score_distribution.excellent}
            </span>
          </div>

          {/* Good */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-[var(--color-text-secondary)]">Good</span>
            <div className="flex-1 h-6 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-amber-500 transition-all duration-500"
                style={{ width: `${distributionPct.good}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-[var(--color-text-secondary)]">
              {statistics.score_distribution.good}
            </span>
          </div>

          {/* Moderate */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-[var(--color-text-secondary)]">Moderate</span>
            <div className="flex-1 h-6 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-amber-500 transition-all duration-500"
                style={{ width: `${distributionPct.moderate}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-[var(--color-text-secondary)]">
              {statistics.score_distribution.moderate}
            </span>
          </div>

          {/* Poor */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-[var(--color-text-secondary)]">Poor</span>
            <div className="flex-1 h-6 bg-[var(--color-surface-sunken)] rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-red-500 transition-all duration-500"
                style={{ width: `${distributionPct.poor}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-[var(--color-text-secondary)]">
              {statistics.score_distribution.poor}
            </span>
          </div>
        </div>

        <div className="mt-3 text-xs text-[var(--color-text-muted)] flex justify-between">
          <span>90-100</span>
          <span>70-89</span>
          <span>50-69</span>
          <span>0-49</span>
        </div>
      </div>

      {/* Alert summary */}
      {Object.keys(statistics.alert_summary).length > 0 && (
        <div className="bg-amber-500/10 dark:bg-amber-500/20 border border-amber-500/30 rounded-lg p-4">
          <h4 className="font-medium text-amber-800 dark:text-amber-400 mb-3">Alert Summary</h4>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
            {Object.entries(statistics.alert_summary).map(([catalog, count]) => (
              <div key={catalog} className="text-center">
                <p className="text-2xl font-semibold text-amber-700 dark:text-amber-400">{count}</p>
                <p className="text-xs text-amber-600 dark:text-amber-500">{catalog}</p>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Export button */}
      <div className="flex justify-end">
        <button
          onClick={() => setIsExportDialogOpen(true)}
          className="px-4 py-2 text-sm font-medium text-white bg-[var(--color-primary)] rounded-lg hover:bg-[var(--color-primary-dark)] transition-colors flex items-center"
        >
          <svg
            className="w-4 h-4 mr-2"
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"
            />
          </svg>
          Export Results
        </button>
      </div>

      {/* Export Dialog */}
      <ExportDialog
        jobId={jobId}
        isOpen={isExportDialogOpen}
        onClose={() => setIsExportDialogOpen(false)}
      />
    </div>
  );
}
