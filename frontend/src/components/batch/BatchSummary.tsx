import type { BatchStatistics } from '../../types/batch';

interface BatchSummaryProps {
  jobId: string;
  statistics: BatchStatistics;
  selectedIndices?: Set<number>;
  onClearSelection?: () => void;
  onFilterByStatus?: (status: 'success' | 'error' | null) => void;
  activeStatusFilter?: 'success' | 'error' | null;
  onAlertClick?: (catalog: string) => void;
  activeAlertFilter?: string | null;
  onIssueClick?: (checkName: string) => void;
  activeIssueFilter?: string | null;
}

/**
 * Summary statistics display with cards and charts.
 */
export function BatchSummary({ jobId: _jobId, statistics, selectedIndices: _selectedIndices, onClearSelection: _onClearSelection, onFilterByStatus, activeStatusFilter, onAlertClick, activeAlertFilter, onIssueClick, activeIssueFilter }: BatchSummaryProps) {
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

  const getPassRateColor = (rate: number | null): string => {
    if (rate === null) return 'text-[var(--color-text-muted)]';
    if (rate >= 80) return 'text-emerald-600 dark:text-emerald-400';
    if (rate >= 50) return 'text-amber-600 dark:text-amber-400';
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
        {/* Total (clears filter) */}
        <button
          type="button"
          onClick={() => onFilterByStatus?.(null)}
          className={`text-left bg-blue-500/10 dark:bg-blue-500/20 rounded-lg p-4 transition-all ${
            onFilterByStatus ? 'cursor-pointer hover:ring-2 hover:ring-blue-400/50' : ''
          } ${activeStatusFilter === null || activeStatusFilter === undefined ? '' : 'opacity-60'}`}
        >
          <p className="text-3xl font-bold text-blue-700 dark:text-blue-400">{statistics.total}</p>
          <p className="text-sm text-blue-600 dark:text-blue-500">Total Molecules</p>
        </button>

        {/* Successful */}
        <button
          type="button"
          onClick={() => onFilterByStatus?.(activeStatusFilter === 'success' ? null : 'success')}
          className={`text-left bg-yellow-500/10 dark:bg-yellow-500/20 rounded-lg p-4 transition-all ${
            onFilterByStatus ? 'cursor-pointer hover:ring-2 hover:ring-yellow-400/50' : ''
          } ${activeStatusFilter === 'success' ? 'ring-2 ring-yellow-500/70' : ''}`}
        >
          <p className="text-3xl font-bold text-amber-700 dark:text-yellow-400">
            {statistics.successful}
          </p>
          <p className="text-sm text-amber-600 dark:text-yellow-500">Successful</p>
        </button>

        {/* Errors */}
        <button
          type="button"
          onClick={() => onFilterByStatus?.(activeStatusFilter === 'error' ? null : 'error')}
          className={`text-left bg-red-500/10 dark:bg-red-500/20 rounded-lg p-4 transition-all ${
            onFilterByStatus ? 'cursor-pointer hover:ring-2 hover:ring-red-400/50' : ''
          } ${activeStatusFilter === 'error' ? 'ring-2 ring-red-500/70' : ''}`}
        >
          <p className="text-3xl font-bold text-red-700 dark:text-red-400">{statistics.errors}</p>
          <p className="text-sm text-red-600 dark:text-red-500">Errors</p>
        </button>

        {/* Processing time */}
        <div className="bg-purple-500/10 dark:bg-purple-500/20 rounded-lg p-4">
          <p className="text-3xl font-bold text-purple-700 dark:text-purple-400">
            {formatTime(statistics.processing_time_seconds)}
          </p>
          <p className="text-sm text-purple-600 dark:text-purple-500">Processing Time</p>
        </div>
      </div>

      {/* Score averages */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
          <p className="text-sm text-[var(--color-text-muted)] mb-1">Avg Validation</p>
          <p
            className={`text-3xl font-bold ${getScoreColor(
              statistics.avg_validation_score
            )}`}
          >
            {statistics.avg_validation_score?.toFixed(1) ?? '-'}
          </p>
        </div>
        <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
          <p className="text-sm text-[var(--color-text-muted)] mb-1">Avg ML-Readiness</p>
          <p
            className={`text-3xl font-bold ${getScoreColor(
              statistics.avg_ml_readiness_score
            )}`}
          >
            {statistics.avg_ml_readiness_score?.toFixed(1) ?? '-'}
          </p>
        </div>
        <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
          <p className="text-sm text-[var(--color-text-muted)] mb-1">Avg QED Score</p>
          <p className="text-3xl font-bold text-purple-600 dark:text-purple-400">
            {statistics.avg_qed_score?.toFixed(2) ?? '-'}
          </p>
        </div>
        <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
          <p className="text-sm text-[var(--color-text-muted)] mb-1">Avg SA Score</p>
          <p className="text-3xl font-bold text-cyan-600 dark:text-cyan-400">
            {statistics.avg_sa_score?.toFixed(1) ?? '-'}
          </p>
        </div>
      </div>

      {/* Pass rates */}
      {(statistics.lipinski_pass_rate !== null || statistics.safety_pass_rate !== null) && (
        <div className="grid grid-cols-2 gap-4">
          {statistics.lipinski_pass_rate !== null && (
            <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
              <p className="text-sm text-[var(--color-text-muted)] mb-1">Lipinski Pass Rate</p>
              <p className={`text-3xl font-bold ${getPassRateColor(statistics.lipinski_pass_rate)}`}>
                {statistics.lipinski_pass_rate?.toFixed(0) ?? '-'}%
              </p>
            </div>
          )}
          {statistics.safety_pass_rate !== null && (
            <div className="bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg p-4">
              <p className="text-sm text-[var(--color-text-muted)] mb-1">Safety Pass Rate</p>
              <p className={`text-3xl font-bold ${getPassRateColor(statistics.safety_pass_rate)}`}>
                {statistics.safety_pass_rate?.toFixed(0) ?? '-'}%
              </p>
            </div>
          )}
        </div>
      )}

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
        <div className="bg-amber-500/10 dark:bg-amber-500/20 border border-amber-500/30 rounded-xl p-5">
          <h4 className="font-medium text-amber-800 dark:text-amber-400 mb-4">Alert Summary</h4>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
            {Object.entries(statistics.alert_summary).map(([catalog, count]) => {
              const isActive = activeAlertFilter === catalog;
              return (
                <button
                  key={catalog}
                  type="button"
                  onClick={() => onAlertClick?.(isActive ? '' : catalog)}
                  className={`group relative text-left rounded-xl p-4 border transition-all duration-200 ${
                    onAlertClick ? 'cursor-pointer' : ''
                  } ${isActive
                    ? 'bg-amber-500/25 border-amber-500/60 shadow-md shadow-amber-500/10 scale-[1.02]'
                    : 'bg-[var(--color-surface-elevated)] border-amber-500/20 hover:border-amber-500/40 hover:shadow-sm hover:scale-[1.01]'
                  }`}
                >
                  <p className="text-3xl font-bold text-amber-700 dark:text-amber-400">{count}</p>
                  <p className="text-xs font-medium text-amber-600 dark:text-amber-500 mt-1">{catalog}</p>
                  {onAlertClick && (
                    <div className={`absolute top-2.5 right-2.5 w-1.5 h-1.5 rounded-full transition-colors ${
                      isActive ? 'bg-amber-500' : 'bg-amber-400/0 group-hover:bg-amber-400/50'
                    }`} />
                  )}
                </button>
              );
            })}
          </div>
        </div>
      )}

      {/* Issue summary */}
      {Object.keys(statistics.issue_summary).length > 0 && (
        <div className="bg-red-500/10 dark:bg-red-500/20 border border-red-500/30 rounded-xl p-5">
          <h4 className="font-medium text-red-800 dark:text-red-400 mb-4">Structural Issues</h4>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
            {Object.entries(statistics.issue_summary)
              .sort(([, a], [, b]) => b - a)
              .map(([checkName, count]) => {
              const isActive = activeIssueFilter === checkName;
              return (
                <button
                  key={checkName}
                  type="button"
                  onClick={() => onIssueClick?.(isActive ? '' : checkName)}
                  className={`group relative text-left rounded-xl p-4 border transition-all duration-200 ${
                    onIssueClick ? 'cursor-pointer' : ''
                  } ${isActive
                    ? 'bg-red-500/25 border-red-500/60 shadow-md shadow-red-500/10 scale-[1.02]'
                    : 'bg-[var(--color-surface-elevated)] border-red-500/20 hover:border-red-500/40 hover:shadow-sm hover:scale-[1.01]'
                  }`}
                >
                  <p className="text-3xl font-bold text-red-700 dark:text-red-400">{count}</p>
                  <p className="text-xs font-medium text-red-600 dark:text-red-500 mt-1">
                    {checkName.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase())}
                  </p>
                  {onIssueClick && (
                    <div className={`absolute top-2.5 right-2.5 w-1.5 h-1.5 rounded-full transition-colors ${
                      isActive ? 'bg-red-500' : 'bg-red-400/0 group-hover:bg-red-400/50'
                    }`} />
                  )}
                </button>
              );
            })}
          </div>
        </div>
      )}

    </div>
  );
}
