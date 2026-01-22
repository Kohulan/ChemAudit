import type { BatchStatistics } from '../../types/batch';

interface BatchSummaryProps {
  statistics: BatchStatistics;
}

/**
 * Summary statistics display with cards and charts.
 */
export function BatchSummary({ statistics }: BatchSummaryProps) {
  const formatTime = (seconds: number | null): string => {
    if (seconds === null) return '-';
    if (seconds < 60) return `${seconds.toFixed(1)}s`;
    const mins = Math.floor(seconds / 60);
    const secs = Math.round(seconds % 60);
    return `${mins}m ${secs}s`;
  };

  const getScoreColor = (score: number | null): string => {
    if (score === null) return 'text-gray-500';
    if (score >= 80) return 'text-green-600';
    if (score >= 50) return 'text-yellow-600';
    return 'text-red-600';
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
        <div className="bg-blue-50 rounded-lg p-4">
          <p className="text-3xl font-bold text-blue-700">{statistics.total}</p>
          <p className="text-sm text-blue-600">Total Molecules</p>
        </div>

        {/* Successful */}
        <div className="bg-green-50 rounded-lg p-4">
          <p className="text-3xl font-bold text-green-700">
            {statistics.successful}
          </p>
          <p className="text-sm text-green-600">Successful</p>
        </div>

        {/* Errors */}
        <div className="bg-red-50 rounded-lg p-4">
          <p className="text-3xl font-bold text-red-700">{statistics.errors}</p>
          <p className="text-sm text-red-600">Errors</p>
        </div>

        {/* Processing time */}
        <div className="bg-purple-50 rounded-lg p-4">
          <p className="text-3xl font-bold text-purple-700">
            {formatTime(statistics.processing_time_seconds)}
          </p>
          <p className="text-sm text-purple-600">Processing Time</p>
        </div>
      </div>

      {/* Score averages */}
      <div className="grid grid-cols-2 gap-4">
        <div className="bg-white border border-gray-200 rounded-lg p-4">
          <p className="text-sm text-gray-500 mb-1">Avg Validation Score</p>
          <p
            className={`text-4xl font-bold ${getScoreColor(
              statistics.avg_validation_score
            )}`}
          >
            {statistics.avg_validation_score?.toFixed(1) ?? '-'}
          </p>
        </div>
        <div className="bg-white border border-gray-200 rounded-lg p-4">
          <p className="text-sm text-gray-500 mb-1">Avg ML-Readiness Score</p>
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
      <div className="bg-white border border-gray-200 rounded-lg p-4">
        <h4 className="font-medium text-gray-900 mb-4">Score Distribution</h4>
        <div className="space-y-3">
          {/* Excellent */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-gray-600">Excellent</span>
            <div className="flex-1 h-6 bg-gray-100 rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-green-500 transition-all duration-500"
                style={{ width: `${distributionPct.excellent}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-gray-700">
              {statistics.score_distribution.excellent}
            </span>
          </div>

          {/* Good */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-gray-600">Good</span>
            <div className="flex-1 h-6 bg-gray-100 rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-blue-500 transition-all duration-500"
                style={{ width: `${distributionPct.good}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-gray-700">
              {statistics.score_distribution.good}
            </span>
          </div>

          {/* Moderate */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-gray-600">Moderate</span>
            <div className="flex-1 h-6 bg-gray-100 rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-yellow-500 transition-all duration-500"
                style={{ width: `${distributionPct.moderate}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-gray-700">
              {statistics.score_distribution.moderate}
            </span>
          </div>

          {/* Poor */}
          <div className="flex items-center">
            <span className="w-20 text-sm text-gray-600">Poor</span>
            <div className="flex-1 h-6 bg-gray-100 rounded-full overflow-hidden mx-3">
              <div
                className="h-full bg-red-500 transition-all duration-500"
                style={{ width: `${distributionPct.poor}%` }}
              />
            </div>
            <span className="w-12 text-sm text-right text-gray-700">
              {statistics.score_distribution.poor}
            </span>
          </div>
        </div>

        <div className="mt-3 text-xs text-gray-500 flex justify-between">
          <span>90-100</span>
          <span>70-89</span>
          <span>50-69</span>
          <span>0-49</span>
        </div>
      </div>

      {/* Alert summary */}
      {Object.keys(statistics.alert_summary).length > 0 && (
        <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
          <h4 className="font-medium text-amber-800 mb-3">Alert Summary</h4>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
            {Object.entries(statistics.alert_summary).map(([catalog, count]) => (
              <div key={catalog} className="text-center">
                <p className="text-2xl font-semibold text-amber-700">{count}</p>
                <p className="text-xs text-amber-600">{catalog}</p>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Export placeholder */}
      <div className="flex justify-end">
        <button
          disabled
          className="px-4 py-2 text-sm text-gray-400 bg-gray-100 rounded-lg cursor-not-allowed"
          title="Coming in Phase 3"
        >
          Export Results (Coming Soon)
        </button>
      </div>
    </div>
  );
}
