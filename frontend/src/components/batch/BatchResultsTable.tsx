import { Fragment, useState } from 'react';
import type { BatchResult, BatchResultsFilters } from '../../types/batch';

interface BatchResultsTableProps {
  results: BatchResult[];
  page: number;
  pageSize: number;
  totalResults: number;
  totalPages: number;
  filters: BatchResultsFilters;
  onPageChange: (page: number) => void;
  onPageSizeChange: (size: number) => void;
  onFiltersChange: (filters: BatchResultsFilters) => void;
  isLoading?: boolean;
}

/**
 * Paginated results table with filtering and sorting.
 */
export function BatchResultsTable({
  results,
  page,
  pageSize,
  totalResults,
  totalPages,
  filters,
  onPageChange,
  onPageSizeChange,
  onFiltersChange,
  isLoading = false,
}: BatchResultsTableProps) {
  const [expandedRow, setExpandedRow] = useState<number | null>(null);
  const [sortField, setSortField] = useState<string>('index');
  const [sortDir, setSortDir] = useState<'asc' | 'desc'>('asc');

  const handleSort = (field: string) => {
    if (sortField === field) {
      setSortDir(sortDir === 'asc' ? 'desc' : 'asc');
    } else {
      setSortField(field);
      setSortDir('asc');
    }
  };

  const getScoreColor = (score: number | null | undefined): string => {
    if (score === null || score === undefined) return 'text-[var(--color-text-muted)]';
    if (score >= 80) return 'text-amber-600 dark:text-yellow-400 bg-yellow-500/10';
    if (score >= 50) return 'text-orange-600 dark:text-orange-400 bg-orange-500/10';
    return 'text-red-600 dark:text-red-400 bg-red-500/10';
  };

  const truncateSmiles = (smiles: string, maxLen: number = 30): string => {
    if (smiles.length <= maxLen) return smiles;
    return smiles.substring(0, maxLen) + '...';
  };

  // Sort results client-side
  const sortedResults = [...results].sort((a, b) => {
    let aVal: any = a[sortField as keyof BatchResult];
    let bVal: any = b[sortField as keyof BatchResult];

    // Handle nested score values
    if (sortField === 'score') {
      aVal = a.validation?.overall_score ?? -1;
      bVal = b.validation?.overall_score ?? -1;
    }

    // Handle nested alerts count
    if (sortField === 'alerts') {
      aVal = a.alerts?.alert_count ?? 0;
      bVal = b.alerts?.alert_count ?? 0;
    }

    if (typeof aVal === 'string') {
      aVal = aVal.toLowerCase();
      bVal = String(bVal).toLowerCase();
    }

    if (aVal < bVal) return sortDir === 'asc' ? -1 : 1;
    if (aVal > bVal) return sortDir === 'asc' ? 1 : -1;
    return 0;
  });

  return (
    <div className="space-y-4">
      {/* Filters */}
      <div className="flex flex-wrap items-center gap-4 p-4 bg-[var(--color-surface-sunken)] rounded-lg">
        {/* Status filter */}
        <div className="flex items-center space-x-2">
          <label className="text-sm text-[var(--color-text-secondary)]">Status:</label>
          <select
            value={filters.status_filter || ''}
            onChange={(e) =>
              onFiltersChange({
                ...filters,
                status_filter: e.target.value as 'success' | 'error' | undefined || undefined,
              })
            }
            className="border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded px-2 py-1 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
          >
            <option value="">All</option>
            <option value="success">Success</option>
            <option value="error">Error</option>
          </select>
        </div>

        {/* Score range filter */}
        <div className="flex items-center space-x-2">
          <label className="text-sm text-[var(--color-text-secondary)]">Score:</label>
          <input
            type="number"
            placeholder="Min"
            min="0"
            max="100"
            value={filters.min_score ?? ''}
            onChange={(e) =>
              onFiltersChange({
                ...filters,
                min_score: e.target.value ? Number(e.target.value) : undefined,
              })
            }
            className="w-16 border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded px-2 py-1 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
          />
          <span className="text-[var(--color-text-muted)]">-</span>
          <input
            type="number"
            placeholder="Max"
            min="0"
            max="100"
            value={filters.max_score ?? ''}
            onChange={(e) =>
              onFiltersChange({
                ...filters,
                max_score: e.target.value ? Number(e.target.value) : undefined,
              })
            }
            className="w-16 border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded px-2 py-1 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
          />
        </div>

        {/* Results info */}
        <div className="ml-auto text-sm text-[var(--color-text-muted)]">
          Showing {results.length} of {totalResults} results
        </div>
      </div>

      {/* Table */}
      <div className="overflow-x-auto">
        <table className="w-full border-collapse">
          <thead>
            <tr className="bg-[var(--color-surface-sunken)]">
              <th
                className="px-4 py-3 text-left text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('index')}
              >
                # {sortField === 'index' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-left text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('name')}
              >
                Name {sortField === 'name' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th className="px-4 py-3 text-left text-xs font-medium text-[var(--color-text-muted)] uppercase">
                SMILES
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('score')}
              >
                Score {sortField === 'score' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                QED
              </th>
              <th className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                Safety
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('status')}
              >
                Status {sortField === 'status' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('alerts')}
              >
                Alerts {sortField === 'alerts' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-[var(--color-border)]">
            {isLoading ? (
              <tr>
                <td colSpan={8} className="px-4 py-8 text-center text-[var(--color-text-muted)]">
                  <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-[var(--color-primary)] mx-auto mb-2" />
                  Loading results...
                </td>
              </tr>
            ) : sortedResults.length === 0 ? (
              <tr>
                <td colSpan={8} className="px-4 py-8 text-center text-[var(--color-text-muted)]">
                  No results match the current filters.
                </td>
              </tr>
            ) : (
              sortedResults.map((result) => (
                <Fragment key={result.index}>
                  <tr
                    className={`
                      hover:bg-[var(--color-surface-sunken)] cursor-pointer transition-colors
                      ${result.status === 'error' ? 'bg-red-500/5' : ''}
                      ${expandedRow === result.index ? 'bg-[var(--color-primary)]/5' : ''}
                    `}
                    onClick={() =>
                      setExpandedRow(expandedRow === result.index ? null : result.index)
                    }
                  >
                    <td className="px-4 py-3 text-sm text-[var(--color-text-muted)]">
                      {result.index + 1}
                    </td>
                    <td className="px-4 py-3 text-sm font-medium text-[var(--color-text-primary)]">
                      {result.name || '-'}
                    </td>
                    <td className="px-4 py-3 text-sm text-[var(--color-text-secondary)] font-mono">
                      {truncateSmiles(result.smiles)}
                    </td>
                    <td className="px-4 py-3 text-center">
                      <span
                        className={`px-2 py-1 rounded text-sm font-medium ${getScoreColor(
                          result.validation?.overall_score
                        )}`}
                      >
                        {result.validation?.overall_score ?? '-'}
                      </span>
                    </td>
                    <td className="px-4 py-3 text-center">
                      {result.scoring?.druglikeness ? (
                        <span className="px-2 py-1 rounded text-xs font-medium bg-purple-500/10 text-purple-600 dark:text-purple-400">
                          {result.scoring.druglikeness.qed_score.toFixed(2)}
                        </span>
                      ) : (
                        <span className="text-[var(--color-text-muted)]">-</span>
                      )}
                    </td>
                    <td className="px-4 py-3 text-center">
                      {result.scoring?.safety_filters ? (
                        <span className={`px-2 py-1 rounded text-xs font-medium ${
                          result.scoring.safety_filters.all_passed
                            ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                            : 'bg-red-500/10 text-red-600 dark:text-red-400'
                        }`}>
                          {result.scoring.safety_filters.all_passed ? '✓' : result.scoring.safety_filters.total_alerts}
                        </span>
                      ) : (
                        <span className="text-[var(--color-text-muted)]">-</span>
                      )}
                    </td>
                    <td className="px-4 py-3 text-center">
                      <span
                        className={`px-2 py-1 rounded text-xs font-medium ${
                          result.status === 'success'
                            ? 'bg-yellow-500/10 text-amber-700 dark:text-yellow-400'
                            : 'bg-red-500/10 text-red-700 dark:text-red-400'
                        }`}
                      >
                        {result.status}
                      </span>
                    </td>
                    <td className="px-4 py-3 text-center">
                      {result.alerts?.alert_count ? (
                        <span className="px-2 py-1 rounded text-xs font-medium bg-amber-500/10 text-amber-700 dark:text-amber-400">
                          {result.alerts.alert_count}
                        </span>
                      ) : (
                        <span className="text-[var(--color-text-muted)]">0</span>
                      )}
                    </td>
                  </tr>

                  {/* Expanded details */}
                  {expandedRow === result.index && (
                    <tr>
                      <td colSpan={8} className="px-4 py-4 bg-[var(--color-surface-sunken)]">
                        <div className="space-y-3">
                          {/* Full SMILES */}
                          <div>
                            <p className="text-xs text-[var(--color-text-muted)] mb-1">Full SMILES:</p>
                            <p className="text-sm font-mono bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] p-2 rounded border border-[var(--color-border)] overflow-x-auto">
                              {result.smiles}
                            </p>
                          </div>

                          {/* Error */}
                          {result.error && (
                            <div className="bg-red-500/10 dark:bg-red-500/20 border border-red-500/30 rounded p-2">
                              <p className="text-sm text-red-600 dark:text-red-400">{result.error}</p>
                            </div>
                          )}

                          {/* Validation issues */}
                          {result.validation?.issues && result.validation.issues.length > 0 && (
                            <div>
                              <p className="text-xs text-[var(--color-text-muted)] mb-1">
                                Validation Issues:
                              </p>
                              <ul className="text-sm space-y-1 text-[var(--color-text-secondary)]">
                                {result.validation.issues.map((issue, i) => (
                                  <li key={i} className="flex items-start">
                                    <span
                                      className={`w-2 h-2 mt-1.5 mr-2 rounded-full ${
                                        issue.severity === 'error'
                                          ? 'bg-red-500'
                                          : 'bg-amber-500'
                                      }`}
                                    />
                                    <span>{issue.message}</span>
                                  </li>
                                ))}
                              </ul>
                            </div>
                          )}

                          {/* Alerts */}
                          {result.alerts?.alerts && result.alerts.alerts.length > 0 && (
                            <div>
                              <p className="text-xs text-[var(--color-text-muted)] mb-1">Alerts:</p>
                              <ul className="text-sm space-y-1">
                                {result.alerts.alerts.map((alert, i) => (
                                  <li key={i} className="text-amber-700 dark:text-amber-400">
                                    {alert.catalog}: {alert.rule_name}
                                  </li>
                                ))}
                              </ul>
                            </div>
                          )}

                          {/* Scoring Section */}
                          {result.scoring && (
                            <div className="space-y-2">
                              <p className="text-xs text-[var(--color-text-muted)] font-medium">Scoring:</p>
                              <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
                                {/* ML-Readiness */}
                                {result.scoring.ml_readiness && (
                                  <div className="bg-[var(--color-surface-elevated)] rounded-lg p-2 border border-[var(--color-border)]">
                                    <p className="text-xs text-[var(--color-text-muted)]">ML-Readiness</p>
                                    <p className="text-lg font-bold text-[var(--color-primary)]">
                                      {result.scoring.ml_readiness.score}
                                    </p>
                                  </div>
                                )}

                                {/* Drug-likeness */}
                                {result.scoring.druglikeness && (
                                  <div className="bg-[var(--color-surface-elevated)] rounded-lg p-2 border border-[var(--color-border)]">
                                    <p className="text-xs text-[var(--color-text-muted)]">QED Score</p>
                                    <p className="text-lg font-bold text-purple-500">
                                      {result.scoring.druglikeness.qed_score.toFixed(2)}
                                    </p>
                                    <p className={`text-xs ${result.scoring.druglikeness.lipinski_passed ? 'text-emerald-500' : 'text-red-500'}`}>
                                      Lipinski: {result.scoring.druglikeness.lipinski_passed ? 'Pass' : 'Fail'}
                                    </p>
                                  </div>
                                )}

                                {/* Safety */}
                                {result.scoring.safety_filters && (
                                  <div className="bg-[var(--color-surface-elevated)] rounded-lg p-2 border border-[var(--color-border)]">
                                    <p className="text-xs text-[var(--color-text-muted)]">Safety</p>
                                    <p className={`text-lg font-bold ${result.scoring.safety_filters.all_passed ? 'text-emerald-500' : 'text-red-500'}`}>
                                      {result.scoring.safety_filters.all_passed ? 'Clear' : `${result.scoring.safety_filters.total_alerts} alerts`}
                                    </p>
                                    <p className="text-xs text-[var(--color-text-muted)]">
                                      PAINS: {result.scoring.safety_filters.pains_passed ? '✓' : '✗'} | Brenk: {result.scoring.safety_filters.brenk_passed ? '✓' : '✗'}
                                    </p>
                                  </div>
                                )}

                                {/* ADMET */}
                                {result.scoring.admet && (
                                  <div className="bg-[var(--color-surface-elevated)] rounded-lg p-2 border border-[var(--color-border)]">
                                    <p className="text-xs text-[var(--color-text-muted)]">ADMET</p>
                                    <p className="text-lg font-bold text-cyan-500">
                                      SA: {result.scoring.admet.sa_score.toFixed(1)}
                                    </p>
                                    <p className="text-xs text-[var(--color-text-muted)]">
                                      {result.scoring.admet.sa_classification} | Fsp3: {result.scoring.admet.fsp3.toFixed(2)}
                                    </p>
                                  </div>
                                )}
                              </div>

                              {/* Interpretation */}
                              {result.scoring.ml_readiness?.interpretation && (
                                <p className="text-sm text-[var(--color-text-secondary)] mt-2">
                                  {result.scoring.ml_readiness.interpretation}
                                </p>
                              )}
                            </div>
                          )}
                        </div>
                      </td>
                    </tr>
                  )}
                </Fragment>
              ))
            )}
          </tbody>
        </table>
      </div>

      {/* Pagination */}
      <div className="flex items-center justify-between px-4 py-3 bg-[var(--color-surface-sunken)] rounded-lg">
        <div className="flex items-center space-x-2">
          <label className="text-sm text-[var(--color-text-secondary)]">Per page:</label>
          <select
            value={pageSize}
            onChange={(e) => onPageSizeChange(Number(e.target.value))}
            className="border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] rounded px-2 py-1 text-sm focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
          >
            <option value={25}>25</option>
            <option value={50}>50</option>
            <option value={100}>100</option>
          </select>
        </div>

        <div className="flex items-center space-x-2">
          <button
            onClick={() => onPageChange(page - 1)}
            disabled={page <= 1}
            className="px-3 py-1 border border-[var(--color-border)] rounded text-sm text-[var(--color-text-secondary)] disabled:opacity-50 disabled:cursor-not-allowed hover:bg-[var(--color-surface-elevated)] transition-colors"
          >
            Previous
          </button>
          <span className="text-sm text-[var(--color-text-secondary)]">
            Page {page} of {totalPages}
          </span>
          <button
            onClick={() => onPageChange(page + 1)}
            disabled={page >= totalPages}
            className="px-3 py-1 border border-[var(--color-border)] rounded text-sm text-[var(--color-text-secondary)] disabled:opacity-50 disabled:cursor-not-allowed hover:bg-[var(--color-surface-elevated)] transition-colors"
          >
            Next
          </button>
        </div>
      </div>
    </div>
  );
}
