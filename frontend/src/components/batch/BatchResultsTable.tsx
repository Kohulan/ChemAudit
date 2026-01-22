import { useState } from 'react';
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
    if (score === null || score === undefined) return 'text-gray-400';
    if (score >= 80) return 'text-green-600 bg-green-50';
    if (score >= 50) return 'text-yellow-600 bg-yellow-50';
    return 'text-red-600 bg-red-50';
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
      <div className="flex flex-wrap items-center gap-4 p-4 bg-gray-50 rounded-lg">
        {/* Status filter */}
        <div className="flex items-center space-x-2">
          <label className="text-sm text-gray-600">Status:</label>
          <select
            value={filters.status_filter || ''}
            onChange={(e) =>
              onFiltersChange({
                ...filters,
                status_filter: e.target.value as 'success' | 'error' | undefined || undefined,
              })
            }
            className="border border-gray-300 rounded px-2 py-1 text-sm"
          >
            <option value="">All</option>
            <option value="success">Success</option>
            <option value="error">Error</option>
          </select>
        </div>

        {/* Score range filter */}
        <div className="flex items-center space-x-2">
          <label className="text-sm text-gray-600">Score:</label>
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
            className="w-16 border border-gray-300 rounded px-2 py-1 text-sm"
          />
          <span className="text-gray-400">-</span>
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
            className="w-16 border border-gray-300 rounded px-2 py-1 text-sm"
          />
        </div>

        {/* Results info */}
        <div className="ml-auto text-sm text-gray-500">
          Showing {results.length} of {totalResults} results
        </div>
      </div>

      {/* Table */}
      <div className="overflow-x-auto">
        <table className="w-full border-collapse">
          <thead>
            <tr className="bg-gray-100">
              <th
                className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase cursor-pointer hover:bg-gray-200"
                onClick={() => handleSort('index')}
              >
                # {sortField === 'index' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase cursor-pointer hover:bg-gray-200"
                onClick={() => handleSort('name')}
              >
                Name {sortField === 'name' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                SMILES
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-gray-500 uppercase cursor-pointer hover:bg-gray-200"
                onClick={() => handleSort('score')}
              >
                Score {sortField === 'score' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-gray-500 uppercase cursor-pointer hover:bg-gray-200"
                onClick={() => handleSort('status')}
              >
                Status {sortField === 'status' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th className="px-4 py-3 text-center text-xs font-medium text-gray-500 uppercase">
                Alerts
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-gray-200">
            {isLoading ? (
              <tr>
                <td colSpan={6} className="px-4 py-8 text-center text-gray-500">
                  <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 mx-auto mb-2" />
                  Loading results...
                </td>
              </tr>
            ) : sortedResults.length === 0 ? (
              <tr>
                <td colSpan={6} className="px-4 py-8 text-center text-gray-500">
                  No results match the current filters.
                </td>
              </tr>
            ) : (
              sortedResults.map((result) => (
                <>
                  <tr
                    key={result.index}
                    className={`
                      hover:bg-gray-50 cursor-pointer
                      ${result.status === 'error' ? 'bg-red-50' : ''}
                      ${expandedRow === result.index ? 'bg-blue-50' : ''}
                    `}
                    onClick={() =>
                      setExpandedRow(expandedRow === result.index ? null : result.index)
                    }
                  >
                    <td className="px-4 py-3 text-sm text-gray-500">
                      {result.index + 1}
                    </td>
                    <td className="px-4 py-3 text-sm font-medium text-gray-900">
                      {result.name || '-'}
                    </td>
                    <td className="px-4 py-3 text-sm text-gray-600 font-mono">
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
                      <span
                        className={`px-2 py-1 rounded text-xs font-medium ${
                          result.status === 'success'
                            ? 'bg-green-100 text-green-800'
                            : 'bg-red-100 text-red-800'
                        }`}
                      >
                        {result.status}
                      </span>
                    </td>
                    <td className="px-4 py-3 text-center">
                      {result.alerts?.alert_count ? (
                        <span className="px-2 py-1 rounded text-xs font-medium bg-amber-100 text-amber-800">
                          {result.alerts.alert_count}
                        </span>
                      ) : (
                        <span className="text-gray-400">0</span>
                      )}
                    </td>
                  </tr>

                  {/* Expanded details */}
                  {expandedRow === result.index && (
                    <tr key={`${result.index}-detail`}>
                      <td colSpan={6} className="px-4 py-4 bg-gray-50">
                        <div className="space-y-3">
                          {/* Full SMILES */}
                          <div>
                            <p className="text-xs text-gray-500 mb-1">Full SMILES:</p>
                            <p className="text-sm font-mono bg-white p-2 rounded border overflow-x-auto">
                              {result.smiles}
                            </p>
                          </div>

                          {/* Error */}
                          {result.error && (
                            <div className="bg-red-50 border border-red-200 rounded p-2">
                              <p className="text-sm text-red-700">{result.error}</p>
                            </div>
                          )}

                          {/* Validation issues */}
                          {result.validation?.issues && result.validation.issues.length > 0 && (
                            <div>
                              <p className="text-xs text-gray-500 mb-1">
                                Validation Issues:
                              </p>
                              <ul className="text-sm space-y-1">
                                {result.validation.issues.map((issue, i) => (
                                  <li key={i} className="flex items-start">
                                    <span
                                      className={`w-2 h-2 mt-1.5 mr-2 rounded-full ${
                                        issue.severity === 'error'
                                          ? 'bg-red-500'
                                          : 'bg-yellow-500'
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
                              <p className="text-xs text-gray-500 mb-1">Alerts:</p>
                              <ul className="text-sm space-y-1">
                                {result.alerts.alerts.map((alert, i) => (
                                  <li key={i} className="text-amber-700">
                                    {alert.catalog}: {alert.rule_name}
                                  </li>
                                ))}
                              </ul>
                            </div>
                          )}

                          {/* ML-Readiness */}
                          {result.scoring?.ml_readiness && (
                            <div>
                              <p className="text-xs text-gray-500 mb-1">
                                ML-Readiness Score: {result.scoring.ml_readiness.score}
                              </p>
                              <p className="text-sm text-gray-600">
                                {result.scoring.ml_readiness.interpretation}
                              </p>
                            </div>
                          )}
                        </div>
                      </td>
                    </tr>
                  )}
                </>
              ))
            )}
          </tbody>
        </table>
      </div>

      {/* Pagination */}
      <div className="flex items-center justify-between px-4 py-3 bg-gray-50 rounded-lg">
        <div className="flex items-center space-x-2">
          <label className="text-sm text-gray-600">Per page:</label>
          <select
            value={pageSize}
            onChange={(e) => onPageSizeChange(Number(e.target.value))}
            className="border border-gray-300 rounded px-2 py-1 text-sm"
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
            className="px-3 py-1 border rounded text-sm disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-100"
          >
            Previous
          </button>
          <span className="text-sm text-gray-600">
            Page {page} of {totalPages}
          </span>
          <button
            onClick={() => onPageChange(page + 1)}
            disabled={page >= totalPages}
            className="px-3 py-1 border rounded text-sm disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-100"
          >
            Next
          </button>
        </div>
      </div>
    </div>
  );
}
