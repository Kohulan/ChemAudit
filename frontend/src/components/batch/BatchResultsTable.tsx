import { Fragment, useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { useNavigate } from 'react-router-dom';
import { AlertTriangle, ShieldAlert, FlaskConical, Atom, ExternalLink, Clipboard, Check, CheckCircle2, XCircle, Beaker, Activity } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { CopyButton } from '../ui/CopyButton';
import { cn } from '../../lib/utils';
import type { BatchResult, BatchResultsFilters, SortField } from '../../types/batch';
import type { RegistrationHashResult } from '../../types/analytics';
import { Tooltip, InfoTooltip } from '../ui/Tooltip';
import { MetricCard } from '../profiler/MetricCard';

interface BatchResultsTableProps {
  results: BatchResult[];
  page: number;
  pageSize: number;
  totalResults: number;
  totalPages: number;
  filters: BatchResultsFilters;
  sortBy: SortField;
  sortDir: 'asc' | 'desc';
  onPageChange: (page: number) => void;
  onPageSizeChange: (size: number) => void;
  onFiltersChange: (filters: BatchResultsFilters) => void;
  onSortChange: (field: SortField, dir: 'asc' | 'desc') => void;
  isLoading?: boolean;
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
  /** When set, auto-expand this molecule row and scroll to it. */
  focusedMoleculeIndex?: number | null;
  /** Called after the focused molecule has been scrolled to (so parent can clear it). */
  onFocusHandled?: () => void;
  /** Registration hash data for the optional RegHash column */
  registrationData?: RegistrationHashResult;
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
  sortBy,
  sortDir,
  onPageChange,
  onPageSizeChange,
  onFiltersChange,
  onSortChange,
  isLoading = false,
  selectedIndices,
  onSelectionChange,
  focusedMoleculeIndex,
  onFocusHandled,
  registrationData,
}: BatchResultsTableProps) {
  const navigate = useNavigate();
  const [expandedRow, setExpandedRow] = useState<number | null>(null);
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([]);
  const headerCheckboxRef = useRef<HTMLInputElement>(null);
  const focusedRowRef = useRef<HTMLTableRowElement>(null);
  const scrollContainerRef = useRef<HTMLDivElement>(null);
  const [scrollContainerWidth, setScrollContainerWidth] = useState(0);

  // Track the scroll container's visible width so expanded rows stay within bounds
  useEffect(() => {
    const el = scrollContainerRef.current;
    if (!el) return;
    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setScrollContainerWidth(entry.contentRect.width);
      }
    });
    observer.observe(el);
    setScrollContainerWidth(el.clientWidth);
    return () => observer.disconnect();
  }, []);

  // When focusedMoleculeIndex is set, expand that row and scroll to it
  useEffect(() => {
    if (focusedMoleculeIndex == null) return;
    if (results.some(r => r.index === focusedMoleculeIndex)) {
      setExpandedRow(focusedMoleculeIndex);
      setHighlightedAtoms([]);
      // Wait for render, then scroll to the specific row
      requestAnimationFrame(() => {
        focusedRowRef.current?.scrollIntoView({ behavior: 'smooth', block: 'center' });
      });
      // Keep highlight visible for 3s before clearing
      const timer = setTimeout(() => onFocusHandled?.(), 3000);
      return () => clearTimeout(timer);
    }
  }, [focusedMoleculeIndex, results, onFocusHandled]);

  const handleSort = (field: SortField) => {
    if (sortBy === field) {
      onSortChange(field, sortDir === 'asc' ? 'desc' : 'asc');
    } else {
      onSortChange(field, 'asc');
    }
  };

  const getScoreColor = (score: number | null | undefined): string => {
    if (score === null || score === undefined) return 'text-[var(--color-text-muted)]';
    if (score >= 80) return 'text-amber-600 dark:text-yellow-400 bg-yellow-500/10';
    if (score >= 50) return 'text-orange-600 dark:text-orange-400 bg-orange-500/10';
    return 'text-red-600 dark:text-red-400 bg-red-500/10';
  };

  // Calculate checkbox states
  const pageIndices = results.map(r => r.index);
  const allPageSelected = pageIndices.length > 0 && pageIndices.every(idx => selectedIndices.has(idx));
  const somePageSelected = pageIndices.some(idx => selectedIndices.has(idx)) && !allPageSelected;

  // Determine if profile column is shown (conditional 10th column)
  const hasProfileColumn = results.some(r => r.scoring?.profile);

  // Detect enrichment data
  const hasProfiling = results.some(r => r.profiling && !r.profiling.error);
  const hasSafetyAssessment = results.some(r => r.safety_assessment && !r.safety_assessment.error);

  // Determine if RegHash column is shown (conditional, only when registration data exists)
  const hasRegHashColumn = !!registrationData;

  const colCount = 9
    + (hasProfileColumn ? 1 : 0)
    + (hasRegHashColumn ? 1 : 0)
    + (hasProfiling ? 3 : 0)
    + (hasSafetyAssessment ? 3 : 0);

  // Build collision hash set for highlighting rows
  const collisionHashes = useMemo(() => {
    if (!registrationData) return new Set<string>();
    const s = new Set<string>();
    for (const g of registrationData.collision_groups) {
      s.add(g.hash);
    }
    return s;
  }, [registrationData]);

  // Build index-to-hash lookup
  const indexToHash = useMemo(() => {
    if (!registrationData) return new Map<number, string>();
    const m = new Map<number, string>();
    for (const mol of registrationData.per_molecule) {
      m.set(mol.index, mol.hash);
    }
    return m;
  }, [registrationData]);

  // Track which hash was recently copied (index -> timeout)
  const [copiedIndex, setCopiedIndex] = useState<number | null>(null);
  const copyTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const handleCopyHash = useCallback((hash: string, molIndex: number, e: React.MouseEvent) => {
    e.stopPropagation();
    navigator.clipboard.writeText(hash);
    setCopiedIndex(molIndex);
    if (copyTimeoutRef.current) clearTimeout(copyTimeoutRef.current);
    copyTimeoutRef.current = setTimeout(() => setCopiedIndex(null), 1500);
  }, []);

  // Update header checkbox indeterminate state
  useEffect(() => {
    if (headerCheckboxRef.current) {
      headerCheckboxRef.current.indeterminate = somePageSelected;
    }
  }, [somePageSelected]);

  // Handle header checkbox click
  const handleHeaderCheckboxClick = () => {
    const newSelection = new Set(selectedIndices);
    if (allPageSelected) {
      // Deselect all on page
      pageIndices.forEach(idx => newSelection.delete(idx));
    } else {
      // Select all on page
      pageIndices.forEach(idx => newSelection.add(idx));
    }
    onSelectionChange(newSelection);
  };

  // Handle individual checkbox click
  const handleRowCheckboxClick = (index: number, e: React.MouseEvent) => {
    e.stopPropagation();
    const newSelection = new Set(selectedIndices);
    if (newSelection.has(index)) {
      newSelection.delete(index);
    } else {
      newSelection.add(index);
    }
    onSelectionChange(newSelection);
  };

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

        {/* Active filter badges */}
        {filters.issue_filter && (
          <button
            type="button"
            onClick={() => onFiltersChange({ ...filters, issue_filter: undefined })}
            className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium bg-red-500/10 text-red-700 dark:text-red-400 border border-red-500/20 hover:bg-red-500/20 transition-colors"
          >
            Issue: {filters.issue_filter.replace(/_/g, ' ')}
            <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5">
              <path d="M18 6L6 18M6 6l12 12" />
            </svg>
          </button>
        )}
        {filters.alert_filter && (
          <button
            type="button"
            onClick={() => onFiltersChange({ ...filters, alert_filter: undefined })}
            className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium bg-amber-500/10 text-amber-700 dark:text-amber-400 border border-amber-500/20 hover:bg-amber-500/20 transition-colors"
          >
            Alert: {filters.alert_filter}
            <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5">
              <path d="M18 6L6 18M6 6l12 12" />
            </svg>
          </button>
        )}
        {filters.min_score !== undefined && filters.max_score !== undefined && (
          <button
            type="button"
            onClick={() => onFiltersChange({ ...filters, min_score: undefined, max_score: undefined })}
            className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium bg-blue-500/10 text-blue-700 dark:text-blue-400 border border-blue-500/20 hover:bg-blue-500/20 transition-colors"
          >
            Score: {filters.min_score}-{filters.max_score}
            <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5">
              <path d="M18 6L6 18M6 6l12 12" />
            </svg>
          </button>
        )}

        {/* Results info */}
        <div className="ml-auto text-sm text-[var(--color-text-muted)]">
          Showing {results.length} of {totalResults} results
        </div>
      </div>

      {/* Table */}
      <div ref={scrollContainerRef} className="overflow-x-auto">
        <table className="w-full border-collapse">
          <thead>
            <tr className="bg-[var(--color-surface-sunken)]">
              <th className="w-10 px-4 py-3 text-center">
                <input
                  ref={headerCheckboxRef}
                  type="checkbox"
                  checked={allPageSelected}
                  onChange={handleHeaderCheckboxClick}
                  className="w-4 h-4 cursor-pointer accent-[var(--color-primary)]"
                  aria-label="Select all on page"
                />
              </th>
              <th
                className="px-4 py-3 text-left text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('index')}
              >
                # {sortBy === 'index' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-left text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('name')}
              >
                Name {sortBy === 'name' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-left text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('smiles')}
              >
                SMILES {sortBy === 'smiles' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('score')}
              >
                <span className="inline-flex items-center gap-1">
                  Score {sortBy === 'score' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
                  <InfoTooltip size="small" asSpan content="Composite validation quality score (0–100)" />
                </span>
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('qed')}
              >
                <span className="inline-flex items-center gap-1">
                  QED {sortBy === 'qed' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
                  <InfoTooltip size="small" asSpan content="Quantitative Estimate of Drug-likeness (0–1)" />
                </span>
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('safety')}
              >
                <span className="inline-flex items-center gap-1">
                  Safety {sortBy === 'safety' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
                  <InfoTooltip size="small" asSpan content="Number of structural alerts flagged" />
                </span>
              </th>
              {results.some(r => r.scoring?.profile) && (
                <th
                  className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                  onClick={() => handleSort('profile_score')}
                >
                  Profile {sortBy === 'profile_score' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
                </th>
              )}
              {hasProfiling && (
                <>
                  <th className="px-3 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                    <span className="inline-flex items-center gap-1">
                      PFI
                      <InfoTooltip size="small" asSpan content="Property Forecast Index — cLogP + aromatic rings" />
                    </span>
                  </th>
                  <th className="px-3 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                    <span className="inline-flex items-center gap-1">
                      Stars
                      <InfoTooltip size="small" asSpan content="Properties outside 95th-percentile drug ranges" />
                    </span>
                  </th>
                  <th className="px-3 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                    <span className="inline-flex items-center gap-1">
                      Abbott
                      <InfoTooltip size="small" asSpan content="Predicted oral bioavailability probability" />
                    </span>
                  </th>
                </>
              )}
              {hasSafetyAssessment && (
                <>
                  <th className="px-3 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                    <span className="inline-flex items-center gap-1">
                      CYP
                      <InfoTooltip size="small" asSpan content="CYP450 metabolic soft-spot count" />
                    </span>
                  </th>
                  <th className="px-3 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                    <span className="inline-flex items-center gap-1">
                      hERG
                      <InfoTooltip size="small" asSpan content="hERG cardiac channel inhibition risk" />
                    </span>
                  </th>
                  <th className="px-3 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase">
                    <span className="inline-flex items-center gap-1">
                      bRo5
                      <InfoTooltip size="small" asSpan content="Beyond Rule of 5 — for molecules with MW > 500" />
                    </span>
                  </th>
                </>
              )}
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('status')}
              >
                <span className="inline-flex items-center gap-1">
                  Status {sortBy === 'status' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
                  <InfoTooltip size="small" asSpan content="Validation pass/fail status" />
                </span>
              </th>
              <th
                className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase cursor-pointer hover:bg-[var(--color-surface-elevated)]"
                onClick={() => handleSort('issues')}
              >
                <span className="inline-flex items-center gap-1">
                  Issues {sortBy === 'issues' && (sortDir === 'asc' ? '\u2191' : '\u2193')}
                  <InfoTooltip size="small" asSpan content="Count of failed validation checks" />
                </span>
              </th>
              {hasRegHashColumn && (
                <th
                  className="px-4 py-3 text-center text-xs font-medium text-[var(--color-text-muted)] uppercase hidden md:table-cell"
                  style={{ width: '120px' }}
                >
                  <span className="inline-flex items-center gap-1">
                    RegHash
                    <InfoTooltip size="small" asSpan content="Registration hash for duplicate detection" />
                  </span>
                </th>
              )}
            </tr>
          </thead>
          <tbody className="divide-y divide-[var(--color-border)]">
            {isLoading ? (
              <tr>
                <td colSpan={colCount} className="px-4 py-8 text-center text-[var(--color-text-muted)]">
                  <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-[var(--color-primary)] mx-auto mb-2" />
                  Loading results...
                </td>
              </tr>
            ) : results.length === 0 ? (
              <tr>
                <td colSpan={colCount} className="px-4 py-8 text-center text-[var(--color-text-muted)]">
                  No results match the current filters.
                </td>
              </tr>
            ) : (
              results.map((result) => (
                <Fragment key={result.index}>
                  <tr
                    ref={result.index === focusedMoleculeIndex ? focusedRowRef : undefined}
                    className={`
                      hover:bg-[var(--color-surface-sunken)] cursor-pointer transition-colors
                      ${result.status === 'error' ? 'bg-red-500/5' : ''}
                      ${expandedRow === result.index ? 'bg-[var(--color-primary)]/5' : ''}
                      ${result.index === focusedMoleculeIndex ? 'ring-2 ring-inset ring-[var(--color-primary)]/40' : ''}
                    `}
                    onClick={() => {
                      setExpandedRow(expandedRow === result.index ? null : result.index);
                      setHighlightedAtoms([]);
                    }}
                  >
                    <td className="px-4 py-3 text-center" onClick={(e) => e.stopPropagation()}>
                      <input
                        type="checkbox"
                        checked={selectedIndices.has(result.index)}
                        onChange={(e) => handleRowCheckboxClick(result.index, e as any)}
                        className="w-4 h-4 cursor-pointer accent-[var(--color-primary)]"
                        aria-label={`Select molecule ${result.index + 1}`}
                      />
                    </td>
                    <td className="px-4 py-3 text-sm text-[var(--color-text-muted)]">
                      {result.index + 1}
                    </td>
                    <td className="px-4 py-3 text-sm font-medium text-[var(--color-text-primary)]">
                      {result.name || '-'}
                    </td>
                    <td className="px-4 py-3 text-sm text-[var(--color-text-secondary)] font-mono">
                      <span className="group/smiles inline-flex items-center gap-1 max-w-[280px]">
                        <span className="truncate" title={result.smiles}>
                          {result.smiles}
                        </span>
                        <CopyButton
                          text={result.smiles}
                          size={13}
                          className="shrink-0 opacity-0 group-hover/smiles:opacity-100 transition-opacity"
                        />
                      </span>
                    </td>
                    <td className="px-4 py-3 text-center">
                      <Tooltip
                        content={
                          result.scoring?.ml_readiness?.interpretation
                            ? `Validation: ${result.validation?.overall_score ?? '-'}/100 — ${result.scoring.ml_readiness.interpretation}`
                            : `Validation Score: ${result.validation?.overall_score ?? '-'}/100`
                        }
                        position="bottom"
                      >
                        <span
                          className={`px-2 py-1 rounded text-sm font-medium ${getScoreColor(
                            result.validation?.overall_score
                          )}`}
                        >
                          {result.validation?.overall_score ?? '-'}
                        </span>
                      </Tooltip>
                    </td>
                    <td className="px-4 py-3 text-center">
                      {result.scoring?.druglikeness ? (
                        <Tooltip
                          content={`QED: ${result.scoring.druglikeness.qed_score.toFixed(2)} — Lipinski: ${result.scoring.druglikeness.lipinski_passed ? 'Pass' : `Fail (${result.scoring.druglikeness.lipinski_violations} violations)`}`}
                          position="bottom"
                        >
                          <span className="px-2 py-1 rounded text-xs font-medium bg-purple-500/10 text-purple-600 dark:text-purple-400">
                            {result.scoring.druglikeness.qed_score.toFixed(2)}
                          </span>
                        </Tooltip>
                      ) : (
                        <span className="text-[var(--color-text-muted)]">-</span>
                      )}
                    </td>
                    <td className="px-4 py-3 text-center">
                      {(() => {
                        const alertCount = result.alerts?.alert_count ?? 0;
                        const allClear = alertCount === 0 && (result.scoring?.safety_filters?.all_passed !== false);

                        // Build catalog breakdown for tooltip
                        const catalogCounts: Record<string, number> = {};
                        result.alerts?.alerts?.forEach(a => {
                          catalogCounts[a.catalog] = (catalogCounts[a.catalog] || 0) + 1;
                        });
                        const breakdown = Object.entries(catalogCounts)
                          .map(([cat, count]) => `${count} ${cat}`)
                          .join(', ');
                        const tooltipText = allClear
                          ? 'All safety screens passed'
                          : `${alertCount} alert${alertCount !== 1 ? 's' : ''}: ${breakdown}`;

                        return (
                          <Tooltip content={tooltipText} position="bottom">
                            {allClear ? (
                              <span className="px-2 py-1 rounded text-xs font-medium bg-emerald-500/10 text-emerald-600 dark:text-emerald-400">
                                ✓
                              </span>
                            ) : (
                              <span className="px-2 py-1 rounded text-xs font-medium bg-red-500/10 text-red-600 dark:text-red-400">
                                {alertCount}
                              </span>
                            )}
                          </Tooltip>
                        );
                      })()}
                    </td>
                    {hasProfileColumn && (
                      <td className="px-4 py-3 text-center">
                        {result.scoring?.profile?.score != null ? (
                          <span className={cn(
                            'inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium',
                            result.scoring.profile.score >= 80
                              ? 'bg-green-500/10 text-green-600 dark:text-green-400'
                              : result.scoring.profile.score >= 50
                              ? 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                              : 'bg-red-500/10 text-red-600 dark:text-red-400'
                          )}>
                            {result.scoring.profile.score.toFixed(1)}
                          </span>
                        ) : (
                          <span className="text-[var(--color-text-muted)]">--</span>
                        )}
                      </td>
                    )}
                    {hasProfiling && (
                      <>
                        <td className="px-3 py-3 text-center text-xs text-[var(--color-text-secondary)]">
                          {result.profiling?.pfi?.pfi != null ? (
                            <Tooltip
                              content={`PFI ${result.profiling.pfi.pfi.toFixed(1)} (${result.profiling.pfi.risk} risk) — cLogP: ${result.profiling.pfi.clogp.toFixed(2)}, Aromatic rings: ${result.profiling.pfi.aromatic_rings}`}
                              position="bottom"
                            >
                              <span className={cn(
                                'px-1.5 py-0.5 rounded',
                                result.profiling.pfi.risk === 'low'
                                  ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                                  : result.profiling.pfi.risk === 'moderate'
                                  ? 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                                  : 'bg-red-500/10 text-red-600 dark:text-red-400'
                              )}>
                                {result.profiling.pfi.pfi.toFixed(1)}
                              </span>
                            </Tooltip>
                          ) : (
                            <span className="text-[var(--color-text-muted)]">&mdash;</span>
                          )}
                        </td>
                        <td className="px-3 py-3 text-center text-xs text-[var(--color-text-secondary)]">
                          {result.profiling?.stars?.stars != null ? (() => {
                            const violated = result.profiling!.stars!.details.filter(d => d.violated);
                            const summary = violated.length > 0
                              ? violated.map(d => `${d.property}: ${d.value.toFixed(1)} (range ${d.range_low}–${d.range_high})`).join('; ')
                              : 'All properties within 95th-percentile drug ranges';
                            return (
                              <Tooltip
                                content={`${result.profiling!.stars!.stars} of ${result.profiling!.stars!.details.length} outside range${violated.length > 0 ? ` — ${summary}` : ''}`}
                                position="bottom"
                                maxWidth={360}
                              >
                                <span>{result.profiling!.stars!.stars}★</span>
                              </Tooltip>
                            );
                          })() : <span className="text-[var(--color-text-muted)]">&mdash;</span>}
                        </td>
                        <td className="px-3 py-3 text-center text-xs text-[var(--color-text-secondary)]">
                          {result.profiling?.abbott?.probability_pct != null ? (
                            <Tooltip
                              content={`${result.profiling.abbott.probability_pct}% predicted bioavailability — TPSA: ${result.profiling.abbott.tpsa.toFixed(1)} Å², Lipinski violations: ${result.profiling.abbott.lipinski_violations}`}
                              position="bottom"
                            >
                              <span>{result.profiling.abbott.probability_pct}%</span>
                            </Tooltip>
                          ) : <span className="text-[var(--color-text-muted)]">&mdash;</span>}
                        </td>
                      </>
                    )}
                    {hasSafetyAssessment && (
                      <>
                        <td className="px-3 py-3 text-center text-xs">
                          {result.safety_assessment?.cyp_softspots != null ? (() => {
                            const sites = result.safety_assessment!.cyp_softspots!;
                            const summary = sites.length === 0
                              ? 'No CYP soft-spots detected'
                              : sites.map(s => s.site_name).join(', ');
                            return (
                              <Tooltip content={summary} position="bottom" maxWidth={320}>
                                <span className={cn(
                                  'px-1.5 py-0.5 rounded',
                                  sites.length === 0
                                    ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                                    : 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                                )}>
                                  {sites.length === 0 ? 'clean' : `${sites.length} sites`}
                                </span>
                              </Tooltip>
                            );
                          })() : <span className="text-[var(--color-text-muted)]">&mdash;</span>}
                        </td>
                        <td className="px-3 py-3 text-center text-xs">
                          {result.safety_assessment?.herg?.herg_risk != null ? (() => {
                            const h = result.safety_assessment!.herg!;
                            const flagsSummary = h.flags.length > 0 ? ` — ${h.flags.join('; ')}` : '';
                            return (
                              <Tooltip
                                content={`hERG risk: ${h.herg_risk} (${h.risk_score}/${h.max_score})${flagsSummary}`}
                                position="bottom"
                                maxWidth={360}
                              >
                                <span className={cn(
                                  'px-1.5 py-0.5 rounded',
                                  h.herg_risk === 'low'
                                    ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                                    : h.herg_risk === 'moderate'
                                    ? 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                                    : 'bg-red-500/10 text-red-600 dark:text-red-400'
                                )}>
                                  {h.herg_risk}
                                </span>
                              </Tooltip>
                            );
                          })() : <span className="text-[var(--color-text-muted)]">&mdash;</span>}
                        </td>
                        <td className="px-3 py-3 text-center text-xs">
                          {result.safety_assessment?.bro5 != null ? (() => {
                            const b = result.safety_assessment!.bro5!;
                            const tooltipText = !b.applicable
                              ? 'Not applicable — MW ≤ 500, standard Ro5 applies'
                              : b.passed
                                ? 'All bRo5 criteria passed'
                                : b.violations.map(v => `${v.property}: ${v.value.toFixed(1)} (threshold: ${v.threshold})`).join('; ');
                            return (
                              <Tooltip content={tooltipText} position="bottom" maxWidth={360}>
                                {!b.applicable ? (
                                  <span className="text-[var(--color-text-muted)]">N/A</span>
                                ) : (
                                  <span className={cn(
                                    'px-1.5 py-0.5 rounded',
                                    b.passed
                                      ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                                      : 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
                                  )}>
                                    {b.passed ? 'pass' : `${b.violations.length} viol`}
                                  </span>
                                )}
                              </Tooltip>
                            );
                          })() : <span className="text-[var(--color-text-muted)]">&mdash;</span>}
                        </td>
                      </>
                    )}
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
                      {(() => {
                        const failedIssues = result.validation?.issues?.filter(i => !i.passed) ?? [];
                        const count = failedIssues.length;
                        if (count === 0) {
                          return <span className="text-[var(--color-text-muted)]">0</span>;
                        }
                        const hasError = failedIssues.some(i => i.severity === 'error');
                        return (
                          <span className={`px-2 py-1 rounded text-xs font-medium ${
                            hasError
                              ? 'bg-red-500/10 text-red-600 dark:text-red-400'
                              : 'bg-amber-500/10 text-amber-700 dark:text-amber-400'
                          }`}>
                            {count}
                          </span>
                        );
                      })()}
                    </td>
                    {hasRegHashColumn && (() => {
                      const hash = indexToHash.get(result.index);
                      const isCollision = hash ? collisionHashes.has(hash) : false;
                      return (
                        <td
                          className={cn(
                            'px-4 py-3 text-center hidden md:table-cell',
                            isCollision && 'bg-amber-100/50 dark:bg-amber-900/20',
                          )}
                          style={{ width: '120px' }}
                        >
                          {hash ? (
                            <span className="inline-flex items-center gap-1">
                              <span className="text-xs font-mono text-[var(--color-text-secondary)]" title={hash}>
                                {hash.slice(0, 12)}
                              </span>
                              <button
                                type="button"
                                onClick={(e) => handleCopyHash(hash, result.index, e)}
                                className="p-0.5 rounded hover:bg-[var(--color-surface-sunken)] transition-colors text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]"
                                aria-label="Copy registration hash to clipboard"
                              >
                                {copiedIndex === result.index ? (
                                  <Check className="w-3 h-3 text-green-500" />
                                ) : (
                                  <Clipboard className="w-3 h-3" />
                                )}
                              </button>
                            </span>
                          ) : (
                            <span className="text-[var(--color-text-muted)]">-</span>
                          )}
                        </td>
                      );
                    })()}
                  </tr>

                  {/* Expanded details — bento grid layout */}
                  {expandedRow === result.index && (
                    <tr>
                      <td colSpan={colCount} className="p-0">
                        <div
                          className="px-3 py-3 bg-[var(--color-surface-sunken)] border-t border-[var(--color-border)]"
                          style={scrollContainerWidth ? { width: scrollContainerWidth, maxWidth: scrollContainerWidth } : undefined}
                        >
                          {/* Error banner — full width */}
                          {result.error && (
                            <div className="mb-2 rounded-xl bg-red-500/10 dark:bg-red-500/15 border border-red-500/20 p-2.5 flex items-center gap-2.5">
                              <div className="w-8 h-8 rounded-lg bg-red-500/20 flex items-center justify-center flex-shrink-0">
                                <AlertTriangle className="w-4 h-4 text-red-500" />
                              </div>
                              <p className="text-sm text-red-600 dark:text-red-400">{result.error}</p>
                            </div>
                          )}

                          {/* ── Expanded detail layout ── */}
                          <div className="grid grid-cols-1 lg:grid-cols-[1fr_1fr_minmax(250px,1fr)] lg:grid-rows-[auto_auto_auto] gap-2">

                                  {/* Validation Issues — col 1, row 1 */}
                                  <div className="card rounded-xl p-3 lg:col-start-1 lg:row-start-1">
                                    <div className="flex items-center gap-2 mb-2">
                                      <div className="w-7 h-7 rounded-lg bg-gradient-to-br from-amber-500/20 to-red-500/10 flex items-center justify-center">
                                        <AlertTriangle className="w-3.5 h-3.5 text-amber-600 dark:text-amber-400" />
                                      </div>
                                      <h4 className="text-xs font-semibold text-[var(--color-text-primary)] tracking-tight uppercase">
                                        Validation Issues
                                        <span className="ml-1.5 text-[var(--color-text-muted)] font-normal">({result.validation?.issues?.length ?? 0})</span>
                                      </h4>
                                    </div>
                                    {result.validation?.issues && result.validation.issues.length > 0 ? (
                                      <ul className="space-y-0.5">
                                        {result.validation.issues.map((issue, i) => (
                                          <li
                                            key={i}
                                            className="flex items-start gap-2 text-sm cursor-default rounded-lg px-2 py-1 hover:bg-[var(--color-surface-sunken)] transition-colors"
                                            onMouseEnter={() => setHighlightedAtoms(issue.affected_atoms ?? [])}
                                            onMouseLeave={() => setHighlightedAtoms([])}
                                          >
                                            <span className={cn("w-2 h-2 mt-[5px] rounded-full flex-shrink-0", issue.severity === 'error' ? 'bg-red-500' : 'bg-amber-500')} />
                                            <span className="flex-1 text-[var(--color-text-secondary)] leading-snug">{issue.message}</span>
                                            {issue.affected_atoms && issue.affected_atoms.length > 0 && (
                                              <Atom className="w-3.5 h-3.5 text-[var(--color-text-muted)] flex-shrink-0 mt-0.5" />
                                            )}
                                          </li>
                                        ))}
                                      </ul>
                                    ) : (
                                      <p className="text-sm text-emerald-600 dark:text-emerald-400 px-2 py-1">No structural issues detected</p>
                                    )}
                                  </div>

                                  {/* Scores — col 2, row 1 */}
                                  {result.scoring && (
                                    <div className="card rounded-xl p-3 lg:col-start-2 lg:row-start-1">
                                      <div className="flex items-center gap-2 mb-2">
                                        <div className="w-7 h-7 rounded-lg bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/10 flex items-center justify-center">
                                          <Activity className="w-3.5 h-3.5 text-[var(--color-primary)]" />
                                        </div>
                                        <h4 className="text-xs font-semibold text-[var(--color-text-primary)] tracking-tight uppercase">Scores</h4>
                                      </div>
                                      <div className="grid grid-cols-2 gap-1.5">
                                        {result.scoring.ml_readiness && (
                                          <div className="rounded-lg bg-gradient-to-br from-[var(--color-primary)]/[0.07] to-transparent border border-[var(--color-border-subtle)] p-2 text-center">
                                            <p className="text-xs uppercase tracking-wider text-[var(--color-text-muted)] font-medium">ML-Ready</p>
                                            <p className="text-xl font-bold text-gradient">{result.scoring.ml_readiness.score}</p>
                                          </div>
                                        )}
                                        {result.scoring.druglikeness && (
                                          <div className="rounded-lg bg-gradient-to-br from-purple-500/[0.07] to-transparent border border-[var(--color-border-subtle)] p-2 text-center">
                                            <p className="text-xs uppercase tracking-wider text-[var(--color-text-muted)] font-medium">QED</p>
                                            <p className="text-xl font-bold text-purple-600 dark:text-purple-400">{result.scoring.druglikeness.qed_score.toFixed(2)}</p>
                                            <p className={cn("text-xs font-medium", result.scoring.druglikeness.lipinski_passed ? 'text-emerald-600 dark:text-emerald-400' : 'text-red-600 dark:text-red-400')}>
                                              Lipinski: {result.scoring.druglikeness.lipinski_passed ? 'Pass' : 'Fail'}
                                            </p>
                                          </div>
                                        )}
                                        {(result.scoring.safety_filters || result.alerts) && (
                                          <div className="rounded-lg bg-gradient-to-br from-emerald-500/[0.07] to-transparent border border-[var(--color-border-subtle)] p-2 text-center">
                                            <p className="text-xs uppercase tracking-wider text-[var(--color-text-muted)] font-medium">Safety</p>
                                            {(() => {
                                              const ac = result.alerts?.alert_count ?? 0;
                                              const ok = ac === 0 && (result.scoring.safety_filters?.all_passed !== false);
                                              return <p className={cn("text-xl font-bold", ok ? 'text-emerald-600 dark:text-emerald-400' : 'text-red-600 dark:text-red-400')}>{ok ? 'Clear' : ac}</p>;
                                            })()}
                                          </div>
                                        )}
                                        {result.scoring.admet && (
                                          <div className="rounded-lg bg-gradient-to-br from-cyan-500/[0.07] to-transparent border border-[var(--color-border-subtle)] p-2 text-center">
                                            <p className="text-xs uppercase tracking-wider text-[var(--color-text-muted)] font-medium">SA</p>
                                            <p className="text-xl font-bold text-cyan-600 dark:text-cyan-400">{result.scoring.admet.sa_score.toFixed(1)}</p>
                                            <p className="text-xs text-[var(--color-text-muted)]">{result.scoring.admet.sa_classification}</p>
                                          </div>
                                        )}
                                      </div>
                                      {result.scoring.ml_readiness?.interpretation && (
                                        <p className="text-sm text-[var(--color-text-secondary)] mt-2 px-1 leading-relaxed">{result.scoring.ml_readiness.interpretation}</p>
                                      )}
                                    </div>
                                  )}

                              {/* Structure — col 3, rows 1-2 */}
                              <div className="card rounded-xl p-3 flex flex-col lg:col-start-3 lg:row-start-1 lg:row-end-3">
                                {result.status === 'success' ? (
                                  <>
                                    <div className="flex items-center gap-2 mb-2">
                                      <div className="w-7 h-7 rounded-lg bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/10 flex items-center justify-center">
                                        <FlaskConical className="w-3.5 h-3.5 text-[var(--color-primary)]" />
                                      </div>
                                      <div className="flex-1 min-w-0">
                                        <div className="flex items-center gap-1.5">
                                          <h4 className="text-xs font-semibold text-[var(--color-text-primary)] tracking-tight uppercase">Structure</h4>
                                          {result.standardization?.standardized_smiles && result.standardization.changed && (
                                            <span className="px-1.5 py-0.5 rounded text-[9px] font-medium bg-emerald-500/15 text-emerald-600 dark:text-emerald-400 uppercase tracking-wide">Std</span>
                                          )}
                                        </div>
                                        {result.name && <p className="text-sm text-[var(--color-text-muted)] truncate">{result.name}</p>}
                                      </div>
                                    </div>
                                    <div className="flex-1 flex items-center justify-center rounded-lg bg-white dark:bg-gray-900/50 border border-[var(--color-border-subtle)]">
                                      <MoleculeViewer smiles={result.standardization?.standardized_smiles || result.smiles} highlightAtoms={highlightedAtoms} width={240} height={162} />
                                    </div>
                                    {highlightedAtoms.length > 0 && (
                                      <div className="mt-1.5 flex justify-center">
                                        <div className="px-2.5 py-0.5 rounded-full bg-[var(--color-primary)]/10 border border-[var(--color-primary)]/20">
                                          <p className="text-xs font-medium text-[var(--color-primary)] whitespace-nowrap">{highlightedAtoms.length} atom{highlightedAtoms.length !== 1 ? 's' : ''} highlighted</p>
                                        </div>
                                      </div>
                                    )}
                                    {result.standardization?.standardized_smiles && result.standardization.changed ? (
                                      <div className="mt-2 space-y-1">
                                        <div className="rounded-lg bg-[var(--color-surface-sunken)] border border-[var(--color-border-subtle)] px-2.5 py-1.5">
                                          <p className="text-[9px] font-medium text-[var(--color-text-muted)] uppercase tracking-wide mb-0.5">Original</p>
                                          <div className="flex items-start gap-1">
                                            <p className="text-xs font-mono text-[var(--color-text-secondary)] break-all flex-1">{result.smiles}</p>
                                            <CopyButton text={result.smiles} size={13} className="shrink-0" />
                                          </div>
                                        </div>
                                        <div className="rounded-lg bg-emerald-500/5 border border-emerald-500/20 px-2.5 py-1.5">
                                          <p className="text-[9px] font-medium text-emerald-600 dark:text-emerald-400 uppercase tracking-wide mb-0.5">Standardized</p>
                                          <div className="flex items-start gap-1">
                                            <p className="text-xs font-mono text-[var(--color-text-secondary)] break-all flex-1">{result.standardization.standardized_smiles}</p>
                                            <CopyButton text={result.standardization.standardized_smiles} size={13} className="shrink-0" />
                                          </div>
                                        </div>
                                      </div>
                                    ) : (
                                      <div className="mt-2 flex items-start gap-1 rounded-lg bg-[var(--color-surface-sunken)] border border-[var(--color-border-subtle)] px-2.5 py-1.5">
                                        <p className="text-xs font-mono text-[var(--color-text-secondary)] break-all flex-1">{result.smiles}</p>
                                        <CopyButton text={result.smiles} size={13} className="shrink-0" />
                                      </div>
                                    )}
                                    <button
                                      onClick={(e) => { e.stopPropagation(); navigate('/', { state: { fromBatch: true, smiles: result.standardization?.standardized_smiles || result.smiles, moleculeName: result.name, moleculeIndex: result.index } }); }}
                                      className={cn('mt-2 w-full flex items-center justify-center gap-2', 'px-3.5 py-2.5 rounded-xl text-sm font-semibold', 'bg-gradient-to-r from-[var(--color-primary)] to-[var(--color-accent)]', 'text-white shadow-md shadow-[var(--color-primary)]/25', 'hover:shadow-lg hover:shadow-[var(--color-primary)]/35 hover:brightness-110', 'active:scale-[0.98]', 'transition-all duration-200 cursor-pointer')}
                                    >
                                      <ExternalLink className="w-4 h-4" />
                                      Full Analysis
                                    </button>
                                  </>
                                ) : (
                                  <div className="flex-1 flex items-center justify-center rounded-lg bg-red-500/5 border border-red-500/10">
                                    <p className="text-sm text-red-500/70">Parse failed</p>
                                  </div>
                                )}
                              </div>

                            {/* Profiling Details — cols 1-2, row 2 */}
                            {result.profiling && !result.profiling.error && (
                              <div className="card rounded-xl p-3 lg:col-start-1 lg:col-end-3 lg:row-start-2">
                                <div className="flex items-center gap-2 mb-2">
                                  <div className="w-7 h-7 rounded-lg bg-gradient-to-br from-blue-500/20 to-cyan-500/10 flex items-center justify-center">
                                    <Beaker className="w-3.5 h-3.5 text-blue-600 dark:text-blue-400" />
                                  </div>
                                  <h4 className="text-xs font-semibold text-[var(--color-text-primary)] tracking-tight uppercase">Profiling Details</h4>
                                </div>
                                <div className="grid grid-cols-1 sm:grid-cols-3 gap-2">
                                  {result.profiling.pfi && (
                                    <MetricCard title="PFI" score={result.profiling.pfi.pfi.toFixed(2)} classification={result.profiling.pfi.risk === 'low' ? 'Low Risk' : result.profiling.pfi.risk === 'moderate' ? 'Moderate' : 'High Risk'} classificationVariant={result.profiling.pfi.risk === 'low' ? 'success' : result.profiling.pfi.risk === 'moderate' ? 'warning' : 'error'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                      <div className="space-y-1 text-xs text-[var(--color-text-secondary)]">
                                        <p className="font-medium text-[var(--color-text-primary)]">cLogP + #Aromatic Rings</p>
                                        <div className="grid grid-cols-2 gap-x-3 gap-y-0.5">
                                          <span className="text-[var(--color-text-muted)]">cLogP</span><span className="tabular-nums text-[var(--color-text-primary)] font-medium">{result.profiling.pfi.clogp.toFixed(2)}</span>
                                          <span className="text-[var(--color-text-muted)]">Aromatic Rings</span><span className="tabular-nums text-[var(--color-text-primary)] font-medium">{result.profiling.pfi.aromatic_rings}</span>
                                        </div>
                                      </div>
                                    </MetricCard>
                                  )}
                                  {result.profiling.stars && (
                                    <MetricCard title="#Stars" score={result.profiling.stars.stars} classification={result.profiling.stars.stars === 0 ? 'Drug-like' : result.profiling.stars.stars <= 2 ? 'Borderline' : 'Outlier'} classificationVariant={result.profiling.stars.stars === 0 ? 'success' : result.profiling.stars.stars <= 2 ? 'warning' : 'error'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                      <div className="space-y-0.5">
                                        {result.profiling.stars.details.map((d) => (
                                          <div key={d.property} className={cn('flex items-center gap-1.5 text-xs', d.violated ? 'text-red-600 dark:text-red-400' : 'text-[var(--color-text-secondary)]')}>
                                            {d.violated ? <XCircle className="w-3.5 h-3.5 shrink-0" /> : <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 shrink-0" />}
                                            <span className="flex-1 truncate">{d.property}</span>
                                            <span className="tabular-nums shrink-0">{d.value.toFixed(1)}</span>
                                            <span className="tabular-nums text-[var(--color-text-muted)] shrink-0 text-[11px]">({d.range_low}-{d.range_high})</span>
                                          </div>
                                        ))}
                                      </div>
                                    </MetricCard>
                                  )}
                                  {result.profiling.abbott && (
                                    <MetricCard title="Abbott" score={`${result.profiling.abbott.probability_pct}%`} classification={result.profiling.abbott.probability_pct >= 85 ? 'High' : result.profiling.abbott.probability_pct >= 56 ? 'Moderate' : 'Low'} classificationVariant={result.profiling.abbott.probability_pct >= 85 ? 'success' : result.profiling.abbott.probability_pct >= 56 ? 'warning' : 'error'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                      <div className="space-y-1 text-xs text-[var(--color-text-secondary)]">
                                        <div className="grid grid-cols-2 gap-x-3 gap-y-0.5">
                                          <span className="text-[var(--color-text-muted)]">TPSA</span><span className="tabular-nums text-[var(--color-text-primary)] font-medium">{result.profiling.abbott.tpsa.toFixed(1)} Å²</span>
                                          <span className="text-[var(--color-text-muted)]">Lipinski viol.</span><span className="tabular-nums text-[var(--color-text-primary)] font-medium">{result.profiling.abbott.lipinski_violations}</span>
                                        </div>
                                      </div>
                                    </MetricCard>
                                  )}
                                </div>
                              </div>
                            )}

                            {/* Safety — cols 1-3, row 3 */}
                            <div className="card rounded-xl p-3 lg:col-span-3 lg:row-start-3">
                              <div className="flex items-center gap-2 mb-2">
                                <div className="w-7 h-7 rounded-lg bg-gradient-to-br from-red-500/20 to-rose-500/10 flex items-center justify-center">
                                  <ShieldAlert className="w-3.5 h-3.5 text-red-600 dark:text-red-400" />
                                </div>
                                <h4 className="text-xs font-semibold text-[var(--color-text-primary)] tracking-tight uppercase">
                                  Safety
                                  <span className="ml-1.5 text-[var(--color-text-muted)] font-normal">({result.alerts?.alert_count ?? 0} alerts)</span>
                                </h4>
                              </div>
                              <div className={cn(
                                "grid gap-2",
                                result.safety_assessment && !result.safety_assessment.error
                                  ? 'grid-cols-2 sm:grid-cols-3 lg:grid-cols-6'
                                  : 'grid-cols-1'
                              )}>
                                {/* Alerts card */}
                                <MetricCard
                                  title="Structural Alerts"
                                  score={result.alerts?.alert_count ?? 0}
                                  classification={(result.alerts?.alert_count ?? 0) === 0 ? 'Clear' : `${result.alerts?.alert_count} flagged`}
                                  classificationVariant={(result.alerts?.alert_count ?? 0) === 0 ? 'success' : 'error'}
                                  defaultExpanded={true}
                                  className="!p-3 !rounded-lg"
                                >
                                  {result.alerts?.alerts && result.alerts.alerts.length > 0 ? (
                                    <div className="space-y-1">
                                      {result.alerts.alerts.map((alert, i) => (
                                        <div key={i} className="flex items-center gap-1.5 text-xs cursor-default hover:bg-[var(--color-surface-sunken)] rounded px-1 py-0.5 transition-colors"
                                          onMouseEnter={() => setHighlightedAtoms(alert.matched_atoms ?? [])} onMouseLeave={() => setHighlightedAtoms([])}>
                                          <span className={cn("px-1 py-0.5 rounded text-[9px] font-semibold shrink-0 uppercase", alert.severity === 'critical' ? 'bg-red-500/20 text-red-600 dark:text-red-400' : 'bg-amber-500/15 text-amber-700 dark:text-amber-400')}>{alert.catalog}</span>
                                          <span className="text-[var(--color-text-secondary)] truncate">{alert.rule_name}</span>
                                          {alert.matched_atoms && alert.matched_atoms.length > 0 && <Atom className="w-3 h-3 text-[var(--color-text-muted)] shrink-0 ml-auto" />}
                                        </div>
                                      ))}
                                      {result.scoring?.safety_filters && (
                                        <div className="pt-1 border-t border-[var(--color-border)] text-xs text-[var(--color-text-muted)] flex flex-wrap gap-x-2 gap-y-0.5">
                                          <span>PAINS: {result.scoring.safety_filters.pains_passed ? '✓' : '✗'}</span>
                                          <span>Brenk: {result.scoring.safety_filters.brenk_passed ? '✓' : '✗'}</span>
                                          {result.scoring.safety_filters.nih_passed !== undefined && <span>NIH: {result.scoring.safety_filters.nih_passed ? '✓' : '✗'}</span>}
                                          {result.scoring.safety_filters.zinc_passed !== undefined && <span>ZINC: {result.scoring.safety_filters.zinc_passed ? '✓' : '✗'}</span>}
                                          {result.scoring.safety_filters.chembl_passed !== undefined && <span>ChEMBL: {result.scoring.safety_filters.chembl_passed ? '✓' : '✗'}</span>}
                                        </div>
                                      )}
                                    </div>
                                  ) : (
                                    <p className="text-xs text-emerald-600 dark:text-emerald-400">All screens passed</p>
                                  )}
                                </MetricCard>
                                {/* Assessment cards inline */}
                                {result.safety_assessment && !result.safety_assessment.error && (
                                  <>
                                    {result.safety_assessment.cyp_softspots && (
                                      <MetricCard title="CYP" score={result.safety_assessment.cyp_softspots.length} classification={result.safety_assessment.cyp_softspots.length > 0 ? `${result.safety_assessment.cyp_softspots.length} site(s)` : 'None'} classificationVariant="default" defaultExpanded={true} className="!p-3 !rounded-lg">
                                        <div className="space-y-1">{result.safety_assessment.cyp_softspots.length === 0 ? <p className="text-xs text-[var(--color-text-muted)]">No soft-spots</p> : result.safety_assessment.cyp_softspots.map((site) => (
                                          <div key={site.site_name} className="text-xs"><span className="font-medium text-[var(--color-text-primary)]">{site.site_name}</span><span className="text-[var(--color-text-muted)]"> — {site.reaction_type}</span></div>
                                        ))}</div>
                                      </MetricCard>
                                    )}
                                    {result.safety_assessment.herg && (
                                      <MetricCard title="hERG" score={`${result.safety_assessment.herg.risk_score}/${result.safety_assessment.herg.max_score}`} classification={result.safety_assessment.herg.herg_risk.charAt(0).toUpperCase() + result.safety_assessment.herg.herg_risk.slice(1)} classificationVariant={result.safety_assessment.herg.herg_risk === 'low' ? 'success' : result.safety_assessment.herg.herg_risk === 'moderate' ? 'warning' : 'error'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                        <div className="space-y-1">{result.safety_assessment.herg.flags.map((flag, i) => <p key={i} className="text-xs text-[var(--color-text-primary)]">&#8226; {flag}</p>)}</div>
                                      </MetricCard>
                                    )}
                                    {result.safety_assessment.bro5 && (
                                      <MetricCard title="bRo5" score={!result.safety_assessment.bro5.applicable ? 'N/A' : result.safety_assessment.bro5.passed ? 'Pass' : `${result.safety_assessment.bro5.violations.length} viol.`} classification={!result.safety_assessment.bro5.applicable ? 'N/A' : result.safety_assessment.bro5.passed ? 'Pass' : `${result.safety_assessment.bro5.violations.length} viol.`} classificationVariant={!result.safety_assessment.bro5.applicable ? 'default' : result.safety_assessment.bro5.passed ? 'success' : 'error'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                        {!result.safety_assessment.bro5.applicable ? <p className="text-xs text-[var(--color-text-muted)]">MW &le; 500</p> : <div className="space-y-0.5">{result.safety_assessment.bro5.violations.map((v, i) => <div key={i} className="flex justify-between text-xs"><span className="text-[var(--color-text-muted)]">{v.property.replace(/_/g, ' ')}</span><span className="text-[var(--color-text-primary)] font-medium tabular-nums">{v.value.toFixed(1)}</span></div>)}</div>}
                                      </MetricCard>
                                    )}
                                    {result.safety_assessment.reos && (
                                      <MetricCard title="REOS" score={result.safety_assessment.reos.n_violations === 0 ? 'Pass' : `${result.safety_assessment.reos.n_violations} viol.`} classification={result.safety_assessment.reos.passed ? 'Pass' : `${result.safety_assessment.reos.n_violations} viol.`} classificationVariant={result.safety_assessment.reos.n_violations === 0 ? 'success' : result.safety_assessment.reos.n_violations === 1 ? 'warning' : 'error'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                        <div className="space-y-0.5">{result.safety_assessment.reos.violations.map((v, i) => <div key={i} className="flex justify-between text-xs"><span className="text-[var(--color-text-muted)]">{v.property.replace(/_/g, ' ')}</span><span className="text-[var(--color-text-primary)] font-medium tabular-nums">{v.value.toFixed(1)}</span></div>)}</div>
                                      </MetricCard>
                                    )}
                                    {result.safety_assessment.complexity && (
                                      <MetricCard title="Complexity" score={result.safety_assessment.complexity.n_outliers === 0 ? 'OK' : `${result.safety_assessment.complexity.n_outliers} outlier(s)`} classification={result.safety_assessment.complexity.within_range ? 'Pass' : `${result.safety_assessment.complexity.n_outliers} outlier(s)`} classificationVariant={result.safety_assessment.complexity.n_outliers === 0 ? 'success' : 'warning'} defaultExpanded={true} className="!p-3 !rounded-lg">
                                        <div className="space-y-0.5">{Object.entries(result.safety_assessment.complexity.properties).map(([key, prop]) => (
                                          <div key={key} className={cn('flex items-center gap-1.5 text-xs', prop.outlier ? 'text-red-600 dark:text-red-400' : 'text-[var(--color-text-secondary)]')}>
                                            {prop.outlier ? <XCircle className="w-3.5 h-3.5 shrink-0" /> : <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 shrink-0" />}
                                            <span className="flex-1 truncate">{key.replace(/_/g, ' ')}</span>
                                            <span className="tabular-nums shrink-0">{prop.value.toFixed(1)}</span>
                                          </div>
                                        ))}</div>
                                      </MetricCard>
                                    )}
                                  </>
                                )}
                              </div>
                            </div>

                          </div>
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
