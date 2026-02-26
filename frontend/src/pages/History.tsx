/**
 * Validation History page - paginated audit trail with filtering
 * by date range, outcome, source, and SMILES search.
 */

import { useState, useEffect, useCallback } from 'react';
import { motion } from 'framer-motion';
import {
  Clock, ChevronLeft, ChevronRight, Search, Filter,
} from 'lucide-react';
import { Badge } from '../components/ui/Badge';
import { ClayButton } from '../components/ui/ClayButton';
import { historyApi } from '../services/api';
import { cn } from '../lib/utils';
import type { AuditEntry, AuditHistoryStats } from '../types/workflow';

export function HistoryPage() {
  const [entries, setEntries] = useState<AuditEntry[]>([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [pageSize] = useState(25);
  const [isLoading, setIsLoading] = useState(true);
  const [stats, setStats] = useState<AuditHistoryStats | null>(null);

  // Filters
  const [dateFrom, setDateFrom] = useState('');
  const [dateTo, setDateTo] = useState('');
  const [outcome, setOutcome] = useState('');
  const [source, setSource] = useState('');
  const [smilesSearch, setSmilesSearch] = useState('');

  const fetchHistory = useCallback(async () => {
    setIsLoading(true);
    try {
      const data = await historyApi.getHistory({
        page,
        page_size: pageSize,
        date_from: dateFrom || undefined,
        date_to: dateTo || undefined,
        outcome: outcome || undefined,
        source: source || undefined,
        smiles_search: smilesSearch || undefined,
      });
      setEntries(data.entries);
      setTotal(data.total);
    } catch {
      // handle error silently
    } finally {
      setIsLoading(false);
    }
  }, [page, pageSize, dateFrom, dateTo, outcome, source, smilesSearch]);

  const fetchStats = useCallback(async () => {
    try {
      const data = await historyApi.getHistoryStats();
      setStats(data);
    } catch {
      // ignore
    }
  }, []);

  useEffect(() => {
    fetchHistory();
  }, [fetchHistory]);

  useEffect(() => {
    fetchStats();
  }, [fetchStats]);

  const totalPages = Math.ceil(total / pageSize);

  const handleFilterApply = () => {
    setPage(1);
    fetchHistory();
  };

  const handleReset = () => {
    setDateFrom('');
    setDateTo('');
    setOutcome('');
    setSource('');
    setSmilesSearch('');
    setPage(1);
  };

  const getOutcomeBadge = (o: string) => {
    switch (o) {
      case 'pass':
        return <Badge variant="success" size="sm">Pass</Badge>;
      case 'warn':
        return <Badge variant="warning" size="sm">Warn</Badge>;
      case 'fail':
        return <Badge variant="error" size="sm">Fail</Badge>;
      default:
        return <Badge variant="default" size="sm">{o}</Badge>;
    }
  };

  const formatDate = (dateStr: string) => {
    const d = new Date(dateStr);
    return d.toLocaleDateString('en-US', {
      month: 'short',
      day: 'numeric',
      year: 'numeric',
      hour: '2-digit',
      minute: '2-digit',
    });
  };

  const truncateSmiles = (smiles: string, maxLen: number = 35): string => {
    if (smiles.length <= maxLen) return smiles;
    return smiles.substring(0, maxLen) + '...';
  };

  return (
    <div className="max-w-7xl mx-auto space-y-6 px-4 sm:px-6">
      {/* Header */}
      <motion.div
        className="text-center pt-4"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: [0.25, 0.46, 0.45, 0.94] }}
      >
        <div className="inline-flex items-center gap-2 px-4 py-1.5 rounded-full bg-[var(--color-primary)]/10 border border-[var(--color-primary)]/20 mb-4">
          <Clock className="w-4 h-4 text-[var(--color-primary)]" />
          <span className="text-sm font-medium text-[var(--color-primary)]">Audit Trail</span>
        </div>
        <h1 className="text-3xl sm:text-4xl font-bold text-gradient tracking-tight font-display">
          Validation History
        </h1>
        <p className="text-[var(--color-text-secondary)] mt-3 text-base sm:text-lg max-w-2xl mx-auto">
          Browse past validations with filtering and search
        </p>
      </motion.div>

      {/* Stats summary */}
      {stats && (
        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.2 }}
          className="grid grid-cols-2 sm:grid-cols-4 gap-4"
        >
          <div className="card p-4 text-center">
            <div className="text-2xl font-bold text-[var(--color-text-primary)]">{stats.total_validations}</div>
            <div className="text-xs text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Total</div>
          </div>
          <div className="card p-4 text-center">
            <div className="text-2xl font-bold text-green-500">{stats.outcome_distribution?.pass ?? 0}</div>
            <div className="text-xs text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Passed</div>
          </div>
          <div className="card p-4 text-center">
            <div className="text-2xl font-bold text-amber-500">{stats.outcome_distribution?.warn ?? 0}</div>
            <div className="text-xs text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Warnings</div>
          </div>
          <div className="card p-4 text-center">
            <div className="text-2xl font-bold text-red-500">{stats.outcome_distribution?.fail ?? 0}</div>
            <div className="text-xs text-[var(--color-text-muted)] uppercase tracking-wider mt-1">Failed</div>
          </div>
        </motion.div>
      )}

      {/* Filters */}
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.3 }}
        className="card p-4"
      >
        <div className="flex items-center gap-2 mb-3">
          <Filter className="w-4 h-4 text-[var(--color-text-muted)]" />
          <span className="text-sm font-medium text-[var(--color-text-secondary)]">Filters</span>
        </div>
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-5 gap-3">
          <div>
            <label className="block text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mb-1">From</label>
            <input
              type="date"
              value={dateFrom}
              onChange={(e) => setDateFrom(e.target.value)}
              className={cn(
                'w-full px-2.5 py-1.5 rounded-lg text-sm',
                'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                'text-[var(--color-text-primary)]',
                'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
              )}
            />
          </div>
          <div>
            <label className="block text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mb-1">To</label>
            <input
              type="date"
              value={dateTo}
              onChange={(e) => setDateTo(e.target.value)}
              className={cn(
                'w-full px-2.5 py-1.5 rounded-lg text-sm',
                'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                'text-[var(--color-text-primary)]',
                'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
              )}
            />
          </div>
          <div>
            <label className="block text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mb-1">Outcome</label>
            <select
              value={outcome}
              onChange={(e) => setOutcome(e.target.value)}
              className={cn(
                'w-full px-2.5 py-1.5 rounded-lg text-sm',
                'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                'text-[var(--color-text-primary)]',
                'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
              )}
            >
              <option value="">All</option>
              <option value="pass">Pass</option>
              <option value="warn">Warn</option>
              <option value="fail">Fail</option>
            </select>
          </div>
          <div>
            <label className="block text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mb-1">Source</label>
            <select
              value={source}
              onChange={(e) => setSource(e.target.value)}
              className={cn(
                'w-full px-2.5 py-1.5 rounded-lg text-sm',
                'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                'text-[var(--color-text-primary)]',
                'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
              )}
            >
              <option value="">All</option>
              <option value="single">Single</option>
              <option value="batch">Batch</option>
            </select>
          </div>
          <div>
            <label className="block text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider mb-1">SMILES Search</label>
            <div className="relative">
              <Search className="absolute left-2.5 top-1/2 -translate-y-1/2 w-3.5 h-3.5 text-[var(--color-text-muted)]" />
              <input
                type="text"
                value={smilesSearch}
                onChange={(e) => setSmilesSearch(e.target.value)}
                placeholder="Search..."
                className={cn(
                  'w-full pl-8 pr-2.5 py-1.5 rounded-lg text-sm',
                  'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                  'text-[var(--color-text-primary)] placeholder:text-[var(--color-text-muted)]',
                  'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
                )}
              />
            </div>
          </div>
        </div>
        <div className="flex items-center gap-2 mt-3">
          <ClayButton size="sm" variant="primary" onClick={handleFilterApply}>
            Apply Filters
          </ClayButton>
          <ClayButton size="sm" onClick={handleReset}>
            Reset
          </ClayButton>
        </div>
      </motion.div>

      {/* Results table */}
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.4 }}
        className="card overflow-hidden"
      >
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-[var(--color-surface-sunken)]">
                <th className="text-left px-4 py-3 font-medium text-[var(--color-text-secondary)]">Date</th>
                <th className="text-left px-4 py-3 font-medium text-[var(--color-text-secondary)]">SMILES</th>
                <th className="text-center px-4 py-3 font-medium text-[var(--color-text-secondary)]">Outcome</th>
                <th className="text-center px-4 py-3 font-medium text-[var(--color-text-secondary)]">Score</th>
                <th className="text-center px-4 py-3 font-medium text-[var(--color-text-secondary)]">Source</th>
                <th className="text-left px-4 py-3 font-medium text-[var(--color-text-secondary)] hidden lg:table-cell">Job ID</th>
              </tr>
            </thead>
            <tbody>
              {isLoading ? (
                <tr>
                  <td colSpan={6} className="px-4 py-8 text-center text-[var(--color-text-muted)]">
                    Loading history...
                  </td>
                </tr>
              ) : entries.length === 0 ? (
                <tr>
                  <td colSpan={6} className="px-4 py-8 text-center text-[var(--color-text-muted)]">
                    No entries match the current filters
                  </td>
                </tr>
              ) : (
                entries.map((entry) => (
                  <tr
                    key={entry.id}
                    className="border-t border-[var(--color-border)] hover:bg-[var(--color-surface-sunken)]/50 transition-colors"
                  >
                    <td className="px-4 py-2.5 text-xs text-[var(--color-text-secondary)] whitespace-nowrap">
                      {formatDate(entry.created_at)}
                    </td>
                    <td className="px-4 py-2.5">
                      <code className="text-xs text-[var(--color-text-primary)] font-mono">
                        {truncateSmiles(entry.smiles)}
                      </code>
                    </td>
                    <td className="px-4 py-2.5 text-center">
                      {getOutcomeBadge(entry.outcome)}
                    </td>
                    <td className="px-4 py-2.5 text-center">
                      <span className={cn(
                        'text-sm font-medium tabular-nums',
                        entry.score !== null && entry.score >= 70
                          ? 'text-green-500'
                          : entry.score !== null && entry.score >= 40
                          ? 'text-amber-500'
                          : 'text-red-500'
                      )}>
                        {entry.score ?? '--'}
                      </span>
                    </td>
                    <td className="px-4 py-2.5 text-center">
                      <Badge variant={entry.source === 'batch' ? 'info' : 'default'} size="sm">
                        {entry.source}
                      </Badge>
                    </td>
                    <td className="px-4 py-2.5 hidden lg:table-cell">
                      {entry.job_id ? (
                        <code className="text-[10px] text-[var(--color-text-muted)] font-mono">
                          {entry.job_id.slice(0, 12)}...
                        </code>
                      ) : (
                        <span className="text-[var(--color-text-muted)]">--</span>
                      )}
                    </td>
                  </tr>
                ))
              )}
            </tbody>
          </table>
        </div>

        {/* Pagination */}
        {totalPages > 1 && (
          <div className="flex items-center justify-between px-4 py-3 border-t border-[var(--color-border)]">
            <span className="text-xs text-[var(--color-text-muted)]">
              Showing {(page - 1) * pageSize + 1}-{Math.min(page * pageSize, total)} of {total}
            </span>
            <div className="flex items-center gap-2">
              <button
                onClick={() => setPage(Math.max(1, page - 1))}
                disabled={page <= 1}
                className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
              >
                <ChevronLeft className="w-4 h-4 text-[var(--color-text-secondary)]" />
              </button>
              <span className="text-sm text-[var(--color-text-primary)] font-medium tabular-nums">
                {page} / {totalPages}
              </span>
              <button
                onClick={() => setPage(Math.min(totalPages, page + 1))}
                disabled={page >= totalPages}
                className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
              >
                <ChevronRight className="w-4 h-4 text-[var(--color-text-secondary)]" />
              </button>
            </div>
          </div>
        )}
      </motion.div>
    </div>
  );
}
