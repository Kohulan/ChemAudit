/**
 * ClusterTable
 *
 * Sortable, paginated cluster table with expandable rows showing
 * cluster member grids. Cluster IDs are colored badges.
 */

import { useState, useMemo, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, ChevronRight } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { ClusterMemberGrid } from './ClusterMemberGrid';
import type { ClusterInfo } from '../../types/analytics';
import type { BatchResult } from '../../types/batch';

interface ClusterTableProps {
  clusters: ClusterInfo[];
  results: BatchResult[];
  expandedClusterId: number | null;
  onToggleExpand: (id: number) => void;
}

const CLUSTER_COLORS = [
  '#c41e3a',
  '#d97706',
  '#059669',
  '#2563eb',
  '#7c3aed',
  '#db2777',
  '#0891b2',
];
const FALLBACK_COLOR = '#6b7280';
const PAGE_SIZE = 20;

type SortKey = 'cluster_id' | 'size';
type SortDir = 'asc' | 'desc';

function getClusterColor(index: number): string {
  return index < CLUSTER_COLORS.length ? CLUSTER_COLORS[index] : FALLBACK_COLOR;
}

export function ClusterTable({
  clusters,
  results,
  expandedClusterId,
  onToggleExpand,
}: ClusterTableProps) {
  const [sortKey, setSortKey] = useState<SortKey>('size');
  const [sortDir, setSortDir] = useState<SortDir>('desc');
  const [currentPage, setCurrentPage] = useState(1);

  // Sort clusters
  const sorted = useMemo(() => {
    const arr = [...clusters];
    arr.sort((a, b) => {
      const aVal = sortKey === 'cluster_id' ? a.cluster_id : a.size;
      const bVal = sortKey === 'cluster_id' ? b.cluster_id : b.size;
      return sortDir === 'asc' ? aVal - bVal : bVal - aVal;
    });
    return arr;
  }, [clusters, sortKey, sortDir]);

  // Paginate
  const totalPages = Math.ceil(sorted.length / PAGE_SIZE);
  const paged = useMemo(
    () => sorted.slice((currentPage - 1) * PAGE_SIZE, currentPage * PAGE_SIZE),
    [sorted, currentPage]
  );

  // Build index->result lookup
  const resultMap = useMemo(() => {
    const map = new Map<number, BatchResult>();
    for (const r of results) {
      map.set(r.index, r);
    }
    return map;
  }, [results]);

  const handleSort = useCallback(
    (key: SortKey) => {
      if (sortKey === key) {
        setSortDir((prev) => (prev === 'asc' ? 'desc' : 'asc'));
      } else {
        setSortKey(key);
        setSortDir(key === 'size' ? 'desc' : 'asc');
      }
      setCurrentPage(1);
    },
    [sortKey]
  );

  const ariaSort = (key: SortKey): 'ascending' | 'descending' | 'none' => {
    if (sortKey !== key) return 'none';
    return sortDir === 'asc' ? 'ascending' : 'descending';
  };

  return (
    <div className="space-y-3">
      <div className="overflow-x-auto rounded-xl border border-[var(--color-border)]">
        <table className="w-full text-sm" role="table">
          <thead>
            <tr className="bg-[var(--color-surface-sunken)]">
              <th
                role="columnheader"
                aria-sort={ariaSort('cluster_id')}
                onClick={() => handleSort('cluster_id')}
                className="px-3 py-2 text-left text-xs font-medium text-[var(--color-text-muted)] cursor-pointer hover:text-[var(--color-text-primary)] select-none"
              >
                Cluster ID {sortKey === 'cluster_id' && (sortDir === 'asc' ? '\u25B2' : '\u25BC')}
              </th>
              <th
                role="columnheader"
                aria-sort={ariaSort('size')}
                onClick={() => handleSort('size')}
                className="px-3 py-2 text-left text-xs font-medium text-[var(--color-text-muted)] cursor-pointer hover:text-[var(--color-text-primary)] select-none"
              >
                Size {sortKey === 'size' && (sortDir === 'asc' ? '\u25B2' : '\u25BC')}
              </th>
              <th className="px-3 py-2 text-left text-xs font-medium text-[var(--color-text-muted)]">
                Representative
              </th>
              <th className="px-3 py-2 text-center text-xs font-medium text-[var(--color-text-muted)]">
                Expand
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-[var(--color-border)]">
            {paged.map((cluster) => {
              const isExpanded = expandedClusterId === cluster.cluster_id;
              const repMol = resultMap.get(cluster.representative_index);
              const repSmiles =
                repMol?.standardization?.standardized_smiles || repMol?.smiles || null;
              const clusterColorIndex =
                clusters.findIndex((c) => c.cluster_id === cluster.cluster_id);
              const color = getClusterColor(clusterColorIndex);

              return (
                <tr key={cluster.cluster_id} className="group">
                  <td className="px-3 py-2">
                    <span className="flex flex-col gap-1">
                      <span
                        className="inline-flex items-center justify-center w-7 h-7 rounded-full text-xs font-bold text-white"
                        style={{ backgroundColor: color }}
                      >
                        {cluster.cluster_id}
                      </span>
                    </span>
                  </td>
                  <td className="px-3 py-2 font-mono tabular-nums text-[var(--color-text-primary)]">
                    {cluster.size}
                  </td>
                  <td className="px-3 py-2">
                    <div className="flex items-center gap-3">
                      <div
                        className="rounded-lg border border-[var(--color-border)] bg-white dark:bg-gray-900/50 overflow-hidden flex-shrink-0"
                        style={{ width: 120, height: 100 }}
                      >
                        <MoleculeViewer smiles={repSmiles} width={120} height={100} />
                      </div>
                      <span
                        className="text-xs font-mono text-[var(--color-text-muted)] truncate max-w-[200px]"
                        title={repSmiles || ''}
                      >
                        {repSmiles || '-'}
                      </span>
                    </div>
                  </td>
                  <td className="px-3 py-2 text-center">
                    <button
                      onClick={() => onToggleExpand(cluster.cluster_id)}
                      aria-expanded={isExpanded}
                      aria-label={`Expand cluster ${cluster.cluster_id} to show ${cluster.size} members`}
                      className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)]"
                    >
                      {isExpanded ? (
                        <ChevronDown className="w-4 h-4" />
                      ) : (
                        <ChevronRight className="w-4 h-4" />
                      )}
                    </button>
                    {/* Expanded member grid */}
                    <AnimatePresence>
                      {isExpanded && (
                        <motion.div
                          initial={{ height: 0, opacity: 0 }}
                          animate={{ height: 'auto', opacity: 1 }}
                          exit={{ height: 0, opacity: 0 }}
                          transition={{ duration: 0.25, ease: 'easeOut' }}
                          className="overflow-hidden"
                        >
                          <div className="pt-2">
                            <ClusterMemberGrid
                              memberIndices={cluster.member_indices}
                              results={results}
                              clusterId={cluster.cluster_id}
                            />
                          </div>
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Pagination */}
      {totalPages > 1 && (
        <div className="flex items-center justify-between text-xs text-[var(--color-text-muted)]">
          <span>
            Page {currentPage} of {totalPages} ({clusters.length} clusters)
          </span>
          <div className="flex gap-2">
            <button
              onClick={() => setCurrentPage((p) => Math.max(1, p - 1))}
              disabled={currentPage <= 1}
              className="px-3 py-1 rounded-lg border border-[var(--color-border)] hover:bg-[var(--color-surface-sunken)] disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
            >
              Previous
            </button>
            <button
              onClick={() => setCurrentPage((p) => Math.min(totalPages, p + 1))}
              disabled={currentPage >= totalPages}
              className="px-3 py-1 rounded-lg border border-[var(--color-border)] hover:bg-[var(--color-surface-sunken)] disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
            >
              Next
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
