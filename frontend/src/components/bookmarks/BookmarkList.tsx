/**
 * Reusable bookmark list with checkboxes, tag filtering, search, and bulk actions.
 */

import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Search, Trash2, Play, ChevronDown, ChevronUp, Tag, ChevronLeft, ChevronRight,
  Eye,
} from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { CopyButton } from '../ui/CopyButton';
import { cn } from '../../lib/utils';
import { getSnapshotIds } from '../../lib/bookmarkStore';
import type { Bookmark } from '../../types/workflow';

interface BookmarkListProps {
  bookmarks: Bookmark[];
  total: number;
  page: number;
  pageSize: number;
  isLoading: boolean;
  selectedIds: Set<number>;
  onToggleSelection: (id: number) => void;
  onSelectAll: () => void;
  onClearSelection: () => void;
  onPageChange: (page: number) => void;
  onSearch: (query: string) => void;
  onTagFilter: (tag: string | undefined) => void;
  onDelete: (id: number) => void;
  onBulkDelete: (ids: number[]) => void;
  onSubmitAsBatch: (ids: number[]) => void;
  onOpenBookmark?: (bookmark: Bookmark) => void;
}

export function BookmarkList({
  bookmarks,
  total,
  page,
  pageSize,
  isLoading,
  selectedIds,
  onToggleSelection,
  onSelectAll,
  onClearSelection,
  onPageChange,
  onSearch,
  onTagFilter,
  onDelete,
  onBulkDelete,
  onSubmitAsBatch,
  onOpenBookmark,
}: BookmarkListProps) {
  const [searchQuery, setSearchQuery] = useState('');
  const [expandedRow, setExpandedRow] = useState<number | null>(null);
  const [activeTag, setActiveTag] = useState<string | undefined>();
  const [snapshotIds, setSnapshotIds] = useState<Set<number>>(new Set());

  const items = bookmarks ?? [];
  const totalPages = Math.ceil(total / pageSize);
  const allSelected = items.length > 0 && items.every((b) => selectedIds.has(b.id));

  // Load snapshot IDs to show indicator badges
  useEffect(() => {
    getSnapshotIds().then(setSnapshotIds).catch(() => {});
  }, [bookmarks]);

  // Collect unique tags across all visible bookmarks
  const allTags = Array.from(new Set(items.flatMap((b) => b.tags)));

  const handleSearchSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onSearch(searchQuery);
  };

  const handleTagClick = (tag: string) => {
    if (activeTag === tag) {
      setActiveTag(undefined);
      onTagFilter(undefined);
    } else {
      setActiveTag(tag);
      onTagFilter(tag);
    }
  };

  const truncateSmiles = (smiles: string, maxLen: number = 40): string => {
    if (smiles.length <= maxLen) return smiles;
    return smiles.substring(0, maxLen) + '...';
  };

  return (
    <div className="relative space-y-4">
      {/* Search and tag filter bar */}
      <div className="flex flex-wrap items-center gap-3">
        <form onSubmit={handleSearchSubmit} className="flex-1 min-w-[200px]">
          <div className="relative">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-[var(--color-text-muted)]" />
            <input
              type="text"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              placeholder="Search SMILES..."
              className={cn(
                'w-full pl-9 pr-3 py-2 rounded-lg text-sm',
                'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                'text-[var(--color-text-primary)] placeholder:text-[var(--color-text-muted)]',
                'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30'
              )}
            />
          </div>
        </form>

        {allTags.length > 0 && (
          <div className="flex flex-wrap items-center gap-1.5">
            <Tag className="w-3.5 h-3.5 text-[var(--color-text-muted)]" />
            {allTags.map((tag) => (
              <button
                key={tag}
                onClick={() => handleTagClick(tag)}
                className={cn(
                  'px-2.5 py-1 text-xs rounded-full transition-all',
                  activeTag === tag
                    ? 'bg-[var(--color-primary)] text-white'
                    : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] hover:bg-[var(--color-surface-elevated)]'
                )}
              >
                {tag}
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Bulk action bar */}
      <AnimatePresence>
        {selectedIds.size > 0 && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: 'auto' }}
            exit={{ opacity: 0, height: 0 }}
            className="flex items-center justify-between p-3 bg-[var(--color-primary)]/5 border border-[var(--color-primary)]/20 rounded-xl"
          >
            <span className="text-sm font-medium text-[var(--color-text-primary)]">
              {selectedIds.size} selected
            </span>
            <div className="flex items-center gap-2">
              <ClayButton
                size="sm"
                variant="primary"
                onClick={() => onSubmitAsBatch(Array.from(selectedIds))}
                leftIcon={<Play className="w-3.5 h-3.5" />}
              >
                Submit as Batch
              </ClayButton>
              <ClayButton
                size="sm"
                onClick={() => onBulkDelete(Array.from(selectedIds))}
                leftIcon={<Trash2 className="w-3.5 h-3.5" />}
              >
                Delete Selected
              </ClayButton>
              <button
                onClick={onClearSelection}
                className="text-xs text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)] ml-2"
              >
                Clear
              </button>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Table */}
      <div className="overflow-x-auto rounded-xl border border-[var(--color-border)]">
        <table className="w-full text-sm">
          <thead>
            <tr className="bg-[var(--color-surface-sunken)]">
              <th className="w-10 px-3 py-2.5">
                <input
                  type="checkbox"
                  checked={allSelected}
                  onChange={() => (allSelected ? onClearSelection() : onSelectAll())}
                  className="rounded border-[var(--color-border-strong)] text-[var(--color-primary)]"
                />
              </th>
              <th className="text-left px-3 py-2.5 font-medium text-[var(--color-text-secondary)]">Molecule</th>
              <th className="text-left px-3 py-2.5 font-medium text-[var(--color-text-secondary)] hidden md:table-cell">InChIKey</th>
              <th className="text-left px-3 py-2.5 font-medium text-[var(--color-text-secondary)] hidden sm:table-cell">Tags</th>
              <th className="text-left px-3 py-2.5 font-medium text-[var(--color-text-secondary)] hidden lg:table-cell">Source</th>
              <th className="text-right px-3 py-2.5 font-medium text-[var(--color-text-secondary)]">Actions</th>
            </tr>
          </thead>
          <tbody>
            {isLoading ? (
              <tr>
                <td colSpan={6} className="px-3 py-8 text-center text-[var(--color-text-muted)]">
                  Loading bookmarks...
                </td>
              </tr>
            ) : items.length === 0 ? (
              <tr>
                <td colSpan={6} className="px-3 py-8 text-center text-[var(--color-text-muted)]">
                  No bookmarks found
                </td>
              </tr>
            ) : (
              items.map((bookmark) => (
                <tr
                  key={bookmark.id}
                  className={cn(
                    'border-t border-[var(--color-border)] transition-colors',
                    selectedIds.has(bookmark.id) && 'bg-[var(--color-primary)]/5',
                    'hover:bg-[var(--color-surface-sunken)]/50'
                  )}
                >
                  <td className="px-3 py-2.5">
                    <input
                      type="checkbox"
                      checked={selectedIds.has(bookmark.id)}
                      onChange={() => onToggleSelection(bookmark.id)}
                      className="rounded border-[var(--color-border-strong)] text-[var(--color-primary)]"
                    />
                  </td>
                  <td className="relative px-3 py-2.5">
                    <button
                      onClick={() => setExpandedRow(expandedRow === bookmark.id ? null : bookmark.id)}
                      className="flex items-center gap-1.5 text-left"
                    >
                      {expandedRow === bookmark.id ? (
                        <ChevronUp className="w-3.5 h-3.5 text-[var(--color-text-muted)] flex-shrink-0" />
                      ) : (
                        <ChevronDown className="w-3.5 h-3.5 text-[var(--color-text-muted)] flex-shrink-0" />
                      )}
                      <div>
                        {bookmark.name && (
                          <div className="font-medium text-[var(--color-text-primary)] text-xs">
                            {bookmark.name}
                          </div>
                        )}
                        <code className="text-xs text-[var(--color-text-secondary)] font-mono">
                          {truncateSmiles(bookmark.smiles)}
                        </code>
                      </div>
                    </button>
                    {/* Snapshot indicator */}
                    {snapshotIds.has(bookmark.id) && (
                      <span
                        className="inline-flex ml-2 px-1.5 py-0.5 text-[9px] font-medium rounded bg-green-500/10 text-green-600 dark:text-green-400 uppercase tracking-wider"
                        title="Full validation results saved"
                      >
                        Results saved
                      </span>
                    )}
                    {/* Expanded row */}
                    <AnimatePresence>
                      {expandedRow === bookmark.id && (
                        <motion.div
                          initial={{ opacity: 0, height: 0 }}
                          animate={{ opacity: 1, height: 'auto' }}
                          exit={{ opacity: 0, height: 0 }}
                          className="mt-2 p-2 bg-[var(--color-surface-sunken)] rounded-lg"
                        >
                          <div className="flex items-center gap-2 mb-1">
                            <code className="text-[10px] font-mono text-[var(--color-text-secondary)] break-all">
                              {bookmark.smiles}
                            </code>
                            <CopyButton text={bookmark.smiles} />
                          </div>
                          {bookmark.notes && (
                            <p className="text-xs text-[var(--color-text-muted)] mt-1">{bookmark.notes}</p>
                          )}
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </td>
                  <td className="px-3 py-2.5 hidden md:table-cell">
                    <code className="text-[10px] text-[var(--color-text-muted)] font-mono">
                      {bookmark.inchikey ? bookmark.inchikey.slice(0, 14) + '...' : '--'}
                    </code>
                  </td>
                  <td className="px-3 py-2.5 hidden sm:table-cell">
                    <div className="flex flex-wrap gap-1">
                      {bookmark.tags.slice(0, 3).map((tag) => (
                        <Badge key={tag} variant="default" size="sm">{tag}</Badge>
                      ))}
                      {bookmark.tags.length > 3 && (
                        <span className="text-[10px] text-[var(--color-text-muted)]">
                          +{bookmark.tags.length - 3}
                        </span>
                      )}
                    </div>
                  </td>
                  <td className="px-3 py-2.5 hidden lg:table-cell">
                    <span className="text-xs text-[var(--color-text-muted)]">{bookmark.source ?? '--'}</span>
                  </td>
                  <td className="px-3 py-2.5 text-right">
                    <div className="flex items-center justify-end gap-1">
                      {onOpenBookmark && (
                        <button
                          onClick={() => onOpenBookmark(bookmark)}
                          className={cn(
                            'inline-flex items-center gap-1 px-2.5 py-1 rounded-lg text-xs font-medium transition-colors cursor-pointer',
                            'bg-[var(--color-primary)]/10 text-[var(--color-primary)]',
                            'hover:bg-[var(--color-primary)]/20'
                          )}
                          title={snapshotIds.has(bookmark.id) ? 'View saved results' : 'View and validate'}
                        >
                          <Eye className="w-3 h-3" />
                          View
                        </button>
                      )}
                      <button
                        onClick={() => onDelete(bookmark.id)}
                        className="p-1.5 rounded-lg hover:bg-red-500/10 text-[var(--color-text-muted)] hover:text-red-500 transition-colors"
                        title="Delete bookmark"
                      >
                        <Trash2 className="w-3.5 h-3.5" />
                      </button>
                    </div>
                  </td>
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>

      {/* Pagination */}
      {totalPages > 1 && (
        <div className="flex items-center justify-between">
          <span className="text-xs text-[var(--color-text-muted)]">
            Showing {(page - 1) * pageSize + 1}-{Math.min(page * pageSize, total)} of {total}
          </span>
          <div className="flex items-center gap-2">
            <button
              onClick={() => onPageChange(page - 1)}
              disabled={page <= 1}
              className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              <ChevronLeft className="w-4 h-4 text-[var(--color-text-secondary)]" />
            </button>
            <span className="text-sm text-[var(--color-text-primary)] font-medium tabular-nums">
              {page} / {totalPages}
            </span>
            <button
              onClick={() => onPageChange(page + 1)}
              disabled={page >= totalPages}
              className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
            >
              <ChevronRight className="w-4 h-4 text-[var(--color-text-secondary)]" />
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
