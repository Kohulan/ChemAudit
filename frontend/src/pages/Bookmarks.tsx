/**
 * Bookmarks page - lists saved molecules with tag filtering, bulk select,
 * and submit-as-batch workflow.
 */

import { useCallback } from 'react';
import { useNavigate } from 'react-router-dom';
import { motion } from 'framer-motion';
import { Star, Beaker } from 'lucide-react';
import { BookmarkList } from '../components/bookmarks/BookmarkList';
import { Badge } from '../components/ui/Badge';
import { useBookmarks } from '../hooks/useBookmarks';
import { deleteSnapshot } from '../lib/bookmarkStore';
import { logger } from '../lib/logger';
import type { Bookmark } from '../types/workflow';

export function BookmarksPage() {
  const navigate = useNavigate();
  const {
    bookmarks,
    total,
    page,
    pageSize,
    isLoading,
    error,
    selectedIds,
    setParams,
    removeBookmark,
    bulkDelete,
    submitAsBatch,
    toggleSelection,
    selectAll,
    clearSelection,
  } = useBookmarks({ page: 1, page_size: 20 });

  const handlePageChange = useCallback(
    (newPage: number) => {
      setParams({ page: newPage, page_size: pageSize });
    },
    [setParams, pageSize]
  );

  const handleSearch = useCallback(
    (query: string) => {
      setParams({ page: 1, page_size: pageSize, search: query || undefined });
    },
    [setParams, pageSize]
  );

  const handleTagFilter = useCallback(
    (tag: string | undefined) => {
      setParams({ page: 1, page_size: pageSize, tag });
    },
    [setParams, pageSize]
  );

  const handleSubmitAsBatch = useCallback(
    async (ids: number[]) => {
      try {
        const result = await submitAsBatch(ids);
        navigate(`/batch?job=${result.job_id}`);
      } catch (err) {
        logger.error('Failed to submit as batch:', err);
      }
    },
    [submitAsBatch, navigate]
  );

  const handleOpenBookmark = useCallback(
    (bookmark: Bookmark) => {
      if (bookmark.source === 'batch' && bookmark.job_id) {
        navigate(`/batch?job=${bookmark.job_id}`);
      } else {
        navigate('/', { state: { bookmarkId: bookmark.id, smiles: bookmark.smiles } });
      }
    },
    [navigate]
  );

  const handleDelete = useCallback(
    async (id: number) => {
      await removeBookmark(id);
      deleteSnapshot(id).catch(() => {});
    },
    [removeBookmark]
  );

  const handleBulkDelete = useCallback(
    async (ids: number[]) => {
      await bulkDelete(ids);
      for (const id of ids) {
        deleteSnapshot(id).catch(() => {});
      }
    },
    [bulkDelete]
  );

  return (
    <div className="max-w-7xl mx-auto space-y-6 px-4 sm:px-6">
      {/* Header */}
      <motion.div
        className="text-center pt-4"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: [0.25, 0.46, 0.45, 0.94] }}
      >
        <div className="inline-flex items-center gap-2 px-4 py-1.5 rounded-full bg-amber-500/10 border border-amber-500/20 mb-4">
          <Star className="w-4 h-4 text-amber-500" />
          <span className="text-sm font-medium text-amber-600 dark:text-amber-400">Molecule Library</span>
        </div>
        <h1 className="text-3xl sm:text-4xl font-bold text-gradient tracking-tight font-display">
          Bookmarks
        </h1>
        <p className="text-[var(--color-text-secondary)] mt-3 text-base sm:text-lg max-w-2xl mx-auto">
          Your saved molecules for quick access and batch processing
        </p>
      </motion.div>

      {/* Quick stats */}
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.2 }}
        className="flex items-center justify-center gap-6"
      >
        <div className="text-center">
          <div className="text-2xl font-bold text-[var(--color-text-primary)]">{total}</div>
          <div className="text-xs text-[var(--color-text-muted)] uppercase tracking-wider">Total</div>
        </div>
      </motion.div>

      {/* Error */}
      {error && (
        <div className="p-4 bg-red-500/10 border border-red-500/30 rounded-xl text-sm text-red-600 dark:text-red-400">
          {error}
        </div>
      )}

      {/* Content */}
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.3 }}
      >
        {total === 0 && !isLoading ? (
          /* Empty state */
          <div className="card p-12 text-center">
            <div className="w-16 h-16 mx-auto mb-4 rounded-2xl bg-amber-500/10 flex items-center justify-center">
              <Beaker className="w-8 h-8 text-amber-500/50" />
            </div>
            <h3 className="text-lg font-semibold text-[var(--color-text-primary)] mb-2 font-display">
              No Bookmarks Yet
            </h3>
            <p className="text-sm text-[var(--color-text-muted)] max-w-md mx-auto">
              Bookmark molecules during validation to build your collection.
              Click the <Badge variant="default" size="sm" icon={<Star className="w-3 h-3" />}>star</Badge> button
              on any validation result to save it here.
            </p>
          </div>
        ) : (
          <div className="card p-6">
            <BookmarkList
              bookmarks={bookmarks}
              total={total}
              page={page}
              pageSize={pageSize}
              isLoading={isLoading}
              selectedIds={selectedIds}
              onToggleSelection={toggleSelection}
              onSelectAll={selectAll}
              onClearSelection={clearSelection}
              onPageChange={handlePageChange}
              onSearch={handleSearch}
              onTagFilter={handleTagFilter}
              onDelete={handleDelete}
              onBulkDelete={handleBulkDelete}
              onSubmitAsBatch={handleSubmitAsBatch}
              onOpenBookmark={handleOpenBookmark}
            />
          </div>
        )}
      </motion.div>
    </div>
  );
}
