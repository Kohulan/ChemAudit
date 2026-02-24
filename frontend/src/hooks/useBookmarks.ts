/**
 * Hook for managing molecule bookmarks with pagination and filtering.
 *
 * Provides CRUD operations, bulk actions, selection state, and
 * submit-as-batch workflow.
 */

import { useState, useEffect, useCallback, useRef } from 'react';
import { bookmarksApi } from '../services/api';
import type { Bookmark, BookmarkCreate, BookmarkUpdate } from '../types/workflow';

export interface BookmarkParams {
  page?: number;
  page_size?: number;
  tag?: string;
  search?: string;
}

export interface UseBookmarksReturn {
  bookmarks: Bookmark[];
  total: number;
  page: number;
  pageSize: number;
  isLoading: boolean;
  error: string | null;
  selectedIds: Set<number>;
  refetch: () => Promise<void>;
  setParams: (params: BookmarkParams) => void;
  addBookmark: (data: BookmarkCreate) => Promise<Bookmark>;
  updateBookmark: (id: number, data: BookmarkUpdate) => Promise<Bookmark>;
  removeBookmark: (id: number) => Promise<void>;
  bulkDelete: (ids: number[]) => Promise<void>;
  submitAsBatch: (ids: number[]) => Promise<{ job_id: string; molecule_count: number }>;
  toggleSelection: (id: number) => void;
  selectAll: () => void;
  clearSelection: () => void;
  isBookmarked: (smiles: string) => boolean;
}

export function useBookmarks(initialParams?: BookmarkParams): UseBookmarksReturn {
  const [bookmarks, setBookmarks] = useState<Bookmark[]>([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(initialParams?.page ?? 1);
  const [pageSize] = useState(initialParams?.page_size ?? 20);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
  const [params, setParamsState] = useState<BookmarkParams>(initialParams ?? {});
  const mountedRef = useRef(true);

  const fetchBookmarks = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    try {
      const result = await bookmarksApi.getBookmarks({
        page: params.page ?? page,
        page_size: params.page_size ?? pageSize,
        tag: params.tag,
        search: params.search,
      });
      if (mountedRef.current) {
        setBookmarks(result.bookmarks);
        setTotal(result.total);
        setPage(result.page);
      }
    } catch (err) {
      if (mountedRef.current) {
        setError(err instanceof Error ? err.message : 'Failed to load bookmarks');
      }
    } finally {
      if (mountedRef.current) {
        setIsLoading(false);
      }
    }
  }, [params, page, pageSize]);

  useEffect(() => {
    mountedRef.current = true;
    fetchBookmarks();
    return () => {
      mountedRef.current = false;
    };
  }, [fetchBookmarks]);

  const setParams = useCallback((newParams: BookmarkParams) => {
    setParamsState(newParams);
    setSelectedIds(new Set());
  }, []);

  const addBookmark = useCallback(
    async (data: BookmarkCreate) => {
      const created = await bookmarksApi.createBookmark(data);
      await fetchBookmarks();
      return created;
    },
    [fetchBookmarks]
  );

  const updateBookmark = useCallback(
    async (id: number, data: BookmarkUpdate) => {
      const updated = await bookmarksApi.updateBookmark(id, data);
      await fetchBookmarks();
      return updated;
    },
    [fetchBookmarks]
  );

  const removeBookmark = useCallback(
    async (id: number) => {
      await bookmarksApi.deleteBookmark(id);
      setSelectedIds((prev) => {
        const next = new Set(prev);
        next.delete(id);
        return next;
      });
      await fetchBookmarks();
    },
    [fetchBookmarks]
  );

  const bulkDelete = useCallback(
    async (ids: number[]) => {
      await bookmarksApi.bulkDeleteBookmarks(ids);
      setSelectedIds(new Set());
      await fetchBookmarks();
    },
    [fetchBookmarks]
  );

  const submitAsBatch = useCallback(
    async (ids: number[]) => {
      return bookmarksApi.submitBookmarksAsBatch(ids);
    },
    []
  );

  const toggleSelection = useCallback((id: number) => {
    setSelectedIds((prev) => {
      const next = new Set(prev);
      if (next.has(id)) {
        next.delete(id);
      } else {
        next.add(id);
      }
      return next;
    });
  }, []);

  const selectAll = useCallback(() => {
    setSelectedIds(new Set(bookmarks.map((b) => b.id)));
  }, [bookmarks]);

  const clearSelection = useCallback(() => {
    setSelectedIds(new Set());
  }, []);

  const isBookmarked = useCallback(
    (smiles: string) => {
      return bookmarks.some((b) => b.smiles === smiles);
    },
    [bookmarks]
  );

  return {
    bookmarks,
    total,
    page,
    pageSize,
    isLoading,
    error,
    selectedIds,
    refetch: fetchBookmarks,
    setParams,
    addBookmark,
    updateBookmark,
    removeBookmark,
    bulkDelete,
    submitAsBatch,
    toggleSelection,
    selectAll,
    clearSelection,
    isBookmarked,
  };
}
