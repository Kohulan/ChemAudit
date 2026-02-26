import { createContext, useContext, useRef, useCallback, type ReactNode } from 'react';
import type {
  BatchPageState,
  BatchResultsResponse,
  BatchResultsFilters,
  SortField,
} from '../types/batch';

export interface BatchCache {
  pageState: BatchPageState;
  jobId: string | null;
  resultsData: BatchResultsResponse | null;
  page: number;
  pageSize: number;
  filters: BatchResultsFilters;
  sortBy: SortField;
  sortDir: 'asc' | 'desc';
  includeAnalytics: boolean;
  selectedProfileId: number | null;
}

interface BatchCacheContextValue {
  getCache: () => BatchCache | null;
  setCache: (cache: BatchCache) => void;
  clearCache: () => void;
}

const BatchCacheContext = createContext<BatchCacheContextValue | undefined>(undefined);

export function BatchCacheProvider({ children }: { children: ReactNode }) {
  const cacheRef = useRef<BatchCache | null>(null);

  const getCache = useCallback(() => cacheRef.current, []);
  const setCache = useCallback((cache: BatchCache) => {
    cacheRef.current = cache;
  }, []);
  const clearCache = useCallback(() => {
    cacheRef.current = null;
  }, []);

  return (
    <BatchCacheContext.Provider value={{ getCache, setCache, clearCache }}>
      {children}
    </BatchCacheContext.Provider>
  );
}

export function useBatchCache() {
  const ctx = useContext(BatchCacheContext);
  if (!ctx) throw new Error('useBatchCache must be used within BatchCacheProvider');
  return ctx;
}
