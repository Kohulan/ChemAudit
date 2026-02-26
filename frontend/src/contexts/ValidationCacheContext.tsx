import { createContext, useContext, useRef, useCallback, type ReactNode } from 'react';
import type { AlertScreenResponse } from '../types/alerts';
import type { ScoringResponse } from '../types/scoring';
import type { StandardizeResponse } from '../types/standardization';
import type { ValidationResponse } from '../types/validation';
import type { PubChemResult, ChEMBLResult, COCONUTResult } from '../types/integrations';

export type TabType = 'validate' | 'deep-validation' | 'scoring-profiles' | 'database' | 'alerts' | 'standardize';

export interface ValidationCache {
  molecule: string;
  activeTab: TabType;
  result: ValidationResponse | null;
  alertResult: AlertScreenResponse | null;
  scoringResult: ScoringResponse | null;
  standardizationResult: StandardizeResponse | null;
  databaseResults: {
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
  } | null;
}

interface ValidationCacheContextValue {
  getCache: () => ValidationCache | null;
  setCache: (cache: ValidationCache) => void;
  clearCache: () => void;
}

const ValidationCacheContext = createContext<ValidationCacheContextValue | undefined>(undefined);

export function ValidationCacheProvider({ children }: { children: ReactNode }) {
  const cacheRef = useRef<ValidationCache | null>(null);

  const getCache = useCallback(() => cacheRef.current, []);
  const setCache = useCallback((cache: ValidationCache) => {
    cacheRef.current = cache;
  }, []);
  const clearCache = useCallback(() => {
    cacheRef.current = null;
  }, []);

  return (
    <ValidationCacheContext.Provider value={{ getCache, setCache, clearCache }}>
      {children}
    </ValidationCacheContext.Provider>
  );
}

export function useValidationCache() {
  const ctx = useContext(ValidationCacheContext);
  if (!ctx) throw new Error('useValidationCache must be used within ValidationCacheProvider');
  return ctx;
}
