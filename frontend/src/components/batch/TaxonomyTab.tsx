/**
 * TaxonomyTab
 *
 * Full taxonomy classification tab with:
 * - Classify trigger button
 * - Horizontal bar chart of top 20 chemotype categories (click-to-filter)
 * - Searchable/paginated category table
 * - Empty, computing, and error states
 */

import { useState, useCallback, useMemo, useEffect, useRef } from 'react';
import { Loader2, AlertTriangle, RotateCcw, Info } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { TaxonomyBarChart } from './TaxonomyBarChart';
import { TaxonomyTable } from './TaxonomyTable';
import type { BatchAnalyticsResponse, TaxonomyResult } from '../../types/analytics';
import type { BatchResult } from '../../types/batch';

interface TaxonomyTabProps {
  /** Full analytics response (includes status and taxonomy data) */
  analyticsData: BatchAnalyticsResponse | null;
  /** Batch results (for context) */
  results: BatchResult[];
  /** Trigger a re-computation of a given analytics type */
  onRetrigger: (type: string, params?: Record<string, string>) => void;
  /** Propagate category filter to batch table */
  onCategoryFilter?: (category: string | null) => void;
  /** Navigate to a molecule in the Detailed Results table */
  onNavigateToMolecule?: (moleculeIndex: number) => void;
  /** Propagate selected molecule indices to parent */
  onSelectionChange?: (indices: Set<number>) => void;
}

export function TaxonomyTab({
  analyticsData,
  results,
  onRetrigger,
  onCategoryFilter,
  onNavigateToMolecule,
  onSelectionChange,
}: TaxonomyTabProps) {
  // Build index→name lookup from batch results for molecule name tooltips
  const nameMap = useMemo(() => {
    const map = new Map<number, string>();
    for (const r of results) {
      if (r.name) map.set(r.index, r.name);
    }
    return map;
  }, [results]);
  const [searchQuery, setSearchQuery] = useState('');
  const [debouncedQuery, setDebouncedQuery] = useState('');
  const [activeCategory, setActiveCategory] = useState<string | null>(null);
  const [, setHasTriggered] = useState(false);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Compute derived state
  const isComputing = analyticsData?.status?.taxonomy?.status === 'computing';
  const taxonomyResult: TaxonomyResult | undefined = analyticsData?.taxonomy;
  const errorMessage = analyticsData?.status?.taxonomy?.error;

  // Debounce search query (300ms)
  useEffect(() => {
    if (debounceRef.current) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => {
      setDebouncedQuery(searchQuery);
    }, 300);
    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, [searchQuery]);

  const handleClassify = useCallback(() => {
    onRetrigger('taxonomy');
    setHasTriggered(true);
  }, [onRetrigger]);

  const handleCategoryClick = useCallback(
    (category: string) => {
      const newCategory = activeCategory === category ? null : category;
      setActiveCategory(newCategory);
      onCategoryFilter?.(newCategory);
      if (onSelectionChange && taxonomyResult) {
        if (newCategory === null) {
          onSelectionChange(new Set());
        } else {
          const indices = new Set<number>();
          for (const mol of taxonomyResult.per_molecule) {
            if (mol.categories.some((c) => c.category === newCategory)) {
              indices.add(mol.index);
            }
          }
          onSelectionChange(indices);
        }
      }
    },
    [activeCategory, onCategoryFilter, onSelectionChange, taxonomyResult],
  );

  // Summary stats
  const summaryText = useMemo(() => {
    if (!taxonomyResult) return null;
    const total = taxonomyResult.total_molecules;
    const classified = taxonomyResult.classified_molecules;
    const catCount = Object.keys(taxonomyResult.category_counts).length;
    return `${classified}/${total} molecules classified into ${catCount} categories`;
  }, [taxonomyResult]);

  return (
    <div className="space-y-4">
      {/* Top row: Classify button */}
      <div className="flex items-center gap-3">
        <ClayButton
          variant="primary"
          onClick={handleClassify}
          disabled={isComputing}
          aria-label="Classify molecules into chemotype categories"
        >
          {isComputing ? (
            <>
              <Loader2 className="w-4 h-4 animate-spin mr-2" />
              Classifying...
            </>
          ) : (
            'Classify'
          )}
        </ClayButton>
        {summaryText && (
          <span className="text-sm text-[var(--color-text-muted)]">{summaryText}</span>
        )}
      </div>

      {/* Method info box */}
      <div className="rounded-2xl border border-[var(--color-primary)]/15 bg-gradient-to-br from-[var(--color-primary)]/[0.03] to-transparent px-5 py-4">
        <div className="flex gap-3">
          <div className="flex-shrink-0 mt-0.5">
            <div className="w-7 h-7 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center">
              <Info className="w-4 h-4 text-[var(--color-primary)]" />
            </div>
          </div>
          <div className="space-y-1.5 text-xs text-[var(--color-text-muted)] leading-relaxed">
            <p className="text-sm font-semibold text-[var(--color-text-primary)]">
              How classification works
            </p>
            <p>
              Each molecule is matched against <strong className="text-[var(--color-text-secondary)]">~50 curated SMARTS rules</strong> using
              RDKit substructure matching. Rules are organized into three groups:{' '}
              <strong className="text-[var(--color-text-secondary)]">Ring Systems</strong> (e.g. indoles, pyridines, benzodiazepines),{' '}
              <strong className="text-[var(--color-text-secondary)]">Functional Groups</strong> (e.g. sulfonamides, amides, phenols), and{' '}
              <strong className="text-[var(--color-text-secondary)]">Pharmacophoric/Drug-class</strong> (e.g. steroids, barbiturates, coumarins).
            </p>
            <p>
              A molecule can match <strong className="text-[var(--color-text-secondary)]">multiple categories</strong> simultaneously, so category
              counts may sum to more than the number of classified molecules. Molecules that match
              no rules are reported as unclassified.
            </p>
          </div>
        </div>
      </div>

      {/* Error state */}
      {errorMessage && (
        <div className="rounded-2xl p-5 bg-amber-500/5 border border-amber-500/20">
          <div className="flex items-center gap-3">
            <AlertTriangle className="w-5 h-5 text-amber-500 flex-shrink-0" />
            <div className="flex-1">
              <p className="text-sm text-amber-700 dark:text-amber-300">{errorMessage}</p>
              <p className="text-xs text-[var(--color-text-muted)] mt-1">
                Try clicking Classify again to retry.
              </p>
            </div>
            <button
              onClick={handleClassify}
              className="p-1.5 rounded-lg bg-amber-500/10 hover:bg-amber-500/20 transition-colors text-amber-600 dark:text-amber-400"
              title="Retry classification"
            >
              <RotateCcw className="w-4 h-4" />
            </button>
          </div>
        </div>
      )}

      {/* Computing state */}
      {isComputing && !taxonomyResult && (
        <div className="flex items-center gap-3 p-5 rounded-2xl bg-[var(--color-surface-sunken)]/50 border border-[var(--color-border)]">
          <Loader2 className="w-5 h-5 text-[var(--color-primary)] animate-spin" />
          <span className="text-sm text-[var(--color-text-primary)]">
            Classifying molecules...
          </span>
        </div>
      )}

      {/* Results state */}
      {taxonomyResult && (
        <div className="space-y-4">
          {/* ChartCard wrapper */}
          <figure
            aria-label="Chemotype distribution bar chart"
            className="rounded-2xl p-5 bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)] border border-[var(--color-border)]"
          >
            <h4 className="text-base font-semibold text-[var(--color-text-primary)] font-display mb-4">
              Chemotype Distribution
            </h4>
            <TaxonomyBarChart
              categoryCounts={taxonomyResult.category_counts}
              activeCategory={activeCategory}
              onCategoryClick={handleCategoryClick}
            />
          </figure>

          {/* Search input */}
          <div className="relative">
            <label htmlFor="taxonomy-search" className="sr-only">
              Search chemotype categories
            </label>
            <input
              id="taxonomy-search"
              type="text"
              placeholder="Search categories..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              className="w-full px-4 py-2.5 rounded-xl border border-[var(--color-border)] bg-[var(--color-surface-elevated)] text-sm text-[var(--color-text-primary)] placeholder-[var(--color-text-muted)] focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)] transition-shadow"
            />
          </div>

          {/* Category table */}
          <TaxonomyTable taxonomyResult={taxonomyResult} searchQuery={debouncedQuery} nameMap={nameMap} onNavigateToMolecule={onNavigateToMolecule} />
        </div>
      )}

      {/* Empty state (no result, not computing, no error) */}
      {!taxonomyResult && !isComputing && !errorMessage && (
        <div className="text-center py-12 space-y-2">
          <p className="text-sm font-medium text-[var(--color-text-primary)]">
            No taxonomy results
          </p>
          <p className="text-xs text-[var(--color-text-muted)]">
            Click &apos;Classify&apos; to categorize molecules into drug-relevant chemotype classes.
          </p>
        </div>
      )}
    </div>
  );
}
