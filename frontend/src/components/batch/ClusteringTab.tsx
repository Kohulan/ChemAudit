/**
 * ClusteringTab
 *
 * Main clustering tab within BatchAnalyticsPanel. Provides Butina clustering
 * with configurable Tanimoto cutoff, summary badges, sortable cluster table,
 * and expandable member grids. Handles empty, computing, and error states.
 */

import { useState, useCallback } from 'react';
import { Loader2, AlertTriangle } from 'lucide-react';
import { ClusterCutoffSelector } from './ClusterCutoffSelector';
import { ClusterSummaryBar } from './ClusterSummaryBar';
import { ClusterTable } from './ClusterTable';
import type { BatchAnalyticsResponse } from '../../types/analytics';
import type { BatchResult } from '../../types/batch';

interface ClusteringTabProps {
  analyticsData: BatchAnalyticsResponse | null;
  results: BatchResult[];
  onRetrigger: (type: string, params?: Record<string, string>) => void;
}

export function ClusteringTab({
  analyticsData,
  results,
  onRetrigger,
}: ClusteringTabProps) {
  const [cutoff, setCutoff] = useState(0.35);
  const [expandedClusterId, setExpandedClusterId] = useState<number | null>(null);
  const [hasTriggered, setHasTriggered] = useState(false);

  const isComputing = analyticsData?.status?.clustering?.status === 'computing';
  const clusteringResult = analyticsData?.clustering ?? null;
  const errorMessage = analyticsData?.status?.clustering?.error ?? null;

  const handleCluster = useCallback(() => {
    onRetrigger('clustering', { cutoff: String(cutoff) });
    setHasTriggered(true);
  }, [onRetrigger, cutoff]);

  const handleToggleExpand = useCallback((id: number) => {
    setExpandedClusterId((prev) => (prev === id ? null : id));
  }, []);

  return (
    <div className="space-y-4">
      {/* Top row: cutoff selector + cluster button */}
      <ClusterCutoffSelector
        cutoff={cutoff}
        onCutoffChange={setCutoff}
        onRecluster={handleCluster}
        isComputing={!!isComputing}
        hasTriggered={hasTriggered}
      />

      {/* Error state */}
      {errorMessage && (
        <div className="rounded-xl p-4 bg-red-500/5 border border-red-500/20">
          <div className="flex items-start gap-3">
            <AlertTriangle className="w-5 h-5 text-red-500 flex-shrink-0 mt-0.5" />
            <div className="text-sm text-red-600 dark:text-red-400">
              {errorMessage.includes('1,000') || errorMessage.includes('1000') ? (
                <p>
                  Clustering is limited to 1,000 molecules. This batch has{' '}
                  {results.length.toLocaleString()} molecules. Filter or subsample before
                  clustering.
                </p>
              ) : (
                <p>{errorMessage}</p>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Computing state: overlay previous results with spinner */}
      {isComputing && (
        <div className="relative">
          {clusteringResult && (
            <div className="opacity-50 pointer-events-none">
              <ClusterSummaryBar
                clusterCount={clusteringResult.cluster_count}
                singletonCount={clusteringResult.singleton_count}
                largestClusterSize={clusteringResult.largest_cluster_size}
              />
              <ClusterTable
                clusters={clusteringResult.clusters}
                results={results}
                smilesMap={clusteringResult.smiles_map}
                expandedClusterId={expandedClusterId}
                onToggleExpand={handleToggleExpand}
              />
            </div>
          )}
          <div className="flex flex-col items-center justify-center py-8 gap-3">
            <Loader2 className="w-6 h-6 text-[var(--color-primary)] animate-spin" />
            <p className="text-sm text-[var(--color-text-muted)]">Clustering molecules...</p>
          </div>
        </div>
      )}

      {/* Results: summary + table */}
      {!isComputing && clusteringResult && (
        <>
          <ClusterSummaryBar
            clusterCount={clusteringResult.cluster_count}
            singletonCount={clusteringResult.singleton_count}
            largestClusterSize={clusteringResult.largest_cluster_size}
          />
          <ClusterTable
            clusters={clusteringResult.clusters}
            results={results}
            smilesMap={clusteringResult.smiles_map}
            expandedClusterId={expandedClusterId}
            onToggleExpand={handleToggleExpand}
          />
        </>
      )}

      {/* Empty state: no clustering, not computing */}
      {!isComputing && !clusteringResult && !errorMessage && (
        <div className="flex flex-col items-center justify-center py-12 text-center">
          <p className="text-sm font-medium text-[var(--color-text-secondary)] mb-1">
            No clustering results
          </p>
          <p className="text-xs text-[var(--color-text-muted)] max-w-xs">
            Click &apos;Cluster&apos; to group molecules by structural similarity using the Butina
            algorithm.
          </p>
        </div>
      )}
    </div>
  );
}
