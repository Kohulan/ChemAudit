import { useState, useCallback } from 'react';
import { SafetySummaryStrip } from './SafetySummaryStrip';
import { AlertDashboard } from './AlertDashboard';
import { SafetyFlagsGrid } from './SafetyFlagsGrid';
import { ComplexityRadar } from './ComplexityRadar';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import type { AlertScreenResponse, SafetyAssessResponse } from '../../types/safety';

interface SafetyAccordionProps {
  /** Current molecule SMILES */
  smiles: string;
  /** Result from /api/v1/alerts/screen (already fetched by parent) */
  alertResult: AlertScreenResponse | null;
  /** Result from /api/v1/safety/assess (CYP/hERG/bRo5/REOS/complexity) */
  safetyResult: SafetyAssessResponse | null;
  /** Whether the safety data is currently loading */
  isLoading: boolean;
  /** Error message if safety screening failed */
  error: string | null;
}

/**
 * SafetyAccordion -- wraps all safety standalone components into a single
 * panel suitable for rendering inside a DrillDownSection accordion on
 * SingleValidation.
 *
 * Layout mirrors Safety page:
 * 1. SafetySummaryStrip -- sticky header with alert count + traffic-light dots
 * 2. MoleculeViewer -- structure display with atom highlighting from alerts
 * 3. AlertDashboard -- dual-view with concern group and catalog source views
 * 4. SafetyFlagsGrid -- 2x2 MetricCard grid (CYP/hERG/bRo5/REOS)
 * 5. ComplexityRadar -- dual-polygon radar chart (molecule vs reference)
 *
 * States:
 * - Loading: skeleton placeholders
 * - Error: error card with message
 * - Empty (no smiles): prompt to enter a molecule
 * - Partial: alertResult present but safetyResult still loading
 */
export function SafetyAccordion({
  smiles,
  alertResult,
  safetyResult,
  isLoading,
  error,
}: SafetyAccordionProps) {
  // Atom highlight state -- shared between AlertDashboard and MoleculeViewer
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([]);
  const [highlightedAlert, setHighlightedAlert] = useState<string | null>(null);

  const handleHighlightAtoms = useCallback(
    (atoms: number[], alertName?: string) => {
      if (alertName !== undefined && alertName === highlightedAlert) {
        // Re-click same alert -> clear
        setHighlightedAtoms([]);
        setHighlightedAlert(null);
      } else {
        setHighlightedAtoms(atoms);
        setHighlightedAlert(alertName ?? null);
      }
    },
    [highlightedAlert]
  );

  // Empty state: no SMILES entered
  if (!smiles && !isLoading && !error) {
    return (
      <div className="text-center py-8 text-text-muted">
        <p className="text-sm">Enter a molecule to assess safety</p>
      </div>
    );
  }

  // Loading state: full skeleton placeholders
  if (isLoading) {
    return (
      <div className="space-y-6">
        {/* Summary strip skeleton */}
        <div className="flex items-center justify-between gap-4 px-5 py-3 rounded-2xl border border-[var(--color-border)]">
          <div className="animate-pulse rounded-xl bg-surface-sunken h-6 w-60" />
          <div className="flex gap-4">
            {[0, 1, 2, 3].map((i) => (
              <div key={i} className="animate-pulse rounded-xl bg-surface-sunken h-4 w-12" />
            ))}
          </div>
          <div className="animate-pulse rounded-xl bg-surface-sunken h-4 w-44" />
        </div>
        {/* MoleculeViewer skeleton */}
        <div className="flex justify-center">
          <div className="animate-pulse rounded-lg bg-surface-sunken h-[300px] w-[400px]" />
        </div>
        {/* Alert cards skeleton */}
        <div className="space-y-3">
          {[0, 1, 2].map((i) => (
            <div key={i} className="animate-pulse rounded-xl bg-surface-sunken h-16 w-full" />
          ))}
        </div>
        {/* Flags grid skeleton */}
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-6">
          {[0, 1, 2, 3].map((i) => (
            <div key={i} className="animate-pulse rounded-xl bg-surface-sunken h-40" />
          ))}
        </div>
        {/* Radar skeleton */}
        <div className="animate-pulse rounded-xl bg-surface-sunken h-[300px] w-full" />
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="p-4 rounded-xl bg-status-error/10 border border-status-error/20 text-status-error text-sm">
        {error}
      </div>
    );
  }

  const hasAnyResult = !!(alertResult || safetyResult);

  // No data yet (smiles set but no results arrived)
  if (!hasAnyResult) {
    return null;
  }

  // Derive SafetySummaryStrip props
  const totalAlerts = alertResult?.total_raw ?? 0;
  const concernGroupCount = Object.keys(alertResult?.concern_groups ?? {}).length;
  const hasCritical = alertResult?.has_critical ?? false;

  // Traffic-light statuses derived from safetyResult
  const cypStatus = 'default' as const;
  const hergStatus: 'success' | 'warning' | 'error' =
    safetyResult?.herg.herg_risk === 'high'
      ? 'error'
      : safetyResult?.herg.herg_risk === 'moderate'
        ? 'warning'
        : 'success';
  const bro5Status: 'success' | 'error' | 'default' =
    !safetyResult?.bro5.applicable
      ? 'default'
      : safetyResult.bro5.passed
        ? 'success'
        : 'error';
  const reosStatus: 'success' | 'warning' | 'error' =
    (safetyResult?.reos.n_violations ?? 0) === 0
      ? 'success'
      : (safetyResult?.reos.n_violations ?? 0) === 1
        ? 'warning'
        : 'error';
  const complexityOutliers = safetyResult?.complexity.n_outliers ?? 0;

  const displaySmiles =
    alertResult?.molecule_info?.smiles ?? safetyResult?.molecule_info?.smiles ?? smiles;

  return (
    <div className="space-y-6">
      {/* 1. Summary strip */}
      <SafetySummaryStrip
        totalAlerts={totalAlerts}
        concernGroupCount={concernGroupCount}
        hasCritical={hasCritical}
        cypStatus={cypStatus}
        hergStatus={hergStatus}
        bro5Status={bro5Status}
        reosStatus={reosStatus}
        complexityOutliers={complexityOutliers}
      />

      {/* 2. MoleculeViewer with atom highlighting */}
      <div className="flex justify-center">
        <MoleculeViewer
          smiles={displaySmiles}
          highlightAtoms={highlightedAtoms}
          width={400}
          height={300}
          className="rounded-lg border border-[var(--color-border)]"
        />
      </div>

      {/* 3. Alert Dashboard */}
      {alertResult && (
        <AlertDashboard
          alertResult={alertResult}
          onHighlightAtoms={(atoms) => handleHighlightAtoms(atoms)}
          highlightedAlert={highlightedAlert}
        />
      )}

      {/* 4. Safety Flags Grid */}
      {safetyResult ? (
        <SafetyFlagsGrid
          cyp={safetyResult.cyp_softspots}
          herg={safetyResult.herg}
          bro5={safetyResult.bro5}
          reos={safetyResult.reos}
          onHighlightAtoms={handleHighlightAtoms}
          highlightedAlert={highlightedAlert}
        />
      ) : (
        // Partial data: safetyResult still loading, show skeleton for flags
        <div className="space-y-4">
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-6">
            {[0, 1, 2, 3].map((i) => (
              <div key={i} className="animate-pulse rounded-xl bg-surface-sunken h-40" />
            ))}
          </div>
        </div>
      )}

      {/* 5. Complexity Radar */}
      {safetyResult ? (
        <ComplexityRadar complexity={safetyResult.complexity} />
      ) : (
        // Partial data: safetyResult still loading, show skeleton for radar
        <div className="animate-pulse rounded-xl bg-surface-sunken h-[300px] w-full" />
      )}
    </div>
  );
}
