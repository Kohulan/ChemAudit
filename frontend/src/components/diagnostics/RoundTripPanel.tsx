import { useState, useEffect, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ArrowRight, AlertTriangle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { Skeleton } from '../ui/Skeleton';
import { Tooltip } from '../ui/Tooltip';
import { LostInfoCard } from './LostInfoCard';
import type { RoundTripResponse } from '../../types/diagnostics';

interface RoundTripPanelProps {
  result: RoundTripResponse | null;
  isLoading: boolean;
  error: string | null;
  currentSmiles: string | null;
  onCheckRoundtrip: (smiles: string, route: string) => void;
  onRetry: () => void;
}

const ROUTES = [
  { value: 'smiles_inchi_smiles', label: 'SMILES \u2192 InChI \u2192 SMILES' },
  { value: 'smiles_mol_smiles', label: 'SMILES \u2192 MOL \u2192 SMILES' },
] as const;

/** Truncate a string to maxLen chars with ellipsis. */
function truncate(str: string | null | undefined, maxLen = 40): string {
  if (!str) return '—';
  return str.length > maxLen ? str.slice(0, maxLen) + '...' : str;
}

/**
 * Round-trip lossiness panel (DIAG-03).
 *
 * Per UI-SPEC (D-08):
 * - Route selector dropdown (SMILES→InChI→SMILES | SMILES→MOL→SMILES)
 * - 3-stage flow visualization: Original → Convert → Roundtrip
 * - Verdict badge: Lossless (warm amber) or Lossy (red)
 * - LostInfoCard shown when lossy=true
 * - Staggered fade-in-up animation: 80ms between boxes
 */
export function RoundTripPanel({
  result,
  isLoading,
  error,
  currentSmiles,
  onCheckRoundtrip,
  onRetry,
}: RoundTripPanelProps) {
  const [selectedRoute, setSelectedRoute] = useState<string>('smiles_inchi_smiles');

  // Auto-trigger on mount if currentSmiles is available and no result yet
  const hasTriggered = useRef(false);
  useEffect(() => {
    if (currentSmiles && !result && !isLoading && !error && !hasTriggered.current) {
      hasTriggered.current = true;
      onCheckRoundtrip(currentSmiles, selectedRoute);
    }
  }, [currentSmiles, result, isLoading, error, onCheckRoundtrip, selectedRoute]);

  const handleRouteChange = (route: string) => {
    setSelectedRoute(route);
    if (currentSmiles) {
      onCheckRoundtrip(currentSmiles, route);
    }
  };

  // Network error state
  if (error) {
    return (
      <motion.div
        initial={{ opacity: 0, y: 8 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.25 }}
      >
        <ClayCard variant="flat" size="sm" className="border border-[var(--color-border)]">
          <div className="flex items-start gap-3">
            <AlertTriangle className="w-4 h-4 text-status-error mt-0.5 shrink-0" />
            <div className="flex-1 min-w-0">
              <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                Round-trip check failed
              </p>
              <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">
                Round-trip check failed. Verify the SMILES is valid and retry.
              </p>
            </div>
            <ClayButton variant="ghost" size="sm" onClick={onRetry} className="shrink-0">
              Retry
            </ClayButton>
          </div>
        </ClayCard>
      </motion.div>
    );
  }

  const stageAnimations = [
    { initial: { opacity: 0, y: 8 }, animate: { opacity: 1, y: 0 }, transition: { duration: 0.4, ease: 'easeOut', delay: 0 } },
    { initial: { opacity: 0, y: 8 }, animate: { opacity: 1, y: 0 }, transition: { duration: 0.4, ease: 'easeOut', delay: 0.08 } },
    { initial: { opacity: 0, y: 8 }, animate: { opacity: 1, y: 0 }, transition: { duration: 0.4, ease: 'easeOut', delay: 0.16 } },
  ];

  return (
    <div className="space-y-4">
      {/* Route selector */}
      <div className="flex items-center gap-3">
        <label
          htmlFor="roundtrip-route"
          className="text-xs font-semibold text-[var(--color-text-primary)] shrink-0"
        >
          Route
        </label>
        <select
          id="roundtrip-route"
          value={selectedRoute}
          onChange={(e) => handleRouteChange(e.target.value)}
          className="text-sm text-[var(--color-text-primary)] bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded-lg px-3 py-1.5 outline-none focus:ring-2 focus:ring-chem-primary-600/30"
        >
          {ROUTES.map((r) => (
            <option key={r.value} value={r.value}>
              {r.label}
            </option>
          ))}
        </select>
      </div>

      {/* Flow visualization */}
      <div
        className="flex items-center gap-2 flex-wrap sm:flex-nowrap"
        aria-label="Round-trip pipeline flow: Original → Convert → Roundtrip"
      >
        {/* Stage 1: Original */}
        <motion.div
          className="flex-1 min-w-0"
          {...stageAnimations[0]}
        >
          <ClayCard variant="default" size="sm" className="clay-card-sm p-3">
            <p className="text-xs font-semibold text-[var(--color-text-muted)] mb-1">Original</p>
            {isLoading ? (
              <Skeleton className="h-4 w-full" />
            ) : (
              <Tooltip content={result?.original_smiles ?? currentSmiles ?? ''}>
                <span className="font-mono text-sm text-[var(--color-text-primary)] break-all">
                  {truncate(result?.original_smiles ?? currentSmiles)}
                </span>
              </Tooltip>
            )}
          </ClayCard>
        </motion.div>

        {/* Arrow */}
        <motion.div {...stageAnimations[0]} className="shrink-0">
          <ArrowRight className="w-4 h-4 text-[var(--color-text-muted)]" />
        </motion.div>

        {/* Stage 2: Convert */}
        <motion.div
          className="flex-1 min-w-0"
          {...stageAnimations[1]}
        >
          <ClayCard variant="default" size="sm" className="clay-card-sm p-3">
            <p className="text-xs font-semibold text-[var(--color-text-muted)] mb-1">Convert</p>
            {isLoading ? (
              <Skeleton className="h-4 w-full" />
            ) : (
              <Tooltip content={result?.intermediate ?? ''}>
                <span className="font-mono text-sm text-[var(--color-text-primary)] break-all">
                  {truncate(result?.intermediate)}
                </span>
              </Tooltip>
            )}
          </ClayCard>
        </motion.div>

        {/* Arrow */}
        <motion.div {...stageAnimations[1]} className="shrink-0">
          <ArrowRight className="w-4 h-4 text-[var(--color-text-muted)]" />
        </motion.div>

        {/* Stage 3: Roundtrip */}
        <motion.div
          className="flex-1 min-w-0"
          {...stageAnimations[2]}
        >
          <ClayCard variant="default" size="sm" className="clay-card-sm p-3">
            <p className="text-xs font-semibold text-[var(--color-text-muted)] mb-1">Roundtrip</p>
            {isLoading ? (
              <Skeleton className="h-4 w-full" />
            ) : (
              <Tooltip content={result?.roundtrip_smiles ?? ''}>
                <span className="font-mono text-sm text-[var(--color-text-primary)] break-all">
                  {truncate(result?.roundtrip_smiles)}
                </span>
              </Tooltip>
            )}
          </ClayCard>
        </motion.div>
      </div>

      {/* Diagnostic error (conversion failure — not a network error) */}
      <AnimatePresence>
        {!isLoading && result?.error && (
          <motion.div
            key="diag-error"
            initial={{ opacity: 0, y: 6 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -6 }}
            transition={{ duration: 0.25 }}
          >
            <ClayCard variant="flat" size="sm" className="border border-[var(--color-border)]">
              <div className="flex items-start gap-2">
                <AlertTriangle className="w-4 h-4 text-status-warning mt-0.5 shrink-0" />
                <p className="text-sm text-[var(--color-text-secondary)]">{result.error}</p>
              </div>
            </ClayCard>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Verdict badge */}
      <AnimatePresence>
        {!isLoading && result && !result.error && (
          <motion.div
            key="verdict"
            initial={{ scale: 0.96, opacity: 0 }}
            animate={{ scale: 1, opacity: 1 }}
            exit={{ scale: 0.96, opacity: 0 }}
            transition={{ duration: 0.3, ease: 'easeOut' }}
            className="flex items-center gap-2"
          >
            {result.lossy === false ? (
              <Badge variant="success">Lossless — no information lost</Badge>
            ) : (
              <Badge variant="error">Lossy — information was lost in conversion</Badge>
            )}
          </motion.div>
        )}
      </AnimatePresence>

      {/* Lost info card (shown when lossy=true) */}
      <AnimatePresence>
        {!isLoading && result?.lossy && result.losses.length > 0 && (
          <LostInfoCard losses={result.losses} />
        )}
      </AnimatePresence>
    </div>
  );
}
