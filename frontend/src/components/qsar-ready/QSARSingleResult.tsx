import { motion } from 'framer-motion';
import { AlertTriangle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { InChIKeyChangeChip } from './InChIKeyChangeChip';
import { PipelineStepTimeline } from './PipelineStepTimeline';
import type { QSARReadyResult } from '../../types/qsar_ready';

// =============================================================================
// Props
// =============================================================================

interface QSARSingleResultProps {
  /** Full QSAR pipeline result for a single molecule. */
  result: QSARReadyResult;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Before/after structure comparison with InChIKey change indicator and step timeline.
 *
 * Layout per UI-SPEC Single Molecule Result Layout:
 * - grid grid-cols-2 gap-6 for before/after structure panels
 * - Centered InChIKeyChangeChip between panels
 * - PipelineStepTimeline below structures
 *
 * Rejected/error status shows error state with rejection_reason.
 * Fade-in-up animation on mount per UI-SPEC animation tokens.
 */
export function QSARSingleResult({ result }: QSARSingleResultProps) {
  const isRejected = result.status === 'rejected' || result.status === 'error';

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4, ease: 'easeOut' }}
      className="space-y-6"
    >
      {/* ── Error / rejected state ── */}
      {isRejected && result.rejection_reason && (
        <div className="flex items-start gap-3 p-4 rounded-xl border border-red-200 bg-red-50 dark:bg-red-900/10 dark:border-red-900/30">
          <AlertTriangle className="w-4 h-4 text-red-500 mt-0.5 shrink-0" />
          <div>
            <p className="text-sm font-semibold text-red-700 dark:text-red-400 font-display">
              {result.status === 'rejected' ? 'Molecule Rejected' : 'Pipeline Error'}
            </p>
            <p className="text-sm text-red-600 dark:text-red-500 mt-0.5">
              {result.rejection_reason}
            </p>
          </div>
        </div>
      )}

      {/* ── Before / After structure grid ── */}
      <div className="grid grid-cols-2 gap-6">
        {/* Before (original) */}
        <ClayCard variant="flat" size="sm" className="p-4 space-y-3">
          <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
            Before (original)
          </p>
          <MoleculeViewer
            smiles={result.original_smiles}
            width={200}
            height={200}
            className="rounded-lg"
          />
          <div className="space-y-1">
            <p
              className="font-mono text-xs text-[var(--color-text-primary)] break-all"
              title={result.original_smiles}
            >
              {result.original_smiles}
            </p>
            {result.original_inchikey && (
              <p
                className="font-mono text-xs text-[var(--color-text-muted)] truncate"
                title={result.original_inchikey}
                aria-label={`InChIKey: ${result.original_inchikey}`}
              >
                {result.original_inchikey}
              </p>
            )}
          </div>
        </ClayCard>

        {/* After (curated) */}
        <ClayCard variant="flat" size="sm" className="p-4 space-y-3">
          <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
            After (curated)
          </p>
          <MoleculeViewer
            smiles={result.curated_smiles}
            width={200}
            height={200}
            className="rounded-lg"
          />
          <div className="space-y-1">
            {result.curated_smiles ? (
              <p
                className="font-mono text-xs text-[var(--color-text-primary)] break-all"
                title={result.curated_smiles}
              >
                {result.curated_smiles}
              </p>
            ) : (
              <p className="font-mono text-xs text-[var(--color-text-muted)]">—</p>
            )}
            {result.standardized_inchikey && (
              <p
                className="font-mono text-xs text-[var(--color-text-muted)] truncate"
                title={result.standardized_inchikey}
                aria-label={`InChIKey: ${result.standardized_inchikey}`}
              >
                {result.standardized_inchikey}
              </p>
            )}
          </div>
        </ClayCard>
      </div>

      {/* ── InChIKey change chip (centered between panels) ── */}
      <div className="flex justify-center">
        <InChIKeyChangeChip
          originalInchikey={result.original_inchikey}
          standardizedInchikey={result.standardized_inchikey}
          inchikeyChanged={result.inchikey_changed}
        />
      </div>

      {/* ── Pipeline step timeline ── */}
      {result.steps.length > 0 && (
        <div className="space-y-3">
          <h3 className="text-sm font-semibold font-display text-[var(--color-text-primary)]">
            Pipeline Steps
          </h3>
          <PipelineStepTimeline steps={result.steps} />
        </div>
      )}
    </motion.div>
  );
}
