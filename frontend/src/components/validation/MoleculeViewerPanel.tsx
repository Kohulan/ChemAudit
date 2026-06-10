import { type ReactNode, type Ref } from 'react';
import { motion } from 'framer-motion';
import { Atom, Download, Lock } from 'lucide-react';

import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { cn } from '../../lib/utils';

function formatMolecularFormula(formula: string): ReactNode {
  return formula
    .split(/(\d+)/)
    .map((part, i) => (/^\d+$/.test(part) ? <sub key={i}>{part}</sub> : <span key={i}>{part}</span>));
}

interface ViewerSummary {
  molecular_formula?: string | null;
  molecular_weight?: number | null;
}

interface ViewerStereo {
  hasStereochemistry: boolean;
  numStereocenters: number;
  hasEZStereo: boolean;
}

interface MoleculeViewerPanelProps {
  /** Ref attached to the preview container (parent uses it to export the SVG). */
  previewRef: Ref<HTMLDivElement>;
  molecule: string;
  canonicalSmiles: string | null | undefined;
  /** Server-side molecule info for the summary line (formula/MW), or null. */
  summary: ViewerSummary | null;
  /** Client-side stereochemistry info, or null. */
  stereo: ViewerStereo | null;
  highlightedAtoms: number[];
  highlightLocked: boolean;
  showCIP: boolean;
  onToggleCIP: () => void;
  onDownload: () => void;
}

/** "Structure Preview" card: RDKit.js viewer + CIP toggle + SVG download + stereo info. */
export function MoleculeViewerPanel({
  previewRef,
  molecule,
  canonicalSmiles,
  summary,
  stereo,
  highlightedAtoms,
  highlightLocked,
  showCIP,
  onToggleCIP,
  onDownload,
}: MoleculeViewerPanelProps) {
  return (
    <div className="card-glow p-4 sm:p-5">
      <div className="flex items-center gap-3 mb-3">
        <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
          <Atom className="w-5 h-5" />
        </div>
        <div className="flex-1">
          <h4 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
            Structure Preview
          </h4>
          <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
            {summary ? (
              <>
                {summary.molecular_formula && (
                  <span className="font-mono text-[var(--color-text-secondary)]">
                    {formatMolecularFormula(summary.molecular_formula)}
                  </span>
                )}
                {summary.molecular_formula && summary.molecular_weight && <span> · </span>}
                {summary.molecular_weight && <span>MW {summary.molecular_weight.toFixed(1)}</span>}
                <span className="text-[var(--color-text-muted)]"> · rendered with RDKit.js</span>
              </>
            ) : molecule ? (
              'Rendered with RDKit.js'
            ) : (
              'Enter a SMILES to preview'
            )}
          </p>
        </div>
        {/* CIP Labels Toggle - Show when molecule has stereochemistry */}
        {stereo?.hasStereochemistry && (
          <button
            onClick={onToggleCIP}
            className={cn(
              'flex items-center gap-1.5 text-xs font-medium px-3 py-1.5 rounded-lg transition-all',
              'border shadow-sm',
              showCIP
                ? 'bg-[var(--color-primary)] text-white border-[var(--color-primary)] shadow-[var(--color-primary)]/20'
                : 'bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] border-[var(--color-border)] hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]'
            )}
            title="Show R/S and E/Z stereochemistry labels"
          >
            <svg className="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="10" />
              <path d="M12 6v6l4 2" />
            </svg>
            {showCIP ? 'Hide CIP' : 'Show CIP'}
          </button>
        )}
        {molecule && (
          <button
            onClick={onDownload}
            className={cn(
              'flex items-center gap-1.5 text-xs font-medium px-3 py-1.5 rounded-lg transition-all',
              'border shadow-sm',
              highlightLocked
                ? 'bg-orange-500/10 text-orange-600 dark:text-orange-400 border-orange-500/30 hover:bg-orange-500/20'
                : 'bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] border-[var(--color-border)] hover:border-[var(--color-primary)] hover:text-[var(--color-primary)]'
            )}
            title={highlightLocked ? 'Download SVG with highlighted atoms' : 'Download structure as SVG'}
          >
            <Download className="w-3.5 h-3.5" />
            {highlightLocked ? 'SVG + Highlights' : 'SVG'}
          </button>
        )}
      </div>
      <div ref={previewRef} className="molecule-preview rounded-xl">
        <MoleculeViewer
          smiles={canonicalSmiles || molecule}
          highlightAtoms={highlightedAtoms}
          width={700}
          height={500}
          showCIP={showCIP}
        />
      </div>
      {highlightedAtoms.length > 0 && (
        <motion.p
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          className={cn(
            'mt-3 text-xs text-center font-medium inline-flex items-center justify-center gap-1.5 w-full',
            highlightLocked ? 'text-orange-500' : 'text-amber-500'
          )}
        >
          {highlightLocked && <Lock className="w-3 h-3" />}
          <span>
            Highlighting atoms: {highlightedAtoms.join(', ')}
            {highlightLocked && ' (locked for download)'}
          </span>
        </motion.p>
      )}
      {/* Stereochemistry info indicator */}
      {stereo?.hasStereochemistry && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          className="mt-3 flex items-center justify-center gap-2"
        >
          <div className="flex items-center gap-1.5 text-xs bg-chem-secondary-500/10 text-chem-secondary-600 dark:text-chem-secondary-400 px-2 py-1 rounded-lg">
            <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5" />
            </svg>
            <span>
              {stereo.numStereocenters > 0 &&
                `${stereo.numStereocenters} stereocenter${stereo.numStereocenters > 1 ? 's' : ''}`}
              {stereo.numStereocenters > 0 && stereo.hasEZStereo && ' + '}
              {stereo.hasEZStereo && 'E/Z bonds'}
            </span>
          </div>
          {showCIP && (
            <span className="text-xs text-[var(--color-text-muted)]">(CIP labels shown)</span>
          )}
        </motion.div>
      )}
    </div>
  );
}
