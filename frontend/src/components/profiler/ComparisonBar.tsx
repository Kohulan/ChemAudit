import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X, Plus, GitCompare } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import type { ProfileResponse } from '../../types/profiler';

/**
 * A molecule that has been pinned for comparison.
 */
export interface PinnedMolecule {
  smiles: string;
  label: string;
  profile: ProfileResponse;
}

interface ComparisonBarProps {
  molecules: PinnedMolecule[];
  onRemove: (smiles: string) => void;
  onCompare: () => void;
  onAddMultiple: (smilesList: string[]) => void;
}

/**
 * ComparisonBar — fixed bottom strip showing pinned molecule chips.
 *
 * Per D-20: appears when >= 1 molecule is pinned. Max 5 molecules enforced.
 * Provides "Compare" CTA (disabled when < 2 pinned) and batch SMILES paste.
 */
export function ComparisonBar({
  molecules,
  onRemove,
  onCompare,
  onAddMultiple,
}: ComparisonBarProps) {
  const [showPasteOverlay, setShowPasteOverlay] = useState(false);
  const [pasteValue, setPasteValue] = useState('');
  const [pasteError, setPasteError] = useState<string | null>(null);

  if (molecules.length === 0) return null;

  const handlePasteSubmit = () => {
    const lines = pasteValue
      .split('\n')
      .map((l) => l.trim())
      .filter(Boolean);

    if (lines.length === 0) {
      setPasteError('Please enter at least one SMILES.');
      return;
    }

    const totalAfter = molecules.length + lines.length;
    if (totalAfter > 5) {
      const allowed = 5 - molecules.length;
      setPasteError(
        `You can add at most ${allowed} more molecule${allowed === 1 ? '' : 's'} (max 5 total).`
      );
      return;
    }

    setPasteError(null);
    onAddMultiple(lines);
    setPasteValue('');
    setShowPasteOverlay(false);
  };

  return (
    <>
      {/* Batch paste overlay */}
      <AnimatePresence>
        {showPasteOverlay && (
          <motion.div
            className="fixed inset-0 z-50 flex items-end justify-center"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
          >
            {/* Backdrop */}
            <div
              className="absolute inset-0 bg-black/40 backdrop-blur-sm"
              onClick={() => {
                setShowPasteOverlay(false);
                setPasteError(null);
              }}
            />
            {/* Paste panel */}
            <motion.div
              className="relative z-10 w-full max-w-lg mx-4 mb-20 clay-card p-5"
              initial={{ opacity: 0, y: 24 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: 24 }}
              transition={{ duration: 0.25, ease: 'easeOut' }}
            >
              <div className="flex items-center justify-between mb-3">
                <h3 className="text-base font-semibold font-display text-text-primary">
                  Add molecules for comparison
                </h3>
                <button
                  className="text-text-muted hover:text-text-primary transition-colors"
                  onClick={() => {
                    setShowPasteOverlay(false);
                    setPasteError(null);
                  }}
                  aria-label="Close add molecules overlay"
                >
                  <X className="w-4 h-4" />
                </button>
              </div>
              <p className="text-xs text-text-secondary mb-3">
                Paste 2–5 SMILES, one per line.
                {molecules.length > 0 && (
                  <span className="ml-1">
                    ({5 - molecules.length} slot{5 - molecules.length === 1 ? '' : 's'} remaining)
                  </span>
                )}
              </p>
              <textarea
                className="w-full h-28 text-xs font-mono bg-surface-sunken border border-border rounded-lg p-3
                           text-text-primary placeholder:text-text-muted resize-none
                           focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30"
                placeholder={`CC(=O)Oc1ccccc1C(=O)O\nCn1c(=O)c2c(ncn2C)n(C)c1=O`}
                value={pasteValue}
                onChange={(e) => {
                  setPasteValue(e.target.value);
                  setPasteError(null);
                }}
              />
              {pasteError && (
                <p className="mt-1 text-xs text-status-error">{pasteError}</p>
              )}
              <div className="flex gap-2 mt-3 justify-end">
                <ClayButton
                  variant="ghost"
                  size="sm"
                  onClick={() => {
                    setShowPasteOverlay(false);
                    setPasteError(null);
                  }}
                >
                  Cancel
                </ClayButton>
                <ClayButton
                  variant="primary"
                  size="sm"
                  onClick={handlePasteSubmit}
                  disabled={!pasteValue.trim()}
                >
                  Add
                </ClayButton>
              </div>
            </motion.div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Fixed bottom bar */}
      <div className="fixed bottom-0 left-0 right-0 z-40 bg-surface-elevated border-t border-border px-4 py-3">
        <div className="max-w-7xl mx-auto flex items-center gap-3 flex-wrap">
          {/* Molecule chips */}
          <div className="flex items-center gap-2 flex-wrap flex-1 min-w-0">
            <AnimatePresence>
              {molecules.map((mol) => (
                <motion.div
                  key={mol.smiles}
                  className="flex items-center gap-1 h-8 px-3 rounded-full
                             bg-surface-sunken border border-border
                             text-xs font-mono text-text-secondary"
                  initial={{ scale: 0, opacity: 0 }}
                  animate={{ scale: 1, opacity: 1 }}
                  exit={{ scale: 0, opacity: 0 }}
                  transition={{ duration: 0.4, ease: 'easeOut' }}
                >
                  <span className="max-w-[120px] truncate" title={mol.smiles}>
                    {mol.label}
                  </span>
                  <button
                    className="ml-1 text-text-muted hover:text-status-error transition-colors flex-shrink-0"
                    onClick={() => onRemove(mol.smiles)}
                    aria-label={`Remove ${mol.label} from comparison`}
                  >
                    <X className="w-3 h-3" />
                  </button>
                </motion.div>
              ))}
            </AnimatePresence>
          </div>

          {/* Action buttons */}
          <div className="flex items-center gap-2 flex-shrink-0">
            {/* Add molecules button — hidden when at max */}
            {molecules.length < 5 && (
              <ClayButton
                variant="ghost"
                size="sm"
                leftIcon={<Plus className="w-3.5 h-3.5" />}
                onClick={() => setShowPasteOverlay(true)}
              >
                Add molecules
              </ClayButton>
            )}
            {molecules.length >= 5 && (
              <span className="text-xs text-text-muted mr-2">Max 5 reached</span>
            )}
            {/* Compare button */}
            <ClayButton
              variant="primary"
              size="sm"
              leftIcon={<GitCompare className="w-3.5 h-3.5" />}
              onClick={onCompare}
              disabled={molecules.length < 2}
            >
              Compare
            </ClayButton>
          </div>
        </div>
      </div>
    </>
  );
}
