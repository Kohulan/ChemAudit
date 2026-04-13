import { motion } from 'framer-motion';
import { ClayButton } from '../ui/ClayButton';
import type { FixSuggestion } from '../../types/diagnostics';

interface FixChipProps {
  /** The fix suggestion from the API. */
  suggestion: FixSuggestion;
  /** Called with the corrected SMILES when the user applies the fix. */
  onApply: (correctedSmiles: string) => void;
  /** When true, shows a spinner and disables the button. */
  isApplying?: boolean;
}

/**
 * Clickable chip that auto-applies a SMILES fix suggestion.
 *
 * Rendered as a ClayButton with variant="outline". If the suggestion has no
 * corrected SMILES it is shown as a disabled informational chip.
 *
 * Keyboard: Enter and Space activate (native button behaviour).
 * Accessibility: aria-label includes the full description.
 * Animation: fade-in-up on mount (0.3s ease-out), staggered by Framer Motion.
 */
export function FixChip({ suggestion, onApply, isApplying = false }: FixChipProps) {
  const hasCorrection = suggestion.corrected_smiles !== null;
  const isDisabled = !hasCorrection || isApplying;

  const handleClick = () => {
    if (suggestion.corrected_smiles) {
      onApply(suggestion.corrected_smiles);
    }
  };

  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3, ease: 'easeOut' }}
    >
      <ClayButton
        variant="outline"
        size="sm"
        aria-label={`Apply fix: ${suggestion.description}`}
        onClick={handleClick}
        disabled={isDisabled}
        loading={isApplying}
        className={[
          'text-xs font-semibold min-h-[32px] min-w-[64px]',
          !hasCorrection ? 'opacity-60 cursor-not-allowed' : '',
        ]
          .filter(Boolean)
          .join(' ')}
      >
        {suggestion.description}
      </ClayButton>
    </motion.div>
  );
}
