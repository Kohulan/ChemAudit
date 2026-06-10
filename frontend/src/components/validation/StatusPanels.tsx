import { type ReactNode } from 'react';
import { AnimatePresence, motion } from 'framer-motion';
import { AlertTriangle, RotateCcw } from 'lucide-react';

import { ClayButton } from '../ui/ClayButton';
import { MoleculeLoader } from '../ui/MoleculeLoader';

/** Animated loading panel shown while any validation/scoring request is in flight. */
export function LoadingPanel({ show, text }: { show: boolean; text: string }) {
  return (
    <AnimatePresence>
      {show && (
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -20 }}
          className="card-glass p-6"
        >
          <MoleculeLoader size="md" text={text} />
        </motion.div>
      )}
    </AnimatePresence>
  );
}

interface ErrorPanelProps {
  show: boolean;
  /** Pre-rendered error detail (parent builds it from the active error). */
  detail: ReactNode;
  onRetry: () => void;
  retryDisabled: boolean;
  onDismiss: () => void;
}

/** Animated error panel with retry/dismiss actions. */
export function ErrorPanel({ show, detail, onRetry, retryDisabled, onDismiss }: ErrorPanelProps) {
  return (
    <AnimatePresence>
      {show && (
        <motion.div
          initial={{ opacity: 0, scale: 0.95 }}
          animate={{ opacity: 1, scale: 1 }}
          exit={{ opacity: 0, scale: 0.95 }}
          className="card p-5 border-red-500/30"
          role="alert"
          aria-live="polite"
        >
          <div className="flex items-start gap-4">
            <div className="w-10 h-10 rounded-xl bg-red-500/10 flex items-center justify-center flex-shrink-0">
              <AlertTriangle className="w-5 h-5 text-red-500" />
            </div>
            <div className="flex-1 min-w-0">
              {detail}
              <div className="mt-3 flex flex-wrap gap-2">
                <ClayButton
                  variant="primary"
                  size="sm"
                  onClick={onRetry}
                  disabled={retryDisabled}
                  leftIcon={<RotateCcw className="w-3.5 h-3.5" />}
                >
                  Try again
                </ClayButton>
                <ClayButton variant="ghost" size="sm" onClick={onDismiss}>
                  Dismiss
                </ClayButton>
              </div>
            </div>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
