/**
 * StereoisomerList Component
 *
 * Collapsible list of enumerated stereoisomer SMILES strings with copy buttons.
 * Shows count summary when collapsed, full list when expanded.
 * Displays cap-exceeded warning when the enumeration limit was hit.
 */
import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, AlertTriangle } from 'lucide-react';
import { CopyButton } from '../ui/CopyButton';
import { cn } from '../../lib/utils';

interface StereoisomerListProps {
  smiles: string[];
  cap: number;
  capExceeded: boolean;
}

/**
 * Renders a collapsible list of enumerated stereoisomer SMILES.
 */
export function StereoisomerList({ smiles, cap, capExceeded }: StereoisomerListProps) {
  const [expanded, setExpanded] = useState(false);

  const count = smiles.length;

  return (
    <div className="mt-3">
      {/* Cap exceeded warning */}
      {capExceeded && (
        <div className="mb-2 flex items-start gap-2 p-2 rounded-lg bg-amber-500/10 border border-amber-500/20">
          <AlertTriangle className="w-3.5 h-3.5 text-amber-500 flex-shrink-0 mt-0.5" />
          <p className="text-xs text-amber-600 dark:text-amber-400">
            Showing count only — {cap} cap exceeded. Too many stereoisomers to enumerate.
          </p>
        </div>
      )}

      {/* Header / Toggle */}
      {!capExceeded && count > 0 && (
        <button
          onClick={() => setExpanded((prev) => !prev)}
          aria-expanded={expanded}
          className={cn(
            'flex items-center gap-2 w-full text-left',
            'text-xs font-medium text-[var(--color-text-secondary)]',
            'hover:text-[var(--color-text-primary)] transition-colors',
            'px-2 py-1.5 rounded-lg',
            'hover:bg-[var(--color-surface-sunken)]'
          )}
        >
          <ChevronDown
            className={cn(
              'w-3.5 h-3.5 text-[var(--color-text-muted)] transition-transform duration-200',
              expanded && 'rotate-180'
            )}
          />
          <span>
            {count} stereoisomer{count !== 1 ? 's' : ''}
          </span>
        </button>
      )}

      {capExceeded && count === 0 && (
        <p className="text-xs text-[var(--color-text-muted)] px-2">
          {cap}+ possible stereoisomers (count only — enumeration skipped)
        </p>
      )}

      {/* Expanded SMILES list */}
      <AnimatePresence>
        {expanded && !capExceeded && count > 0 && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: 'auto' }}
            exit={{ opacity: 0, height: 0 }}
            transition={{ duration: 0.2 }}
            className="overflow-hidden"
          >
            <div className="mt-1 space-y-1 max-h-48 overflow-y-auto pr-1">
              {smiles.map((smi, idx) => (
                <div
                  key={idx}
                  className={cn(
                    'flex items-center gap-2',
                    'px-2 py-1 rounded-lg',
                    'bg-[var(--color-surface-sunken)]',
                    'border border-[var(--color-border)]'
                  )}
                >
                  <span className="text-[10px] font-medium text-[var(--color-text-muted)] w-5 flex-shrink-0 text-right">
                    {idx + 1}.
                  </span>
                  <code className="flex-1 text-xs text-[var(--color-text-secondary)] font-mono break-all min-w-0">
                    {smi}
                  </code>
                  <CopyButton text={smi} size={12} />
                </div>
              ))}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
