import { type ReactNode, useId, useState, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import { cn } from '../../lib/utils';

interface DrillDownSectionProps {
  /** Section heading text */
  title: string;
  /** Optional icon rendered before the title */
  icon?: ReactNode;
  /** Summary badge shown on the collapsed header (hidden when expanded) */
  summary?: ReactNode;
  /** Show a loading skeleton in place of the summary badge */
  summaryLoading?: boolean;
  /** Whether the section starts expanded */
  defaultOpen?: boolean;
  /** Called when the section is toggled, receives the new open state */
  onToggle?: (isOpen: boolean) => void;
  /** Section body content */
  children: ReactNode;
  /** Additional classes on the outer wrapper */
  className?: string;
}

/**
 * Shared accordion component for drill-down sections.
 *
 * Used by Profiler, Safety, and Diagnostics accordions in SingleValidation.
 * Collapses by default; expands with a 300ms framer-motion animation.
 */
export function DrillDownSection({
  title,
  icon,
  summary,
  summaryLoading = false,
  defaultOpen = false,
  onToggle,
  children,
  className,
}: DrillDownSectionProps) {
  const [isOpen, setIsOpen] = useState(defaultOpen);
  const panelId = useId();

  const toggle = useCallback(() => {
    const next = !isOpen;
    setIsOpen(next);
    onToggle?.(next);
  }, [isOpen, onToggle]);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent<HTMLButtonElement>) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        toggle();
      }
    },
    [toggle],
  );

  return (
    <div
      className={cn(
        'rounded-xl border border-border/50 bg-surface-secondary overflow-hidden',
        className,
      )}
    >
      {/* Header button */}
      <button
        type="button"
        className={cn(
          'flex w-full items-center gap-3 px-4 py-3 text-left',
          'cursor-pointer select-none',
          'transition-colors duration-150 hover:bg-surface-secondary/80',
          'focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-primary/50',
        )}
        aria-expanded={isOpen}
        aria-controls={panelId}
        onClick={toggle}
        onKeyDown={handleKeyDown}
      >
        {/* Chevron */}
        <motion.div
          animate={{ rotate: isOpen ? 180 : 0 }}
          transition={{ duration: 0.3 }}
          className="flex-shrink-0 text-text-tertiary"
        >
          <ChevronDown className="h-4 w-4" />
        </motion.div>

        {/* Icon (optional) */}
        {icon && <span className="flex-shrink-0">{icon}</span>}

        {/* Title */}
        <span className="flex-1 text-sm font-medium text-text-primary">{title}</span>

        {/* Summary badge (collapsed only) */}
        {!isOpen && summary && !summaryLoading && (
          <span className="flex-shrink-0 text-xs text-text-tertiary">{summary}</span>
        )}

        {/* Summary loading skeleton */}
        {!isOpen && summaryLoading && (
          <span className="h-4 w-16 animate-pulse rounded bg-border/50" />
        )}
      </button>

      {/* Expandable panel */}
      <AnimatePresence initial={false}>
        {isOpen && (
          <motion.div
            id={panelId}
            role="region"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.3, ease: 'easeInOut' }}
            className="overflow-hidden"
          >
            <div className="px-4 pb-4">{children}</div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
