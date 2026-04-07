import { useState, useId } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { Badge } from '../ui/Badge';
import { cn } from '../../lib/utils';

interface MetricCardProps {
  title: string;
  score: string | number;
  classification?: string;
  classificationVariant?: 'success' | 'warning' | 'error' | 'default';
  unit?: string;
  children?: React.ReactNode;
  defaultExpanded?: boolean;
  className?: string;
}

/**
 * Reusable expandable metric card for profiler panels.
 *
 * Shows score + classification badge in the header.
 * Click header to toggle the expandable detail section.
 * Conforms to UI-SPEC D-05 interaction and accessibility contracts.
 */
export function MetricCard({
  title,
  score,
  classification,
  classificationVariant = 'default',
  unit,
  children,
  defaultExpanded = false,
  className,
}: MetricCardProps) {
  const [expanded, setExpanded] = useState(defaultExpanded);
  const panelId = useId();

  const hasChildren = !!children;

  return (
    <ClayCard size="md" className={cn('overflow-hidden', className)}>
      {/* Header — always visible */}
      <button
        type="button"
        onClick={() => hasChildren && setExpanded(prev => !prev)}
        aria-expanded={hasChildren ? expanded : undefined}
        aria-controls={hasChildren ? panelId : undefined}
        className={cn(
          'w-full text-left',
          hasChildren && 'cursor-pointer'
        )}
      >
        <div className="flex items-start justify-between gap-3">
          {/* Left: title + score */}
          <div className="flex-1 min-w-0">
            <p className="text-sm font-semibold text-text-secondary font-display truncate">
              {title}
            </p>
            <div className="flex items-baseline gap-1.5 mt-1">
              <span className="text-2xl font-semibold tabular-nums text-text-primary font-display">
                {typeof score === 'number' && !isNaN(score) ? score : score}
              </span>
              {unit && (
                <span className="text-xs text-text-muted">{unit}</span>
              )}
            </div>
          </div>

          {/* Right: classification badge + chevron */}
          <div className="flex items-center gap-2 shrink-0 pt-0.5">
            {classification && (
              <Badge variant={classificationVariant} size="sm">
                {classification}
              </Badge>
            )}
            {hasChildren && (
              <motion.div
                animate={{ rotate: expanded ? 180 : 0 }}
                transition={{ duration: 0.2, ease: 'easeOut' }}
                className="text-text-muted"
              >
                <ChevronDown className="w-4 h-4" />
              </motion.div>
            )}
          </div>
        </div>
      </button>

      {/* Expandable detail panel */}
      {hasChildren && (
        <AnimatePresence initial={false}>
          {expanded && (
            <motion.div
              id={panelId}
              role="region"
              key="content"
              initial={{ height: 0, opacity: 0 }}
              animate={{ height: 'auto', opacity: 1 }}
              exit={{ height: 0, opacity: 0 }}
              transition={{ duration: 0.5, ease: 'easeOut' }}
              className="overflow-hidden"
            >
              <div className="pt-4 border-t border-[var(--color-border)] mt-4">
                {children}
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      )}
    </ClayCard>
  );
}
