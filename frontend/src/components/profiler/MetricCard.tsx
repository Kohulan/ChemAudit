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
        {/* Title row — full width, wraps freely */}
        <div className="flex items-start justify-between gap-2 mb-2">
          <p className="text-xs font-medium text-text-muted uppercase tracking-wider leading-snug flex-1">
            {title}
          </p>
          {hasChildren && (
            <motion.div
              animate={{ rotate: expanded ? 180 : 0 }}
              transition={{ duration: 0.2, ease: 'easeOut' }}
              className="text-text-muted shrink-0 mt-0.5"
            >
              <ChevronDown className="w-3.5 h-3.5" />
            </motion.div>
          )}
        </div>

        {/* Score + badge row */}
        <div className="flex items-baseline justify-between gap-3">
          <div className="flex items-baseline gap-1.5">
            <span className="text-2xl font-semibold tabular-nums text-text-primary font-display">
              {score}
            </span>
            {unit && (
              <span className="text-xs text-text-muted">{unit}</span>
            )}
          </div>
          {classification && (
            <Badge variant={classificationVariant} size="sm">
              {classification}
            </Badge>
          )}
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
