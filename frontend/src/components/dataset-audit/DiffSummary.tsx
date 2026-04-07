import { motion } from 'framer-motion';

// =============================================================================
// Types
// =============================================================================

interface DiffSummaryProps {
  /** Number of molecules added in comparison dataset. */
  addedCount: number;
  /** Number of molecules removed from comparison dataset. */
  removedCount: number;
  /** Number of molecules modified between datasets. */
  modifiedCount: number;
  /** Number of molecules unchanged between datasets. */
  unchangedCount: number;
}

// =============================================================================
// Badge config
// =============================================================================

const BADGES = [
  {
    key: 'added',
    prefix: '+',
    label: 'Added',
    classes: 'bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400',
  },
  {
    key: 'removed',
    prefix: '-',
    label: 'Removed',
    classes: 'bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-400',
  },
  {
    key: 'modified',
    prefix: '~',
    label: 'Modified',
    classes: 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400',
  },
  {
    key: 'unchanged',
    prefix: '',
    label: 'Unchanged',
    classes: 'bg-stone-100 text-stone-600 dark:bg-stone-800/30 dark:text-stone-400',
  },
] as const;

// =============================================================================
// Component
// =============================================================================

/**
 * Row of 4 inline count badges summarizing dataset diff results.
 *
 * Per UI-SPEC D-12:
 * - Added (green): "+{N} Added"
 * - Removed (red): "-{N} Removed"
 * - Modified (amber): "~{N} Modified"
 * - Unchanged (stone): "{N} Unchanged"
 *
 * Entry animation: opacity 0->1, scale 0.9->1, stagger 0.1s, 0.3s ease-out.
 */
export function DiffSummary({
  addedCount,
  removedCount,
  modifiedCount,
  unchangedCount,
}: DiffSummaryProps) {
  const counts: Record<string, number> = {
    added: addedCount,
    removed: removedCount,
    modified: modifiedCount,
    unchanged: unchangedCount,
  };

  return (
    <div className="flex flex-wrap gap-3">
      {BADGES.map((badge, index) => (
        <motion.span
          key={badge.key}
          initial={{ opacity: 0, scale: 0.9 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{
            duration: 0.3,
            ease: 'easeOut',
            delay: index * 0.1,
          }}
          className={`inline-flex items-center px-3 py-1.5 text-sm font-semibold rounded-full ${badge.classes}`}
        >
          {badge.prefix}
          {counts[badge.key]} {badge.label}
        </motion.span>
      ))}
    </div>
  );
}
