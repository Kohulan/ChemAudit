/**
 * ClusterSummaryBar
 *
 * Displays three metric badges summarizing clustering results:
 * cluster count, singleton count, and largest cluster size.
 */

import { motion } from 'framer-motion';

interface ClusterSummaryBarProps {
  clusterCount: number;
  singletonCount: number;
  largestClusterSize: number;
}

const badgeVariants = {
  hidden: { opacity: 0, scale: 0.9 },
  visible: (i: number) => ({
    opacity: 1,
    scale: 1,
    transition: { delay: i * 0.1, duration: 0.3, ease: 'easeOut' },
  }),
};

export function ClusterSummaryBar({
  clusterCount,
  singletonCount,
  largestClusterSize,
}: ClusterSummaryBarProps) {
  const badges = [
    { label: 'Clusters', value: clusterCount },
    { label: 'Singletons', value: singletonCount },
    { label: 'Largest Cluster', value: largestClusterSize },
  ];

  return (
    <div className="space-y-2">
      <div className="flex items-center gap-3 flex-wrap">
        {badges.map((badge, i) => (
          <motion.span
            key={badge.label}
            custom={i}
            initial="hidden"
            animate="visible"
            variants={badgeVariants}
            className="inline-flex flex-col items-center bg-[var(--color-surface-elevated)] rounded-xl px-4 py-2 border border-[var(--color-border)]"
          >
            <span className="text-sm font-semibold text-[var(--color-text-primary)]">
              {badge.value}
            </span>
            <span className="text-xs text-[var(--color-text-muted)]">
              {badge.label}
            </span>
          </motion.span>
        ))}
      </div>
      <p className="text-xs text-[var(--color-text-muted)]">
        {clusterCount} clusters found, {singletonCount} singletons, largest cluster has{' '}
        {largestClusterSize} molecules
      </p>
    </div>
  );
}
