/**
 * ClusterMemberGrid
 *
 * Flex-wrap grid of mini 2D structure thumbnails for cluster members.
 * Each thumbnail is 80x60px rendered via MoleculeViewer.
 */

import { motion } from 'framer-motion';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import type { BatchResult } from '../../types/batch';

interface ClusterMemberGridProps {
  memberIndices: number[];
  results: BatchResult[];
  clusterId?: number;
}

const thumbVariants = {
  hidden: { opacity: 0 },
  visible: (i: number) => ({
    opacity: 1,
    transition: { delay: i * 0.03, duration: 0.2, ease: 'easeOut' },
  }),
};

export function ClusterMemberGrid({
  memberIndices,
  results,
  clusterId,
}: ClusterMemberGridProps) {
  // Build index->result lookup for efficient access
  const resultMap = new Map<number, BatchResult>();
  for (const r of results) {
    resultMap.set(r.index, r);
  }

  return (
    <div
      className="flex flex-wrap gap-2 p-2"
      aria-label={`${memberIndices.length} molecules in cluster ${clusterId ?? '?'}`}
    >
      {memberIndices.map((idx, i) => {
        const mol = resultMap.get(idx);
        const smiles = mol?.standardization?.standardized_smiles || mol?.smiles || null;
        return (
          <motion.div
            key={idx}
            custom={i}
            initial="hidden"
            animate="visible"
            variants={thumbVariants}
            className="rounded-lg border border-[var(--color-border)] bg-white dark:bg-gray-900/50 overflow-hidden flex-shrink-0"
            title={smiles || `Molecule ${idx}`}
            style={{ width: 80, height: 60 }}
          >
            <MoleculeViewer smiles={smiles} width={80} height={60} />
          </motion.div>
        );
      })}
    </div>
  );
}
