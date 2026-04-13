/**
 * ClusterMemberGrid
 *
 * Flex-wrap grid of mini 2D structure thumbnails for cluster members.
 * Each thumbnail is 80x60px rendered via MoleculeViewer.
 */

import { useState, useCallback, useRef } from 'react';
import { motion } from 'framer-motion';
import { Copy, Check } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import type { BatchResult } from '../../types/batch';

interface ClusterMemberGridProps {
  memberIndices: number[];
  results: BatchResult[];
  /** Index→SMILES map from clustering response (covers all molecules, not just current page) */
  smilesMap?: Record<string, string>;
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
  smilesMap,
  clusterId,
}: ClusterMemberGridProps) {
  const [copiedIdx, setCopiedIdx] = useState<number | null>(null);
  const copyTimeoutRef = useRef<ReturnType<typeof setTimeout>>();

  const handleCopy = useCallback((smiles: string | null, idx: number) => {
    if (!smiles) return;
    navigator.clipboard.writeText(smiles);
    setCopiedIdx(idx);
    if (copyTimeoutRef.current) clearTimeout(copyTimeoutRef.current);
    copyTimeoutRef.current = setTimeout(() => setCopiedIdx(null), 1500);
  }, []);

  // Build index->result lookup for efficient access (fallback for paginated results)
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
        const smiles =
          smilesMap?.[String(idx)] ??
          (() => {
            const mol = resultMap.get(idx);
            return mol?.standardization?.standardized_smiles || mol?.smiles || null;
          })();
        return (
          <motion.button
            key={idx}
            type="button"
            custom={i}
            initial="hidden"
            animate="visible"
            variants={thumbVariants}
            className="relative rounded-lg border border-[var(--color-border)] bg-white dark:bg-gray-900/50 overflow-hidden flex-shrink-0 cursor-pointer hover:ring-2 hover:ring-[var(--color-primary)]/30 transition-shadow group"
            title={smiles ? `${smiles}\nClick to copy` : `Molecule ${idx}`}
            style={{ width: 80, height: 60 }}
            onClick={() => handleCopy(smiles, idx)}
            disabled={!smiles}
          >
            <MoleculeViewer smiles={smiles} width={80} height={60} />
            {smiles && (
              <span className="absolute inset-0 flex items-center justify-center bg-black/50 opacity-0 group-hover:opacity-100 transition-opacity rounded-lg">
                {copiedIdx === idx ? (
                  <Check className="w-4 h-4 text-green-400" />
                ) : (
                  <Copy className="w-3 h-3 text-white" />
                )}
              </span>
            )}
          </motion.button>
        );
      })}
    </div>
  );
}
