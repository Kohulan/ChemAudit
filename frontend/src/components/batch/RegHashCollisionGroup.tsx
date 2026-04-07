/**
 * RegHashCollisionGroup
 *
 * Expandable card showing molecules that share the same registration hash.
 * Displays truncated hash in header, expands to show full hash and molecule list.
 * Uses Framer Motion for smooth expand/collapse animation.
 */

import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import type { RegistrationHashCollisionGroup as CollisionGroupType } from '../../types/analytics';
import type { BatchResult } from '../../types/batch';

interface RegHashCollisionGroupProps {
  /** Collision group data (shared hash, molecule indices, count) */
  group: CollisionGroupType;
  /** Batch results for looking up molecule SMILES */
  results: BatchResult[];
}

export function RegHashCollisionGroup({ group, results }: RegHashCollisionGroupProps) {
  const [isExpanded, setIsExpanded] = useState(false);

  const truncatedHash = group.hash.slice(0, 12);

  // Build molecule details from batch results
  const molecules = group.molecule_indices.map((idx) => {
    const result = results.find((r) => r.index === idx);
    return {
      index: idx,
      smiles: result?.smiles ?? 'N/A',
    };
  });

  return (
    <div
      className="rounded-xl border border-amber-500/20 bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400 overflow-hidden"
    >
      {/* Header (always visible) */}
      <button
        type="button"
        className="w-full flex items-center justify-between px-4 py-3 text-left hover:bg-amber-200/50 dark:hover:bg-amber-800/20 transition-colors"
        onClick={() => setIsExpanded(!isExpanded)}
        aria-expanded={isExpanded}
        aria-label={`Hash collision group: ${truncatedHash}, ${group.count} molecules`}
      >
        <div className="flex items-center gap-3">
          <span className="text-xs font-mono bg-amber-200/60 dark:bg-amber-800/40 px-2 py-0.5 rounded">
            {truncatedHash}...
          </span>
          <span className="text-sm font-medium">
            {group.count} molecules
          </span>
        </div>
        <motion.span
          animate={{ rotate: isExpanded ? 180 : 0 }}
          transition={{ duration: 0.2 }}
        >
          <ChevronDown className="w-4 h-4" />
        </motion.span>
      </button>

      {/* Expandable content */}
      <AnimatePresence initial={false}>
        {isExpanded && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{
              height: { duration: 0.25, ease: 'easeOut' },
              opacity: { duration: 0.2 },
            }}
            className="overflow-hidden"
          >
            <div className="px-4 pb-3 space-y-2 border-t border-amber-500/20">
              {/* Full hash */}
              <div className="pt-2">
                <span className="text-[10px] uppercase font-medium text-amber-600/70 dark:text-amber-500/70">
                  Full Hash
                </span>
                <p className="text-xs font-mono break-all mt-0.5">{group.hash}</p>
              </div>

              {/* Molecule list */}
              <div className="space-y-1">
                {molecules.map((mol) => (
                  <div
                    key={mol.index}
                    className="flex items-center gap-2 text-xs py-1 px-2 rounded-lg bg-amber-200/30 dark:bg-amber-800/20"
                  >
                    <span className="text-amber-600 dark:text-amber-500 font-medium flex-shrink-0">
                      (row {mol.index + 1})
                    </span>
                    <span className="font-mono text-amber-800 dark:text-amber-300 truncate" title={mol.smiles}>
                      {mol.smiles}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
