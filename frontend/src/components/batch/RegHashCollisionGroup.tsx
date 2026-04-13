/**
 * RegHashCollisionGroup
 *
 * Expandable card showing molecules that share the same registration hash.
 * Displays truncated hash in header, expands to show full hash and molecule list.
 * Clicking a molecule scrolls to it in the results table.
 * "Compare" button selects all group molecules and triggers the compare workflow.
 * Uses Framer Motion for smooth expand/collapse animation.
 */

import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, GitCompare, ArrowRight } from 'lucide-react';
import type { RegistrationHashCollisionGroup as CollisionGroupType, RegistrationHashMolecule } from '../../types/analytics';
import type { BatchResult } from '../../types/batch';

interface RegHashCollisionGroupProps {
  /** Collision group data (shared hash, molecule indices, count) */
  group: CollisionGroupType;
  /** Batch results for looking up molecule names (current page only) */
  results: BatchResult[];
  /** Full per-molecule registration data (all molecules, always available) */
  perMolecule: RegistrationHashMolecule[];
  /** Navigate to a molecule in the results table */
  onNavigateToMolecule?: (moleculeIndex: number) => void;
  /** Select specific molecule indices and trigger comparison */
  onCompareGroup?: (indices: number[]) => void;
}

export function RegHashCollisionGroup({ group, results, perMolecule, onNavigateToMolecule, onCompareGroup }: RegHashCollisionGroupProps) {
  const [isExpanded, setIsExpanded] = useState(false);

  const truncatedHash = group.hash.slice(0, 12);

  // Build molecule details: SMILES from perMolecule (always complete), name from results (current page, best-effort)
  const molecules = group.molecule_indices.map((idx) => {
    const regMol = perMolecule.find((m) => m.index === idx);
    const batchResult = results.find((r) => r.index === idx);
    return {
      index: idx,
      name: batchResult?.name || null,
      smiles: regMol?.smiles ?? batchResult?.smiles ?? '',
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
                  <button
                    key={mol.index}
                    type="button"
                    onClick={(e) => {
                      e.stopPropagation();
                      onNavigateToMolecule?.(mol.index);
                    }}
                    className="w-full flex items-center gap-2 text-xs py-1.5 px-2 rounded-lg bg-amber-200/30 dark:bg-amber-800/20 hover:bg-amber-300/40 dark:hover:bg-amber-700/30 transition-colors text-left group/mol cursor-pointer"
                  >
                    <span className="text-amber-600 dark:text-amber-500 font-semibold flex-shrink-0">
                      #{mol.index + 1}
                    </span>
                    {mol.name && (
                      <span className="text-amber-800 dark:text-amber-300 font-medium truncate max-w-[200px]" title={mol.name}>
                        {mol.name}
                      </span>
                    )}
                    <span className="font-mono text-amber-700/70 dark:text-amber-400/70 truncate flex-1" title={mol.smiles}>
                      {mol.smiles}
                    </span>
                    <ArrowRight className="w-3 h-3 opacity-0 group-hover/mol:opacity-100 transition-opacity flex-shrink-0 text-amber-600 dark:text-amber-400" />
                  </button>
                ))}
              </div>

              {/* Compare button */}
              {group.count >= 2 && onCompareGroup && (
                <button
                  type="button"
                  onClick={(e) => {
                    e.stopPropagation();
                    onCompareGroup(group.molecule_indices);
                  }}
                  className="w-full flex items-center justify-center gap-2 px-3 py-2 rounded-lg text-xs font-semibold bg-amber-600 text-white hover:bg-amber-700 dark:bg-amber-500 dark:text-amber-950 dark:hover:bg-amber-400 transition-colors"
                >
                  <GitCompare className="w-3.5 h-3.5" />
                  Compare {group.count} Molecules
                </button>
              )}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
