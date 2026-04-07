import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, Box, AlertCircle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { PMITernaryPlot } from './PMITernaryPlot';
import type { Shape3DResult } from '../../types/profiler';
import { cn } from '../../lib/utils';

interface Shape3DPanelProps {
  smiles: string;
  compute3DShape: (smiles: string) => Promise<Shape3DResult | null>;
}

const shapeClassLabels: Record<string, { label: string; variant: 'success' | 'warning' | 'default' }> = {
  rod: { label: 'Rod', variant: 'warning' },
  disc: { label: 'Disc', variant: 'default' },
  sphere: { label: 'Sphere', variant: 'success' },
};

/**
 * Shape & 3D Descriptors panel.
 *
 * Starts collapsed per D-26 (lazy compute). On expand, shows "Compute 3D Shape"
 * button. API call fires only on button click. Renders PMI ternary SVG plot
 * and descriptor values on success. Shows graceful message on conformer failure.
 */
export function Shape3DPanel({ smiles, compute3DShape }: Shape3DPanelProps) {
  const [expanded, setExpanded] = useState(false);
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<Shape3DResult | null>(null);
  const [computed, setComputed] = useState(false);
  const panelId = 'shape-3d-panel';

  const handleCompute = async () => {
    setLoading(true);
    try {
      const data = await compute3DShape(smiles);
      setResult(data);
    } finally {
      setLoading(false);
      setComputed(true);
    }
  };

  const shapeInfo = result ? (shapeClassLabels[result.shape_class] ?? { label: result.shape_class, variant: 'default' as const }) : null;

  return (
    <ClayCard size="md" className="overflow-hidden">
      {/* Section header — collapsible toggle */}
      <button
        type="button"
        onClick={() => setExpanded((prev) => !prev)}
        aria-expanded={expanded}
        aria-controls={panelId}
        className="w-full text-left cursor-pointer"
      >
        <div className="flex items-center justify-between gap-3">
          <div className="flex items-center gap-2">
            <Box className="w-4 h-4 text-text-muted flex-shrink-0" />
            <h3 className="text-lg font-semibold text-text-primary font-display">
              Shape &amp; 3D Descriptors
            </h3>
          </div>
          <motion.div
            animate={{ rotate: expanded ? 180 : 0 }}
            transition={{ duration: 0.2, ease: 'easeOut' }}
            className="text-text-muted flex-shrink-0"
          >
            <ChevronDown className="w-5 h-5" />
          </motion.div>
        </div>
      </button>

      {/* Expandable content */}
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
              {/* Not yet computed — show compute button */}
              {!computed && !loading && (
                <div className="flex flex-col items-center gap-3 py-4">
                  <p className="text-sm text-text-secondary text-center">
                    3D shape descriptors are computed on demand via ETKDGv3 conformer generation.
                  </p>
                  <ClayButton
                    variant="default"
                    size="md"
                    onClick={handleCompute}
                  >
                    Compute 3D Shape
                  </ClayButton>
                </div>
              )}

              {/* Loading state */}
              {loading && (
                <div className="flex flex-col items-center gap-3 py-6">
                  <motion.div
                    animate={{ rotate: 360 }}
                    transition={{ duration: 1.5, repeat: Infinity, ease: 'linear' }}
                    className="w-8 h-8"
                  >
                    <svg viewBox="0 0 24 24" fill="none" className="w-8 h-8">
                      <circle className="opacity-25" cx="12" cy="12" r="10" stroke="var(--color-primary)" strokeWidth="3" />
                      <path
                        className="opacity-75"
                        fill="var(--color-primary)"
                        d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                      />
                    </svg>
                  </motion.div>
                  <p className="text-sm text-text-muted">Computing 3D conformer...</p>
                </div>
              )}

              {/* Computed — conformer failed */}
              {computed && !loading && result && result['3d_conformer_failed'] && (
                <div className="flex items-start gap-2 py-2 text-text-muted">
                  <AlertCircle className="w-4 h-4 mt-0.5 flex-shrink-0" />
                  <p className="text-sm">
                    3D shape descriptors unavailable — conformer generation failed after 3 attempts.
                    All other metrics are shown.
                  </p>
                </div>
              )}

              {/* Computed — null result (API error) */}
              {computed && !loading && result === null && (
                <div className="flex items-start gap-2 py-2 text-text-muted">
                  <AlertCircle className="w-4 h-4 mt-0.5 flex-shrink-0" />
                  <p className="text-sm">
                    3D shape descriptors unavailable — conformer generation failed after 3 attempts.
                    All other metrics are shown.
                  </p>
                </div>
              )}

              {/* Computed — success */}
              {computed && !loading && result && !result['3d_conformer_failed'] && (
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6 items-start">
                  {/* PMI ternary plot */}
                  <div className="flex flex-col items-center gap-2">
                    <PMITernaryPlot
                      npr1={result.npr1}
                      npr2={result.npr2}
                      shapeClass={result.shape_class}
                    />
                    {shapeInfo && (
                      <div className="flex items-center gap-2">
                        <span className="text-sm text-text-secondary">Shape class:</span>
                        <Badge variant={shapeInfo.variant} size="sm">
                          {shapeInfo.label}
                        </Badge>
                      </div>
                    )}
                  </div>

                  {/* Descriptor values */}
                  <div className="space-y-4">
                    <div>
                      <p className="text-xs font-semibold text-text-muted uppercase tracking-wide mb-2">
                        Principal Moments of Inertia
                      </p>
                      <div className="grid grid-cols-3 gap-2">
                        {(['pmi1', 'pmi2', 'pmi3'] as const).map((key) => (
                          <div key={key} className={cn('text-center p-2 rounded-lg bg-surface-sunken')}>
                            <p className="text-xs text-text-muted uppercase">{key.toUpperCase()}</p>
                            <p className="text-sm font-semibold tabular-nums text-text-primary mt-0.5">
                              {result[key].toFixed(3)}
                            </p>
                          </div>
                        ))}
                      </div>
                    </div>

                    <div>
                      <p className="text-xs font-semibold text-text-muted uppercase tracking-wide mb-2">
                        Normalized PMI Ratios &amp; PBF
                      </p>
                      <div className="grid grid-cols-3 gap-2">
                        {[
                          { key: 'npr1', label: 'NPR1' },
                          { key: 'npr2', label: 'NPR2' },
                          { key: 'pbf', label: 'PBF' },
                        ].map(({ key, label }) => (
                          <div key={key} className={cn('text-center p-2 rounded-lg bg-surface-sunken')}>
                            <p className="text-xs text-text-muted uppercase">{label}</p>
                            <p className="text-sm font-semibold tabular-nums text-text-primary mt-0.5">
                              {result[key as keyof Shape3DResult] !== undefined
                                ? Number(result[key as keyof Shape3DResult]).toFixed(3)
                                : '—'}
                            </p>
                          </div>
                        ))}
                      </div>
                    </div>

                    {/* Recompute button */}
                    <div className="pt-1">
                      <ClayButton
                        variant="ghost"
                        size="sm"
                        onClick={handleCompute}
                      >
                        Recompute
                      </ClayButton>
                    </div>
                  </div>
                </div>
              )}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </ClayCard>
  );
}
