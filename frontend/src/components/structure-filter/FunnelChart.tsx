import { useCallback } from 'react';
import { motion } from 'framer-motion';
import { STAGE_COLORS, type StageResult } from '../../types/structure_filter';
import { STAGE_CRITERIA } from './StructureFilterInput';

// =============================================================================
// Animation variants (D-05 locked — Pitfall 8: must use variants pattern)
// =============================================================================

const containerVariants = {
  hidden: {},
  visible: { transition: { staggerChildren: 0.18, delayChildren: 0.1 } },
};

const stageVariants = {
  hidden: { opacity: 0, scaleX: 0, y: -12 },
  visible: {
    opacity: 1,
    scaleX: 1,
    y: 0,
    transition: {
      duration: 0.5,
      ease: [0.34, 1.56, 0.64, 1] as [number, number, number, number], // spring-like overshoot
      opacity: { duration: 0.3 },
    },
  },
};

// =============================================================================
// Constants
// =============================================================================

/** Height of each trapezoid in SVG units. */
const STAGE_HEIGHT = 56;
/** Vertical gap between trapezoids. */
const STAGE_GAP = 16;
/** Total vertical slot per stage. */
const STAGE_SLOT = STAGE_HEIGHT + STAGE_GAP;
/** SVG canvas width. */
const SVG_WIDTH = 400;
/** Maximum trapezoid width (full-width stage). */
const MAX_TRAP_WIDTH = 360;
/** Minimum trapezoid width even at 0% pass rate per UI-SPEC. */
const MIN_TRAP_WIDTH = 40;
/** Trapezoid side angle offset (px) for the slanted sides — creates trapezoid shape. */
const TRAP_SKEW = 8;

// =============================================================================
// Props
// =============================================================================

interface FunnelChartProps {
  /** Stage-level statistics from the filter result. */
  stages: StageResult[];
  /** Total input count for width percentage calculation. */
  inputCount: number;
  /** Name of the currently selected stage (null = no selection). */
  selectedStage: string | null;
  /** Called when a stage is clicked. Pass null to deselect. */
  onStageClick: (stageName: string | null) => void;
}

// =============================================================================
// Helper: build SVG polygon points string for a trapezoid
// =============================================================================

function trapezoidPoints(
  x: number,
  y: number,
  topWidth: number,
  bottomWidth: number,
  height: number,
): string {
  const topLeft = x;
  const topRight = x + topWidth;
  const bottomLeft = x + (topWidth - bottomWidth) / 2;
  const bottomRight = bottomLeft + bottomWidth;
  return `${topLeft},${y} ${topRight},${y} ${bottomRight},${y + height} ${bottomLeft},${y + height}`;
}

// =============================================================================
// Component
// =============================================================================

/**
 * Custom SVG funnel chart for generative chemistry filter results.
 *
 * Renders colored trapezoids (one per funnel stage) with Framer Motion stagger
 * animation. Each trapezoid width represents the fraction of molecules that
 * passed that stage.
 *
 * Design decisions (locked):
 * - D-05: Framer Motion stagger via containerVariants/stageVariants pattern
 * - D-06: Stage colors from STAGE_COLORS in types/structure_filter.ts
 * - D-07: Click to select/deselect stage (re-click clears selection)
 * - D-08: Count and percentage labels below each trapezoid
 *
 * Accessibility: each stage group has role="button", aria-label, tabIndex=0,
 * and onKeyDown for Enter/Space activation.
 */
export function FunnelChart({
  stages,
  inputCount,
  selectedStage,
  onStageClick,
}: FunnelChartProps) {
  const numStages = stages.length;
  const svgHeight = numStages * STAGE_SLOT;

  const handleStageKeyDown = useCallback(
    (e: React.KeyboardEvent, stageName: string) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        // D-07: re-click same stage clears selection
        onStageClick(stageName === selectedStage ? null : stageName);
      }
    },
    [onStageClick, selectedStage],
  );

  return (
    <div className="w-full overflow-x-auto">
      <motion.svg
        viewBox={`0 0 ${SVG_WIDTH} ${svgHeight}`}
        width="100%"
        variants={containerVariants}
        initial="hidden"
        animate="visible"
        aria-label="Funnel chart showing molecules passing each filter stage"
        role="img"
        style={{ display: 'block', maxWidth: '100%' }}
      >
        {stages.map((stage, index) => {
          const y = index * STAGE_SLOT;
          const passRate =
            inputCount > 0
              ? stage.passed_count / inputCount
              : 0;

          // Trapezoid width: proportional to pass rate, min 40px per UI-SPEC (Pitfall 7 guard above)
          const topWidth = Math.max(passRate * MAX_TRAP_WIDTH, MIN_TRAP_WIDTH);
          // Bottom is slightly narrower to create the funnel taper
          const bottomWidth = Math.max(topWidth - TRAP_SKEW * 2, MIN_TRAP_WIDTH - TRAP_SKEW);
          const x = (SVG_WIDTH - topWidth) / 2;

          const fillColor = STAGE_COLORS[stage.stage_name] ?? '#94a3b8';
          const isSelected = selectedStage === stage.stage_name;

          const pct =
            inputCount > 0
              ? Math.round((stage.passed_count / inputCount) * 100)
              : 0;

          const criterion = STAGE_CRITERIA[stage.stage_name] ?? 'applies filters';
          const ariaLabel = `${stage.stage_name}: ${stage.passed_count} passed (${pct}%) — ${criterion}`;

          return (
            <motion.g
              key={stage.stage_name}
              variants={stageVariants}
              role="button"
              aria-label={ariaLabel}
              tabIndex={0}
              onClick={() =>
                // D-07: re-click same stage clears selection
                onStageClick(stage.stage_name === selectedStage ? null : stage.stage_name)
              }
              onKeyDown={(e) => handleStageKeyDown(e, stage.stage_name)}
              style={{ cursor: 'pointer', outline: 'none' }}
              whileHover={{ scale: 1.02 }}
              transition={{ duration: 0.15 }}
              filter={
                isSelected
                  ? 'drop-shadow(0 0 8px rgba(255,255,255,0.6))'
                  : undefined
              }
            >
              {/* Native SVG tooltip fallback */}
              <title>{`${stage.stage_name}: ${criterion}. ${stage.passed_count} of ${stage.input_count} passed (${pct}%).`}</title>

              {/* Trapezoid shape */}
              <polygon
                points={trapezoidPoints(x, y, topWidth, bottomWidth, STAGE_HEIGHT)}
                fill={fillColor}
                stroke={isSelected ? 'white' : 'none'}
                strokeWidth={isSelected ? 2 : 0}
                opacity={stage.enabled ? 1 : 0.4}
                rx="4"
              />

              {/* Stage name label — centered inside trapezoid */}
              <text
                x={SVG_WIDTH / 2}
                y={y + STAGE_HEIGHT / 2}
                textAnchor="middle"
                dominantBaseline="middle"
                fontSize="12"
                fontWeight="600"
                fill="white"
                style={{ filter: 'drop-shadow(0 1px 1px rgba(0,0,0,0.3))' }}
                pointerEvents="none"
              >
                {stage.stage_name}
              </text>

              {/* Count + percentage badge below trapezoid */}
              <motion.text
                x={SVG_WIDTH / 2}
                y={y + STAGE_HEIGHT + 10}
                textAnchor="middle"
                dominantBaseline="auto"
                fontSize="11"
                fill="var(--color-text-secondary, #6b7280)"
                pointerEvents="none"
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ delay: 0.15 + index * 0.18 + 0.3, duration: 0.3 }}
              >
                {stage.passed_count.toLocaleString()} passed ({pct}%)
              </motion.text>
            </motion.g>
          );
        })}
      </motion.svg>
    </div>
  );
}
