import { motion } from 'framer-motion';
import { useState, useEffect, useRef } from 'react';
import {
  RadialBarChart,
  RadialBar,
  ResponsiveContainer,
  PolarAngleAxis,
} from 'recharts';
import { ClayCard } from '../ui/ClayCard';
import { cn } from '../../lib/utils';
import { useThemeContext } from '../../contexts/ThemeContext';
import type { SubScoreDetail } from '../../types/dataset_intelligence';

// =============================================================================
// Types
// =============================================================================

interface SubScoreCardsProps {
  /** Array of 5 sub-score details from health audit. */
  subScores: SubScoreDetail[];
  /** Callback when a sub-score card is clicked. */
  onCardClick: (name: string) => void;
}

// =============================================================================
// Sub-score display names
// =============================================================================

const SUB_SCORE_LABELS: Record<string, string> = {
  parsability: 'Parsability',
  stereo: 'Stereo Completeness',
  uniqueness: 'Uniqueness',
  alerts: 'Alert Prevalence',
  std_consistency: 'Std. Consistency',
};

// =============================================================================
// Color thresholds (0-1 range mapped to same colors as HealthScoreGauge)
// =============================================================================

function getSubScoreColor(score: number) {
  if (score >= 0.8) {
    return {
      fill: '#b45309', // amber-700
      text: 'text-amber-600 dark:text-yellow-400',
      gradientId: 'subScoreExcellent',
    };
  }
  if (score >= 0.5) {
    return {
      fill: '#d97706', // amber-600
      text: 'text-amber-600 dark:text-amber-400',
      gradientId: 'subScoreFair',
    };
  }
  return {
    fill: '#dc2626', // red-600
    text: 'text-red-600 dark:text-red-400',
    gradientId: 'subScorePoor',
  };
}

// =============================================================================
// Mini gauge component
// =============================================================================

function MiniGauge({ score }: { score: number }) {
  const { isDark } = useThemeContext();
  const containerRef = useRef<HTMLDivElement>(null);
  const [isReady, setIsReady] = useState(false);

  useEffect(() => {
    const checkDimensions = () => {
      if (containerRef.current) {
        const { offsetWidth, offsetHeight } = containerRef.current;
        if (offsetWidth > 0 && offsetHeight > 0) {
          setIsReady(true);
        }
      }
    };
    checkDimensions();
    const timer = setTimeout(checkDimensions, 50);
    return () => clearTimeout(timer);
  }, []);

  const percent = Math.max(0, Math.min(100, score * 100));
  const color = getSubScoreColor(score);
  const backgroundFill = isDark ? '#374151' : '#e5e7eb';

  const data = [
    {
      name: 'SubScore',
      value: percent,
      fill: `url(#${color.gradientId})`,
    },
  ];

  return (
    <div ref={containerRef} className="relative" style={{ width: 80, height: 80 }}>
      {/* SVG Gradient Definitions */}
      <svg width="0" height="0" className="absolute">
        <defs>
          <linearGradient id="subScoreExcellent" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor={isDark ? '#fde68a' : '#fcd34d'} />
            <stop offset="100%" stopColor={isDark ? '#fbbf24' : '#b45309'} />
          </linearGradient>
          <linearGradient id="subScoreFair" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor={isDark ? '#fcd34d' : '#fbbf24'} />
            <stop offset="100%" stopColor={isDark ? '#f59e0b' : '#d97706'} />
          </linearGradient>
          <linearGradient id="subScorePoor" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor={isDark ? '#fca5a5' : '#f87171'} />
            <stop offset="100%" stopColor={isDark ? '#ef4444' : '#dc2626'} />
          </linearGradient>
        </defs>
      </svg>

      {isReady && (
        <ResponsiveContainer width="100%" height="100%" initialDimension={{ width: 1, height: 1 }}>
          <RadialBarChart
            innerRadius="60%"
            outerRadius="100%"
            data={data}
            startAngle={90}
            endAngle={-270}
            barSize={8}
          >
            <PolarAngleAxis
              type="number"
              domain={[0, 100]}
              angleAxisId={0}
              tick={false}
            />
            <RadialBar
              background={{ fill: backgroundFill }}
              dataKey="value"
              cornerRadius={8}
              animationDuration={1000}
              animationEasing="ease-out"
            />
          </RadialBarChart>
        </ResponsiveContainer>
      )}

      {/* Center text */}
      <div className="absolute inset-0 flex items-center justify-center">
        <span className={cn('text-sm font-semibold', getSubScoreColor(score).text)}>
          {Math.round(score * 100)}
        </span>
      </div>
    </div>
  );
}

// =============================================================================
// Component
// =============================================================================

/**
 * 5 mini gauge cards showing individual sub-scores.
 *
 * Layout:
 * - Mobile: horizontal scroll with flex overflow-x-auto
 * - Desktop: 5-column grid
 *
 * Each card contains an 80px mini RadialBarChart gauge,
 * sub-score name, numeric value, and weight percentage.
 *
 * Entry animation: staggered fade-in with y offset (0.08s per card, 0.3s ease-out).
 */
export function SubScoreCards({ subScores, onCardClick }: SubScoreCardsProps) {
  return (
    <div className="flex overflow-x-auto gap-4 pb-2 lg:grid lg:grid-cols-5 lg:gap-4 lg:overflow-visible lg:pb-0">
      {subScores.map((sub, idx) => {
        const displayName = SUB_SCORE_LABELS[sub.name] ?? sub.name;
        const scorePercent = Math.round(sub.score * 100);
        const weightPercent = Math.round(sub.weight * 100);

        return (
          <motion.div
            key={sub.name}
            initial={{ opacity: 0, y: 8 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: idx * 0.08, duration: 0.3, ease: 'easeOut' }}
            className="flex-shrink-0 min-w-[140px] lg:min-w-0"
          >
            <ClayCard
              variant="default"
              size="sm"
              hover
              className="flex flex-col items-center cursor-pointer"
              role="button"
              tabIndex={0}
              aria-label={`${displayName}: ${scorePercent} out of 100, weight ${weightPercent}%`}
              onClick={() => onCardClick(sub.name)}
              onKeyDown={(e: React.KeyboardEvent) => {
                if (e.key === 'Enter' || e.key === ' ') {
                  e.preventDefault();
                  onCardClick(sub.name);
                }
              }}
            >
              {/* Mini gauge */}
              <MiniGauge score={sub.score} />

              {/* Sub-score name */}
              <span className="text-xs text-[var(--color-text-primary)] font-medium mt-2 text-center leading-tight">
                {displayName}
              </span>

              {/* Numeric value */}
              <span className={cn('text-sm font-semibold mt-1', getSubScoreColor(sub.score).text)}>
                {scorePercent}
              </span>

              {/* Weight percentage */}
              <span className="text-xs text-[var(--color-text-secondary)] mt-0.5">
                {weightPercent}%
              </span>
            </ClayCard>
          </motion.div>
        );
      })}
    </div>
  );
}
