import { useState, useEffect, useRef } from 'react';
import {
  RadialBarChart,
  RadialBar,
  ResponsiveContainer,
  PolarAngleAxis,
} from 'recharts';
import { cn } from '../../lib/utils';
import { useThemeContext } from '../../contexts/ThemeContext';

// =============================================================================
// Types
// =============================================================================

interface HealthScoreGaugeProps {
  /** Composite health score 0-100. */
  score: number;
  /** Optional CSS class name. */
  className?: string;
}

// =============================================================================
// Color thresholds (per UI-SPEC)
// =============================================================================

function getScoreColor(score: number) {
  if (score >= 80) {
    return {
      fill: '#b45309', // amber-700
      text: 'text-amber-600 dark:text-yellow-400',
      bg: 'bg-yellow-500/10 dark:bg-yellow-400/15',
      label: 'Excellent' as const,
      gradientId: 'healthGaugeExcellent',
    };
  }
  if (score >= 50) {
    return {
      fill: '#d97706', // amber-600
      text: 'text-amber-600 dark:text-amber-400',
      bg: 'bg-amber-500/10 dark:bg-amber-400/15',
      label: 'Fair' as const,
      gradientId: 'healthGaugeFair',
    };
  }
  return {
    fill: '#dc2626', // red-600
    text: 'text-red-600 dark:text-red-400',
    bg: 'bg-red-500/10 dark:bg-red-400/15',
    label: 'Poor' as const,
    gradientId: 'healthGaugePoor',
  };
}

/** Contextual message below the gauge based on score label. */
function getContextMessage(label: 'Excellent' | 'Fair' | 'Poor'): string {
  switch (label) {
    case 'Excellent':
      return 'This dataset is well-curated and ready for most ML workflows.';
    case 'Fair':
      return 'Some quality issues detected. Review sub-scores for targeted curation.';
    case 'Poor':
      return 'Significant quality issues found. Curation is strongly recommended before use.';
  }
}

// =============================================================================
// Component
// =============================================================================

/**
 * Large 200px radial gauge showing the composite 0-100 health score.
 *
 * Color thresholds per UI-SPEC:
 * - score >= 80 -> amber-700 (#b45309) "Excellent"
 * - score >= 50 -> amber-600 (#d97706) "Fair"
 * - score <  50 -> red-600   (#dc2626) "Poor"
 *
 * Reuses the RadialBarChart pattern from ScoreGauge.tsx.
 */
export function HealthScoreGauge({ score, className }: HealthScoreGaugeProps) {
  const { isDark } = useThemeContext();
  const containerRef = useRef<HTMLDivElement>(null);
  const [isReady, setIsReady] = useState(false);

  // Delay chart render until container has valid dimensions
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

  const clamped = Math.max(0, Math.min(100, score));
  const color = getScoreColor(clamped);
  const backgroundFill = isDark ? '#374151' : '#e5e7eb';

  const data = [
    {
      name: 'Dataset Health',
      value: clamped,
      fill: `url(#${color.gradientId})`,
    },
  ];

  return (
    <div
      className={cn('flex flex-col items-center', className)}
      aria-label={`Dataset health score: ${Math.round(clamped)} out of 100, rated ${color.label}`}
    >
      {/* Gauge */}
      <div ref={containerRef} className="relative" style={{ width: 200, height: 200 }}>
        {/* SVG Gradient Definitions */}
        <svg width="0" height="0" className="absolute">
          <defs>
            <linearGradient id="healthGaugeExcellent" x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" stopColor={isDark ? '#fde68a' : '#fcd34d'} />
              <stop offset="100%" stopColor={isDark ? '#fbbf24' : '#b45309'} />
            </linearGradient>
            <linearGradient id="healthGaugeFair" x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" stopColor={isDark ? '#fcd34d' : '#fbbf24'} />
              <stop offset="100%" stopColor={isDark ? '#f59e0b' : '#d97706'} />
            </linearGradient>
            <linearGradient id="healthGaugePoor" x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" stopColor={isDark ? '#fca5a5' : '#f87171'} />
              <stop offset="100%" stopColor={isDark ? '#ef4444' : '#dc2626'} />
            </linearGradient>
          </defs>
        </svg>

        {isReady && (
          <ResponsiveContainer width="100%" height="100%" initialDimension={{ width: 1, height: 1 }}>
            <RadialBarChart
              innerRadius="65%"
              outerRadius="100%"
              data={data}
              startAngle={90}
              endAngle={-270}
              barSize={14}
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
                cornerRadius={10}
                animationDuration={1000}
                animationEasing="ease-out"
              />
            </RadialBarChart>
          </ResponsiveContainer>
        )}

        {/* Center text */}
        <div className="absolute inset-0 flex flex-col items-center justify-center">
          <span className={cn('text-2xl font-semibold font-display', color.text)}>
            {Math.round(clamped)}
          </span>
          <span className="text-xs text-text-muted uppercase tracking-wide">Score</span>
        </div>
      </div>

      {/* Heading */}
      <h2 className="text-base font-semibold font-display mt-3 text-[var(--color-text-primary)]">
        Dataset Health Score
      </h2>

      {/* Label badge */}
      <div className={cn('mt-2 px-4 py-1.5 rounded-full', color.bg)}>
        <span className={cn('text-sm font-medium', color.text)}>{color.label}</span>
      </div>

      {/* Context message */}
      <p className="mt-2 text-xs text-[var(--color-text-secondary)] max-w-[280px] text-center leading-relaxed">
        {getContextMessage(color.label)}
      </p>
    </div>
  );
}
