import { useState, useEffect, useRef } from 'react';
import { RadialBarChart, RadialBar, ResponsiveContainer, PolarAngleAxis } from 'recharts';
import { CalculationTooltip } from '../ui/Tooltip';
import { cn } from '../../lib/utils';
import { useThemeContext } from '../../contexts/ThemeContext';

interface ScoreChartProps {
  /** Score value (0-100) */
  score: number;
  /** Label for the score */
  label: string;
  /** Size of the chart */
  size?: number;
  /** Calculation method for tooltip */
  calculation?: string;
  /** Interpretation for tooltip */
  interpretation?: string;
  /** Show as compact version */
  compact?: boolean;
  /** Color variant: 'warm' (amber/yellow) or 'cool' (blue/cyan) */
  variant?: 'warm' | 'cool';
}

type ScoreColorConfig = {
  fill: string;
  text: string;
  bg: string;
  label: string;
  gradientId: string;
  startColor: string;
  endColor: string;
};

/**
 * Get color configuration based on score value
 */
function getScoreColor(score: number, isDark: boolean, variant: 'warm' | 'cool' = 'warm'): ScoreColorConfig {
  if (variant === 'cool') return getCoolColor(score, isDark);
  return getWarmColor(score, isDark);
}

function getWarmColor(score: number, isDark: boolean): ScoreColorConfig {
  if (score >= 80) {
    return {
      fill: isDark ? '#fcd34d' : '#b45309',
      text: 'text-amber-600 dark:text-yellow-400',
      bg: 'bg-yellow-500/10 dark:bg-yellow-400/15',
      label: 'Excellent',
      gradientId: 'scoreGradientExcellent',
      startColor: isDark ? '#fde68a' : '#fcd34d',
      endColor: isDark ? '#fbbf24' : '#b45309',
    };
  }
  if (score >= 50) {
    return {
      fill: isDark ? '#fbbf24' : '#d97706',
      text: 'text-amber-600 dark:text-amber-400',
      bg: 'bg-amber-500/10 dark:bg-amber-400/15',
      label: 'Fair',
      gradientId: 'scoreGradientFair',
      startColor: isDark ? '#fcd34d' : '#fbbf24',
      endColor: isDark ? '#f59e0b' : '#d97706',
    };
  }
  return {
    fill: isDark ? '#f87171' : '#dc2626',
    text: 'text-red-600 dark:text-red-400',
    bg: 'bg-red-500/10 dark:bg-red-400/15',
    label: 'Poor',
    gradientId: 'scoreGradientPoor',
    startColor: isDark ? '#fca5a5' : '#f87171',
    endColor: isDark ? '#ef4444' : '#dc2626',
  };
}

function getCoolColor(score: number, isDark: boolean): ScoreColorConfig {
  // Matches the app's primary "Laboratory Crimson" (rose-500 → red-600 → rose-700)
  if (score >= 80) {
    return {
      fill: isDark ? '#e11d48' : '#be123c',
      text: 'text-rose-600 dark:text-rose-400',
      bg: 'bg-rose-500/10 dark:bg-rose-400/15',
      label: 'Excellent',
      gradientId: 'scoreGradientCoolExcellent',
      startColor: isDark ? '#f43f5e' : '#f43f5e',  // rose-500
      endColor: isDark ? '#be123c' : '#9f1239',     // rose-700 → rose-800
    };
  }
  if (score >= 50) {
    return {
      fill: isDark ? '#dc2626' : '#b91c1c',
      text: 'text-red-600 dark:text-red-400',
      bg: 'bg-red-500/10 dark:bg-red-400/15',
      label: 'Fair',
      gradientId: 'scoreGradientCoolFair',
      startColor: isDark ? '#ef4444' : '#ef4444',   // red-500
      endColor: isDark ? '#991b1b' : '#7f1d1d',     // red-800 → red-900
    };
  }
  return {
    fill: isDark ? '#78716c' : '#57534e',
    text: 'text-stone-600 dark:text-stone-400',
    bg: 'bg-stone-500/10 dark:bg-stone-400/15',
    label: 'Poor',
    gradientId: 'scoreGradientCoolPoor',
    startColor: isDark ? '#a8a29e' : '#78716c',
    endColor: isDark ? '#57534e' : '#44403c',
  };
}

/**
 * Radial score chart with gradient fill and interactive tooltip
 */
export function ScoreChart({
  score,
  label,
  size = 160,
  calculation,
  interpretation,
  compact = false,
  variant = 'warm',
}: ScoreChartProps) {
  const { isDark } = useThemeContext();
  const containerRef = useRef<HTMLDivElement>(null);
  const [isReady, setIsReady] = useState(false);
  const [animated, setAnimated] = useState(false);
  const clampedScore = Math.max(0, Math.min(100, score));
  const color = getScoreColor(clampedScore, isDark, variant);
  const backgroundFill = isDark ? '#374151' : '#e5e7eb';

  // Trigger mount animation for compact gauge
  useEffect(() => {
    if (compact) {
      const timer = setTimeout(() => setAnimated(true), 50);
      return () => clearTimeout(timer);
    }
  }, [compact]);

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
    
    // Check immediately and after a short delay for animations
    checkDimensions();
    const timer = setTimeout(checkDimensions, 50);
    return () => clearTimeout(timer);
  }, []);

  const data = [
    {
      name: label,
      value: clampedScore,
      fill: `url(#${color.gradientId})`,
    },
  ];

  // Custom SVG gauge for compact mode — glow, gradient arc, tick marks, mount animation
  const renderCompactGauge = () => {
    const cx = size / 2;
    const cy = size / 2;
    const strokeWidth = 7;
    const radius = (size / 2) - strokeWidth - 4;
    const circumference = 2 * Math.PI * radius;
    const targetArc = circumference * (clampedScore / 100);
    const arcLength = animated ? targetArc : 0;
    const tickCount = 20;
    const tickRadius = radius + strokeWidth / 2 + 3;
    const glowId = `glow-${color.gradientId}`;
    const gradId = `grad-${color.gradientId}`;

    return (
      <div className="relative" style={{ width: size, height: size }}>
        <svg width={size} height={size} className="-rotate-90">
          <defs>
            <linearGradient id={gradId} x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stopColor={color.startColor} />
              <stop offset="100%" stopColor={color.endColor} />
            </linearGradient>
            <filter id={glowId} x="-30%" y="-30%" width="160%" height="160%">
              <feGaussianBlur stdDeviation="3" result="blur" />
              <feFlood floodColor={color.endColor} floodOpacity="0.4" result="color" />
              <feComposite in="color" in2="blur" operator="in" result="glow" />
              <feMerge>
                <feMergeNode in="glow" />
                <feMergeNode in="SourceGraphic" />
              </feMerge>
            </filter>
          </defs>

          {/* Tick marks */}
          {Array.from({ length: tickCount }, (_, i) => {
            const angle = (i / tickCount) * 360;
            const rad = (angle * Math.PI) / 180;
            const filled = animated && i / tickCount <= clampedScore / 100;
            const x1 = cx + (tickRadius - 2) * Math.cos(rad);
            const y1 = cy + (tickRadius - 2) * Math.sin(rad);
            const x2 = cx + (tickRadius + 2) * Math.cos(rad);
            const y2 = cy + (tickRadius + 2) * Math.sin(rad);
            return (
              <line
                key={i}
                x1={x1} y1={y1} x2={x2} y2={y2}
                stroke={filled ? color.endColor : (isDark ? '#374151' : '#d1d5db')}
                strokeWidth={1.5}
                strokeLinecap="round"
                opacity={filled ? 0.7 : 0.3}
                style={{ transition: `all 800ms ease-out ${i * 30}ms` }}
              />
            );
          })}

          {/* Track */}
          <circle
            cx={cx} cy={cy} r={radius}
            fill="none"
            stroke={backgroundFill}
            strokeWidth={strokeWidth}
            strokeLinecap="round"
          />

          {/* Score arc with glow — animated from 0 on mount */}
          <circle
            cx={cx} cy={cy} r={radius}
            fill="none"
            stroke={`url(#${gradId})`}
            strokeWidth={strokeWidth}
            strokeLinecap="round"
            strokeDasharray={`${arcLength} ${circumference}`}
            filter={`url(#${glowId})`}
            style={{ transition: 'stroke-dasharray 1s cubic-bezier(0.4, 0, 0.2, 1)' }}
          />
        </svg>

        {/* Center content */}
        <div className="absolute inset-0 flex flex-col items-center justify-center rotate-0">
          <span className={cn('text-2xl font-bold tabular-nums', color.text)}>
            {Math.round(clampedScore)}
          </span>
        </div>
      </div>
    );
  };

  const chartContent = compact ? renderCompactGauge() : (
    <div ref={containerRef} className="relative p-2" style={{ width: size, height: size }}>
      {/* SVG Gradient Definitions */}
      <svg width="0" height="0" className="absolute">
        <defs>
          <linearGradient id="scoreGradientExcellent" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor={color.gradientId === 'scoreGradientExcellent' ? color.startColor : '#fcd34d'} />
            <stop offset="100%" stopColor={color.gradientId === 'scoreGradientExcellent' ? color.endColor : '#b45309'} />
          </linearGradient>
          <linearGradient id="scoreGradientFair" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor={color.gradientId === 'scoreGradientFair' ? color.startColor : '#fbbf24'} />
            <stop offset="100%" stopColor={color.gradientId === 'scoreGradientFair' ? color.endColor : '#d97706'} />
          </linearGradient>
          <linearGradient id="scoreGradientPoor" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor={color.gradientId === 'scoreGradientPoor' ? color.startColor : '#f87171'} />
            <stop offset="100%" stopColor={color.gradientId === 'scoreGradientPoor' ? color.endColor : '#dc2626'} />
          </linearGradient>
        </defs>
      </svg>

      {/* Chart - only render when container has valid dimensions */}
      {isReady && (
        <ResponsiveContainer width="100%" height="100%" initialDimension={{ width: 1, height: 1 }}>
          <RadialBarChart
            innerRadius="70%"
            outerRadius="100%"
            data={data}
            startAngle={90}
            endAngle={-270}
            barSize={12}
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

      {/* Center content */}
      <div className="absolute inset-0 flex flex-col items-center justify-center">
        <span className={cn('text-4xl font-bold', color.text)}>
          {Math.round(clampedScore)}
        </span>
        <span className="text-xs text-text-muted uppercase tracking-wider mt-1">
          Score
        </span>
      </div>
    </div>
  );

  // If we have calculation/interpretation, wrap in tooltip
  if (calculation && interpretation) {
    return (
      <div className="flex flex-col items-center">
        <CalculationTooltip
          calculation={calculation}
          interpretation={interpretation}
          title={label}
          value={`${Math.round(clampedScore)}/100`}
          position="top"
        >
          {chartContent}
        </CalculationTooltip>
        {!compact && (
          <div className={cn('mt-2 px-4 py-1.5 rounded-full', color.bg)}>
            <span className={cn('text-sm font-medium', color.text)}>{color.label}</span>
          </div>
        )}
        {!compact && (
          <p className="text-sm text-text-secondary mt-2 text-center">{label}</p>
        )}
      </div>
    );
  }

  return (
    <div className="flex flex-col items-center">
      {chartContent}
      {!compact && (
        <>
          <div className={cn('mt-2 px-4 py-1.5 rounded-full', color.bg)}>
            <span className={cn('text-sm font-medium', color.text)}>{color.label}</span>
          </div>
          <p className="text-sm text-text-secondary mt-2 text-center">{label}</p>
        </>
      )}
    </div>
  );
}

/**
 * Compact score indicator for tables/lists
 */
export function ScoreIndicator({ score, size = 40 }: { score: number; size?: number }) {
  const { isDark } = useThemeContext();
  const color = getScoreColor(score, isDark);

  return (
    <div
      className={cn('rounded-full flex items-center justify-center', color.bg)}
      style={{ width: size, height: size }}
    >
      <span className={cn('text-sm font-bold', color.text)}>{Math.round(score)}</span>
    </div>
  );
}

/**
 * Score breakdown bar component
 */
interface ScoreBreakdownBarProps {
  label: string;
  score: number;
  maxScore: number;
  detail?: string;
  calculation?: string;
  interpretation?: string;
}

export function ScoreBreakdownBar({
  label,
  score,
  maxScore,
  detail,
  calculation,
  interpretation,
}: ScoreBreakdownBarProps) {
  const { isDark } = useThemeContext();
  const percentage = maxScore > 0 ? (score / maxScore) * 100 : 0;
  const color = getScoreColor(percentage, isDark);

  const content = (
    <div className="space-y-1.5">
      <div className="flex justify-between text-sm">
        <span className="text-text-primary font-medium">{label}</span>
        <span className={color.text}>
          {score.toFixed(0)}/{maxScore}
          {detail && (
            <span className="text-text-muted ml-2 text-xs">({detail})</span>
          )}
        </span>
      </div>
      <div className="h-2.5 bg-gray-200 dark:bg-gray-700 rounded-full overflow-hidden">
        <div
          className="h-full rounded-full transition-all duration-700 ease-out"
          style={{
            width: `${percentage}%`,
            background: `linear-gradient(90deg, ${color.fill}dd, ${color.fill})`,
          }}
        />
      </div>
    </div>
  );

  if (calculation && interpretation) {
    return (
      <CalculationTooltip
        calculation={calculation}
        interpretation={interpretation}
        title={label}
        value={`${score.toFixed(0)}/${maxScore}`}
        position="right"
      >
        {content}
      </CalculationTooltip>
    );
  }

  return content;
}

export default ScoreChart;
