import { RadialBarChart, RadialBar, ResponsiveContainer, PolarAngleAxis } from 'recharts';
import { CalculationTooltip } from '../ui/Tooltip';

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
}

/**
 * Get color configuration based on score value
 */
function getScoreColor(score: number): {
  fill: string;
  text: string;
  bg: string;
  label: string;
  gradientId: string;
} {
  if (score >= 80) {
    return {
      fill: '#059669', // emerald-600
      text: 'text-score-excellent',
      bg: 'bg-score-excellent/10',
      label: 'Excellent',
      gradientId: 'scoreGradientExcellent',
    };
  }
  if (score >= 50) {
    return {
      fill: '#d97706', // amber-600
      text: 'text-score-fair',
      bg: 'bg-score-fair/10',
      label: 'Fair',
      gradientId: 'scoreGradientFair',
    };
  }
  return {
    fill: '#dc2626', // red-600
    text: 'text-score-poor',
    bg: 'bg-score-poor/10',
    label: 'Poor',
    gradientId: 'scoreGradientPoor',
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
}: ScoreChartProps) {
  const clampedScore = Math.max(0, Math.min(100, score));
  const color = getScoreColor(clampedScore);

  const data = [
    {
      name: label,
      value: clampedScore,
      fill: `url(#${color.gradientId})`,
    },
  ];

  const chartContent = (
    <div className={`relative ${compact ? '' : 'p-2'}`} style={{ width: size, height: size }}>
      {/* SVG Gradient Definitions */}
      <svg width="0" height="0" className="absolute">
        <defs>
          <linearGradient id="scoreGradientExcellent" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#10b981" />
            <stop offset="100%" stopColor="#059669" />
          </linearGradient>
          <linearGradient id="scoreGradientFair" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#fbbf24" />
            <stop offset="100%" stopColor="#d97706" />
          </linearGradient>
          <linearGradient id="scoreGradientPoor" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#f87171" />
            <stop offset="100%" stopColor="#dc2626" />
          </linearGradient>
        </defs>
      </svg>

      {/* Chart */}
      <ResponsiveContainer width="100%" height="100%">
        <RadialBarChart
          innerRadius="70%"
          outerRadius="100%"
          data={data}
          startAngle={90}
          endAngle={-270}
          barSize={compact ? 8 : 12}
        >
          <PolarAngleAxis
            type="number"
            domain={[0, 100]}
            angleAxisId={0}
            tick={false}
          />
          <RadialBar
            background={{ fill: '#e5e7eb' }}
            dataKey="value"
            cornerRadius={10}
            animationDuration={1000}
            animationEasing="ease-out"
          />
        </RadialBarChart>
      </ResponsiveContainer>

      {/* Center content */}
      <div className="absolute inset-0 flex flex-col items-center justify-center">
        <span className={`${compact ? 'text-2xl' : 'text-4xl'} font-bold ${color.text}`}>
          {Math.round(clampedScore)}
        </span>
        {!compact && (
          <span className="text-xs text-chem-dark/50 uppercase tracking-wider mt-1">
            Score
          </span>
        )}
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
          <div className={`mt-2 px-4 py-1.5 rounded-full ${color.bg}`}>
            <span className={`text-sm font-medium ${color.text}`}>{color.label}</span>
          </div>
        )}
        {!compact && (
          <p className="text-sm text-chem-dark/60 mt-2 text-center">{label}</p>
        )}
      </div>
    );
  }

  return (
    <div className="flex flex-col items-center">
      {chartContent}
      {!compact && (
        <>
          <div className={`mt-2 px-4 py-1.5 rounded-full ${color.bg}`}>
            <span className={`text-sm font-medium ${color.text}`}>{color.label}</span>
          </div>
          <p className="text-sm text-chem-dark/60 mt-2 text-center">{label}</p>
        </>
      )}
    </div>
  );
}

/**
 * Compact score indicator for tables/lists
 */
export function ScoreIndicator({ score, size = 40 }: { score: number; size?: number }) {
  const color = getScoreColor(score);

  return (
    <div
      className={`rounded-full flex items-center justify-center ${color.bg}`}
      style={{ width: size, height: size }}
    >
      <span className={`text-sm font-bold ${color.text}`}>{Math.round(score)}</span>
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
  const percentage = maxScore > 0 ? (score / maxScore) * 100 : 0;
  const color = getScoreColor(percentage);

  const content = (
    <div className="space-y-1.5">
      <div className="flex justify-between text-sm">
        <span className="text-chem-dark font-medium">{label}</span>
        <span className={color.text}>
          {score.toFixed(0)}/{maxScore}
          {detail && (
            <span className="text-chem-dark/40 ml-2 text-xs">({detail})</span>
          )}
        </span>
      </div>
      <div className="h-2.5 bg-chem-dark/10 rounded-full overflow-hidden">
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
