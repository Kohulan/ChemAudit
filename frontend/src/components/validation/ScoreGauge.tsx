import { RadialBarChart, RadialBar, ResponsiveContainer, PolarAngleAxis } from 'recharts';
import { CalculationTooltip } from '../ui/Tooltip';

interface ScoreGaugeProps {
  score: number;
  size?: number;
  className?: string;
  showCalculation?: boolean;
}

/**
 * Circular score gauge using recharts RadialBarChart.
 * Displays validation score with color-coded feedback and optional calculation tooltip.
 */
export function ScoreGauge({ score, size = 140, className = '', showCalculation = true }: ScoreGaugeProps) {
  // Clamp score to 0-100
  const clampedScore = Math.max(0, Math.min(100, score));

  // Determine color based on score
  const getColor = (score: number) => {
    if (score >= 80) return {
      fill: '#059669',
      text: 'text-score-excellent',
      bg: 'bg-score-excellent/10',
      label: 'Excellent',
      gradientId: 'gaugeExcellent',
    };
    if (score >= 50) return {
      fill: '#d97706',
      text: 'text-score-fair',
      bg: 'bg-score-fair/10',
      label: 'Fair',
      gradientId: 'gaugeFair',
    };
    return {
      fill: '#dc2626',
      text: 'text-score-poor',
      bg: 'bg-score-poor/10',
      label: 'Poor',
      gradientId: 'gaugePoor',
    };
  };

  const color = getColor(clampedScore);

  const data = [
    {
      name: 'Validation Score',
      value: clampedScore,
      fill: `url(#${color.gradientId})`,
    },
  ];

  const calculation = `Score = 100 - (CRITICAL * 50 + ERROR * 20 + WARNING * 5)
Clamped to range 0-100`;

  const interpretation = clampedScore >= 80
    ? 'This molecule passes all critical validation checks and has minimal issues. It is suitable for most applications.'
    : clampedScore >= 50
    ? 'This molecule has some validation issues that may need attention. Review the warnings and errors below.'
    : 'This molecule has significant validation problems. Critical issues must be resolved before use.';

  const chartContent = (
    <div className={`relative ${className}`} style={{ width: size, height: size }}>
      {/* SVG Gradient Definitions */}
      <svg width="0" height="0" className="absolute">
        <defs>
          <linearGradient id="gaugeExcellent" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#10b981" />
            <stop offset="100%" stopColor="#059669" />
          </linearGradient>
          <linearGradient id="gaugeFair" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#fbbf24" />
            <stop offset="100%" stopColor="#d97706" />
          </linearGradient>
          <linearGradient id="gaugePoor" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#f87171" />
            <stop offset="100%" stopColor="#dc2626" />
          </linearGradient>
        </defs>
      </svg>

      {/* Chart */}
      <ResponsiveContainer width="100%" height="100%">
        <RadialBarChart
          innerRadius="65%"
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
            background={{ fill: '#e5e7eb' }}
            dataKey="value"
            cornerRadius={10}
            animationDuration={1000}
            animationEasing="ease-out"
          />
        </RadialBarChart>
      </ResponsiveContainer>

      {/* Score text in center */}
      <div className="absolute inset-0 flex flex-col items-center justify-center">
        <span className={`text-3xl font-bold ${color.text}`}>
          {Math.round(clampedScore)}
        </span>
        <span className="text-xs text-chem-dark/50 uppercase tracking-wide">Score</span>
      </div>
    </div>
  );

  if (showCalculation) {
    return (
      <div className="flex flex-col items-center">
        <CalculationTooltip
          calculation={calculation}
          interpretation={interpretation}
          title="Validation Score"
          value={`${Math.round(clampedScore)}/100`}
          position="top"
        >
          {chartContent}
        </CalculationTooltip>
        <div className={`mt-3 px-4 py-1.5 rounded-full ${color.bg}`}>
          <span className={`text-sm font-medium ${color.text}`}>{color.label}</span>
        </div>
      </div>
    );
  }

  return (
    <div className="flex flex-col items-center">
      {chartContent}
      <div className={`mt-3 px-4 py-1.5 rounded-full ${color.bg}`}>
        <span className={`text-sm font-medium ${color.text}`}>{color.label}</span>
      </div>
    </div>
  );
}
