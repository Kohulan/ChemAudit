import {
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
  ResponsiveContainer,
} from 'recharts';
import { ClayCard } from '../ui/ClayCard';
import { Badge } from '../ui/Badge';
import type { ComplexityResult } from '../../types/safety';

interface ComplexityRadarProps {
  complexity: ComplexityResult;
}

interface RadarDataPoint {
  property: string;
  molecule: number;
  reference: number;
  isOutlier: boolean;
}

type SVGTextAnchor = 'start' | 'middle' | 'end' | 'inherit';

/**
 * Custom tick renderer for PolarAngleAxis.
 * Outlier axes are rendered in status-error red with bold weight.
 * In-range axes use text-secondary color.
 */
function CustomTick(props: {
  x?: number | string;
  y?: number | string;
  payload?: { value: string };
  textAnchor?: string;
  data?: RadarDataPoint[];
}) {
  const { x = 0, y = 0, payload, textAnchor, data = [] } = props;
  if (!payload) return null;
  const point = data.find((d) => d.property === payload.value);
  const isOutlier = point?.isOutlier ?? false;

  // Narrow textAnchor to the SVG-allowed set
  const validAnchors: SVGTextAnchor[] = ['start', 'middle', 'end', 'inherit'];
  const anchor: SVGTextAnchor | undefined = validAnchors.includes(textAnchor as SVGTextAnchor)
    ? (textAnchor as SVGTextAnchor)
    : undefined;

  return (
    <text
      x={x}
      y={y}
      textAnchor={anchor}
      fill={isOutlier ? '#ef4444' : 'var(--color-text-secondary, #5c5650)'}
      fontWeight={isOutlier ? 600 : 400}
      fontSize={12}
    >
      {payload.value}
    </text>
  );
}

/**
 * Complexity Percentile Radar Chart (D-17).
 *
 * Renders a Recharts RadarChart with dual polygons:
 * - Green polygon (#bbf7d0) for the 5th-95th percentile reference zone (locked D-17)
 * - Blue polygon (#93c5fd) for molecule values (locked D-17)
 *
 * Normalizes molecule values to [0, 1] relative to [p5, p95] with slight
 * overshoot allowed (clamped to [0, 1.5]) for visual clarity.
 * Outlier axis labels are colored red (#ef4444).
 */
export function ComplexityRadar({ complexity }: ComplexityRadarProps) {
  const radarData: RadarDataPoint[] = Object.entries(complexity.properties).map(
    ([name, prop]) => ({
      property: name,
      molecule: Math.min(
        Math.max((prop.value - prop.p5) / (prop.p95 - prop.p5), 0),
        1.5
      ),
      reference: 1, // p95 boundary always 1
      isOutlier: prop.outlier,
    })
  );

  return (
    <ClayCard size="md">
      <div className="space-y-4">
        <div className="flex items-center justify-between gap-3 flex-wrap">
          <h2 className="text-lg font-semibold font-display">Complexity Percentile</h2>
          {complexity.n_outliers > 0 ? (
            <Badge variant="error">{complexity.n_outliers} complexity outlier(s)</Badge>
          ) : (
            <Badge variant="success">Complexity within normal range</Badge>
          )}
        </div>

        <div
          aria-label={`Complexity percentile radar chart showing ${complexity.n_outliers} outlier axes`}
        >
          <ResponsiveContainer width="100%" height={300}>
            <RadarChart data={radarData}>
              <PolarGrid />
              <PolarAngleAxis
                dataKey="property"
                tick={(tickProps) => (
                  <CustomTick {...tickProps} data={radarData} />
                )}
              />
              <PolarRadiusAxis
                domain={[0, 1.5]}
                tick={false}
                axisLine={false}
              />
              {/* Green reference zone polygon: 5th-95th percentile range per D-17 */}
              <Radar
                dataKey="reference"
                stroke="none"
                fill="#bbf7d0"
                fillOpacity={0.3}
                name="Reference Zone"
              />
              {/* Blue molecule polygon per D-17 */}
              <Radar
                dataKey="molecule"
                stroke="#93c5fd"
                fill="#93c5fd"
                fillOpacity={0.5}
                name="Molecule"
              />
            </RadarChart>
          </ResponsiveContainer>
        </div>

        <div className="flex items-center gap-6 text-xs text-text-muted justify-center">
          <span className="flex items-center gap-1.5">
            <span
              className="inline-block w-3 h-3 rounded-sm"
              style={{ background: '#bbf7d0', opacity: 0.7 }}
            />
            Reference zone (p5–p95)
          </span>
          <span className="flex items-center gap-1.5">
            <span
              className="inline-block w-3 h-3 rounded-sm"
              style={{ background: '#93c5fd', opacity: 0.7 }}
            />
            Molecule
          </span>
        </div>
      </div>
    </ClayCard>
  );
}
