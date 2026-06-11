/**
 * Theme-aware chart color scales.
 *
 * Single source of truth for score-quality colors used in Recharts/SVG
 * fills, where Tailwind classes and CSS variables don't reach gradient
 * stops reliably. Values mirror the warm hierarchy in tailwind.config.js
 * (score-excellent/good/fair/poor): gold = excellent ... red = poor.
 * Success is NEVER green in this app (see PRODUCT.md, "Warm precision").
 */
export type ScoreBucket = 'excellent' | 'good' | 'moderate' | 'poor';

const SCORE_FILL: Record<ScoreBucket, { light: string; dark: string }> = {
  excellent: { light: '#b45309', dark: '#fcd34d' }, // amber-700 / amber-300
  good:      { light: '#d97706', dark: '#fbbf24' }, // amber-600 / amber-400
  moderate:  { light: '#ea580c', dark: '#fb923c' }, // orange-600 / orange-400
  poor:      { light: '#dc2626', dark: '#f87171' }, // red-600 / red-400
};

export function scoreFill(bucket: ScoreBucket, isDark: boolean): string {
  return SCORE_FILL[bucket][isDark ? 'dark' : 'light'];
}

/** Standard score → bucket mapping used by histograms and gauges. */
export function scoreBucket(score: number): ScoreBucket {
  if (score >= 80) return 'excellent';
  if (score >= 50) return 'good';
  if (score >= 20) return 'moderate';
  return 'poor';
}

/** Neutral track behind gauges/radials: warm stone, not cool gray. */
export function chartTrackFill(isDark: boolean): string {
  return isDark ? '#292524' : '#e7e5e4'; // stone-800 / stone-200
}

/**
 * Translucent reference-zone tint for radar/area charts.
 * Derived from the "good" fill at 12% alpha so the amber hue lives in
 * one place (used e.g. for the 5th-95th percentile band in ComplexityRadar).
 */
export function referenceZoneFill(isDark: boolean): string {
  const hex = SCORE_FILL.good[isDark ? 'dark' : 'light'];
  const r = parseInt(hex.slice(1, 3), 16);
  const g = parseInt(hex.slice(3, 5), 16);
  const b = parseInt(hex.slice(5, 7), 16);
  return `rgba(${r}, ${g}, ${b}, 0.12)`;
}
