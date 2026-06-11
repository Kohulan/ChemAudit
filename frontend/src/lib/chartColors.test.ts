import { describe, it, expect } from 'vitest';
import { scoreFill, scoreBucket, chartTrackFill } from './chartColors';

describe('chartColors', () => {
  it('maps scores to the warm bucket hierarchy', () => {
    expect(scoreBucket(85)).toBe('excellent');
    expect(scoreBucket(80)).toBe('excellent');
    expect(scoreBucket(50)).toBe('good');
    expect(scoreBucket(20)).toBe('moderate');
    expect(scoreBucket(5)).toBe('poor');
  });

  it('returns theme-specific fills (gold for excellent, never green)', () => {
    expect(scoreFill('excellent', false)).toBe('#b45309');
    expect(scoreFill('excellent', true)).toBe('#fcd34d');
    expect(scoreFill('poor', false)).toBe('#dc2626');
  });

  it('uses warm stone for the gauge track', () => {
    expect(chartTrackFill(false)).toBe('#e7e5e4');
    expect(chartTrackFill(true)).toBe('#292524');
  });
});
