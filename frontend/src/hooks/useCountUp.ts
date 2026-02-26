import { useState, useEffect } from 'react';

/**
 * Animates a numeric value from 0 to the target using cubic ease-out.
 *
 * @param target  - Final value to animate towards
 * @param decimals - Decimal places to preserve (0 for integers)
 * @param duration - Animation duration in ms
 * @param delay    - Delay before animation starts in ms
 */
export function useCountUp(
  target: number,
  decimals = 2,
  duration = 800,
  delay = 200,
): number {
  const [value, setValue] = useState(0);

  useEffect(() => {
    const timeout = setTimeout(() => {
      const start = performance.now();
      const tick = (now: number) => {
        const p = Math.min((now - start) / duration, 1);
        const eased = 1 - Math.pow(1 - p, 3);
        setValue(target * eased);
        if (p < 1) requestAnimationFrame(tick);
      };
      requestAnimationFrame(tick);
    }, delay);
    return () => clearTimeout(timeout);
  }, [target, duration, delay]);

  return decimals > 0 ? +value.toFixed(decimals) : Math.round(value);
}
