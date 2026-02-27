/**
 * Debug-gated logger. Only outputs in development mode.
 */
const DEBUG_MODE = import.meta.env.DEV;

export const logger = {
  error(...args: unknown[]): void {
    if (DEBUG_MODE) console.error(...args);
  },
  warn(...args: unknown[]): void {
    if (DEBUG_MODE) console.warn(...args);
  },
  info(...args: unknown[]): void {
    if (DEBUG_MODE) console.info(...args);
  },
  debug(...args: unknown[]): void {
    if (DEBUG_MODE) console.debug(...args);
  },
};
