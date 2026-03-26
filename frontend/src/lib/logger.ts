/**
 * Logger utility.
 * error/warn always log (including production).
 * info/debug are gated to development mode only.
 */
const DEBUG_MODE = import.meta.env.DEV;

export const logger = {
  error(...args: unknown[]): void {
    console.error(...args);
  },
  warn(...args: unknown[]): void {
    console.warn(...args);
  },
  info(...args: unknown[]): void {
    if (DEBUG_MODE) console.info(...args);
  },
  debug(...args: unknown[]): void {
    if (DEBUG_MODE) console.debug(...args);
  },
};
