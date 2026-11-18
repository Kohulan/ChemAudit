/**
 * CSV cell encoding for client-side downloads.
 *
 * Molecule names and SMILES come from uploaded files. Text starting with
 * = + - @ (or tab / carriage return, or "{=") runs as a formula when the CSV is
 * opened in a spreadsheet, so it gets a leading single quote. Numbers such as
 * "-1.5" are left alone. Mirrors backend/app/services/export/formula_guard.py.
 */

const FORMULA_START = /^(?:[=+\-@\t\r]|\{=)/;

/** Encode one value as a safe, correctly quoted CSV cell. */
export function csvCell(value: string | number | boolean | null | undefined): string {
  if (value === null || value === undefined) return '';
  let text = String(value);
  if (typeof value === 'string' && FORMULA_START.test(text) && Number.isNaN(Number(text))) {
    text = `'${text}`;
  }
  return /[",\n\r]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}
