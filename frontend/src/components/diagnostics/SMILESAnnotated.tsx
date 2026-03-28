import { cn } from '../../lib/utils';
import { Badge } from '../ui/Badge';

interface SMILESAnnotatedProps {
  /** The full SMILES string to display. */
  smiles: string;
  /** Zero-indexed character position of the parse error, if any. */
  errorPosition?: number | null;
  /** Character indices of atoms with RDKit warnings (valid SMILES). */
  warningPositions?: number[];
}

/**
 * Renders a SMILES string with character-level position highlighting.
 *
 * - Invalid SMILES: highlights the error character with a red background and
 *   renders a compiler-style caret (^) on the line below pointing at the error.
 * - Valid SMILES with warnings: underlines flagged characters.
 * - Valid SMILES with no issues: renders the string plainly with a "Valid" badge.
 *
 * Uses JetBrains Mono (font-mono) as specified in the Phase 9 UI-SPEC.
 */
export function SMILESAnnotated({
  smiles,
  errorPosition,
  warningPositions = [],
}: SMILESAnnotatedProps) {
  const hasError = errorPosition !== null && errorPosition !== undefined;
  const hasWarnings = warningPositions.length > 0;
  const warningSet = new Set(warningPositions);

  // ── Valid SMILES — no errors ──
  if (!hasError) {
    return (
      <div className="space-y-1">
        <div className="flex items-start gap-3">
          <pre className="font-mono text-sm text-[var(--color-text-primary)] whitespace-pre-wrap break-all flex-1 m-0">
            {smiles.split('').map((char, i) => {
              if (hasWarnings && warningSet.has(i)) {
                return (
                  <span
                    key={i}
                    className="border-b border-status-warning"
                    title={`Warning at position ${i}`}
                  >
                    {char}
                  </span>
                );
              }
              return char;
            })}
          </pre>
          {!hasWarnings && (
            <Badge variant="success" size="md" className="shrink-0 mt-0.5">
              Valid
            </Badge>
          )}
        </div>
      </div>
    );
  }

  // ── Invalid SMILES — error position highlighting + caret ──
  const prefix = smiles.slice(0, errorPosition);
  const errorChar = smiles[errorPosition] ?? ' ';
  const suffix = smiles.slice(errorPosition + 1);

  return (
    <div className="space-y-1">
      {/* SMILES string with highlighted error character */}
      <pre className="font-mono text-sm text-[var(--color-text-primary)] whitespace-pre-wrap break-all m-0">
        <span>{prefix}</span>
        <span
          className="bg-[rgba(239,68,68,0.2)] rounded-sm"
          title={`Parse error at position ${errorPosition}`}
        >
          {errorChar}
        </span>
        <span>{suffix}</span>
      </pre>

      {/* Caret line — aligns ^ under the error character */}
      <pre
        className={cn(
          'font-mono text-sm text-status-error m-0 select-none',
          'leading-none',
        )}
        aria-hidden="true"
      >
        {/* Invisible prefix to align the caret */}
        <span className="invisible">{prefix}</span>
        <span>^</span>
      </pre>
    </div>
  );
}
