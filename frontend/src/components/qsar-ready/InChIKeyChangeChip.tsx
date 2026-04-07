import { cn } from '../../lib/utils';

// =============================================================================
// Props
// =============================================================================

interface InChIKeyChangeChipProps {
  /** The original InChIKey before curation. */
  originalInchikey: string | null;
  /** The standardized InChIKey after curation. */
  standardizedInchikey: string | null;
  /** Whether the InChIKey changed during curation (authoritative flag from backend). */
  inchikeyChanged: boolean;
}

// =============================================================================
// Helpers
// =============================================================================

/** Truncate an InChIKey to first 14 chars + "…" for display. */
function truncateInchiKey(key: string | null): string {
  if (!key) return '—';
  if (key.length <= 14) return key;
  return key.slice(0, 14) + '…';
}

// =============================================================================
// Component
// =============================================================================

/**
 * Inline chip showing InChIKey change status between original and curated molecule.
 *
 * - Changed: amber accent badge "InChIKey changed" per UI-SPEC copywriting.
 * - Unchanged: muted badge "InChIKey preserved" per UI-SPEC copywriting.
 *
 * Both InChIKey values are shown in font-mono text-xs truncated with title tooltip.
 * Follows FixChip pattern from diagnostics per UI-SPEC Component Inventory.
 */
export function InChIKeyChangeChip({
  originalInchikey,
  standardizedInchikey,
  inchikeyChanged,
}: InChIKeyChangeChipProps) {
  return (
    <div className="flex flex-col items-center gap-2">
      {/* Status badge */}
      <span
        className={cn(
          'inline-flex items-center gap-1.5 px-3 py-1.5 rounded-full text-xs font-semibold',
          inchikeyChanged
            ? 'bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-400'
            : 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]',
        )}
      >
        {inchikeyChanged ? '⇄ InChIKey changed' : '✓ InChIKey preserved'}
      </span>

      {/* InChIKey values */}
      <div className="flex flex-col items-center gap-1">
        {originalInchikey && (
          <span
            className="font-mono text-xs text-[var(--color-text-muted)]"
            title={originalInchikey}
            aria-label={`InChIKey: ${originalInchikey}`}
          >
            {truncateInchiKey(originalInchikey)}
          </span>
        )}
        {inchikeyChanged && standardizedInchikey && (
          <>
            <span className="text-[var(--color-text-muted)] text-xs">↓</span>
            <span
              className="font-mono text-xs text-[var(--color-text-muted)]"
              title={standardizedInchikey}
              aria-label={`InChIKey: ${standardizedInchikey}`}
            >
              {truncateInchiKey(standardizedInchikey)}
            </span>
          </>
        )}
      </div>
    </div>
  );
}
