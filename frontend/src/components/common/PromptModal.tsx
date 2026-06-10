import { useEffect, useId, useRef, useState } from 'react';

export interface PromptModalProps {
  isOpen: boolean;
  title: string;
  placeholder?: string;
  initialValue?: string;
  confirmLabel?: string;
  cancelLabel?: string;
  /** Called with the trimmed, non-empty value. */
  onConfirm: (value: string) => void;
  onCancel: () => void;
}

/**
 * Accessible single-line text-input dialog — a testable replacement for the
 * native blocking `window.prompt`. Submits on Enter, cancels on Escape, and
 * disables confirm while the trimmed value is empty.
 */
export function PromptModal({
  isOpen,
  title,
  placeholder,
  initialValue = '',
  confirmLabel = 'Save',
  cancelLabel = 'Cancel',
  onConfirm,
  onCancel,
}: PromptModalProps) {
  const titleId = useId();
  const [value, setValue] = useState(initialValue);
  const inputRef = useRef<HTMLInputElement>(null);

  // Reset to the initial value and focus the input each time the dialog opens.
  useEffect(() => {
    if (!isOpen) return;
    setValue(initialValue);
    inputRef.current?.focus();
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') onCancel();
    };
    document.addEventListener('keydown', onKey);
    return () => document.removeEventListener('keydown', onKey);
  }, [isOpen, initialValue, onCancel]);

  if (!isOpen) return null;

  const trimmed = value.trim();

  const submit = () => {
    if (trimmed) onConfirm(trimmed);
  };

  return (
    <div
      className="fixed inset-0 z-[60] flex items-center justify-center bg-[var(--color-text-primary)]/50 p-4"
      onClick={onCancel}
    >
      <div
        role="dialog"
        aria-modal="true"
        aria-labelledby={titleId}
        className="relative w-full max-w-sm rounded-2xl border border-[var(--color-border)] bg-[var(--color-surface-elevated)] p-6 shadow-2xl"
        onClick={(e) => e.stopPropagation()}
      >
        <h2
          id={titleId}
          className="font-display text-lg font-semibold text-[var(--color-text-primary)]"
        >
          {title}
        </h2>
        <input
          ref={inputRef}
          type="text"
          value={value}
          placeholder={placeholder}
          onChange={(e) => setValue(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === 'Enter') {
              e.preventDefault();
              submit();
            }
          }}
          className="mt-4 w-full rounded-xl border border-[var(--color-border)] bg-[var(--color-surface-elevated)] px-3 py-2 text-sm text-[var(--color-text-primary)] focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
        />
        <div className="mt-5 flex justify-end gap-2">
          <button
            type="button"
            onClick={onCancel}
            className="rounded-xl px-3 py-2 text-sm text-[var(--color-text-secondary)] transition-colors hover:bg-[var(--color-surface-sunken)] focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]"
          >
            {cancelLabel}
          </button>
          <button
            type="button"
            onClick={submit}
            disabled={!trimmed}
            className="rounded-xl bg-chem-primary-600 px-3 py-2 text-sm font-medium text-white transition-colors hover:bg-chem-primary-700 focus:outline-none focus:ring-2 focus:ring-chem-primary-600 focus:ring-offset-1 disabled:cursor-not-allowed disabled:opacity-50"
          >
            {confirmLabel}
          </button>
        </div>
      </div>
    </div>
  );
}

export default PromptModal;
