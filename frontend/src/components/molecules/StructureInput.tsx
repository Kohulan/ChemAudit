import { useHotkeys } from 'react-hotkeys-hook';
import { cn } from '../../lib/utils';

interface StructureInputProps {
  value: string;
  onChange: (value: string) => void;
  onSubmit?: () => void;
  disabled?: boolean;
  placeholder?: string;
}

export function StructureInput({
  value,
  onChange,
  onSubmit,
  disabled = false,
  placeholder = 'Enter SMILES, InChI, or paste MOL block...'
}: StructureInputProps) {
  // Detect platform for keyboard hint
  const isMac = typeof navigator !== 'undefined' && /Mac|iPhone|iPad|iPod/.test(navigator.platform);
  const modifierKey = isMac ? 'Cmd' : 'Ctrl';

  // Use react-hotkeys-hook for Ctrl+Enter / Cmd+Enter
  useHotkeys(
    'ctrl+enter, meta+enter',
    (e) => {
      e.preventDefault();
      if (onSubmit && !disabled) {
        onSubmit();
      }
    },
    {
      enableOnFormTags: ['TEXTAREA'],
      enabled: !disabled,
    },
    [onSubmit, disabled]
  );

  const handleChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
    onChange(e.target.value);
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey && onSubmit) {
      e.preventDefault();
      onSubmit();
    }
  };

  return (
    <div className="space-y-2">
      <textarea
        value={value}
        onChange={handleChange}
        onKeyDown={handleKeyDown}
        disabled={disabled}
        placeholder={placeholder}
        className={cn(
          'w-full h-24 px-4 py-3 rounded-xl font-mono text-sm resize-none',
          'bg-[var(--color-surface-elevated)] dark:bg-[var(--color-surface-sunken)]',
          'border border-[var(--color-border-strong)]',
          'text-[var(--color-text-primary)]',
          'placeholder:text-[var(--color-text-muted)]',
          'focus:outline-none focus:border-[var(--color-primary)] focus:ring-2 focus:ring-[var(--glow-primary)]',
          'disabled:opacity-50 disabled:cursor-not-allowed',
          'transition-all duration-200'
        )}
        spellCheck={false}
      />
      <div className="flex items-center justify-between text-xs text-[var(--color-text-muted)]">
        <span>Press {modifierKey}+Enter to validate, Shift+Enter for new line</span>
        <span>{value.length} chars</span>
      </div>
    </div>
  );
}
