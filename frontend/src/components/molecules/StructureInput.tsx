import { useState, useCallback } from 'react';

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
  const [format, setFormat] = useState<string>('auto');

  const detectFormat = useCallback((input: string) => {
    const trimmed = input.trim();
    if (trimmed.startsWith('InChI=')) return 'inchi';
    if (trimmed.includes('M  END') || trimmed.includes('V2000')) return 'mol';
    return 'smiles';
  }, []);

  const handleChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
    const newValue = e.target.value;
    onChange(newValue);
    setFormat(detectFormat(newValue));
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey && onSubmit) {
      e.preventDefault();
      onSubmit();
    }
  };

  return (
    <div className="space-y-2">
      <div className="relative">
        <textarea
          value={value}
          onChange={handleChange}
          onKeyDown={handleKeyDown}
          disabled={disabled}
          placeholder={placeholder}
          className="w-full h-24 px-4 py-3 border border-gray-300 rounded-lg font-mono text-sm
                     focus:ring-2 focus:ring-blue-500 focus:border-transparent
                     disabled:bg-gray-100 disabled:cursor-not-allowed resize-none"
          spellCheck={false}
        />
        <div className="absolute bottom-2 right-2">
          <span className="px-2 py-1 text-xs bg-gray-100 text-gray-500 rounded">
            {format.toUpperCase()}
          </span>
        </div>
      </div>
      <div className="flex items-center justify-between text-xs text-gray-500">
        <span>Press Enter to validate, Shift+Enter for new line</span>
        <span>{value.length} characters</span>
      </div>
    </div>
  );
}
