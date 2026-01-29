import { useState, useRef, useEffect, ReactNode } from 'react';
import { cn } from '../../lib/utils';

interface TooltipProps {
  children: ReactNode;
  content: ReactNode;
  title?: string;
  position?: 'top' | 'bottom' | 'left' | 'right';
  delay?: number;
  maxWidth?: number;
  disabled?: boolean;
  className?: string;
}

const positionClasses = {
  top: 'bottom-full left-1/2 -translate-x-1/2 mb-2',
  bottom: 'top-full left-1/2 -translate-x-1/2 mt-2',
  left: 'right-full top-1/2 -translate-y-1/2 mr-2',
  right: 'left-full top-1/2 -translate-y-1/2 ml-2',
};

const arrowClasses = {
  top: 'top-full left-1/2 -translate-x-1/2 border-t-zinc-800 border-x-transparent border-b-transparent',
  bottom: 'bottom-full left-1/2 -translate-x-1/2 border-b-zinc-800 border-x-transparent border-t-transparent',
  left: 'left-full top-1/2 -translate-y-1/2 border-l-zinc-800 border-y-transparent border-r-transparent',
  right: 'right-full top-1/2 -translate-y-1/2 border-r-zinc-800 border-y-transparent border-l-transparent',
};

/**
 * Simple tooltip with relative positioning
 */
export function Tooltip({
  children,
  content,
  title,
  position = 'top',
  delay = 150,
  maxWidth = 280,
  disabled = false,
  className = '',
}: TooltipProps) {
  const [isVisible, setIsVisible] = useState(false);
  const timeoutRef = useRef<ReturnType<typeof setTimeout>>();

  // Cleanup timeout on unmount to prevent memory leaks
  useEffect(() => {
    return () => {
      if (timeoutRef.current) clearTimeout(timeoutRef.current);
    };
  }, []);

  const show = () => {
    if (disabled) return;
    timeoutRef.current = setTimeout(() => setIsVisible(true), delay);
  };

  const hide = () => {
    if (timeoutRef.current) clearTimeout(timeoutRef.current);
    setIsVisible(false);
  };

  return (
    <div
      className="relative inline-flex"
      onMouseEnter={show}
      onMouseLeave={hide}
      onFocus={show}
      onBlur={hide}
    >
      {children}

      {isVisible && (
        <div
          className={cn(
            'absolute z-[9999] px-3 py-2 rounded-lg shadow-lg',
            'bg-zinc-800 text-white text-sm',
            'whitespace-normal',
            positionClasses[position],
            className
          )}
          style={{ maxWidth, minWidth: 120 }}
          role="tooltip"
        >
          <div
            className={cn(
              'absolute w-0 h-0 border-[6px]',
              arrowClasses[position]
            )}
          />

          {title && (
            <div className="font-semibold text-white mb-1.5 pb-1.5 border-b border-white/20">
              {title}
            </div>
          )}
          <div className="text-zinc-200">{content}</div>
        </div>
      )}
    </div>
  );
}

/**
 * Info icon with tooltip
 */
export function InfoTooltip({
  content,
  title,
  position = 'right',
}: {
  content: ReactNode;
  title?: string;
  position?: 'top' | 'bottom' | 'left' | 'right';
}) {
  return (
    <Tooltip content={content} title={title} position={position} maxWidth={260}>
      <button
        type="button"
        className={cn(
          'inline-flex items-center justify-center w-4 h-4 rounded-full',
          'bg-[var(--color-primary)]/15 text-[var(--color-primary)]',
          'hover:bg-[var(--color-primary)]/25',
          'transition-colors cursor-help'
        )}
      >
        <svg className="w-2.5 h-2.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="3">
          <circle cx="12" cy="12" r="10" />
          <path d="M12 16v-4M12 8h.01" />
        </svg>
      </button>
    </Tooltip>
  );
}

/**
 * Calculation tooltip for explaining scores
 */
export function CalculationTooltip({
  children,
  calculation,
  interpretation,
  title,
  value,
  position = 'top',
}: {
  children: ReactNode;
  calculation: string;
  interpretation: string;
  title?: string;
  value?: string | number;
  position?: 'top' | 'bottom' | 'left' | 'right';
}) {
  const content = (
    <div className="space-y-2">
      {value !== undefined && (
        <div className="flex items-center gap-2 pb-2 border-b border-white/20">
          <span className="text-zinc-400 text-xs">Value:</span>
          <span className="font-mono font-semibold text-red-400">{value}</span>
        </div>
      )}
      <div>
        <span className="text-zinc-400 text-xs block mb-1">Calculation:</span>
        <code className="text-xs bg-black/30 px-2 py-1 rounded block">{calculation}</code>
      </div>
      <div>
        <span className="text-zinc-400 text-xs block mb-1">Meaning:</span>
        <span className="text-zinc-200">{interpretation}</span>
      </div>
    </div>
  );

  return (
    <Tooltip content={content} title={title} position={position} maxWidth={300}>
      <span className="cursor-help border-b border-dashed border-[var(--color-text-muted)] hover:border-[var(--color-primary)] transition-colors inline-flex items-center gap-1">
        {children}
        <svg className="w-3 h-3 text-[var(--color-text-muted)]" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="10" />
          <path d="M9.09 9a3 3 0 015.83 1c0 2-3 3-3 3M12 17h.01" />
        </svg>
      </span>
    </Tooltip>
  );
}

export default Tooltip;
