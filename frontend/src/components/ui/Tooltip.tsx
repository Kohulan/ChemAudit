import { useState, useRef, useEffect, ReactNode } from 'react';

interface TooltipProps {
  /** The content that triggers the tooltip */
  children: ReactNode;
  /** The tooltip content */
  content: ReactNode;
  /** Tooltip title (optional) */
  title?: string;
  /** Position of tooltip relative to trigger */
  position?: 'top' | 'bottom' | 'left' | 'right';
  /** Delay before showing tooltip (ms) */
  delay?: number;
  /** Maximum width of tooltip */
  maxWidth?: number;
  /** Whether tooltip is disabled */
  disabled?: boolean;
  /** Additional className for tooltip content */
  className?: string;
}

interface CalculationTooltipProps {
  /** The content that triggers the tooltip */
  children: ReactNode;
  /** Calculation formula or method */
  calculation: string;
  /** Interpretation/meaning of the value */
  interpretation: string;
  /** Optional title */
  title?: string;
  /** Value being explained */
  value?: string | number;
  /** Position of tooltip */
  position?: 'top' | 'bottom' | 'left' | 'right';
}

/**
 * Enhanced tooltip component with positioning and animations
 */
export function Tooltip({
  children,
  content,
  title,
  position = 'top',
  delay = 200,
  maxWidth = 280,
  disabled = false,
  className = '',
}: TooltipProps) {
  const [isVisible, setIsVisible] = useState(false);
  const [coords, setCoords] = useState({ x: 0, y: 0 });
  const triggerRef = useRef<HTMLDivElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);
  const timeoutRef = useRef<ReturnType<typeof setTimeout>>();

  const showTooltip = () => {
    if (disabled) return;
    timeoutRef.current = setTimeout(() => {
      setIsVisible(true);
    }, delay);
  };

  const hideTooltip = () => {
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
    }
    setIsVisible(false);
  };

  useEffect(() => {
    if (isVisible && triggerRef.current && tooltipRef.current) {
      const triggerRect = triggerRef.current.getBoundingClientRect();
      const tooltipRect = tooltipRef.current.getBoundingClientRect();

      let x = 0;
      let y = 0;

      switch (position) {
        case 'top':
          x = triggerRect.left + triggerRect.width / 2 - tooltipRect.width / 2;
          y = triggerRect.top - tooltipRect.height - 8;
          break;
        case 'bottom':
          x = triggerRect.left + triggerRect.width / 2 - tooltipRect.width / 2;
          y = triggerRect.bottom + 8;
          break;
        case 'left':
          x = triggerRect.left - tooltipRect.width - 8;
          y = triggerRect.top + triggerRect.height / 2 - tooltipRect.height / 2;
          break;
        case 'right':
          x = triggerRect.right + 8;
          y = triggerRect.top + triggerRect.height / 2 - tooltipRect.height / 2;
          break;
      }

      // Keep tooltip within viewport
      const padding = 10;
      x = Math.max(padding, Math.min(x, window.innerWidth - tooltipRect.width - padding));
      y = Math.max(padding, Math.min(y, window.innerHeight - tooltipRect.height - padding));

      setCoords({ x, y });
    }
  }, [isVisible, position]);

  useEffect(() => {
    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, []);

  return (
    <div className="inline-block">
      <div
        ref={triggerRef}
        onMouseEnter={showTooltip}
        onMouseLeave={hideTooltip}
        onFocus={showTooltip}
        onBlur={hideTooltip}
      >
        {children}
      </div>

      {isVisible && (
        <div
          ref={tooltipRef}
          className={`fixed z-[100] px-4 py-3 bg-chem-dark text-white rounded-xl shadow-xl animate-fade-in ${className}`}
          style={{
            left: coords.x,
            top: coords.y,
            maxWidth,
          }}
          role="tooltip"
        >
          {/* Arrow indicator */}
          <div
            className={`absolute w-2 h-2 bg-chem-dark transform rotate-45 ${
              position === 'top'
                ? 'bottom-[-4px] left-1/2 -translate-x-1/2'
                : position === 'bottom'
                ? 'top-[-4px] left-1/2 -translate-x-1/2'
                : position === 'left'
                ? 'right-[-4px] top-1/2 -translate-y-1/2'
                : 'left-[-4px] top-1/2 -translate-y-1/2'
            }`}
          />

          {title && (
            <div className="font-semibold text-sm mb-2 text-white/90 border-b border-white/10 pb-2">
              {title}
            </div>
          )}
          <div className="text-sm text-white/80">{content}</div>
        </div>
      )}
    </div>
  );
}

/**
 * Specialized tooltip for score calculations with formula and interpretation
 */
export function CalculationTooltip({
  children,
  calculation,
  interpretation,
  title,
  value,
  position = 'top',
}: CalculationTooltipProps) {
  const content = (
    <div className="space-y-3">
      {/* Value display */}
      {value !== undefined && (
        <div className="flex items-center gap-2 pb-2 border-b border-white/10">
          <span className="text-white/60 text-xs uppercase tracking-wide">Value:</span>
          <span className="font-mono font-semibold text-chem-primary-300">{value}</span>
        </div>
      )}

      {/* Calculation method */}
      <div>
        <div className="flex items-center gap-2 mb-1">
          <svg className="w-4 h-4 text-chem-secondary-300" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <path d="M4 4h16v16H4z" />
            <path d="M8 8h2v2H8zM14 8h2v2h-2zM8 14h8v2H8z" />
          </svg>
          <span className="text-xs uppercase tracking-wide text-white/60">Calculation</span>
        </div>
        <p className="text-sm text-white/90 font-mono bg-white/5 rounded-lg px-3 py-2">
          {calculation}
        </p>
      </div>

      {/* Interpretation */}
      <div>
        <div className="flex items-center gap-2 mb-1">
          <svg className="w-4 h-4 text-chem-accent-300" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <circle cx="12" cy="12" r="10" />
            <path d="M12 16v-4M12 8h.01" />
          </svg>
          <span className="text-xs uppercase tracking-wide text-white/60">What This Means</span>
        </div>
        <p className="text-sm text-white/90">{interpretation}</p>
      </div>
    </div>
  );

  return (
    <Tooltip content={content} title={title} position={position} maxWidth={320}>
      <span className="cursor-help inline-flex items-center gap-1 border-b border-dashed border-chem-dark/30 hover:border-chem-primary transition-colors">
        {children}
        <svg className="w-3.5 h-3.5 text-chem-dark/40" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="10" />
          <path d="M9.09 9a3 3 0 015.83 1c0 2-3 3-3 3M12 17h.01" />
        </svg>
      </span>
    </Tooltip>
  );
}

/**
 * Info icon tooltip for quick explanations
 */
export function InfoTooltip({
  content,
  title,
  position = 'top',
}: {
  content: ReactNode;
  title?: string;
  position?: 'top' | 'bottom' | 'left' | 'right';
}) {
  return (
    <Tooltip content={content} title={title} position={position}>
      <button
        type="button"
        className="inline-flex items-center justify-center w-5 h-5 rounded-full bg-chem-primary/10 text-chem-primary hover:bg-chem-primary/20 transition-colors focus:outline-none focus:ring-2 focus:ring-chem-primary/30"
      >
        <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5">
          <circle cx="12" cy="12" r="10" />
          <path d="M12 16v-4M12 8h.01" />
        </svg>
      </button>
    </Tooltip>
  );
}

export default Tooltip;
