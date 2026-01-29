import { forwardRef } from 'react';
import { motion, type HTMLMotionProps } from 'framer-motion';
import { cva, type VariantProps } from 'class-variance-authority';
import { cn } from '../../lib/utils';

const buttonVariants = cva(
  'inline-flex items-center justify-center gap-2 font-medium transition-all duration-200 disabled:opacity-50 disabled:cursor-not-allowed focus-ring',
  {
    variants: {
      variant: {
        default: 'btn btn-default',
        primary: 'btn btn-primary',
        accent: 'btn btn-accent',
        ghost: 'btn btn-ghost',
        outline: [
          'btn bg-transparent',
          'border border-[var(--color-primary)]',
          'text-[var(--color-primary)]',
          'hover:bg-[var(--color-primary)]/10',
        ].join(' '),
        danger: [
          'btn',
          'bg-gradient-to-r from-red-500 to-red-600',
          'text-white border-none',
          'shadow-[0_2px_8px_rgba(239,68,68,0.25)]',
          'hover:shadow-[0_4px_16px_rgba(239,68,68,0.35)]',
        ].join(' '),
      },
      size: {
        sm: 'px-3 py-1.5 text-xs rounded-lg',
        md: 'px-4 py-2.5 text-sm rounded-xl',
        lg: 'px-6 py-3 text-base rounded-xl',
        icon: 'p-2.5 rounded-lg',
      },
    },
    defaultVariants: {
      variant: 'default',
      size: 'md',
    },
  }
);

interface ClayButtonProps
  extends Omit<HTMLMotionProps<'button'>, 'children'>,
    VariantProps<typeof buttonVariants> {
  children: React.ReactNode;
  loading?: boolean;
  leftIcon?: React.ReactNode;
  rightIcon?: React.ReactNode;
}

/**
 * Premium button with gradient effects and smooth animations
 */
export const ClayButton = forwardRef<HTMLButtonElement, ClayButtonProps>(
  ({
    children,
    variant,
    size,
    loading = false,
    leftIcon,
    rightIcon,
    className,
    disabled,
    ...motionProps
  }, ref) => {
    const isDisabled = disabled || loading;

    return (
      <motion.button
        ref={ref}
        className={cn(buttonVariants({ variant, size }), className)}
        disabled={isDisabled}
        whileHover={!isDisabled ? { scale: 1.02, y: -1 } : undefined}
        whileTap={!isDisabled ? { scale: 0.98 } : undefined}
        transition={{ duration: 0.15 }}
        {...motionProps}
      >
        {loading ? (
          <svg
            className="w-4 h-4 animate-spin"
            viewBox="0 0 24 24"
            fill="none"
          >
            <circle
              className="opacity-25"
              cx="12"
              cy="12"
              r="10"
              stroke="currentColor"
              strokeWidth="3"
            />
            <path
              className="opacity-75"
              fill="currentColor"
              d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
            />
          </svg>
        ) : leftIcon ? (
          <span className="flex-shrink-0">{leftIcon}</span>
        ) : null}

        <span>{children}</span>

        {!loading && rightIcon && (
          <span className="flex-shrink-0">{rightIcon}</span>
        )}
      </motion.button>
    );
  }
);

ClayButton.displayName = 'ClayButton';
