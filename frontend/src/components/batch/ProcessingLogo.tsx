import { motion } from 'framer-motion';
import { cn } from '../../lib/utils';

interface ProcessingLogoProps {
  progress?: number;
  status?: 'pending' | 'processing' | 'complete' | 'failed' | 'cancelled';
  className?: string;
}

/**
 * Stunning animated logo for batch processing visualization.
 * Uses Laboratory Crimson (#c41e3a) and Warm Amber (#d97706) color scheme.
 * Features orbital particles, pulsing glow, and progress-reactive animations.
 */
export function ProcessingLogo({ progress = 0, status = 'processing', className }: ProcessingLogoProps) {
  const isActive = status === 'processing' || status === 'pending';
  const isComplete = status === 'complete';
  const isFailed = status === 'failed' || status === 'cancelled';

  return (
    <div className={cn('relative flex items-center justify-center', className)}>
      {/* Outer glow rings - using primary crimson and accent amber */}
      <motion.div
        className="absolute inset-0 rounded-full"
        style={{
          background: isComplete
            ? 'radial-gradient(circle, rgba(234,179,8,0.3) 0%, transparent 70%)'
            : isFailed
            ? 'radial-gradient(circle, rgba(239,68,68,0.3) 0%, transparent 70%)'
            : 'radial-gradient(circle, rgba(196,30,58,0.4) 0%, rgba(217,119,6,0.2) 50%, transparent 70%)',
        }}
        animate={isActive ? {
          scale: [1, 1.2, 1],
          opacity: [0.5, 0.8, 0.5],
        } : {}}
        transition={{
          duration: 2,
          repeat: Infinity,
          ease: 'easeInOut',
        }}
      />

      {/* Secondary pulse ring - crimson */}
      <motion.div
        className="absolute w-40 h-40 rounded-full border-2"
        style={{
          borderColor: isComplete
            ? 'rgba(234,179,8,0.4)'
            : isFailed
            ? 'rgba(239,68,68,0.4)'
            : 'rgba(196,30,58,0.3)',
        }}
        animate={isActive ? {
          scale: [1, 1.5, 1],
          opacity: [0.6, 0, 0.6],
        } : {}}
        transition={{
          duration: 2.5,
          repeat: Infinity,
          ease: 'easeOut',
        }}
      />

      {/* Tertiary pulse ring - amber accent */}
      <motion.div
        className="absolute w-48 h-48 rounded-full border"
        style={{
          borderColor: isComplete
            ? 'rgba(234,179,8,0.2)'
            : isFailed
            ? 'rgba(239,68,68,0.2)'
            : 'rgba(217,119,6,0.2)',
        }}
        animate={isActive ? {
          scale: [1.1, 1.8, 1.1],
          opacity: [0.4, 0, 0.4],
        } : {}}
        transition={{
          duration: 3,
          repeat: Infinity,
          ease: 'easeOut',
          delay: 0.5,
        }}
      />

      {/* Rotating orbital ring 1 - primary crimson */}
      <motion.div
        className="absolute w-36 h-36"
        animate={{ rotate: 360 }}
        transition={{
          duration: 8,
          repeat: Infinity,
          ease: 'linear',
        }}
      >
        <div className="absolute w-full h-full rounded-full border border-dashed border-[var(--color-primary)]/30" />
        {/* Orbital particle 1 - crimson */}
        <motion.div
          className="absolute w-3 h-3 rounded-full bg-gradient-to-r from-[#c41e3a] to-[#e11d48]"
          style={{
            top: 0,
            left: '50%',
            marginLeft: -6,
            boxShadow: '0 0 10px rgba(196,30,58,0.6)',
          }}
          animate={{
            scale: [1, 1.3, 1],
            boxShadow: [
              '0 0 10px rgba(196,30,58,0.5)',
              '0 0 20px rgba(196,30,58,0.8)',
              '0 0 10px rgba(196,30,58,0.5)',
            ],
          }}
          transition={{
            duration: 1.5,
            repeat: Infinity,
            ease: 'easeInOut',
          }}
        />
        {/* Orbital particle 2 - amber */}
        <motion.div
          className="absolute w-2 h-2 rounded-full bg-gradient-to-r from-[#d97706] to-[#f59e0b]"
          style={{
            bottom: 0,
            left: '50%',
            marginLeft: -4,
            boxShadow: '0 0 8px rgba(217,119,6,0.5)',
          }}
          animate={{
            scale: [1.2, 0.8, 1.2],
          }}
          transition={{
            duration: 2,
            repeat: Infinity,
            ease: 'easeInOut',
          }}
        />
      </motion.div>

      {/* Rotating orbital ring 2 - Counter rotation */}
      <motion.div
        className="absolute w-44 h-44"
        animate={{ rotate: -360 }}
        transition={{
          duration: 12,
          repeat: Infinity,
          ease: 'linear',
        }}
      >
        <div className="absolute w-full h-full rounded-full border border-dotted border-[var(--color-accent)]/20" />
        {/* Orbital particle 3 - warm rose */}
        <motion.div
          className="absolute w-2.5 h-2.5 rounded-full bg-gradient-to-r from-[#e11d48] to-[#f43f5e]"
          style={{
            top: '50%',
            right: 0,
            marginTop: -5,
            boxShadow: '0 0 8px rgba(225,29,72,0.5)',
          }}
          animate={{
            scale: [1, 1.4, 1],
          }}
          transition={{
            duration: 1.8,
            repeat: Infinity,
            ease: 'easeInOut',
            delay: 0.3,
          }}
        />
        {/* Orbital particle 4 - golden amber */}
        <motion.div
          className="absolute w-2 h-2 rounded-full bg-gradient-to-r from-[#f59e0b] to-[#fbbf24]"
          style={{
            top: '50%',
            left: 0,
            marginTop: -4,
            boxShadow: '0 0 6px rgba(245,158,11,0.5)',
          }}
        />
      </motion.div>

      {/* Rotating orbital ring 3 - Tilted */}
      <motion.div
        className="absolute w-32 h-32"
        style={{ transform: 'rotateX(60deg) rotateZ(0deg)' }}
        animate={{ rotateZ: 360 }}
        transition={{
          duration: 6,
          repeat: Infinity,
          ease: 'linear',
        }}
      >
        <div className="absolute w-full h-full rounded-full border border-[var(--color-primary)]/20" />
        <motion.div
          className="absolute w-2 h-2 rounded-full bg-[var(--color-accent)]"
          style={{
            top: 0,
            left: '50%',
            marginLeft: -4,
            boxShadow: '0 0 8px rgba(217,119,6,0.6)',
          }}
        />
      </motion.div>

      {/* Progress ring - using primary to accent gradient */}
      <svg className="absolute w-32 h-32" viewBox="0 0 100 100">
        {/* Background ring */}
        <circle
          cx="50"
          cy="50"
          r="46"
          fill="none"
          stroke="currentColor"
          strokeWidth="3"
          className="text-gray-200 dark:text-gray-700"
        />
        {/* Progress arc - crimson to amber gradient */}
        <motion.circle
          cx="50"
          cy="50"
          r="46"
          fill="none"
          stroke="url(#progressGradient)"
          strokeWidth="4"
          strokeLinecap="round"
          strokeDasharray={`${2 * Math.PI * 46}`}
          strokeDashoffset={`${2 * Math.PI * 46 * (1 - progress / 100)}`}
          style={{
            transformOrigin: 'center',
            transform: 'rotate(-90deg)',
          }}
          animate={{
            strokeDashoffset: `${2 * Math.PI * 46 * (1 - progress / 100)}`,
          }}
          transition={{ duration: 0.5, ease: 'easeOut' }}
        />
        <defs>
          <linearGradient id="progressGradient" x1="0%" y1="0%" x2="100%" y2="100%">
            <stop offset="0%" stopColor="#c41e3a" />
            <stop offset="50%" stopColor="#e11d48" />
            <stop offset="100%" stopColor="#d97706" />
          </linearGradient>
        </defs>
      </svg>

      {/* Logo container with glow - crimson/amber glow */}
      <motion.div
        className="relative w-24 h-24 rounded-full overflow-hidden"
        style={{
          boxShadow: isComplete
            ? '0 0 40px rgba(234,179,8,0.6), 0 0 80px rgba(234,179,8,0.3)'
            : isFailed
            ? '0 0 40px rgba(239,68,68,0.5), 0 0 60px rgba(239,68,68,0.2)'
            : '0 0 30px rgba(196,30,58,0.5), 0 0 60px rgba(217,119,6,0.3), 0 0 90px rgba(196,30,58,0.2)',
        }}
        animate={isActive ? {
          boxShadow: [
            '0 0 30px rgba(196,30,58,0.5), 0 0 60px rgba(217,119,6,0.3), 0 0 90px rgba(196,30,58,0.2)',
            '0 0 50px rgba(217,119,6,0.6), 0 0 80px rgba(196,30,58,0.4), 0 0 120px rgba(217,119,6,0.3)',
            '0 0 30px rgba(196,30,58,0.5), 0 0 60px rgba(217,119,6,0.3), 0 0 90px rgba(196,30,58,0.2)',
          ],
        } : {}}
        transition={{
          duration: 2,
          repeat: Infinity,
          ease: 'easeInOut',
        }}
      >
        {/* Inner glow overlay */}
        <div className="absolute inset-0 bg-gradient-to-br from-white/20 to-transparent pointer-events-none z-10" />

        {/* Logo image */}
        <motion.img
          src="/logo.png"
          alt="ChemVault"
          className="w-full h-full object-cover"
          animate={isActive ? {
            scale: [1, 1.05, 1],
          } : isComplete ? {
            scale: [1, 1.1, 1],
          } : {}}
          transition={{
            duration: isComplete ? 0.5 : 3,
            repeat: isComplete ? 0 : Infinity,
            ease: 'easeInOut',
          }}
        />

        {/* Shimmer effect - warm amber tint */}
        {isActive && (
          <motion.div
            className="absolute inset-0 bg-gradient-to-r from-transparent via-amber-100/40 to-transparent"
            style={{ transform: 'skewX(-20deg)' }}
            animate={{
              x: ['-200%', '200%'],
            }}
            transition={{
              duration: 2,
              repeat: Infinity,
              ease: 'easeInOut',
              repeatDelay: 1,
            }}
          />
        )}
      </motion.div>

      {/* Floating particles - crimson and amber */}
      {isActive && (
        <>
          {[...Array(6)].map((_, i) => (
            <motion.div
              key={i}
              className="absolute w-1.5 h-1.5 rounded-full"
              style={{
                background: `linear-gradient(135deg, ${
                  ['#c41e3a', '#d97706', '#e11d48', '#f59e0b', '#9d1830', '#b45309'][i]
                }, ${
                  ['#e11d48', '#f59e0b', '#f43f5e', '#fbbf24', '#c41e3a', '#d97706'][i]
                })`,
                left: '50%',
                top: '50%',
              }}
              animate={{
                x: [0, (Math.random() - 0.5) * 120],
                y: [0, (Math.random() - 0.5) * 120],
                opacity: [0, 1, 0],
                scale: [0, 1.5, 0],
              }}
              transition={{
                duration: 2 + Math.random() * 2,
                repeat: Infinity,
                delay: i * 0.4,
                ease: 'easeOut',
              }}
            />
          ))}
        </>
      )}

      {/* Success burst effect - golden amber */}
      {isComplete && (
        <motion.div
          className="absolute inset-0 flex items-center justify-center"
          initial={{ scale: 0.8, opacity: 0 }}
          animate={{ scale: 1, opacity: 1 }}
        >
          {[...Array(12)].map((_, i) => (
            <motion.div
              key={i}
              className="absolute w-2 h-2 rounded-full bg-gradient-to-r from-yellow-400 to-amber-500"
              style={{
                transformOrigin: 'center',
              }}
              initial={{
                x: 0,
                y: 0,
                opacity: 1,
              }}
              animate={{
                x: Math.cos((i * Math.PI * 2) / 12) * 100,
                y: Math.sin((i * Math.PI * 2) / 12) * 100,
                opacity: 0,
                scale: [1, 0],
              }}
              transition={{
                duration: 1,
                ease: 'easeOut',
              }}
            />
          ))}
        </motion.div>
      )}
    </div>
  );
}
