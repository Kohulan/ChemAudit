import { motion, AnimatePresence } from 'framer-motion';
import { useEffect, useState } from 'react';
import { useThemeContext } from '../../contexts/ThemeContext';

interface SplashScreenProps {
  isVisible: boolean;
  onComplete: () => void;
}

/**
 * Theme-aware splash screen with elegant logo reveal
 * Uses app's core colors: crimson/rose primary + amber accent
 */
export function SplashScreen({ isVisible, onComplete }: SplashScreenProps) {
  const { resolvedTheme } = useThemeContext();
  const isDark = resolvedTheme === 'dark';
  const [phase, setPhase] = useState<'enter' | 'reveal' | 'exit'>('enter');

  useEffect(() => {
    if (isVisible) {
      setPhase('enter');
      const revealTimer = setTimeout(() => setPhase('reveal'), 100);
      const exitTimer = setTimeout(() => setPhase('exit'), 1000);
      const completeTimer = setTimeout(() => {
        onComplete();
        setPhase('enter');
      }, 1400);

      return () => {
        clearTimeout(revealTimer);
        clearTimeout(exitTimer);
        clearTimeout(completeTimer);
      };
    }
  }, [isVisible, onComplete]);

  // Theme-aware colors
  const colors = isDark
    ? {
        bg: 'from-slate-950 via-slate-900 to-slate-950',
        glow: 'rgba(248,113,113,0.2)', // primary-rgb for dark
        glowAccent: 'rgba(251,191,36,0.15)',
        primary: '#f87171',
        primaryDark: '#ef4444',
        accent: '#fbbf24',
        accentDark: '#f59e0b',
        text: '#ffffff',
        textMuted: 'rgba(255,255,255,0.6)',
        cardBg: 'rgba(30,30,40,0.95)',
        cardBorder: 'rgba(255,255,255,0.1)',
      }
    : {
        bg: 'from-stone-100 via-amber-50/30 to-rose-50/20',
        glow: 'rgba(196,30,58,0.12)', // primary-rgb for light
        glowAccent: 'rgba(217,119,6,0.1)',
        primary: '#c41e3a',
        primaryDark: '#9d1830',
        accent: '#d97706',
        accentDark: '#b45309',
        text: '#1a1a1a',
        textMuted: 'rgba(0,0,0,0.5)',
        cardBg: 'rgba(255,255,255,0.95)',
        cardBorder: 'rgba(196,30,58,0.15)',
      };

  return (
    <AnimatePresence>
      {isVisible && (
        <motion.div
          className="fixed inset-0 z-[100] flex items-center justify-center overflow-hidden"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          transition={{ duration: 0.15 }}
        >
          {/* Background */}
          <div className={`absolute inset-0 bg-gradient-to-br ${colors.bg}`} />

          {/* Ambient glow - primary */}
          <motion.div
            className="absolute w-[500px] h-[500px] rounded-full"
            style={{
              background: `radial-gradient(circle, ${colors.glow} 0%, transparent 70%)`,
            }}
            initial={{ scale: 0.6, opacity: 0 }}
            animate={{
              scale: phase === 'exit' ? 1.8 : 1,
              opacity: phase === 'exit' ? 0 : 1,
            }}
            transition={{ duration: 0.6, ease: 'easeOut' }}
          />

          {/* Ambient glow - accent (offset) */}
          <motion.div
            className="absolute w-[400px] h-[400px] rounded-full translate-x-20 translate-y-10"
            style={{
              background: `radial-gradient(circle, ${colors.glowAccent} 0%, transparent 70%)`,
            }}
            initial={{ scale: 0.6, opacity: 0 }}
            animate={{
              scale: phase === 'exit' ? 1.5 : 1,
              opacity: phase === 'exit' ? 0 : 0.8,
            }}
            transition={{ duration: 0.6, ease: 'easeOut', delay: 0.1 }}
          />

          {/* Center content */}
          <div className="relative flex flex-col items-center">
            {/* Rotating ring behind logo */}
            <motion.div
              className="absolute w-36 h-36"
              initial={{ opacity: 0, scale: 0.5, rotate: -90 }}
              animate={{
                opacity: phase === 'exit' ? 0 : 0.6,
                scale: phase === 'exit' ? 1.5 : 1,
                rotate: phase === 'exit' ? 90 : 0,
              }}
              transition={{
                duration: phase === 'exit' ? 0.3 : 0.5,
                ease: [0.22, 1, 0.36, 1],
              }}
            >
              <svg viewBox="0 0 144 144" className="w-full h-full">
                <defs>
                  <linearGradient id="ringGrad" x1="0%" y1="0%" x2="100%" y2="100%">
                    <stop offset="0%" stopColor={colors.primary} stopOpacity="0.8" />
                    <stop offset="50%" stopColor={colors.accent} stopOpacity="0.6" />
                    <stop offset="100%" stopColor={colors.primary} stopOpacity="0.8" />
                  </linearGradient>
                </defs>
                <circle
                  cx="72"
                  cy="72"
                  r="68"
                  fill="none"
                  stroke="url(#ringGrad)"
                  strokeWidth="2"
                  strokeDasharray="12 8"
                  strokeLinecap="round"
                />
              </svg>
            </motion.div>

            {/* Logo container */}
            <motion.div
              className="relative w-20 h-20 rounded-2xl flex items-center justify-center overflow-hidden z-10"
              style={{
                background: colors.cardBg,
                boxShadow: isDark
                  ? `0 0 0 1px ${colors.cardBorder}, 0 25px 50px -12px rgba(0,0,0,0.5), 0 0 40px ${colors.glow}`
                  : `0 0 0 1px ${colors.cardBorder}, 0 25px 50px -12px rgba(0,0,0,0.15), 0 0 40px ${colors.glow}`,
              }}
              initial={{ scale: 0, rotate: -90 }}
              animate={{
                scale: phase === 'exit' ? [1, 1.1, 0] : [0, 1.08, 1],
                rotate: phase === 'exit' ? 90 : [-90, 5, 0],
              }}
              transition={{
                duration: phase === 'exit' ? 0.3 : 0.5,
                ease: [0.22, 1, 0.36, 1],
              }}
            >
              <img
                src="/logo.png"
                alt="ChemVault"
                className="w-14 h-14 object-contain"
              />
            </motion.div>

            {/* Brand name */}
            <motion.div
              className="mt-5 text-center"
              initial={{ opacity: 0, y: 15 }}
              animate={{
                opacity: phase === 'exit' ? 0 : 1,
                y: phase === 'exit' ? -10 : 0,
              }}
              transition={{
                delay: phase === 'exit' ? 0 : 0.25,
                duration: 0.35,
                ease: [0.22, 1, 0.36, 1],
              }}
            >
              <h1 className="text-2xl font-bold tracking-tight">
                <span style={{ color: colors.primary }} className="font-extrabold">
                  Chem
                </span>
                <span style={{ color: colors.text }} className="font-semibold">
                  Vault
                </span>
              </h1>

              {/* Animated line */}
              <motion.div
                className="h-[2px] mt-2 mx-auto rounded-full"
                style={{
                  background: `linear-gradient(90deg, transparent, ${colors.primary}, ${colors.accent}, transparent)`,
                }}
                initial={{ width: 0, opacity: 0 }}
                animate={{
                  width: phase === 'exit' ? 0 : 100,
                  opacity: phase === 'exit' ? 0 : 1,
                }}
                transition={{
                  delay: phase === 'exit' ? 0 : 0.35,
                  duration: 0.3,
                  ease: [0.22, 1, 0.36, 1],
                }}
              />

              {/* Tagline */}
              <motion.p
                className="text-xs mt-2 tracking-widest uppercase font-medium"
                style={{ color: colors.textMuted }}
                initial={{ opacity: 0 }}
                animate={{ opacity: phase === 'exit' ? 0 : 1 }}
                transition={{
                  delay: phase === 'exit' ? 0 : 0.4,
                  duration: 0.25,
                }}
              >
                Structure Validation
              </motion.p>
            </motion.div>

            {/* Orbiting dots */}
            <motion.div
              className="absolute w-32 h-32 pointer-events-none"
              initial={{ opacity: 0 }}
              animate={{
                opacity: phase === 'exit' ? 0 : 1,
                rotate: 360,
              }}
              transition={{
                opacity: { duration: 0.3 },
                rotate: { duration: 4, repeat: Infinity, ease: 'linear' },
              }}
            >
              {[0, 180].map((angle) => (
                <motion.div
                  key={angle}
                  className="absolute w-2 h-2 rounded-full"
                  style={{
                    background: angle === 0 ? colors.primary : colors.accent,
                    boxShadow: `0 0 8px ${angle === 0 ? colors.primary : colors.accent}`,
                    top: '50%',
                    left: '50%',
                    transform: `rotate(${angle}deg) translateX(60px) translateY(-50%)`,
                  }}
                />
              ))}
            </motion.div>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
