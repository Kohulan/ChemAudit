import { AnimatePresence, motion, useReducedMotion } from 'framer-motion';
import { useEffect, useMemo, useState } from 'react';
import { useThemeContext } from '../../contexts/ThemeContext';

interface SplashScreenProps {
  isVisible: boolean;
  onComplete: () => void;
}

// ─── Caffeine molecule (C₈H₁₀N₄O₂) ─── viewBox "0 0 260 180"
const N1 = [75, 55] as const;
const C6 = [112, 40] as const;
const C5 = [148, 58] as const;
const C4 = [148, 103] as const;
const N3 = [112, 120] as const;
const C2 = [75, 103] as const;
const N7 = [182, 42] as const;
const C8 = [202, 76] as const;
const N9 = [186, 112] as const;
const O6 = [112, 6] as const;
const O2 = [42, 120] as const;
const Me1 = [42, 35] as const;
const Me3 = [112, 158] as const;
const Me7 = [216, 18] as const;

const RING6_PATH = `M ${N1} L ${C6} L ${C5} L ${C4} L ${N3} L ${C2} Z`;
const RING5_PATH = `M ${C5} L ${N7} L ${C8} L ${N9} L ${C4}`;
const SUBS_PATH = [
  `M ${C6} L ${O6}`, `M ${C2} L ${O2}`,
  `M ${N1} L ${Me1}`, `M ${N3} L ${Me3}`, `M ${N7} L ${Me7}`,
].join(' ');
const DOUBLES_PATH = [
  'M 143 63 L 143 98', 'M 197 74 L 181 110',
  'M 117 38 L 117 10', 'M 73 99 L 40 116',
].join(' ');

const HETEROATOMS: { x: number; y: number; el: 'N' | 'O' }[] = [
  { x: N1[0], y: N1[1], el: 'N' },
  { x: N3[0], y: N3[1], el: 'N' },
  { x: N7[0], y: N7[1], el: 'N' },
  { x: N9[0], y: N9[1], el: 'N' },
  { x: O6[0], y: O6[1], el: 'O' },
  { x: O2[0], y: O2[1], el: 'O' },
];

const CARBONS = [C6, C5, C4, C2, C8] as const;

const METHYLS = [
  { x: Me1[0] - 8, y: Me1[1] + 1 },
  { x: Me3[0], y: Me3[1] + 8 },
  { x: Me7[0] + 2, y: Me7[1] - 2 },
];

/**
 * Cinematic splash screen — large caffeine molecule fades in over 1.5s,
 * then booms out explosively as the brand logo reveals through the blast.
 */
export function SplashScreen({ isVisible, onComplete }: SplashScreenProps) {
  const { isDark } = useThemeContext();
  const prefersReducedMotion = useReducedMotion();
  const [exiting, setExiting] = useState(false);

  // Timeline: 0-1.5s molecule draws → 1.5s BOOM → 1.55s logo → 2.3s exit → 2.7s done
  useEffect(() => {
    if (!isVisible) return;
    setExiting(false);
    const t1 = setTimeout(() => setExiting(true), 2300);
    const t2 = setTimeout(() => { onComplete(); setExiting(false); }, 2700);
    return () => { clearTimeout(t1); clearTimeout(t2); };
  }, [isVisible, onComplete]);

  const c = isDark
    ? {
        bg: '#0c0a09',
        primary: '#f87171', primarySolid: '#dc2626',
        primaryGlow: 'rgba(248,113,113,0.30)',
        accent: '#fbbf24', accentGlow: 'rgba(251,191,36,0.20)',
        bond: 'rgba(248,113,113,0.50)', bondDouble: 'rgba(248,113,113,0.70)',
        nitrogen: '#60a5fa', oxygen: '#f87171',
        carbon: 'rgba(168,162,158,0.30)', methyl: 'rgba(168,162,158,0.40)',
        text: '#fafaf9', textMuted: 'rgba(250,250,249,0.40)',
        particle: 'rgba(248,113,113,0.40)', particleAlt: 'rgba(251,191,36,0.35)',
        flash: 'rgba(248,113,113,0.35)',
      }
    : {
        bg: '#faf9f7',
        primary: '#c41e3a', primarySolid: '#c41e3a',
        primaryGlow: 'rgba(196,30,58,0.18)',
        accent: '#d97706', accentGlow: 'rgba(217,119,6,0.12)',
        bond: 'rgba(196,30,58,0.45)', bondDouble: 'rgba(196,30,58,0.60)',
        nitrogen: '#2563eb', oxygen: '#dc2626',
        carbon: 'rgba(92,86,80,0.25)', methyl: 'rgba(92,86,80,0.35)',
        text: '#1a1815', textMuted: 'rgba(26,24,21,0.38)',
        particle: 'rgba(196,30,58,0.35)', particleAlt: 'rgba(217,119,6,0.30)',
        flash: 'rgba(196,30,58,0.25)',
      };

  const electrons = useMemo(
    () => Array.from({ length: 6 }, (_, i) => ({
      id: i,
      radius: 140 + Math.random() * 40,
      size: 2.5 + Math.random() * 2,
      duration: 5 + Math.random() * 5,
      delay: Math.random() * 2,
      startAngle: Math.random() * 360,
      ccw: i % 2 === 0,
    })),
    [],
  );

  const particles = useMemo(
    () => Array.from({ length: 10 }, (_, i) => ({
      id: i, x: 8 + Math.random() * 84, y: 8 + Math.random() * 84,
      size: 2 + Math.random() * 2.5, dur: 5 + Math.random() * 5, delay: Math.random() * 2,
    })),
    [],
  );

  const drawEase = [0.22, 1, 0.36, 1] as const;

  // Molecule keyframe timing (1.9s total):
  //   0→0.79 (0-1.5s):    fade in slowly
  //   0.79→0.815 (1.5-1.55s): brief hold at peak
  //   0.815→1 (1.55-1.9s):  BOOM — explode outward
  const molTimes = [0, 0.79, 0.815, 1];

  return (
    <AnimatePresence>
      {isVisible && (
        <motion.div
          className="fixed inset-0 z-[100] flex items-center justify-center overflow-hidden"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          transition={{ duration: 0.3 }}
        >
          {/* ─── Background ─── */}
          <div className="absolute inset-0" style={{ backgroundColor: c.bg }} />

          {/* Gradient mesh */}
          <div className="absolute inset-0 overflow-hidden">
            {[
              { w: 700, x: '50%', y: '50%', color: c.primaryGlow, blur: 90, dur: 4, delay: 0 },
              { w: 380, x: '38%', y: '38%', color: c.accentGlow, blur: 65, dur: 5.5, delay: 0.8 },
              { w: 300, x: '62%', y: '58%', color: c.primaryGlow, blur: 75, dur: 6.5, delay: 1.4 },
            ].map((b, i) => (
              <motion.div
                key={i}
                className="absolute rounded-full"
                style={{
                  width: b.w, height: b.w,
                  left: b.x, top: b.y,
                  transform: 'translate(-50%,-50%)',
                  background: `radial-gradient(circle, ${b.color} 0%, transparent 65%)`,
                  filter: `blur(${b.blur}px)`,
                }}
                animate={{
                  scale: exiting ? [1, 1.8] : [0.85, 1.08, 0.85],
                  opacity: exiting ? [0.7, 0] : [0.4, 0.85, 0.4],
                }}
                transition={{
                  duration: exiting ? 0.4 : b.dur,
                  repeat: exiting ? 0 : Infinity,
                  ease: 'easeInOut', delay: b.delay,
                }}
              />
            ))}
          </div>

          {/* Background particles */}
          {!prefersReducedMotion && particles.map(p => (
            <motion.div
              key={`p-${p.id}`}
              className="absolute rounded-full"
              style={{
                width: p.size, height: p.size,
                left: `${p.x}%`, top: `${p.y}%`,
                background: p.id % 3 === 0 ? c.particleAlt : c.particle,
              }}
              initial={{ opacity: 0, scale: 0 }}
              animate={{
                opacity: exiting ? 0 : [0, 0.65, 0],
                scale: exiting ? 0 : [0, 1, 0],
                y: exiting ? 0 : [0, -35 - Math.random() * 30],
              }}
              transition={{
                duration: p.dur, repeat: exiting ? 0 : Infinity,
                delay: p.delay, ease: 'easeOut',
              }}
            />
          ))}

          {/* ═══════ CAFFEINE MOLECULE — 200% scale (720×500) ═══════
              Fades in over 1.5s, then BOOMS outward explosively */}
          <motion.div
            className="absolute"
            style={{
              width: 720, height: 500,
              top: 'calc(50% - 250px)',
              left: 'calc(50% - 360px)',
            }}
            initial={{ opacity: 0, scale: 0.3, y: 15 }}
            animate={
              exiting
                ? { opacity: 0 }
                : {
                    opacity: [0, 0.92, 0.95, 0],
                    scale: [0.3, 1, 1.02, 3.5],
                    y: [15, 0, 0, 0],
                    filter: ['blur(3px)', 'blur(0px)', 'blur(0px)', 'blur(14px)'],
                  }
            }
            transition={
              exiting
                ? { duration: 0.15 }
                : { duration: 1.9, times: molTimes, ease: drawEase }
            }
          >
            <svg viewBox="0 0 260 180" className="w-full h-full" fill="none">
              <defs>
                <filter id="molGlow" x="-50%" y="-50%" width="200%" height="200%">
                  <feGaussianBlur stdDeviation="4" result="blur" />
                  <feMerge>
                    <feMergeNode in="blur" />
                    <feMergeNode in="SourceGraphic" />
                  </feMerge>
                </filter>
              </defs>

              {/* ── Bond glow layer ── */}
              {!prefersReducedMotion && (
                <g filter="url(#molGlow)" opacity="0.5">
                  <motion.path
                    d={RING6_PATH} stroke={c.bond} strokeWidth="3.5" strokeLinejoin="round"
                    initial={{ pathLength: 0 }} animate={{ pathLength: exiting ? 0 : 1 }}
                    transition={{ duration: 1.0, ease: drawEase, delay: 0.1 }}
                  />
                  <motion.path
                    d={RING5_PATH} stroke={c.bond} strokeWidth="3.5" strokeLinejoin="round"
                    initial={{ pathLength: 0 }} animate={{ pathLength: exiting ? 0 : 1 }}
                    transition={{ duration: 0.85, ease: drawEase, delay: 0.3 }}
                  />
                </g>
              )}

              {/* ── Crisp bond strokes ── */}
              <motion.path
                d={RING6_PATH} stroke={c.bond} strokeWidth="1.8"
                strokeLinejoin="round" strokeLinecap="round"
                initial={{ pathLength: 0 }} animate={{ pathLength: exiting ? 0 : 1 }}
                transition={{ duration: prefersReducedMotion ? 0 : 1.0, ease: drawEase, delay: 0.1 }}
              />
              <motion.path
                d={RING5_PATH} stroke={c.bond} strokeWidth="1.8"
                strokeLinejoin="round" strokeLinecap="round"
                initial={{ pathLength: 0 }} animate={{ pathLength: exiting ? 0 : 1 }}
                transition={{ duration: prefersReducedMotion ? 0 : 0.85, ease: drawEase, delay: 0.3 }}
              />
              <motion.path
                d={SUBS_PATH} stroke={c.bond} strokeWidth="1.5" strokeLinecap="round"
                initial={{ pathLength: 0 }} animate={{ pathLength: exiting ? 0 : 1 }}
                transition={{ duration: prefersReducedMotion ? 0 : 0.65, ease: drawEase, delay: 0.5 }}
              />
              <motion.path
                d={DOUBLES_PATH} stroke={c.bondDouble} strokeWidth="1.2" strokeLinecap="round"
                initial={{ pathLength: 0 }} animate={{ pathLength: exiting ? 0 : 1 }}
                transition={{ duration: prefersReducedMotion ? 0 : 0.5, ease: drawEase, delay: 0.65 }}
              />

              {/* ── Heteroatoms ── */}
              {HETEROATOMS.map((a, i) => {
                const col = a.el === 'N' ? c.nitrogen : c.oxygen;
                return (
                  <motion.g key={`ha-${i}`}>
                    <motion.circle
                      cx={a.x} cy={a.y} r="9" fill={c.bg}
                      initial={{ scale: 0, opacity: 0 }}
                      animate={{ scale: exiting ? 0 : 1, opacity: exiting ? 0 : 1 }}
                      transition={{ duration: 0.25, delay: exiting ? 0 : 0.7 + i * 0.04 }}
                    />
                    {!prefersReducedMotion && (
                      <motion.circle
                        cx={a.x} cy={a.y} r="10" fill="none" stroke={col} strokeWidth="1"
                        initial={{ scale: 0, opacity: 0 }}
                        animate={{
                          scale: exiting ? 1.5 : [1, 1.3, 1],
                          opacity: exiting ? 0 : [0.3, 0.6, 0.3],
                        }}
                        transition={{
                          duration: exiting ? 0.2 : 2.5,
                          repeat: exiting ? 0 : Infinity,
                          delay: exiting ? 0 : 0.75 + i * 0.04,
                          ease: 'easeInOut',
                        }}
                      />
                    )}
                    <motion.text
                      x={a.x} y={a.y} textAnchor="middle" dominantBaseline="central"
                      fontSize="11" fontWeight="700" fontFamily="'Source Sans 3', sans-serif"
                      fill={col}
                      initial={{ opacity: 0, scale: 0.5 }}
                      animate={{ opacity: exiting ? 0 : 0.9, scale: exiting ? 0 : 1 }}
                      transition={{ duration: 0.25, delay: exiting ? 0 : 0.75 + i * 0.04, ease: drawEase }}
                    >
                      {a.el}
                    </motion.text>
                  </motion.g>
                );
              })}

              {/* ── Carbon dots ── */}
              {CARBONS.map(([cx, cy], i) => (
                <motion.circle
                  key={`c-${i}`} cx={cx} cy={cy} r="2" fill={c.carbon}
                  initial={{ scale: 0 }}
                  animate={{ scale: exiting ? 0 : 1 }}
                  transition={{ duration: 0.15, delay: exiting ? 0 : 0.8 + i * 0.03 }}
                />
              ))}

              {/* ── Methyl labels ── */}
              {METHYLS.map((m, i) => (
                <motion.text
                  key={`me-${i}`} x={m.x} y={m.y}
                  textAnchor="middle" dominantBaseline="central"
                  fontSize="8" fontWeight="500" fontFamily="'Source Sans 3', sans-serif"
                  fill={c.methyl}
                  initial={{ opacity: 0 }}
                  animate={{ opacity: exiting ? 0 : 0.65 }}
                  transition={{ duration: 0.15, delay: exiting ? 0 : 0.9 + i * 0.04 }}
                >
                  CH₃
                </motion.text>
              ))}
            </svg>
          </motion.div>

          {/* ═══════ BOOM FLASH — radial pulse at explosion moment ═══════ */}
          {!prefersReducedMotion && (
            <motion.div
              className="absolute inset-0 pointer-events-none"
              style={{
                background: `radial-gradient(circle, ${c.flash} 0%, transparent 65%)`,
              }}
              initial={{ opacity: 0 }}
              animate={
                exiting
                  ? { opacity: 0 }
                  : {
                      opacity: [0, 0, 0, 0.6, 0],
                      scale: [1, 1, 1, 1.8, 2.5],
                    }
              }
              transition={
                exiting
                  ? { duration: 0.15 }
                  : { duration: 1.9, times: [0, 0.79, 0.815, 0.88, 1] }
              }
            />
          )}

          {/* ═══════ ORBITAL RINGS ═══════ */}
          {!prefersReducedMotion && [
            { size: 240, rx: 68, ry: 12, dur: 14, dash: '5 12', color: c.primary },
            { size: 200, rx: -60, ry: -22, dur: 11, dash: '3 16', color: c.accent },
          ].map((ring, idx) => (
            <motion.div
              key={`ring-${idx}`}
              className="absolute"
              style={{
                width: ring.size, height: ring.size,
                top: `calc(50% - ${ring.size / 2 + 25}px)`,
                left: `calc(50% - ${ring.size / 2}px)`,
                transform: `rotateX(${ring.rx}deg) rotateY(${ring.ry}deg)`,
                transformStyle: 'preserve-3d',
              }}
              initial={{ opacity: 0 }}
              animate={{ opacity: exiting ? 0 : [0.12, 0.30, 0.12] }}
              transition={{
                duration: exiting ? 0.25 : 3.5,
                repeat: exiting ? 0 : Infinity,
                ease: 'easeInOut', delay: idx * 0.5,
              }}
            >
              <motion.svg
                viewBox={`0 0 ${ring.size} ${ring.size}`} className="w-full h-full"
                animate={{ rotate: 360 }}
                transition={{ duration: ring.dur, repeat: Infinity, ease: 'linear' }}
              >
                <circle
                  cx={ring.size / 2} cy={ring.size / 2} r={ring.size / 2 - 2}
                  fill="none" stroke={ring.color} strokeWidth="0.8"
                  strokeDasharray={ring.dash} strokeLinecap="round" opacity="0.5"
                />
              </motion.svg>
            </motion.div>
          ))}

          {/* ═══════ ORBITING ELECTRONS ═══════ */}
          {!prefersReducedMotion && electrons.map(e => (
            <motion.div
              key={`e-${e.id}`}
              className="absolute rounded-full"
              style={{
                width: e.size, height: e.size,
                top: `calc(50% - ${e.size / 2 + 25}px)`,
                left: `calc(50% - ${e.size / 2}px)`,
                background: e.id % 3 === 0 ? c.accent : c.primary,
                boxShadow: `0 0 ${e.size * 3}px ${e.id % 3 === 0 ? c.accentGlow : c.primaryGlow}`,
              }}
              initial={{ opacity: 0 }}
              animate={
                exiting
                  ? { opacity: 0, scale: 0 }
                  : {
                      opacity: [0, 0.85, 0.85, 0],
                      x: Array.from({ length: 37 }, (_, i) => {
                        const a = ((e.startAngle + (e.ccw ? -1 : 1) * i * 10) * Math.PI) / 180;
                        return Math.cos(a) * e.radius;
                      }),
                      y: Array.from({ length: 37 }, (_, i) => {
                        const a = ((e.startAngle + (e.ccw ? -1 : 1) * i * 10) * Math.PI) / 180;
                        return Math.sin(a) * e.radius;
                      }),
                    }
              }
              transition={{
                duration: exiting ? 0.2 : e.duration,
                repeat: exiting ? 0 : Infinity,
                delay: e.delay, ease: 'linear',
              }}
            />
          ))}

          {/* ═══════ CENTER CONTENT — appears after boom (~1.55s) ═══════ */}
          <div className="relative z-10 flex flex-col items-center">
            {/* Logo — 3D vertical flip: starts edge-on (rotateY 90°), slowly reveals */}
            <motion.div
              className="relative"
              style={{ perspective: 600 }}
            >
              <motion.div
                className="relative"
                initial={{ rotateY: 90, opacity: 0, scale: 0.9 }}
                animate={{
                  rotateY: exiting ? [0, 90] : [90, 0],
                  opacity: exiting ? [1, 0] : [0, 1],
                  scale: exiting ? [1, 0.8] : [0.9, 1],
                }}
                transition={{
                  duration: exiting ? 0.3 : 0.7,
                  ease: [0.25, 0.1, 0.25, 1],
                  delay: exiting ? 0 : 1.5,
                }}
              >
              <motion.div
                className="absolute inset-0 rounded-full"
                style={{
                  background: `radial-gradient(circle, ${c.primaryGlow} 0%, transparent 60%)`,
                  filter: 'blur(18px)',
                }}
                animate={{
                  opacity: exiting ? 0 : [0.35, 0.75, 0.35],
                  scale: exiting ? 2.5 : [1.8, 2.3, 1.8],
                }}
                transition={{
                  duration: exiting ? 0.3 : 2.8,
                  repeat: exiting ? 0 : Infinity, ease: 'easeInOut',
                }}
              />
              <img
                src="/logo.png"
                alt="ChemAudit"
                className="relative w-[110px] h-[110px] object-contain"
                style={{
                  filter: isDark
                    ? `drop-shadow(0 0 20px ${c.primaryGlow}) drop-shadow(0 0 6px ${c.primaryGlow})`
                    : `drop-shadow(0 3px 10px rgba(0,0,0,0.10)) drop-shadow(0 0 14px ${c.primaryGlow})`,
                }}
              />
              </motion.div>
            </motion.div>

            {/* Brand name — letters assemble from scattered positions */}
            <div className="mt-6 text-center">
              <h1 className="text-[35px] font-bold tracking-tight font-display flex justify-center">
                {[
                  { ch: 'C', primary: true, x: -40, y: -30, r: -25 },
                  { ch: 'h', primary: true, x: 20, y: -45, r: 15 },
                  { ch: 'e', primary: true, x: -30, y: 35, r: -20 },
                  { ch: 'm', primary: true, x: 45, y: 25, r: 30 },
                  { ch: 'A', primary: false, x: -35, y: -40, r: 20 },
                  { ch: 'u', primary: false, x: 30, y: 40, r: -15 },
                  { ch: 'd', primary: false, x: -25, y: -35, r: 25 },
                  { ch: 'i', primary: false, x: 40, y: -25, r: -30 },
                  { ch: 't', primary: false, x: -20, y: 30, r: 18 },
                ].map((l, i) => (
                  <motion.span
                    key={i}
                    className={l.primary ? 'font-extrabold' : 'font-semibold'}
                    style={{ color: l.primary ? c.primarySolid : c.text, display: 'inline-block' }}
                    initial={{ opacity: 0, x: l.x, y: l.y, rotate: l.r, scale: 0.3 }}
                    animate={{
                      opacity: exiting ? 0 : 1,
                      x: exiting ? l.x * 2 : 0,
                      y: exiting ? l.y * 2 : 0,
                      rotate: exiting ? l.r * 2 : 0,
                      scale: exiting ? 0 : 1,
                    }}
                    transition={{
                      duration: exiting ? 0.2 : 0.5,
                      ease: drawEase,
                      delay: exiting ? i * 0.02 : 1.7 + i * 0.04,
                    }}
                  >
                    {l.ch}
                  </motion.span>
                ))}
              </h1>

              <motion.div
                className="h-[1.5px] mt-2.5 mx-auto rounded-full"
                style={{
                  background: `linear-gradient(90deg, transparent, ${c.primarySolid}, ${c.accent}, transparent)`,
                }}
                initial={{ width: 0, opacity: 0 }}
                animate={{ width: exiting ? 0 : 156, opacity: exiting ? 0 : 1 }}
                transition={{ delay: exiting ? 0 : 1.85, duration: 0.4, ease: drawEase }}
              />

              <motion.p
                className="text-[13px] mt-3 tracking-[0.22em] uppercase font-medium"
                style={{ color: c.textMuted }}
                initial={{ opacity: 0, y: 5 }}
                animate={{ opacity: exiting ? 0 : 1, y: exiting ? -6 : 0 }}
                transition={{ delay: exiting ? 0 : 1.92, duration: 0.3 }}
              >
                Chemical Structure Validation
              </motion.p>
            </div>

            {/* Sweep loader */}
            <motion.div
              className="mt-7 overflow-hidden rounded-full"
              style={{ width: 112, height: 1.5, backgroundColor: c.carbon }}
              initial={{ opacity: 0 }}
              animate={{ opacity: exiting ? 0 : 1 }}
              transition={{ delay: exiting ? 0 : 2.0, duration: 0.25 }}
            >
              <motion.div
                className="h-full rounded-full"
                style={{
                  background: `linear-gradient(90deg, transparent, ${c.primarySolid}, ${c.accent}, transparent)`,
                  width: '35%',
                }}
                animate={{ x: exiting ? 100 : ['-35%', '135%'] }}
                transition={{
                  duration: exiting ? 0.15 : 1.1,
                  repeat: exiting ? 0 : Infinity, ease: 'easeInOut',
                }}
              />
            </motion.div>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
