import { motion } from 'framer-motion';
import { Header } from './Header';
import { CinematicFooter } from '../ui/motion-footer';
import { cn } from '../../lib/utils';

interface LayoutProps {
  children: React.ReactNode;
}

// Refined ambient orb configurations
const orbConfigs = [
  {
    position: '-top-40 -right-40',
    size: 'w-[600px] h-[600px]',
    color: 'bg-[var(--color-primary)]',
    blur: 'blur-[120px]',
    scale: [1, 1.08, 1],
    opacity: [0.06, 0.1, 0.06],
    duration: 10,
    delay: 0,
  },
  {
    position: 'top-1/4 -left-40',
    size: 'w-[500px] h-[500px]',
    color: 'bg-[var(--color-accent)]',
    blur: 'blur-[100px]',
    scale: [1, 1.12, 1],
    opacity: [0.04, 0.08, 0.04],
    duration: 12,
    delay: 3,
  },
  {
    position: '-bottom-32 right-1/4',
    size: 'w-[450px] h-[450px]',
    color: 'bg-[var(--color-secondary)]',
    blur: 'blur-[90px]',
    scale: [1, 1.06, 1],
    opacity: [0.03, 0.06, 0.03],
    duration: 14,
    delay: 6,
  },
];

/**
 * Layout wrapper with ambient gradient background and cinematic footer.
 * The inner content div sits above the fixed footer via z-index stacking.
 */
export function Layout({ children }: LayoutProps) {
  return (
    <div className="relative w-full bg-[var(--color-surface)] min-h-screen">
      <div className="relative z-10 w-full bg-[var(--color-surface)]">
        {/* Ambient gradient orbs */}
        <div className="fixed inset-0 overflow-hidden pointer-events-none z-0">
          {orbConfigs.map((orb, i) => (
            <motion.div
              key={i}
              className={cn('absolute rounded-full', orb.position, orb.size, orb.color, orb.blur)}
              animate={{ scale: orb.scale, opacity: orb.opacity }}
              transition={{ duration: orb.duration, repeat: Infinity, ease: 'easeInOut', delay: orb.delay }}
            />
          ))}
        </div>

        {/* Grid pattern overlay */}
        <div
          className="fixed inset-0 pointer-events-none opacity-[0.015] dark:opacity-[0.025] z-0"
          style={{
            backgroundImage: `
              linear-gradient(var(--color-text-primary) 1px, transparent 1px),
              linear-gradient(90deg, var(--color-text-primary) 1px, transparent 1px)
            `,
            backgroundSize: '80px 80px',
          }}
        />

        <Header />

        <motion.main
          className="relative flex-1 py-10 z-0"
          initial={{ opacity: 0, y: 12 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
        >
          {children}
        </motion.main>
      </div>

      <CinematicFooter />
    </div>
  );
}
