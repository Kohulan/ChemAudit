import { motion } from 'framer-motion';
import { Header } from './Header';
import { cn } from '../../lib/utils';

interface LayoutProps {
  children: React.ReactNode;
}

/**
 * Premium layout wrapper with ambient gradient background
 */
export function Layout({ children }: LayoutProps) {
  return (
    <div className="min-h-screen bg-[var(--color-surface)]">
      {/* Ambient gradient orbs - subtle and premium */}
      <div className="fixed inset-0 overflow-hidden pointer-events-none">
        {/* Primary orb - top right */}
        <motion.div
          className={cn(
            'absolute -top-32 -right-32 w-[500px] h-[500px] rounded-full',
            'bg-[var(--color-primary)]/8 dark:bg-[var(--color-primary)]/6',
            'blur-[100px]'
          )}
          animate={{
            scale: [1, 1.1, 1],
            opacity: [0.5, 0.7, 0.5],
          }}
          transition={{
            duration: 8,
            repeat: Infinity,
            ease: 'easeInOut',
          }}
        />
        {/* Accent orb - left */}
        <motion.div
          className={cn(
            'absolute top-1/3 -left-32 w-[400px] h-[400px] rounded-full',
            'bg-[var(--color-accent)]/6 dark:bg-[var(--color-accent)]/4',
            'blur-[100px]'
          )}
          animate={{
            scale: [1, 1.15, 1],
            opacity: [0.4, 0.6, 0.4],
          }}
          transition={{
            duration: 10,
            repeat: Infinity,
            ease: 'easeInOut',
            delay: 2,
          }}
        />
        {/* Bottom orb */}
        <motion.div
          className={cn(
            'absolute -bottom-20 right-1/4 w-[350px] h-[350px] rounded-full',
            'bg-[var(--color-primary)]/5 dark:bg-[var(--color-accent)]/4',
            'blur-[80px]'
          )}
          animate={{
            scale: [1, 1.08, 1],
            opacity: [0.3, 0.5, 0.3],
          }}
          transition={{
            duration: 12,
            repeat: Infinity,
            ease: 'easeInOut',
            delay: 4,
          }}
        />
      </div>

      {/* Subtle grid pattern */}
      <div
        className="fixed inset-0 pointer-events-none opacity-[0.02] dark:opacity-[0.03]"
        style={{
          backgroundImage: `
            linear-gradient(var(--color-primary) 1px, transparent 1px),
            linear-gradient(90deg, var(--color-primary) 1px, transparent 1px)
          `,
          backgroundSize: '60px 60px',
        }}
      />

      {/* Content */}
      <div className="relative">
        <Header />

        <motion.main
          className="py-8"
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.4, ease: [0.25, 0.46, 0.45, 0.94] }}
        >
          {children}
        </motion.main>

        {/* Footer */}
        <footer className="border-t border-[var(--color-border)] mt-16">
          <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
            <div className="flex flex-col sm:flex-row items-center justify-between gap-4">
              <div className="flex items-center gap-3">
                <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[var(--color-primary)] to-[var(--color-accent)] flex items-center justify-center">
                  <svg className="w-4 h-4 text-white" viewBox="0 0 24 24" fill="currentColor">
                    <circle cx="12" cy="12" r="3" />
                    <circle cx="12" cy="4" r="1.5" />
                    <circle cx="18.9" cy="8" r="1.5" />
                    <circle cx="18.9" cy="16" r="1.5" />
                    <circle cx="12" cy="20" r="1.5" />
                    <circle cx="5.1" cy="16" r="1.5" />
                    <circle cx="5.1" cy="8" r="1.5" />
                  </svg>
                </div>
                <p className="text-sm text-[var(--color-text-muted)]">
                  ChemVault - Chemical Structure Validation
                </p>
              </div>
              <div className="flex items-center gap-6">
                <FooterLink href="http://localhost:8000/docs">
                  API Docs
                </FooterLink>
                <FooterLink href="https://github.com">
                  GitHub
                </FooterLink>
              </div>
            </div>
          </div>
        </footer>
      </div>
    </div>
  );
}

interface FooterLinkProps {
  href: string;
  children: React.ReactNode;
}

function FooterLink({ href, children }: FooterLinkProps) {
  return (
    <a
      href={href}
      target="_blank"
      rel="noopener noreferrer"
      className={cn(
        'text-sm transition-all duration-200',
        'text-[var(--color-text-muted)]',
        'hover:text-[var(--color-primary)]'
      )}
    >
      {children}
    </a>
  );
}
