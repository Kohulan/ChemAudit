import { useRef } from 'react';
import { motion, useScroll, useTransform, useInView } from 'framer-motion';
import {
  Building2,
  Mail,
  Globe,
  Github,
  Code2,
  Database,
  Server,
  Layout,
  TestTube,
  Heart,
  ExternalLink,
  MapPin,
  Sparkles,
  Zap,
  Shield,
  BarChart3,
} from 'lucide-react';
import { cn } from '../lib/utils';

// Floating molecule configurations for background
const moleculeConfigs = [
  { x: '10%', y: '15%', size: 60, delay: 0, duration: 20, opacity: 0.15 },
  { x: '85%', y: '25%', size: 80, delay: 2, duration: 25, opacity: 0.12 },
  { x: '75%', y: '60%', size: 50, delay: 4, duration: 18, opacity: 0.18 },
  { x: '15%', y: '70%', size: 70, delay: 1, duration: 22, opacity: 0.14 },
  { x: '50%', y: '85%', size: 45, delay: 3, duration: 16, opacity: 0.16 },
  { x: '90%', y: '80%', size: 55, delay: 5, duration: 19, opacity: 0.13 },
  { x: '5%', y: '45%', size: 65, delay: 2.5, duration: 21, opacity: 0.15 },
];

// Glowing orb configurations
const glowOrbConfigs = [
  { x: '20%', y: '30%', size: 300, color: 'var(--color-primary)', opacity: 0.08, blur: 100, duration: 8 },
  { x: '70%', y: '50%', size: 250, color: 'var(--color-accent)', opacity: 0.06, blur: 80, duration: 10 },
  { x: '40%', y: '80%', size: 200, color: 'var(--color-primary)', opacity: 0.05, blur: 60, duration: 12 },
];

/**
 * Stunning About page with claymorphism, animations, and visual effects
 */
export function AboutPage() {
  const containerRef = useRef<HTMLDivElement>(null);
  const { scrollYProgress } = useScroll({ target: containerRef });
  const heroY = useTransform(scrollYProgress, [0, 1], [0, -100]);

  return (
    <div ref={containerRef} className="relative min-h-screen overflow-hidden">
      {/* Animated gradient mesh background */}
      <div className="fixed inset-0 pointer-events-none overflow-hidden">
        {/* Aurora gradient effect */}
        <motion.div
          className="absolute inset-0"
          style={{
            background: `
              radial-gradient(ellipse 80% 50% at 20% 40%, rgba(var(--color-primary-rgb, 220, 38, 38), 0.08) 0%, transparent 50%),
              radial-gradient(ellipse 60% 40% at 80% 60%, rgba(var(--color-accent-rgb, 234, 88, 12), 0.06) 0%, transparent 50%),
              radial-gradient(ellipse 50% 30% at 50% 90%, rgba(var(--color-primary-rgb, 220, 38, 38), 0.05) 0%, transparent 50%)
            `,
          }}
          animate={{
            opacity: [0.8, 1, 0.8],
          }}
          transition={{
            duration: 8,
            repeat: Infinity,
            ease: 'easeInOut',
          }}
        />

        {/* Animated glowing orbs */}
        {glowOrbConfigs.map((orb, i) => (
          <motion.div
            key={i}
            className="absolute rounded-full"
            style={{
              left: orb.x,
              top: orb.y,
              width: orb.size,
              height: orb.size,
              background: orb.color,
              opacity: orb.opacity,
              filter: `blur(${orb.blur}px)`,
              transform: 'translate(-50%, -50%)',
            }}
            animate={{
              scale: [1, 1.2, 1],
              opacity: [orb.opacity, orb.opacity * 1.5, orb.opacity],
            }}
            transition={{
              duration: orb.duration,
              repeat: Infinity,
              ease: 'easeInOut',
            }}
          />
        ))}

        {/* Floating molecules */}
        {moleculeConfigs.map((mol, i) => (
          <FloatingMolecule key={i} {...mol} />
        ))}
      </div>

      <div className="relative max-w-6xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Hero Section with parallax */}
        <motion.div style={{ y: heroY }} className="mb-16">
          <HeroSection />
        </motion.div>

        {/* Main Content Grid */}
        <div className="grid grid-cols-1 lg:grid-cols-12 gap-6 mb-16">
          {/* What is ChemVault - Large card spanning 8 columns */}
          <AnimatedCard className="lg:col-span-8" delay={0.1}>
            <WhatIsChemVault />
          </AnimatedCard>

          {/* Research Group - Tall card spanning 4 columns, 2 rows */}
          <AnimatedCard className="lg:col-span-4 lg:row-span-2" delay={0.2}>
            <ResearchGroup />
          </AnimatedCard>

          {/* Tech Stack - Wide card */}
          <AnimatedCard className="lg:col-span-8" delay={0.3}>
            <TechStack />
          </AnimatedCard>
        </div>

        {/* Bottom row */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 mb-16">
          <AnimatedCard delay={0.4}>
            <SourceCode />
          </AnimatedCard>
          <AnimatedCard delay={0.5}>
            <Contact />
          </AnimatedCard>
          <AnimatedCard delay={0.6} className="md:col-span-2 lg:col-span-1">
            <QuickStats />
          </AnimatedCard>
        </div>

        {/* Acknowledgments - Full width */}
        <AnimatedCard delay={0.7} className="mb-12">
          <Acknowledgments />
        </AnimatedCard>

        {/* License Footer */}
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ delay: 0.8, duration: 0.5 }}
          className="text-center text-sm text-[var(--color-text-muted)]"
        >
          <p>
            ChemVault is open-source software released under the{' '}
            <a
              href="https://opensource.org/licenses/MIT"
              className="text-[var(--color-primary)] hover:underline transition-colors"
            >
              MIT License
            </a>
          </p>
        </motion.div>
      </div>
    </div>
  );
}

// ============================================================================
// FLOATING MOLECULE COMPONENT
// ============================================================================

interface FloatingMoleculeProps {
  x: string;
  y: string;
  size: number;
  delay: number;
  duration: number;
  opacity: number;
}

function FloatingMolecule({ x, y, size, delay, duration, opacity }: FloatingMoleculeProps) {
  return (
    <motion.div
      className="absolute pointer-events-none"
      style={{ left: x, top: y }}
      initial={{ opacity: 0, scale: 0 }}
      animate={{
        opacity: [0, opacity, opacity, 0],
        scale: [0.8, 1, 1, 0.8],
        y: [0, -30, -60, -90],
        rotate: [0, 90, 180, 270],
      }}
      transition={{
        duration,
        delay,
        repeat: Infinity,
        ease: 'easeInOut',
      }}
    >
      <svg
        width={size}
        height={size}
        viewBox="0 0 100 100"
        fill="none"
        className="text-[var(--color-primary)]"
      >
        {/* Benzene-like ring */}
        <circle cx="50" cy="50" r="25" stroke="currentColor" strokeWidth="2" opacity="0.5" />
        {/* Atoms */}
        <circle cx="50" cy="25" r="6" fill="currentColor" opacity="0.6" />
        <circle cx="71.65" cy="37.5" r="6" fill="currentColor" opacity="0.6" />
        <circle cx="71.65" cy="62.5" r="6" fill="currentColor" opacity="0.6" />
        <circle cx="50" cy="75" r="6" fill="currentColor" opacity="0.6" />
        <circle cx="28.35" cy="62.5" r="6" fill="currentColor" opacity="0.6" />
        <circle cx="28.35" cy="37.5" r="6" fill="currentColor" opacity="0.6" />
        {/* Bonds */}
        <line x1="50" y1="31" x2="50" y2="69" stroke="currentColor" strokeWidth="1.5" opacity="0.3" />
      </svg>
    </motion.div>
  );
}

// ============================================================================
// ANIMATED CARD WRAPPER
// ============================================================================

interface AnimatedCardProps {
  children: React.ReactNode;
  className?: string;
  delay?: number;
}

function AnimatedCard({ children, className, delay = 0 }: AnimatedCardProps) {
  const ref = useRef<HTMLDivElement>(null);
  const isInView = useInView(ref, { once: true, margin: '-50px' });

  return (
    <motion.div
      ref={ref}
      initial={{ opacity: 0, y: 40, scale: 0.95 }}
      animate={isInView ? { opacity: 1, y: 0, scale: 1 } : {}}
      transition={{
        duration: 0.6,
        delay,
        ease: [0.25, 0.46, 0.45, 0.94],
      }}
      className={className}
    >
      <ClayCard>{children}</ClayCard>
    </motion.div>
  );
}

// ============================================================================
// CLAYMORPHISM CARD
// ============================================================================

function ClayCard({ children, className }: { children: React.ReactNode; className?: string }) {
  return (
    <motion.div
      className={cn(
        'relative h-full p-6 rounded-3xl overflow-hidden',
        // Claymorphism base
        'bg-[var(--color-surface-elevated)]',
        // Soft outer shadow for depth
        'shadow-[0_8px_32px_rgba(0,0,0,0.08),0_2px_8px_rgba(0,0,0,0.04)]',
        'dark:shadow-[0_8px_32px_rgba(0,0,0,0.3),0_2px_8px_rgba(0,0,0,0.2)]',
        // Inner highlight for 3D effect
        'before:absolute before:inset-0 before:rounded-3xl',
        'before:bg-gradient-to-br before:from-white/50 before:via-transparent before:to-transparent',
        'before:dark:from-white/10 before:dark:via-transparent before:dark:to-transparent',
        'before:pointer-events-none',
        // Border for definition
        'border border-[var(--color-border)]/50',
        // Hover effects
        'transition-all duration-300',
        'hover:shadow-[0_12px_48px_rgba(var(--color-primary-rgb,220,38,38),0.12),0_4px_16px_rgba(0,0,0,0.06)]',
        'hover:dark:shadow-[0_12px_48px_rgba(var(--color-primary-rgb,220,38,38),0.2),0_4px_16px_rgba(0,0,0,0.3)]',
        'hover:border-[var(--color-primary)]/30',
        'hover:-translate-y-1',
        className
      )}
      whileHover={{ scale: 1.01 }}
      transition={{ duration: 0.2 }}
    >
      <div className="relative z-10">{children}</div>
    </motion.div>
  );
}

// ============================================================================
// HERO SECTION
// ============================================================================

function HeroSection() {
  return (
    <motion.div
      initial={{ opacity: 0, y: 30 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.8, ease: [0.25, 0.46, 0.45, 0.94] }}
      className="text-center"
    >
      {/* Logo with glow effect */}
      <motion.div
        className="relative inline-block mb-8"
        whileHover={{ scale: 1.05 }}
        transition={{ duration: 0.3 }}
      >
        {/* Glow behind logo */}
        <motion.div
          className="absolute inset-0 rounded-3xl blur-2xl"
          style={{
            background: 'linear-gradient(135deg, var(--color-primary), var(--color-accent))',
            opacity: 0.3,
          }}
          animate={{
            opacity: [0.2, 0.4, 0.2],
            scale: [1, 1.1, 1],
          }}
          transition={{
            duration: 4,
            repeat: Infinity,
            ease: 'easeInOut',
          }}
        />
        <div className="relative w-28 h-28 rounded-3xl overflow-hidden shadow-2xl border-2 border-white/20">
          <img src="/logo.png" alt="ChemVault Logo" className="w-full h-full object-contain" />
        </div>
      </motion.div>

      {/* Title with gradient */}
      <motion.h1
        className="text-5xl md:text-6xl font-bold mb-6 font-display"
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.2, duration: 0.6 }}
      >
        <span className="text-[var(--color-text-primary)]">
          About <span className="font-extrabold text-[var(--color-primary)]">Chem</span>Vault
        </span>
      </motion.h1>

      {/* Subtitle */}
      <motion.p
        className="text-lg md:text-xl text-[var(--color-text-secondary)] max-w-2xl mx-auto leading-relaxed"
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.3, duration: 0.6 }}
      >
        A comprehensive web-based chemical structure validation and quality assessment platform
        for cheminformatics workflows, drug discovery, and ML dataset curation.
      </motion.p>

      {/* Decorative line */}
      <motion.div
        className="mt-8 mx-auto w-24 h-1 rounded-full bg-gradient-to-r from-[var(--color-primary)] to-[var(--color-accent)]"
        initial={{ scaleX: 0 }}
        animate={{ scaleX: 1 }}
        transition={{ delay: 0.5, duration: 0.8, ease: 'easeOut' }}
      />
    </motion.div>
  );
}

// ============================================================================
// SECTION HEADER
// ============================================================================

function SectionHeader({ icon, title, isChemVault }: { icon: React.ReactNode; title: string; isChemVault?: boolean }) {
  return (
    <div className="flex items-center gap-3 mb-5">
      <motion.div
        className={cn(
          'p-2.5 rounded-2xl',
          'bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/10',
          'text-[var(--color-primary)]',
          'shadow-inner'
        )}
        whileHover={{ rotate: [0, -10, 10, 0] }}
        transition={{ duration: 0.5 }}
      >
        {icon}
      </motion.div>
      <h2 className="text-xl font-bold text-[var(--color-text-primary)] font-display">
        {isChemVault ? (
          <>What is <span className="font-extrabold text-[var(--color-primary)]">Chem</span>Vault?</>
        ) : (
          title
        )}
      </h2>
    </div>
  );
}

// ============================================================================
// WHAT IS CHEMVAULT
// ============================================================================

function WhatIsChemVault() {
  const features = [
    { icon: <Shield className="w-4 h-4" />, title: '15+ Validation Checks', desc: 'Comprehensive analysis' },
    { icon: <Zap className="w-4 h-4" />, title: 'ChEMBL Standardization', desc: 'Trusted pipeline' },
    { icon: <BarChart3 className="w-4 h-4" />, title: 'ML-Readiness Scoring', desc: 'Dataset quality' },
    { icon: <Sparkles className="w-4 h-4" />, title: 'Batch Processing', desc: 'Millions of molecules' },
  ];

  return (
    <>
      <SectionHeader icon={<TestTube className="w-5 h-5" />} title="What is ChemVault?" isChemVault />
      <div className="space-y-4 text-[var(--color-text-secondary)]">
        <p className="leading-relaxed">
          ChemVault is an open-source platform designed to validate, standardize, and assess
          the quality of chemical structures. It provides researchers and data scientists with
          powerful tools to ensure their molecular data is accurate and ready for downstream
          applications.
        </p>
        <p className="leading-relaxed">
          Whether you're preparing datasets for machine learning, curating compound libraries,
          or validating structures for publication, ChemVault offers a comprehensive suite of
          validation checks, standardization pipelines, and quality scoring metrics.
        </p>

        {/* Feature grid with staggered animation */}
        <div className="grid grid-cols-2 gap-3 pt-4">
          {features.map((feature, i) => (
            <motion.div
              key={feature.title}
              initial={{ opacity: 0, x: -20 }}
              animate={{ opacity: 1, x: 0 }}
              transition={{ delay: 0.4 + i * 0.1, duration: 0.4 }}
              className={cn(
                'p-4 rounded-2xl',
                'bg-[var(--color-surface-sunken)]',
                'border border-[var(--color-border)]/30',
                'hover:border-[var(--color-primary)]/30',
                'hover:bg-[var(--color-primary)]/5',
                'transition-all duration-300',
                'group'
              )}
            >
              <div className="flex items-center gap-2 mb-1">
                <span className="text-[var(--color-primary)] group-hover:scale-110 transition-transform">
                  {feature.icon}
                </span>
                <span className="font-semibold text-[var(--color-text-primary)] text-sm">
                  {feature.title}
                </span>
              </div>
              <p className="text-xs text-[var(--color-text-muted)]">{feature.desc}</p>
            </motion.div>
          ))}
        </div>
      </div>
    </>
  );
}

// ============================================================================
// RESEARCH GROUP
// ============================================================================

function ResearchGroup() {
  return (
    <div className="h-full flex flex-col">
      <SectionHeader icon={<Building2 className="w-5 h-5" />} title="Research Group" />

      {/* Group Logo */}
      <motion.a
        href="http://cheminf.uni-jena.de/"
        target="_blank"
        rel="noopener noreferrer"
        className="block mb-4 group"
        whileHover={{ scale: 1.02 }}
      >
        <div className="relative overflow-hidden rounded-2xl bg-white dark:bg-white/10 p-4">
          <img
            src="/cheminf-logo.png"
            alt="Cheminformatics and Computational Metabolomics"
            className="h-20 mx-auto object-contain group-hover:scale-105 transition-transform duration-300"
          />
        </div>
      </motion.a>

      <div className="text-center mb-4">
        <h3 className="font-semibold text-[var(--color-text-primary)]">
          Cheminformatics and Computational Metabolomics
        </h3>
        <p className="text-sm text-[var(--color-text-muted)]">Prof. Dr. Christoph Steinbeck</p>
      </div>

      <p className="text-sm text-[var(--color-text-secondary)] mb-4 flex-1">
        Research focus on chemical structure annotation, deep learning for chemical
        information mining, and development of open-source cheminformatics tools.
      </p>

      <div className="space-y-3 pt-4 border-t border-[var(--color-border)]/50">
        <div className="flex items-start gap-2 text-sm text-[var(--color-text-muted)]">
          <MapPin className="w-4 h-4 mt-0.5 flex-shrink-0 text-[var(--color-primary)]" />
          <span>
            Friedrich Schiller University Jena<br />
            Institute for Inorganic and Analytical Chemistry<br />
            Lessingstr 8, 07743 Jena, Germany
          </span>
        </div>

        <div className="flex flex-col gap-2">
          <ExternalLinkButton href="http://cheminf.uni-jena.de/" icon={<Globe className="w-4 h-4" />}>
            Group Website
          </ExternalLinkButton>
          <ExternalLinkButton href="https://www.uni-jena.de" icon={<Building2 className="w-4 h-4" />}>
            University Website
          </ExternalLinkButton>
        </div>
      </div>
    </div>
  );
}

// ============================================================================
// TECH STACK
// ============================================================================

function TechStack() {
  const categories = [
    {
      icon: <Layout className="w-5 h-5" />,
      title: 'Frontend',
      items: ['React 18', 'TypeScript', 'Vite', 'Tailwind CSS', 'Framer Motion'],
      color: 'from-blue-500/20 to-cyan-500/10',
    },
    {
      icon: <Server className="w-5 h-5" />,
      title: 'Backend',
      items: ['Python 3.11+', 'FastAPI', 'Celery', 'Redis'],
      color: 'from-amber-500/20 to-yellow-500/10',
    },
    {
      icon: <TestTube className="w-5 h-5" />,
      title: 'Chemistry',
      items: ['RDKit', 'RDKit.js', 'MolVS', 'ChEMBL Pipeline'],
      color: 'from-purple-500/20 to-pink-500/10',
    },
    {
      icon: <Database className="w-5 h-5" />,
      title: 'Infrastructure',
      items: ['PostgreSQL', 'Docker', 'Nginx', 'Prometheus'],
      color: 'from-orange-500/20 to-amber-500/10',
    },
  ];

  return (
    <>
      <SectionHeader icon={<Code2 className="w-5 h-5" />} title="Technology Stack" />
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        {categories.map((cat, i) => (
          <motion.div
            key={cat.title}
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.3 + i * 0.1, duration: 0.4 }}
            className={cn(
              'p-4 rounded-2xl',
              'bg-gradient-to-br',
              cat.color,
              'border border-[var(--color-border)]/20'
            )}
          >
            <div className="flex items-center gap-2 mb-3">
              <span className="text-[var(--color-primary)]">{cat.icon}</span>
              <span className="font-semibold text-sm text-[var(--color-text-primary)]">{cat.title}</span>
            </div>
            <ul className="space-y-1">
              {cat.items.map((item) => (
                <li key={item} className="text-xs text-[var(--color-text-muted)]">
                  {item}
                </li>
              ))}
            </ul>
          </motion.div>
        ))}
      </div>
    </>
  );
}

// ============================================================================
// SOURCE CODE
// ============================================================================

function SourceCode() {
  return (
    <>
      <SectionHeader icon={<Github className="w-5 h-5" />} title="Source Code" />
      <p className="text-sm text-[var(--color-text-secondary)] mb-4">
        ChemVault is open-source software. Contributions, bug reports, and feature
        requests are welcome!
      </p>
      <motion.a
        href="https://github.com/Kohulan/ChemVault"
        target="_blank"
        rel="noopener noreferrer"
        className={cn(
          'flex items-center justify-center gap-2 w-full py-3 px-4 rounded-2xl',
          'bg-gradient-to-r from-[var(--color-surface-sunken)] to-[var(--color-surface-sunken)]',
          'border border-[var(--color-border)]/50',
          'text-[var(--color-text-primary)] font-medium',
          'hover:from-[var(--color-primary)]/10 hover:to-[var(--color-accent)]/10',
          'hover:border-[var(--color-primary)]/30',
          'transition-all duration-300'
        )}
        whileHover={{ scale: 1.02 }}
        whileTap={{ scale: 0.98 }}
      >
        <Github className="w-5 h-5" />
        View on GitHub
        <ExternalLink className="w-4 h-4 opacity-50" />
      </motion.a>
    </>
  );
}

// ============================================================================
// CONTACT
// ============================================================================

function Contact() {
  return (
    <>
      <SectionHeader icon={<Mail className="w-5 h-5" />} title="Contact" />
      <p className="text-sm text-[var(--color-text-secondary)] mb-4">
        Have questions, suggestions, or want to contribute? We'd love to hear from you!
      </p>
      <motion.a
        href="mailto:kohulan.rajan@uni-jena.de"
        className={cn(
          'flex items-center gap-3 p-3 rounded-2xl',
          'bg-[var(--color-surface-sunken)]',
          'border border-[var(--color-border)]/30',
          'hover:border-[var(--color-primary)]/30',
          'hover:bg-[var(--color-primary)]/5',
          'transition-all duration-300',
          'group'
        )}
        whileHover={{ x: 4 }}
      >
        <div className="p-2 rounded-xl bg-[var(--color-primary)]/10 text-[var(--color-primary)]">
          <Mail className="w-4 h-4" />
        </div>
        <span className="text-sm text-[var(--color-text-secondary)] group-hover:text-[var(--color-primary)] transition-colors">
          kohulan.rajan@uni-jena.de
        </span>
      </motion.a>
    </>
  );
}

// ============================================================================
// QUICK STATS
// ============================================================================

function QuickStats() {
  const stats = [
    { label: 'Validation Checks', value: '15+' },
    { label: 'Export Formats', value: '5' },
    { label: 'Open Source', value: 'MIT' },
  ];

  return (
    <>
      <SectionHeader icon={<Sparkles className="w-5 h-5" />} title="At a Glance" />
      <div className="grid grid-cols-3 gap-3">
        {stats.map((stat, i) => (
          <motion.div
            key={stat.label}
            initial={{ opacity: 0, scale: 0.8 }}
            animate={{ opacity: 1, scale: 1 }}
            transition={{ delay: 0.5 + i * 0.1, duration: 0.4 }}
            className="text-center p-3 rounded-2xl bg-[var(--color-surface-sunken)]"
          >
            <div className="text-2xl font-bold text-[var(--color-primary)]">{stat.value}</div>
            <div className="text-xs text-[var(--color-text-muted)]">{stat.label}</div>
          </motion.div>
        ))}
      </div>
    </>
  );
}

// ============================================================================
// ACKNOWLEDGMENTS
// ============================================================================

function Acknowledgments() {
  const acknowledgments = [
    { name: 'RDKit', description: 'Cheminformatics toolkit', href: 'https://www.rdkit.org/' },
    { name: 'ChEMBL', description: 'Bioactivity database', href: 'https://www.ebi.ac.uk/chembl/' },
    { name: 'PubChem', description: 'Chemical database', href: 'https://pubchem.ncbi.nlm.nih.gov/' },
    { name: 'COCONUT', description: 'Natural products DB', href: 'https://coconut.naturalproducts.net/' },
  ];

  return (
    <>
      <SectionHeader icon={<Heart className="w-5 h-5 text-red-500" />} title="Acknowledgments" />
      <p className="text-[var(--color-text-secondary)] mb-6">
        ChemVault is built upon the shoulders of giants. We gratefully acknowledge the
        following open-source projects and communities:
      </p>
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
        {acknowledgments.map((ack, i) => (
          <motion.a
            key={ack.name}
            href={ack.href}
            target="_blank"
            rel="noopener noreferrer"
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.6 + i * 0.1, duration: 0.4 }}
            className={cn(
              'p-4 rounded-2xl transition-all duration-300',
              'bg-[var(--color-surface-sunken)]',
              'border border-[var(--color-border)]/30',
              'hover:border-[var(--color-primary)]/30',
              'hover:bg-[var(--color-primary)]/5',
              'hover:-translate-y-1',
              'group'
            )}
            whileHover={{ scale: 1.02 }}
          >
            <div className="font-semibold text-sm text-[var(--color-text-primary)] group-hover:text-[var(--color-primary)] flex items-center gap-1 mb-1">
              {ack.name}
              <ExternalLink className="w-3 h-3 opacity-0 group-hover:opacity-100 transition-opacity" />
            </div>
            <div className="text-xs text-[var(--color-text-muted)]">{ack.description}</div>
          </motion.a>
        ))}
      </div>
      <p className="text-sm text-[var(--color-text-muted)] pt-4 border-t border-[var(--color-border)]/50">
        Special thanks to the open-source cheminformatics community for their continuous
        contributions to making chemical data more accessible and usable.
      </p>
    </>
  );
}

// ============================================================================
// EXTERNAL LINK BUTTON
// ============================================================================

interface ExternalLinkButtonProps {
  href: string;
  icon: React.ReactNode;
  children: React.ReactNode;
}

function ExternalLinkButton({ href, icon, children }: ExternalLinkButtonProps) {
  return (
    <motion.a
      href={href}
      target="_blank"
      rel="noopener noreferrer"
      className={cn(
        'inline-flex items-center gap-2 text-sm',
        'text-[var(--color-primary)] hover:text-[var(--color-accent)]',
        'transition-colors'
      )}
      whileHover={{ x: 4 }}
    >
      {icon}
      {children}
      <ExternalLink className="w-3 h-3" />
    </motion.a>
  );
}

export default AboutPage;
