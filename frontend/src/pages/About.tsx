import { useRef, useState } from 'react';
import { Link } from 'react-router-dom';
import { motion, useInView, useReducedMotion } from 'framer-motion';
import {
  Mail,
  Globe,
  Github,
  Coffee,
  ArrowRight,
  ExternalLink,
  MapPin,
  User,
} from 'lucide-react';
import { cn } from '../lib/utils';

// Static ambient glow. No floating molecules, no perpetual animation:
// the About page is read, not watched.
const glowOrbConfigs = [
  { x: '20%', y: '30%', size: 300, color: 'var(--color-primary)', opacity: 0.08, blur: 100 },
  { x: '70%', y: '50%', size: 250, color: 'var(--color-accent)', opacity: 0.06, blur: 80 },
];

/**
 * About page. One claim per section, every number stated once,
 * methods traceable to their citations.
 */
export function AboutPage() {
  return (
    <div className="relative min-h-screen overflow-hidden">
      {/* Ambient background */}
      <div className="fixed inset-0 pointer-events-none overflow-hidden">
        <div
          className="absolute inset-0"
          style={{
            background: `
              radial-gradient(ellipse 80% 50% at 20% 40%, rgba(var(--color-primary-rgb, 220, 38, 38), 0.08) 0%, transparent 50%),
              radial-gradient(ellipse 60% 40% at 80% 60%, rgba(var(--color-accent-rgb, 234, 88, 12), 0.06) 0%, transparent 50%),
              radial-gradient(ellipse 50% 30% at 50% 90%, rgba(var(--color-primary-rgb, 220, 38, 38), 0.05) 0%, transparent 50%)
            `,
          }}
        />
        {glowOrbConfigs.map((orb, i) => (
          <div
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
          />
        ))}
      </div>

      <div className="relative max-w-5xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="mb-12">
          <HeroSection />
        </div>

        {/* Alternating chrome: card, open, card, open, warm card. Six identical
            containers read as monotony; the rhythm is the design. */}
        <div className="space-y-6">
          <AnimatedCard>
            <WhatIsChemAudit />
          </AnimatedCard>

          <AnimatedCard variant="open">
            <Capabilities />
          </AnimatedCard>

          <AnimatedCard>
            <ScientificReferences />
          </AnimatedCard>

          <AnimatedCard variant="open">
            <BuiltInTheOpen />
          </AnimatedCard>

          <AnimatedCard tone="warm">
            <BuiltInJena />
          </AnimatedCard>
        </div>

        <LicenseFooter />
      </div>
    </div>
  );
}

// ============================================================================
// ANIMATED CARD WRAPPER
// ============================================================================

interface AnimatedCardProps {
  children: React.ReactNode;
  className?: string;
  /** 'card' wraps in claymorphism chrome; 'open' breathes directly on the page surface. */
  variant?: 'card' | 'open';
  /** 'warm' washes the card with the brand gradient tint. */
  tone?: 'default' | 'warm';
}

function AnimatedCard({ children, className, variant = 'card', tone = 'default' }: AnimatedCardProps) {
  const ref = useRef<HTMLDivElement>(null);
  const isInView = useInView(ref, { once: true, margin: '-50px' });

  return (
    <motion.div
      ref={ref}
      initial={{ opacity: 0, y: 24 }}
      animate={isInView ? { opacity: 1, y: 0 } : {}}
      transition={{ duration: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
      className={className}
    >
      {variant === 'open' ? (
        <section className="px-1 sm:px-2 py-10 sm:py-14">{children}</section>
      ) : (
        <ClayCard tone={tone}>{children}</ClayCard>
      )}
    </motion.div>
  );
}

// ============================================================================
// CLAYMORPHISM CARD
// Static surface: these sections are read, not clicked, so no hover lift.
// ============================================================================

function ClayCard({
  children,
  className,
  tone = 'default',
}: {
  children: React.ReactNode;
  className?: string;
  tone?: 'default' | 'warm';
}) {
  return (
    <div
      style={tone === 'warm' ? { backgroundImage: 'var(--gradient-primary-subtle)' } : undefined}
      className={cn(
        'relative h-full p-6 sm:p-8 rounded-3xl overflow-hidden',
        'bg-[var(--color-surface-elevated)]',
        'shadow-[0_8px_32px_rgba(0,0,0,0.08),0_2px_8px_rgba(0,0,0,0.04)]',
        'dark:shadow-[0_8px_32px_rgba(0,0,0,0.3),0_2px_8px_rgba(0,0,0,0.2)]',
        'before:absolute before:inset-0 before:rounded-3xl',
        'before:bg-gradient-to-br before:from-white/50 before:via-transparent before:to-transparent',
        'before:dark:from-white/10 before:dark:via-transparent before:dark:to-transparent',
        'before:pointer-events-none',
        'border border-[var(--color-border)]/50',
        className
      )}
    >
      <div className="relative z-10">{children}</div>
    </div>
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
      <div className="relative inline-block mb-8">
        {/* Static glow behind logo */}
        <div
          className="absolute inset-0 rounded-3xl blur-2xl"
          style={{
            background: 'linear-gradient(135deg, var(--color-primary), var(--color-accent))',
            opacity: 0.25,
          }}
        />
        <div className="relative w-28 h-28 rounded-3xl overflow-hidden shadow-2xl border-2 border-white/20">
          <img src="/logo-512.png" alt="ChemAudit Logo" className="w-full h-full object-contain" />
        </div>
      </div>

      <h1 className="text-5xl md:text-7xl font-bold mb-6 font-display tracking-tight text-[var(--color-text-primary)]">
        About <span className="font-extrabold text-[var(--color-primary)]">Chem</span>Audit
      </h1>

      <p className="text-lg md:text-xl text-[var(--color-text-secondary)] max-w-2xl mx-auto leading-relaxed">
        Validation, standardization, and quality scoring for chemical structures:
        built for cheminformatics workflows, drug discovery, and ML dataset curation.
      </p>

      <ReactionScheme />
    </motion.div>
  );
}

// ============================================================================
// REACTION SCHEME
// The product stated in the audience's native notation: reactant, conditions
// over the arrow, product. The arrow draws itself once on load.
// ============================================================================

function ReactionScheme() {
  const reduceMotion = useReducedMotion();

  return (
    <div className="mt-12 flex flex-col sm:flex-row items-center justify-center gap-4 sm:gap-6">
      <Link
        to={`/?smiles=${encodeURIComponent('CC(=O)Oc1ccccc1C(=O)O')}`}
        title="Run this structure through the validator"
        className={cn(
          'font-mono text-xs sm:text-sm px-3 py-1.5 rounded-lg',
          'bg-[var(--color-surface-elevated)] border border-[var(--color-border)]',
          'text-[var(--color-text-secondary)] shadow-sm',
          'hover:text-[var(--color-primary)] hover:border-[var(--color-primary)]/40 transition-colors'
        )}
      >
        CC(=O)Oc1ccccc1C(=O)O
      </Link>

      <div className="relative flex flex-col items-center rotate-90 sm:rotate-0 my-2 sm:my-0" aria-hidden="true">
        <span className="-rotate-90 sm:rotate-0 font-display text-sm font-semibold text-[var(--color-primary)] mb-1">
          ChemAudit
        </span>
        <svg width="120" height="12" viewBox="0 0 120 12" fill="none" className="overflow-visible">
          <motion.path
            d="M2 6 H 112 M 106 1.5 L 114 6 L 106 10.5"
            stroke="var(--color-text-secondary)"
            strokeWidth="1.5"
            strokeLinecap="round"
            strokeLinejoin="round"
            initial={reduceMotion ? false : { pathLength: 0 }}
            animate={{ pathLength: 1 }}
            transition={{ duration: 0.9, delay: 0.5, ease: [0.25, 1, 0.5, 1] }}
          />
        </svg>
      </div>

      <span
        className={cn(
          'text-xs sm:text-sm font-medium px-3 py-1.5 rounded-lg',
          'bg-[rgba(251,191,36,0.12)] dark:bg-[rgba(251,191,36,0.15)]',
          'text-[#b45309] dark:text-[#fcd34d]',
          'border border-[rgba(251,191,36,0.3)]'
        )}
      >
        valid &middot; standardized &middot; scored
      </span>
    </div>
  );
}

// ============================================================================
// SECTION HEADER
// ============================================================================

// Sections are numbered like a paper: this audience navigates by section
// number. The mono index replaces generic icon tiles.
function SectionHeader({ index, title }: { index: string; title: React.ReactNode }) {
  return (
    <div className="flex items-center gap-3 mb-6">
      <div
        className={cn(
          'w-10 h-10 rounded-2xl flex items-center justify-center',
          'bg-gradient-to-br from-[var(--color-primary)]/20 to-[var(--color-accent)]/10',
          'shadow-inner'
        )}
        aria-hidden="true"
      >
        <span className="font-mono text-sm font-semibold text-[var(--color-primary)]">{index}</span>
      </div>
      <h2 className="text-xl sm:text-2xl font-bold text-[var(--color-text-primary)] font-display tracking-tight">
        {title}
      </h2>
    </div>
  );
}

// ============================================================================
// WHAT IS CHEMAUDIT
// The numbers live here, once. No stat tiles restating them below.
// ============================================================================

function WhatIsChemAudit() {
  return (
    <>
      <SectionHeader
        index="01"
        title={
          <>
            What is <span className="font-extrabold text-[var(--color-primary)]">Chem</span>Audit?
          </>
        }
      />
      <div className="grid grid-cols-1 lg:grid-cols-5 gap-8 lg:gap-12">
        <div className="lg:col-span-3 space-y-4 text-[var(--color-text-secondary)]">
          <p className="leading-relaxed">
            ChemAudit is an open-source platform that checks chemical structures before they reach
            your model, your library, or your paper. Every molecule passes through{' '}
            <strong className="font-semibold text-[var(--color-text-primary)]">27 validation checks</strong>,
            the ChEMBL standardization pipeline, and screening against more than{' '}
            <strong className="font-semibold text-[var(--color-text-primary)]">1,500 structural-alert patterns</strong>.
          </p>
          <p className="leading-relaxed">
            It scales from{' '}
            <Link
              to={`/?smiles=${encodeURIComponent('CC(=O)Oc1ccccc1C(=O)O')}`}
              title="Try it: this link validates aspirin"
              className="underline decoration-dotted decoration-[var(--color-text-muted)] underline-offset-2 hover:text-[var(--color-primary)] hover:decoration-[var(--color-primary)] transition-colors"
            >
              a single pasted SMILES
            </Link>{' '}
            to batch jobs of millions of molecules, computing{' '}
            <strong className="font-semibold text-[var(--color-text-primary)]">451 descriptors</strong> and
            seven fingerprint types along the way, with results exportable in five formats. The
            verdicts are traceable: every score links back to the published method that produced it.
          </p>
        </div>

        <div className="lg:col-span-2">
          <h3 className="text-xs font-semibold uppercase tracking-widest text-[var(--color-text-muted)] pb-2 mb-4 border-b border-[var(--color-border)]">
            Start here
          </h3>
          <ul className="space-y-4">
            {START_HERE.map((entry) => (
              <li key={entry.title}>
                {entry.external ? (
                  <a href={entry.href} target="_blank" rel="noopener noreferrer" className="group block">
                    <StartHereEntry entry={entry} />
                  </a>
                ) : (
                  <Link to={entry.href} className="group block">
                    <StartHereEntry entry={entry} />
                  </Link>
                )}
              </li>
            ))}
          </ul>
        </div>
      </div>
    </>
  );
}

interface StartHereItem {
  title: string;
  desc: string;
  href: string;
  external?: boolean;
}

const START_HERE: StartHereItem[] = [
  {
    title: 'Validate a molecule',
    desc: 'Paste a SMILES, get a verdict in seconds.',
    href: '/',
  },
  {
    title: 'Process a batch',
    desc: 'Upload an SDF or CSV and watch progress live.',
    href: '/batch',
  },
  {
    title: 'Audit a dataset',
    desc: 'Health-check an entire collection before training on it.',
    href: '/dataset-audit',
  },
  {
    title: 'Read the docs',
    desc: 'Guides, API reference, and deployment.',
    href: 'https://kohulan.github.io/ChemAudit/',
    external: true,
  },
];

function StartHereEntry({ entry }: { entry: StartHereItem }) {
  return (
    <>
      <span className="flex items-center gap-1.5 text-sm font-semibold text-[var(--color-text-primary)] group-hover:text-[var(--color-primary)] transition-colors">
        {entry.title}
        <ArrowRight
          className="w-3.5 h-3.5 opacity-0 -translate-x-1 group-hover:opacity-100 group-hover:translate-x-0 transition-all duration-200 ease-out"
          aria-hidden="true"
        />
      </span>
      <span className="block text-sm text-[var(--color-text-secondary)] leading-relaxed">
        {entry.desc}
      </span>
    </>
  );
}

// ============================================================================
// CAPABILITIES
// One inventory for the whole platform. Typographic list, not a card wall:
// the references section below documents the methods behind each line.
// ============================================================================

interface Capability {
  title: string;
  desc: string;
}

const SCORING_CAPABILITIES: Capability[] = [
  {
    title: 'Drug-likeness',
    desc: "Lipinski, QED, Veber, Rule of Three, Ghose, and Muegge filters in one pass.",
  },
  {
    title: 'ADMET predictions',
    desc: 'Solubility (ESOL), CNS MPO, synthetic accessibility, and the Pfizer 3/75, GSK 4/400, and Golden Triangle rules.',
  },
  {
    title: 'Safety filters',
    desc: 'Structural alerts across PAINS (480 patterns), Brenk, NIH, ZINC, and seven ChEMBL filter sets.',
  },
  {
    title: 'Aggregator likelihood',
    desc: 'Flags colloidal aggregation risk before it produces false positives in your screen.',
  },
  {
    title: 'ML-readiness',
    desc: 'Descriptor and fingerprint coverage scored for QSAR/QSPR suitability, with dataset-level quality metrics.',
  },
  {
    title: 'Natural-product likeness',
    desc: 'Positions each molecule on the synthetic-to-natural axis for library triage.',
  },
];

const PIPELINE_CAPABILITIES: Capability[] = [
  {
    title: 'QSAR-ready pipeline',
    desc: 'Configurable standardization with salt stripping, neutralization, tautomer canonicalization, and duplicate removal.',
  },
  {
    title: 'Structure filter',
    desc: 'Sequential property and SMARTS funnels that narrow generative-chemistry output to viable candidates.',
  },
  {
    title: 'Dataset intelligence',
    desc: 'Health scoring, contradictory-label detection, dataset diffing, and curation reports.',
  },
  {
    title: 'Identifier resolution',
    desc: 'SMILES, InChI, InChIKey, CAS, ChEBI, UNII, PubChem CID, ChEMBL ID, and compound names, cross-linked through UniChem.',
  },
  {
    title: 'Batch analytics',
    desc: 'Butina clustering, chemotype taxonomy, t-SNE chemical space, scaffold analysis, outlier detection, and registration hashing.',
  },
  {
    title: 'Diagnostics & provenance',
    desc: 'SMILES diagnostics, InChI layer diffs, round-trip validation, and a full audit trail of every transformation applied.',
  },
];

function CapabilityGroup({ label, items }: { label: string; items: Capability[] }) {
  return (
    <div>
      <h3 className="text-xs font-semibold uppercase tracking-widest text-[var(--color-text-muted)] pb-2 mb-4 border-b border-[var(--color-border)]">
        {label}
      </h3>
      <dl className="space-y-4">
        {items.map((item, i) => (
          <div key={item.title} className="grid grid-cols-[2rem_1fr]">
            <dt className="contents">
              <span
                className="font-mono text-xs text-[var(--color-accent-dark)] dark:text-[var(--color-accent)] pt-0.5"
                aria-hidden="true"
              >
                {String(i + 1).padStart(2, '0')}
              </span>
              <span className="font-semibold text-sm text-[var(--color-text-primary)]">
                {item.title}
              </span>
            </dt>
            <dd className="col-start-2 text-sm text-[var(--color-text-secondary)] leading-relaxed">
              {item.desc}
            </dd>
          </div>
        ))}
      </dl>
    </div>
  );
}

function Capabilities() {
  return (
    <>
      <SectionHeader index="02" title="What's inside" />
      <p className="text-[var(--color-text-secondary)] mb-8 max-w-prose">
        Twelve modules: six for scoring and screening, six for dataset preparation and analysis.
        The scoring side implements industry rules from Pfizer, GSK, Lilly, and BMS alongside the
        academic literature; each method is cited in the references below.
      </p>
      <div className="grid grid-cols-1 md:grid-cols-2 gap-x-12 gap-y-8">
        <CapabilityGroup label="Score & screen" items={SCORING_CAPABILITIES} />
        <CapabilityGroup label="Prepare & analyze" items={PIPELINE_CAPABILITIES} />
      </div>
    </>
  );
}

// ============================================================================
// SCIENTIFIC REFERENCES
// ============================================================================

interface ReferenceCategory {
  title: string;
  references: {
    method: string;
    citation: string;
    doi?: string;
  }[];
}

function ScientificReferences() {
  const [expandedCategory, setExpandedCategory] = useState<string | null>(null);

  const categories: ReferenceCategory[] = [
    {
      title: 'Drug-likeness Rules',
      references: [
        {
          method: "Lipinski's Rule of Five",
          citation: 'Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Adv Drug Deliv Rev. 2001;46(1-3):3-26.',
          doi: '10.1016/S0169-409X(00)00129-0',
        },
        {
          method: 'QED (Quantitative Estimate of Drug-likeness)',
          citation: 'Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying the chemical beauty of drugs. Nat Chem. 2012;4(2):90-98.',
          doi: '10.1038/nchem.1243',
        },
        {
          method: 'Veber Rules',
          citation: 'Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD. Molecular properties that influence the oral bioavailability of drug candidates. J Med Chem. 2002;45(12):2615-2623.',
          doi: '10.1021/jm020017n',
        },
        {
          method: 'Rule of Three (Ro3)',
          citation: "Congreve M, Carr R, Murray C, Jhoti H. A 'rule of three' for fragment-based lead discovery? Drug Discov Today. 2003;8(19):876-877.",
          doi: '10.1016/S1359-6446(03)02831-9',
        },
        {
          method: 'Ghose Filter',
          citation: 'Ghose AK, Viswanadhan VN, Wendoloski JJ. A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. J Comb Chem. 1999;1(1):55-68.',
          doi: '10.1021/cc9800071',
        },
        {
          method: 'Muegge Filter',
          citation: 'Muegge I, Heald SL, Brittelli D. Simple selection criteria for drug-like chemical matter. J Med Chem. 2001;44(12):1841-1846.',
          doi: '10.1021/jm015507e',
        },
      ],
    },
    {
      title: 'ADMET Predictions',
      references: [
        {
          method: 'Synthetic Accessibility (SA) Score',
          citation: 'Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. J Cheminform. 2009;1:8.',
          doi: '10.1186/1758-2946-1-8',
        },
        {
          method: 'ESOL Solubility',
          citation: 'Delaney JS. ESOL: estimating aqueous solubility directly from molecular structure. J Chem Inf Comput Sci. 2004;44(3):1000-1005.',
          doi: '10.1021/ci034243x',
        },
        {
          method: 'Fsp3 (3D Complexity)',
          citation: 'Lovering F, Bikker J, Humblet C. Escape from flatland: increasing saturation as an approach to improving clinical success. J Med Chem. 2009;52(21):6752-6756.',
          doi: '10.1021/jm901241e',
        },
        {
          method: 'CNS MPO Score',
          citation: 'Wager TT, Hou X, Verhoest PR, Villalobos A. Moving beyond rules: the development of a central nervous system multiparameter optimization (CNS MPO) approach. ACS Chem Neurosci. 2010;1(6):435-449.',
          doi: '10.1021/cn100008c',
        },
        {
          method: 'Pfizer 3/75 Rule',
          citation: 'Hughes JD, Blagg J, Price DA, et al. Physiochemical drug properties associated with in vivo toxicological outcomes. Bioorg Med Chem Lett. 2008;18(17):4872-4875.',
          doi: '10.1016/j.bmcl.2008.07.071',
        },
        {
          method: 'GSK 4/400 Rule',
          citation: 'Gleeson MP. Generation of a set of simple, interpretable ADMET rules of thumb. J Med Chem. 2008;51(4):817-834.',
          doi: '10.1021/jm701122q',
        },
        {
          method: 'Golden Triangle',
          citation: 'Johnson TW, Dress KR, Edwards M. Using the Golden Triangle to optimize clearance and oral absorption. Bioorg Med Chem Lett. 2009;19(19):5560-5564.',
          doi: '10.1016/j.bmcl.2009.08.045',
        },
      ],
    },
    {
      title: 'Safety Filters & Structural Alerts',
      references: [
        {
          method: 'PAINS (Pan-Assay Interference Compounds)',
          citation: 'Baell JB, Holloway GA. New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays. J Med Chem. 2010;53(7):2719-2740.',
          doi: '10.1021/jm901137j',
        },
        {
          method: 'Brenk / Dundee NTD Screening Alerts',
          citation: 'Brenk R, Schipani A, James D, et al. Lessons learnt from assembling screening libraries for drug discovery for neglected diseases. ChemMedChem. 2008;3(3):435-444.',
          doi: '10.1002/cmdc.200700139',
        },
        {
          method: 'NIH MLPCN Exclusion Filters',
          citation: 'Jadhav A, Ferreira RS, Klumpp C, et al. Quantitative analyses of aggregation, autofluorescence, and reactivity artifacts in a screen for inhibitors of a thiol protease. J Med Chem. 2010;53(1):37-51.',
          doi: '10.1021/jm901070c',
        },
        {
          method: 'ZINC Druglike Filters',
          citation: 'Irwin JJ, Shoichet BK. ZINC — a free database of commercially available compounds for virtual screening. J Chem Inf Model. 2005;45(1):177-182.',
          doi: '10.1021/ci049714+',
        },
        {
          method: 'BMS HTS Desirability Filters',
          citation: 'Pearce BC, Sofia MJ, Good AC, Drexler DM, Stock DA. An empirical process for the design of high-throughput screening deck filters. J Chem Inf Model. 2006;46(3):1060-1068.',
          doi: '10.1021/ci050504m',
        },
        {
          method: 'Glaxo Hard Filters',
          citation: 'Hann M, Hudson B, Lewell X, Lifely R, Miller L, Ramsden N. Strategic pooling of compounds for high-throughput screening. J Chem Inf Comput Sci. 1999;39(5):897-902.',
          doi: '10.1021/ci990423o',
        },
        {
          method: 'Lilly MedChem Rules (LINT)',
          citation: 'Bruns RF, Watson IA. Rules for identifying potentially reactive or promiscuous compounds. J Med Chem. 2012;55(22):9763-9772.',
          doi: '10.1021/jm301008n',
        },
        {
          method: 'SureChEMBL Non-chemical Filters',
          citation: 'Papadatos G, Davies M, Dedber N, et al. SureChEMBL: a large-scale, chemically annotated patent document database. Nucleic Acids Res. 2016;44(D1):D1220-D1228.',
          doi: '10.1093/nar/gkv1253',
        },
        {
          method: 'Phantom PAINS: Context for Alerts',
          citation: 'Jasial S, Hu Y, Bajorath J. How frequently are pan-assay interference compounds active? Large-scale analysis of screening data reveals diverse activity profiles, low global hit frequency, and many consistently inactive compounds. J Med Chem. 2017;60(9):3879-3886.',
          doi: '10.1021/acs.jmedchem.7b00154',
        },
      ],
    },
    {
      title: 'Scoring & Analysis',
      references: [
        {
          method: 'NP-likeness Score',
          citation: 'Ertl P, Roggo S, Schuffenhauer A. Natural product-likeness score and its application for prioritization of compound libraries. J Chem Inf Model. 2008;48(1):68-74.',
          doi: '10.1021/ci700286x',
        },
        {
          method: 'Murcko Scaffold',
          citation: 'Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. J Med Chem. 1996;39(15):2887-2893.',
          doi: '10.1021/jm9602928',
        },
        {
          method: 'Aggregator Detection',
          citation: 'McGovern SL, Caselli E, Grigorieff N, Shoichet BK. A common mechanism underlying promiscuous inhibitors from virtual and high-throughput screening. J Med Chem. 2002;45(8):1712-1722.',
          doi: '10.1021/jm010533y',
        },
        {
          method: 'Ligand Efficiency (LE)',
          citation: 'Hopkins AL, Keserü GM, Leeson PD, Rees DC, Reynolds CH. The role of ligand efficiency metrics in drug discovery. Nat Rev Drug Discov. 2014;13(2):105-121.',
          doi: '10.1038/nrd4163',
        },
      ],
    },
    {
      title: 'Bioavailability & Permeation',
      references: [
        {
          method: 'SwissADME / Bioavailability Radar',
          citation: 'Daina A, Michielin O, Zoete V. SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. Sci Rep. 2017;7:42717.',
          doi: '10.1038/srep42717',
        },
        {
          method: 'BOILED-Egg (GI Absorption / BBB Permeation)',
          citation: 'Daina A, Zoete V. A BOILED-Egg to predict gastrointestinal absorption and brain penetration of small molecules. ChemMedChem. 2016;11(11):1117-1121.',
          doi: '10.1002/cmdc.201600182',
        },
        {
          method: 'Wildman-Crippen LogP (WLOGP)',
          citation: 'Wildman SA, Crippen GM. Prediction of physicochemical parameters by atomic contributions. J Chem Inf Comput Sci. 1999;39(5):868-873.',
          doi: '10.1021/ci990307l',
        },
      ],
    },
    {
      title: 'Clustering & Analytics',
      references: [
        {
          method: 'Butina Clustering',
          citation: 'Butina D. Unsupervised data base clustering based on Daylight\'s fingerprint and Tanimoto similarity: a fast and automated way to cluster small and large data sets. J Chem Inf Comput Sci. 1999;39(4):747-750.',
          doi: '10.1021/ci9803381',
        },
        {
          method: 'Morgan Fingerprints (ECFP)',
          citation: 'Rogers D, Hahn M. Extended-connectivity fingerprints. J Chem Inf Model. 2010;50(5):742-754.',
          doi: '10.1021/ci100050t',
        },
        {
          method: 'Tanimoto Similarity',
          citation: 'Bajusz D, Rácz A, Héberger K. Why is Tanimoto index an appropriate choice for fingerprint-based similarity calculations? J Cheminform. 2015;7:20.',
          doi: '10.1186/s13321-015-0069-3',
        },
        {
          method: 't-SNE Visualization',
          citation: 'van der Maaten L, Hinton G. Visualizing data using t-SNE. J Mach Learn Res. 2008;9:2579-2605.',
          doi: undefined,
        },
      ],
    },
    {
      title: 'Software & Pipelines',
      references: [
        {
          method: 'RDKit',
          citation: 'Landrum G. RDKit: Open-Source Cheminformatics Software.',
          doi: undefined,
        },
        {
          method: 'ChEMBL Structure Pipeline',
          citation: 'Bento AP, Hersey A, Félix E, et al. An open source chemical structure curation pipeline using RDKit. J Cheminform. 2020;12:51.',
          doi: '10.1186/s13321-020-00456-1',
        },
        {
          method: 'MolVS Standardizer',
          citation: 'Swain M. MolVS: Molecule Validation and Standardization.',
          doi: undefined,
        },
      ],
    },
  ];

  return (
    <>
      <SectionHeader index="03" title="Methods & Scientific References" />
      <p className="text-[var(--color-text-secondary)] mb-6 max-w-prose">
        Every algorithm and scoring function in ChemAudit comes from the published literature.
        These are the primary references.
      </p>

      <div className="space-y-3">
        {categories.map((category) => {
          const panelId = `refs-${category.title.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`;
          return (
          <div
            key={category.title}
            className="rounded-xl overflow-hidden border border-[var(--color-border)]/50 bg-[var(--color-surface-sunken)]/50"
          >
            <button
              onClick={() => setExpandedCategory(expandedCategory === category.title ? null : category.title)}
              aria-expanded={expandedCategory === category.title}
              aria-controls={panelId}
              className="w-full flex items-center justify-between p-4 text-left hover:bg-[var(--color-primary)]/5 transition-colors"
            >
              <div>
                <span className="font-semibold text-sm text-[var(--color-text-primary)]">
                  {category.title}
                </span>
                <span className="ml-2 text-xs text-[var(--color-text-muted)]">
                  ({category.references.length} references)
                </span>
              </div>
              <motion.div
                animate={{ rotate: expandedCategory === category.title ? 180 : 0 }}
                transition={{ duration: 0.2 }}
              >
                <svg
                  className="w-5 h-5 text-[var(--color-text-muted)]"
                  fill="none"
                  viewBox="0 0 24 24"
                  stroke="currentColor"
                  aria-hidden="true"
                >
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
                </svg>
              </motion.div>
            </button>

            {expandedCategory === category.title && (
              <motion.div
                id={panelId}
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ duration: 0.2 }}
                className="px-4 pb-4"
              >
                <div className="space-y-3">
                  {category.references.map((ref) => (
                    <div
                      key={ref.method}
                      className="p-3 rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)]/30"
                    >
                      <div className="flex items-start justify-between gap-2 mb-1">
                        <h3 className="font-medium text-sm text-[var(--color-text-primary)]">
                          {ref.method}
                        </h3>
                        {ref.doi && (
                          <a
                            href={`https://doi.org/${ref.doi}`}
                            target="_blank"
                            rel="noopener noreferrer"
                            className={cn(
                              'inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs',
                              'bg-[var(--color-primary)]/10 text-[var(--color-primary)]',
                              'hover:bg-[var(--color-primary)]/20 transition-colors'
                            )}
                          >
                            DOI
                            <ExternalLink className="w-2.5 h-2.5" aria-hidden="true" />
                          </a>
                        )}
                      </div>
                      <p className="text-xs text-[var(--color-text-secondary)] leading-relaxed">
                        {ref.citation}
                      </p>
                    </div>
                  ))}
                </div>
              </motion.div>
            )}
          </div>
          );
        })}
      </div>

      <p className="mt-5 text-xs text-[var(--color-text-muted)]">
        Throughout the app, hover the{' '}
        <span className="inline-flex items-center justify-center w-4 h-4 rounded-full border border-[var(--color-text-muted)]/40 bg-[var(--color-surface-sunken)] text-[9px] font-semibold mx-0.5">
          i
        </span>{' '}
        icons next to scores for the citation behind each method.
      </p>
    </>
  );
}

// ============================================================================
// BUILT IN THE OPEN (tech stack + acknowledgments)
// ============================================================================

const TECH_STACK: { label: string; items: string }[] = [
  { label: 'Frontend', items: 'React 18, TypeScript, Vite, Tailwind CSS, Framer Motion, Recharts, RDKit.js' },
  { label: 'Backend', items: 'Python 3.11+, FastAPI, Celery, Redis, Pandas, asyncpg' },
  { label: 'Chemistry', items: 'RDKit, MolVS, ChEMBL Structure Pipeline, openTSNE' },
  { label: 'Infrastructure', items: 'PostgreSQL, Docker, Nginx, Prometheus' },
];

const ACKNOWLEDGMENTS = [
  {
    name: 'RDKit',
    description: 'Open-source cheminformatics toolkit powering molecular operations',
    href: 'https://www.rdkit.org/',
    logo: 'https://www.rdkit.org/Images/logo.png',
  },
  {
    name: 'ChEMBL',
    description: 'Bioactivity database for drug discovery from EMBL-EBI',
    href: 'https://www.ebi.ac.uk/chembl/',
    logo: 'https://cfde-gene-pages.cloud/logos/chEMBL_logo.png',
  },
  {
    name: 'PubChem',
    description: "World's largest collection of freely accessible chemical information",
    href: 'https://pubchem.ncbi.nlm.nih.gov/',
    logo: 'https://upload.wikimedia.org/wikipedia/commons/thumb/b/b6/PubChem_logo.svg/1280px-PubChem_logo.svg.png',
  },
  {
    name: 'COCONUT',
    description: 'Collection of Open Natural Products database',
    href: 'https://coconut.naturalproducts.net/',
    logo: 'https://raw.githubusercontent.com/Steinbeck-Lab/coconut/main/public/img/logo.svg',
  },
  {
    name: 'ChEBI',
    description: 'Chemical Entities of Biological Interest ontology from EMBL-EBI',
    href: 'https://www.ebi.ac.uk/chebi/',
    logo: 'https://www.ebi.ac.uk/chebi/chebi_logo.svg',
  },
  {
    name: 'UniChem',
    description: 'Cross-reference mapping between chemistry databases from EMBL-EBI',
    href: 'https://www.ebi.ac.uk/unichem/',
    logo: 'https://www.ebi.ac.uk/unichem/_nuxt/img/unichem_logo.de3b448.png',
  },
  {
    name: 'Wikidata',
    description: 'Free and open knowledge base for structured chemical data',
    href: 'https://www.wikidata.org/',
    logo: 'https://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Wikidata-logo-en.svg/1920px-Wikidata-logo-en.svg.png',
  },
  {
    name: 'SureChEMBL',
    description: 'Patent chemistry database for compound-patent literature linkage',
    href: 'https://www.surechembl.org/',
    logo: 'https://www.surechembl.org/img/icons/chembl_logo_pink.png',
  },
];

function BuiltInTheOpen() {
  return (
    <>
      <SectionHeader index="04" title="Built in the open" />

      {/* Tech stack as quiet definition rows */}
      <div className="mb-8 divide-y divide-[var(--color-border)]/50">
        {TECH_STACK.map((row) => (
          <div key={row.label} className="flex flex-col sm:flex-row sm:items-baseline gap-1 sm:gap-4 py-2.5">
            <span className="w-32 shrink-0 text-sm font-semibold text-[var(--color-text-primary)]">
              {row.label}
            </span>
            <span className="text-sm text-[var(--color-text-secondary)]">{row.items}</span>
          </div>
        ))}
      </div>

      <p className="text-[var(--color-text-secondary)] mb-5 max-w-prose">
        ChemAudit stands on open scientific resources. We gratefully acknowledge:
      </p>

      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        {ACKNOWLEDGMENTS.map((ack) => (
          <a
            key={ack.name}
            href={ack.href}
            target="_blank"
            rel="noopener noreferrer"
            className={cn(
              'flex flex-col p-4 rounded-2xl',
              'bg-[var(--color-surface-sunken)]',
              'border border-[var(--color-border)]/30',
              'hover:border-[var(--color-primary)]/40',
              'transition-colors duration-200',
              'group'
            )}
          >
            <div className="relative w-full h-14 rounded-xl mb-3 flex items-center justify-center bg-white dark:bg-white/95 overflow-hidden">
              <img
                src={ack.logo}
                alt={`${ack.name} logo`}
                className="max-h-10 max-w-[85%] object-contain"
                onError={(e) => {
                  // Neutral text fallback when a remote logo fails to load
                  const target = e.target as HTMLImageElement;
                  const parent = target.parentElement as HTMLDivElement;
                  target.style.display = 'none';
                  parent.classList.remove('bg-white', 'dark:bg-white/95');
                  parent.classList.add('bg-chem-dark-800');
                  const span = document.createElement('span');
                  span.className = 'text-base font-bold text-white';
                  span.textContent = ack.name;
                  parent.appendChild(span);
                }}
              />
            </div>
            <div className="flex items-center gap-2 mb-1">
              <h3 className="font-semibold text-sm text-[var(--color-text-primary)] group-hover:text-[var(--color-primary)] transition-colors">
                {ack.name}
              </h3>
              <ExternalLink
                className="w-3 h-3 text-[var(--color-text-muted)] opacity-0 group-hover:opacity-100 transition-opacity"
                aria-hidden="true"
              />
            </div>
            <p className="text-xs text-[var(--color-text-muted)] leading-relaxed">{ack.description}</p>
          </a>
        ))}
      </div>
    </>
  );
}

// ============================================================================
// BUILT IN JENA (research group, developer, contact)
// ============================================================================

const CONTACT_LINKS = [
  { icon: <Github className="w-4 h-4" />, label: 'Source on GitHub', href: 'https://github.com/Kohulan/ChemAudit' },
  { icon: <Github className="w-4 h-4" />, label: 'Steinbeck Lab', href: 'https://github.com/Steinbeck-Lab' },
  { icon: <Mail className="w-4 h-4" />, label: 'kohulan.rajan@uni-jena.de', href: 'mailto:kohulan.rajan@uni-jena.de' },
  { icon: <Globe className="w-4 h-4" />, label: 'cheminf.uni-jena.de', href: 'http://cheminf.uni-jena.de/' },
];

function BuiltInJena() {
  const mapUrl = 'https://www.google.com/maps/place/Lessingstra%C3%9Fe+8,+07743+Jena,+Germany';

  return (
    <>
      <SectionHeader index="05" title="Built in Jena" />

      <div className="grid grid-cols-1 md:grid-cols-2 gap-x-12 gap-y-8 mb-7">
        {/* Research group */}
        <div>
          <a
            href="http://cheminf.uni-jena.de/"
            target="_blank"
            rel="noopener noreferrer"
            className="inline-block mb-4"
          >
            <span className="inline-flex items-center rounded-xl bg-white dark:bg-white/95 px-4 py-2.5 shadow-sm hover:shadow-md transition-shadow">
              <img
                src="/cheminf-logo.png"
                alt="Natural Products Cheminformatics research group"
                className="h-12 object-contain"
              />
            </span>
          </a>
          <p className="text-sm text-[var(--color-text-secondary)] max-w-prose leading-relaxed">
            The group works on chemical structure annotation, deep learning for chemical
            information mining, and open-source cheminformatics tools.
          </p>
        </div>

        {/* Developer + contact */}
        <div>
          <div className="flex items-center gap-4 mb-4">
            <div className="w-16 h-16 rounded-full overflow-hidden border-2 border-[var(--color-surface-elevated)] shadow-md bg-[var(--color-surface-sunken)] shrink-0">
              <img
                src="https://github.com/Kohulan.png"
                alt="Kohulan Rajan"
                className="w-full h-full object-cover"
              />
            </div>
            <div>
              <div className="flex items-center gap-2">
                <User className="w-4 h-4 text-[var(--color-primary)]" aria-hidden="true" />
                <h3 className="text-base font-bold text-[var(--color-text-primary)]">Kohulan Rajan</h3>
              </div>
              <p className="text-sm text-[var(--color-text-secondary)]">
                Senior Researcher: AI, deep learning, and cheminformatics
              </p>
              <div className="flex items-center gap-4 mt-1.5">
                <a
                  href="https://kohulanr.com"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="inline-flex items-center gap-1.5 text-xs text-[var(--color-text-muted)] hover:text-[var(--color-primary)] transition-colors"
                >
                  <Globe className="w-3.5 h-3.5" aria-hidden="true" />
                  kohulanr.com
                </a>
                <a
                  href="https://github.com/Kohulan"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="inline-flex items-center gap-1.5 text-xs text-[var(--color-text-muted)] hover:text-[var(--color-primary)] transition-colors"
                >
                  <Github className="w-3.5 h-3.5" aria-hidden="true" />
                  github.com/Kohulan
                </a>
              </div>
            </div>
          </div>

          <p className="text-sm text-[var(--color-text-secondary)] mb-3 max-w-prose">
            ChemAudit is open source. Contributions, bug reports, and feature requests are welcome.
          </p>

          <ul className="grid grid-cols-1 sm:grid-cols-2 gap-x-6 gap-y-2.5">
            {CONTACT_LINKS.map((link) => (
              <li key={link.label}>
                <a
                  href={link.href}
                  target={link.href.startsWith('mailto') ? undefined : '_blank'}
                  rel={link.href.startsWith('mailto') ? undefined : 'noopener noreferrer'}
                  className={cn(
                    'inline-flex items-center gap-2.5 text-sm',
                    'text-[var(--color-text-secondary)]',
                    'hover:text-[var(--color-primary)] transition-colors'
                  )}
                >
                  <span className="text-[var(--color-primary)]">{link.icon}</span>
                  {link.label}
                </a>
              </li>
            ))}
          </ul>
        </div>
      </div>

      {/* Map anchors the section's base at full width; the address rides on it.
          Top-left so the OpenStreetMap attribution stays unobstructed. */}
      <a
        href={mapUrl}
        target="_blank"
        rel="noopener noreferrer"
        className="relative block h-40 md:h-48 rounded-2xl overflow-hidden border border-[var(--color-border)]/50 hover:border-[var(--color-primary)]/30 transition-colors"
      >
        <iframe
          src="https://www.openstreetmap.org/export/embed.html?bbox=11.5825%2C50.9245%2C11.5955%2C50.9305&layer=mapnik&marker=50.9275%2C11.589"
          className="absolute inset-0 w-full h-full border-0 pointer-events-none"
          title="Location Map"
        />
        <span
          className={cn(
            'absolute top-3 left-3 z-10 inline-flex items-start gap-2',
            'rounded-lg px-3 py-2 text-xs leading-snug',
            'bg-[var(--color-surface-elevated)]/95 text-[var(--color-text-secondary)]',
            'border border-[var(--color-border)] shadow-md'
          )}
        >
          <MapPin className="w-3.5 h-3.5 mt-0.5 shrink-0 text-[var(--color-primary)]" aria-hidden="true" />
          <span>
            Friedrich Schiller University Jena
            <br />
            Lessingstr 8, 07743 Jena, Germany
          </span>
        </span>
      </a>
    </>
  );
}

// ============================================================================
// LICENSE FOOTER
// The coffee icon hides this page's one easter egg: caffeine, with a SMILES
// that really runs through the validator. Discovery reward for the curious;
// screen-reader users get it through the button label.
// ============================================================================

const CAFFEINE_SMILES = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C';

function LicenseFooter() {
  const [brewing, setBrewing] = useState(false);

  return (
    <div className="text-center py-10">
      <div className="flex items-center justify-center gap-1.5 mb-2">
        <span className="text-sm font-medium text-[var(--color-text-secondary)]">Made with</span>
        <button
          onClick={() => setBrewing((v) => !v)}
          aria-expanded={brewing}
          aria-controls="caffeine-note"
          aria-label="About the coffee in this line"
          className="group p-1 -m-1 rounded-lg cursor-pointer hover:bg-[var(--color-accent)]/10 transition-colors"
        >
          <Coffee
            className="w-4 h-4 text-[var(--color-accent)] transition-transform duration-200 ease-out group-hover:-rotate-12"
            aria-hidden="true"
          />
        </button>
        <span className="text-sm font-medium text-[var(--color-text-secondary)]">
          for the chemistry community
        </span>
      </div>

      {brewing && (
        <motion.div
          id="caffeine-note"
          initial={{ opacity: 0, y: 6 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.25, ease: [0.25, 1, 0.5, 1] }}
          className="mx-auto mb-4 max-w-sm"
        >
          <p className="text-xs text-[var(--color-text-secondary)] mb-2">
            Caffeine (C&#8328;H&#8321;&#8320;N&#8324;O&#8322;). The project's other dependency.
          </p>
          <Link
            to={`/?smiles=${encodeURIComponent(CAFFEINE_SMILES)}`}
            title="Run it through the validator"
            className={cn(
              'inline-block font-mono text-[11px] px-2.5 py-1 rounded-lg',
              'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
              'text-[var(--color-text-secondary)]',
              'hover:text-[var(--color-primary)] hover:border-[var(--color-primary)]/40',
              'transition-colors'
            )}
          >
            {CAFFEINE_SMILES}
          </Link>
        </motion.div>
      )}

      <p className="text-xs text-[var(--color-text-muted)]">
        ChemAudit is open-source software released under the{' '}
        <a
          href="https://opensource.org/licenses/MIT"
          target="_blank"
          rel="noopener noreferrer"
          className="text-[var(--color-primary)] hover:underline transition-colors"
        >
          MIT License
        </a>
      </p>
    </div>
  );
}

// ============================================================================
// EXPORTS
// ============================================================================

export default AboutPage;
