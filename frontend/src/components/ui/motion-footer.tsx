import { forwardRef, useEffect, useRef, type ReactNode, type ElementType, type ButtonHTMLAttributes, type AnchorHTMLAttributes } from 'react';
import { Link } from 'react-router-dom';
import { gsap } from 'gsap';
import { ScrollTrigger } from 'gsap/ScrollTrigger';
import { cn } from '../../lib/utils';

gsap.registerPlugin(ScrollTrigger);

const FOOTER_CSS = `
@import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700;800;900&display=swap');

.cinematic-footer-wrapper {
  font-family: 'Plus Jakarta Sans', sans-serif;
  -webkit-font-smoothing: antialiased;

  --pill-bg-1: color-mix(in oklch, var(--color-text-primary) 3%, transparent);
  --pill-bg-2: color-mix(in oklch, var(--color-text-primary) 1%, transparent);
  --pill-shadow: color-mix(in oklch, var(--color-surface) 50%, transparent);
  --pill-highlight: color-mix(in oklch, var(--color-text-primary) 10%, transparent);
  --pill-inset-shadow: color-mix(in oklch, var(--color-surface) 80%, transparent);
  --pill-border: color-mix(in oklch, var(--color-text-primary) 8%, transparent);

  --pill-bg-1-hover: color-mix(in oklch, var(--color-text-primary) 8%, transparent);
  --pill-bg-2-hover: color-mix(in oklch, var(--color-text-primary) 2%, transparent);
  --pill-border-hover: color-mix(in oklch, var(--color-text-primary) 20%, transparent);
  --pill-shadow-hover: color-mix(in oklch, var(--color-surface) 70%, transparent);
  --pill-highlight-hover: color-mix(in oklch, var(--color-text-primary) 20%, transparent);
}

@keyframes footer-breathe {
  0% { transform: translate(-50%, -50%) scale(1); opacity: 0.6; }
  100% { transform: translate(-50%, -50%) scale(1.1); opacity: 1; }
}

@keyframes footer-scroll-marquee {
  from { transform: translateX(0); }
  to { transform: translateX(-50%); }
}

@keyframes footer-cup-wobble {
  0%, 100% { transform: rotate(-0.5deg) translateY(0); }
  50% { transform: rotate(0.5deg) translateY(0); }
}

@keyframes footer-cup-hover {
  0% { transform: rotate(0deg) translateY(-2px) scale(1.15); }
  20% { transform: rotate(-1.5deg) translateY(-3px) scale(1.18); }
  40% { transform: rotate(1deg) translateY(-4px) scale(1.2); }
  60% { transform: rotate(-0.5deg) translateY(-3px) scale(1.18); }
  80% { transform: rotate(1.5deg) translateY(-2.5px) scale(1.16); }
  100% { transform: rotate(0deg) translateY(-2px) scale(1.15); }
}

@keyframes footer-steam {
  0% { opacity: 0; transform: translateY(0) translateX(0) scaleY(1); }
  30% { opacity: 0.5; }
  60% { opacity: 0.25; }
  100% { opacity: 0; transform: translateY(-10px) translateX(1px) scaleY(1); }
}

@keyframes footer-steam-hover {
  0% { opacity: 0; transform: translateY(0) translateX(0) scaleY(1); }
  15% { opacity: 0.7; }
  40% { opacity: 0.5; transform: translateY(-8px) translateX(0.5px) scaleY(1.2); }
  70% { opacity: 0.25; transform: translateY(-14px) translateX(1px) scaleY(1.35); }
  100% { opacity: 0; transform: translateY(-20px) translateX(1.5px) scaleY(1.5); }
}

.animate-footer-breathe {
  animation: footer-breathe 8s ease-in-out infinite alternate;
}

.animate-footer-scroll-marquee {
  animation: footer-scroll-marquee 40s linear infinite;
}

.footer-coffee-icon {
  animation: footer-cup-wobble 3s ease-in-out infinite;
  transition: filter 0.6s cubic-bezier(0.4, 0, 0.2, 1), color 0.6s ease, transform 0.6s cubic-bezier(0.4, 0, 0.2, 1);
  filter: drop-shadow(0 0 0px transparent);
}

.footer-coffee-group:hover .footer-coffee-icon {
  animation: footer-cup-hover 2.5s cubic-bezier(0.4, 0, 0.2, 1) infinite;
  filter: drop-shadow(0 0 6px color-mix(in oklch, var(--color-primary) 60%, transparent)) drop-shadow(0 2px 10px color-mix(in oklch, var(--color-primary) 25%, transparent));
  color: var(--color-primary-light, var(--color-primary));
}

.footer-steam-wisp {
  animation: footer-steam 2s ease-out infinite;
  transition: filter 0.6s ease, animation 0.6s ease;
}

.footer-steam-wisp:nth-child(2) {
  animation-delay: 0.5s;
  animation-duration: 2.2s;
}

.footer-steam-wisp:nth-child(3) {
  animation-delay: 0.3s;
  animation-duration: 1.8s;
}

.footer-coffee-group:hover .footer-steam-wisp {
  animation-name: footer-steam-hover;
  animation-duration: 1.8s;
  filter: drop-shadow(0 0 3px color-mix(in oklch, var(--color-primary) 20%, transparent));
}

.footer-coffee-group:hover .footer-steam-wisp:nth-child(1) {
  animation-delay: 0s;
  animation-duration: 1.8s;
}

.footer-coffee-group:hover .footer-steam-wisp:nth-child(2) {
  animation-delay: 0.4s;
  animation-duration: 1.6s;
}

.footer-coffee-group:hover .footer-steam-wisp:nth-child(3) {
  animation-delay: 0.2s;
  animation-duration: 1.7s;
}

.footer-bg-grid {
  background-size: 60px 60px;
  background-image:
    linear-gradient(to right, color-mix(in oklch, var(--color-text-primary) 3%, transparent) 1px, transparent 1px),
    linear-gradient(to bottom, color-mix(in oklch, var(--color-text-primary) 3%, transparent) 1px, transparent 1px);
  mask-image: linear-gradient(to bottom, transparent, black 30%, black 70%, transparent);
  -webkit-mask-image: linear-gradient(to bottom, transparent, black 30%, black 70%, transparent);
}

.footer-aurora {
  background: radial-gradient(
    circle at 50% 50%,
    color-mix(in oklch, var(--color-primary) 15%, transparent) 0%,
    color-mix(in oklch, var(--color-secondary) 15%, transparent) 40%,
    transparent 70%
  );
}

.footer-glass-pill {
  background: linear-gradient(145deg, var(--pill-bg-1) 0%, var(--pill-bg-2) 100%);
  box-shadow:
      0 10px 30px -10px var(--pill-shadow),
      inset 0 1px 1px var(--pill-highlight),
      inset 0 -1px 2px var(--pill-inset-shadow);
  border: 1px solid var(--pill-border);
  backdrop-filter: blur(16px);
  -webkit-backdrop-filter: blur(16px);
  transition: all 0.4s cubic-bezier(0.16, 1, 0.3, 1);
}

.footer-glass-pill:hover {
  background: linear-gradient(145deg, var(--pill-bg-1-hover) 0%, var(--pill-bg-2-hover) 100%);
  border-color: var(--pill-border-hover);
  box-shadow:
      0 20px 40px -10px var(--pill-shadow-hover),
      inset 0 1px 1px var(--pill-highlight-hover);
  color: var(--color-text-primary);
}

.footer-giant-bg-text {
  font-size: 20vw;
  line-height: 0.75;
  font-weight: 900;
  letter-spacing: -0.05em;
  color: transparent;
  -webkit-text-stroke: 1px color-mix(in oklch, var(--color-text-primary) 5%, transparent);
  background: linear-gradient(180deg, color-mix(in oklch, var(--color-text-primary) 10%, transparent) 0%, transparent 60%);
  -webkit-background-clip: text;
  background-clip: text;
}

.footer-text-glow {
  background: linear-gradient(180deg, var(--color-text-primary) 0%, color-mix(in oklch, var(--color-text-primary) 40%, transparent) 100%);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
  filter: drop-shadow(0px 0px 20px color-mix(in oklch, var(--color-text-primary) 15%, transparent));
}
`;

export type MagneticButtonProps = {
  as?: ElementType;
  to?: string;
  className?: string;
  children?: ReactNode;
} & ButtonHTMLAttributes<HTMLButtonElement> &
  AnchorHTMLAttributes<HTMLAnchorElement>;

const MagneticButton = forwardRef<HTMLElement, MagneticButtonProps>(
  function MagneticButton({ className, children, as: Component = 'button', ...props }, forwardedRef) {
    const localRef = useRef<HTMLElement>(null);

    useEffect(() => {
      const element = localRef.current;
      if (!element) return;

      const ctx = gsap.context(() => {
        const handleMouseMove = (e: MouseEvent) => {
          const rect = element.getBoundingClientRect();
          const cx = rect.width / 2;
          const cy = rect.height / 2;
          const x = e.clientX - rect.left - cx;
          const y = e.clientY - rect.top - cy;

          gsap.to(element, {
            x: x * 0.4,
            y: y * 0.4,
            rotationX: -y * 0.15,
            rotationY: x * 0.15,
            scale: 1.05,
            ease: 'power2.out',
            duration: 0.4,
          });
        };

        const handleMouseLeave = () => {
          gsap.to(element, {
            x: 0,
            y: 0,
            rotationX: 0,
            rotationY: 0,
            scale: 1,
            ease: 'elastic.out(1, 0.3)',
            duration: 1.2,
          });
        };

        element.addEventListener('mousemove', handleMouseMove as EventListener);
        element.addEventListener('mouseleave', handleMouseLeave);
      }, element);

      return () => ctx.revert();
    }, []);

    return (
      <Component
        ref={(node: HTMLElement) => {
          (localRef as React.MutableRefObject<HTMLElement | null>).current = node;
          if (typeof forwardedRef === 'function') {
            forwardedRef(node);
          } else if (forwardedRef) {
            (forwardedRef as React.MutableRefObject<HTMLElement | null>).current = node;
          }
        }}
        className={cn('cursor-pointer', className)}
        {...props}
      >
        {children}
      </Component>
    );
  }
);
MagneticButton.displayName = 'MagneticButton';

function MarqueeItem() {
  return (
    <div className="flex items-center space-x-12 px-6">
      <span>Structure Validation</span>
      <span style={{ color: 'var(--color-primary)', opacity: 0.6 }}>&#x2726;</span>
      <span>SMILES Standardization</span>
      <span style={{ color: 'var(--color-secondary)', opacity: 0.6 }}>&#x2726;</span>
      <span>ML-Readiness Scoring</span>
      <span style={{ color: 'var(--color-primary)', opacity: 0.6 }}>&#x2726;</span>
      <span>Batch Processing</span>
      <span style={{ color: 'var(--color-secondary)', opacity: 0.6 }}>&#x2726;</span>
      <span>Open Source</span>
      <span style={{ color: 'var(--color-primary)', opacity: 0.6 }}>&#x2726;</span>
    </div>
  );
}

const steamOffsets = [-2, 0, 2];

export function CinematicFooter() {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const giantTextRef = useRef<HTMLDivElement>(null);
  const headingRef = useRef<HTMLHeadingElement>(null);
  const linksRef = useRef<HTMLDivElement>(null);
  const styleRef = useRef<HTMLStyleElement>(null);

  useEffect(() => {
    if (styleRef.current) {
      styleRef.current.textContent = FOOTER_CSS;
    }
  }, []);

  useEffect(() => {
    if (!wrapperRef.current) return;

    const ctx = gsap.context(() => {
      gsap.fromTo(
        giantTextRef.current,
        { y: '10vh', scale: 0.8, opacity: 0 },
        {
          y: '0vh',
          scale: 1,
          opacity: 1,
          ease: 'power1.out',
          scrollTrigger: {
            trigger: wrapperRef.current,
            start: 'top 80%',
            end: 'bottom bottom',
            scrub: 1,
          },
        }
      );

      gsap.fromTo(
        [headingRef.current, linksRef.current],
        { y: 50, opacity: 0 },
        {
          y: 0,
          opacity: 1,
          stagger: 0.15,
          ease: 'power3.out',
          scrollTrigger: {
            trigger: wrapperRef.current,
            start: 'top 40%',
            end: 'bottom bottom',
            scrub: 1,
          },
        }
      );
    }, wrapperRef);

    return () => ctx.revert();
  }, []);

  function scrollToTop() {
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }

  return (
    <>
      <style ref={styleRef} />

      <div
        ref={wrapperRef}
        className="relative h-screen w-full"
        style={{ clipPath: 'polygon(0% 0, 100% 0%, 100% 100%, 0 100%)' }}
      >
        <footer
          className="fixed bottom-0 left-0 flex h-screen w-full flex-col justify-between overflow-hidden cinematic-footer-wrapper"
          style={{
            backgroundColor: 'var(--color-surface)',
            color: 'var(--color-text-primary)',
          }}
        >
          <div className="footer-aurora absolute left-1/2 top-1/2 h-[60vh] w-[80vw] -translate-x-1/2 -translate-y-1/2 animate-footer-breathe rounded-[50%] blur-[80px] pointer-events-none z-0" />
          <div className="footer-bg-grid absolute inset-0 z-0 pointer-events-none" />

          <div
            ref={giantTextRef}
            className="footer-giant-bg-text absolute -bottom-[5vh] left-1/2 -translate-x-1/2 whitespace-nowrap z-0 pointer-events-none select-none"
          >
            CHEMAUDIT
          </div>

          <div
            className="absolute top-12 left-0 w-full overflow-hidden py-4 z-10 -rotate-2 scale-110 shadow-2xl backdrop-blur-md"
            style={{
              borderTop: '1px solid var(--color-border)',
              borderBottom: '1px solid var(--color-border)',
              backgroundColor: 'color-mix(in oklch, var(--color-surface) 60%, transparent)',
            }}
          >
            <div
              className="flex w-max animate-footer-scroll-marquee text-xs md:text-sm font-bold tracking-[0.3em] uppercase"
              style={{ color: 'var(--color-text-muted)' }}
            >
              <MarqueeItem />
              <MarqueeItem />
            </div>
          </div>

          <div className="relative z-10 flex flex-1 flex-col items-center justify-center px-6 mt-20 w-full max-w-5xl mx-auto">
            <h2
              ref={headingRef}
              className="text-5xl md:text-8xl font-black footer-text-glow tracking-tighter mb-12 text-center"
              style={{ fontFamily: "'Outfit', sans-serif" }}
            >
              Every structure,<br />audited.
            </h2>

            <div ref={linksRef} className="flex flex-col items-center gap-6 w-full">
              <div className="flex flex-wrap justify-center gap-4 w-full">
                <MagneticButton
                  as="a"
                  href="https://github.com/Kohulan/ChemAudit"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="footer-glass-pill px-10 py-5 rounded-full font-bold text-sm md:text-base flex items-center gap-3 group"
                  style={{ color: 'var(--color-text-primary)' }}
                >
                  <svg
                    className="w-6 h-6 transition-colors"
                    style={{ color: 'var(--color-text-muted)' }}
                    viewBox="0 0 24 24"
                    fill="currentColor"
                  >
                    <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z" />
                  </svg>
                  Star on GitHub
                </MagneticButton>

                <MagneticButton
                  as={Link}
                  to="/"
                  onClick={() => window.scrollTo({ top: 0 })}
                  className="footer-glass-pill px-10 py-5 rounded-full font-bold text-sm md:text-base flex items-center gap-3 group"
                  style={{ color: 'var(--color-text-primary)' }}
                >
                  <svg
                    className="w-6 h-6 transition-colors"
                    style={{ color: 'var(--color-text-muted)' }}
                    viewBox="0 0 24 24"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  >
                    <path d="M10 2v7.527a2 2 0 0 1-.211.896L4.72 20.55a1 1 0 0 0 .9 1.45h12.76a1 1 0 0 0 .9-1.45l-5.069-10.127A2 2 0 0 1 14 9.527V2" />
                    <path d="M8.5 2h7" />
                    <path d="M7 16.5h10" />
                  </svg>
                  Try Validation
                </MagneticButton>
              </div>

              <div className="flex flex-wrap justify-center gap-3 md:gap-6 w-full mt-2">
                <MagneticButton
                  as={Link}
                  to="/privacy"
                  onClick={() => window.scrollTo({ top: 0 })}
                  className="footer-glass-pill px-6 py-3 rounded-full font-medium text-xs md:text-sm"
                  style={{ color: 'var(--color-text-muted)' }}
                >
                  Privacy Policy
                </MagneticButton>
                <MagneticButton
                  as={Link}
                  to="/about"
                  onClick={() => window.scrollTo({ top: 0 })}
                  className="footer-glass-pill px-6 py-3 rounded-full font-medium text-xs md:text-sm"
                  style={{ color: 'var(--color-text-muted)' }}
                >
                  About
                </MagneticButton>
                <MagneticButton
                  as="a"
                  href="https://github.com/Kohulan/ChemAudit/issues"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="footer-glass-pill px-6 py-3 rounded-full font-medium text-xs md:text-sm"
                  style={{ color: 'var(--color-text-muted)' }}
                >
                  Support
                </MagneticButton>
              </div>
            </div>
          </div>

          <div className="relative z-20 w-full pb-8 px-6 md:px-12 flex flex-col md:flex-row items-center justify-between gap-6">
            <div
              className="text-[10px] md:text-xs font-semibold tracking-widest uppercase order-2 md:order-1"
              style={{ color: 'var(--color-text-muted)' }}
            >
              &copy; {new Date().getFullYear()} ChemAudit. All rights reserved.
            </div>

            <div className="footer-glass-pill px-6 py-3 rounded-full flex items-center gap-2 order-1 md:order-2 cursor-default">
              <span
                className="text-[10px] md:text-xs font-bold uppercase tracking-widest"
                style={{ color: 'var(--color-text-muted)' }}
              >
                Made with
              </span>
              <div className="footer-coffee-group relative flex items-center justify-center cursor-pointer">
                <div className="absolute -top-2.5 left-1/2 -translate-x-1/2 flex gap-[1px]">
                  {steamOffsets.map((offset, i) => (
                    <svg
                      key={i}
                      width="5"
                      height="10"
                      viewBox="0 0 5 10"
                      className="footer-steam-wisp"
                      style={{ color: 'var(--color-text-muted)', marginLeft: offset }}
                    >
                      <path
                        d="M2.5 10C2.5 10 1 7 1 5C1 3 2.5 1 2.5 0C2.5 1 4 3 4 5C4 7 2.5 10 2.5 10Z"
                        fill="currentColor"
                        fillOpacity="0.6"
                      />
                    </svg>
                  ))}
                </div>
                <svg
                  className="w-6 h-6 footer-coffee-icon"
                  style={{ color: 'var(--color-primary)' }}
                  viewBox="0 0 24 24"
                  fill="none"
                  stroke="currentColor"
                  strokeWidth="2"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                >
                  <path d="M17 8h1a4 4 0 1 1 0 8h-1" />
                  <path d="M3 8h14v9a4 4 0 0 1-4 4H7a4 4 0 0 1-4-4Z" />
                </svg>
              </div>
              <span
                className="text-[10px] md:text-xs font-bold uppercase tracking-widest"
                style={{ color: 'var(--color-text-muted)' }}
              >
                by
              </span>
              <a
                href="https://kohulanr.com"
                target="_blank"
                rel="noopener noreferrer"
                className="font-black text-xs md:text-sm tracking-normal ml-1 hover:opacity-80 transition-opacity"
                style={{ color: 'var(--color-text-primary)' }}
              >
                Kohulan.R
              </a>
              <span
                className="text-[10px] md:text-xs font-bold uppercase tracking-widest"
                style={{ color: 'var(--color-text-muted)' }}
              >
                at
              </span>
              <a
                href="https://www.uni-jena.de"
                target="_blank"
                rel="noopener noreferrer"
                className="font-semibold text-xs md:text-sm tracking-normal hover:opacity-80 transition-opacity"
                style={{ color: 'var(--color-text-secondary)' }}
              >
                Friedrich Schiller University Jena
              </a>
            </div>

            <MagneticButton
              onClick={scrollToTop}
              className="w-12 h-12 rounded-full footer-glass-pill flex items-center justify-center group order-3"
              style={{ color: 'var(--color-text-muted)' }}
            >
              <svg
                className="w-5 h-5 transform group-hover:-translate-y-1.5 transition-transform duration-300"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth="2"
                  d="M5 10l7-7m0 0l7 7m-7-7v18"
                />
              </svg>
            </MagneticButton>
          </div>
        </footer>
      </div>
    </>
  );
}
