import { motion } from 'framer-motion';

interface PMITernaryPlotProps {
  npr1: number;
  npr2: number;
  shapeClass: 'rod' | 'disc' | 'sphere';
}

/**
 * Converts NPR1/NPR2 barycentric coordinates to SVG pixel coordinates.
 *
 * Triangle vertices:
 *   Rod    (150, 10)  — top
 *   Disc   (10, 250)  — bottom-left
 *   Sphere (290, 250) — bottom-right
 *
 * Barycentric weights:
 *   rodW   = 1 - npr2
 *   discW  = npr2 * (1 - npr1)
 *   sphW   = npr2 * npr1
 */
function nprToSVG(npr1: number, npr2: number): { x: number; y: number } {
  const rod = { x: 150, y: 10 };
  const disc = { x: 10, y: 250 };
  const sphere = { x: 290, y: 250 };
  const rodW = 1 - npr2;
  const discW = npr2 * (1 - npr1);
  const sphW = npr2 * npr1;
  return {
    x: rodW * rod.x + discW * disc.x + sphW * sphere.x,
    y: rodW * rod.y + discW * disc.y + sphW * sphere.y,
  };
}

/**
 * Custom SVG PMI ternary plot.
 *
 * Renders an equilateral triangle with rod/disc/sphere vertices
 * and plots the molecule's NPR1/NPR2 position as an animated point.
 *
 * Conforms to UI-SPEC Animation Contract and Accessibility Contract.
 */
export function PMITernaryPlot({ npr1, npr2, shapeClass }: PMITernaryPlotProps) {
  const point = nprToSVG(npr1, npr2);

  return (
    <svg
      viewBox="0 0 300 270"
      role="img"
      aria-label={`PMI ternary plot showing ${shapeClass} classification`}
      className="w-full max-w-[280px] mx-auto"
    >
      <title>{`PMI ternary plot -- ${shapeClass} region, position: NPR1=${npr1.toFixed(3)}, NPR2=${npr2.toFixed(3)}`}</title>

      {/* Region fills — subtle background shading for each vertex region */}
      {/* Rod region (top) */}
      <path
        d="M150,10 L80,130 L220,130 Z"
        fill="var(--color-primary, #c41e3a)"
        opacity="0.05"
      />
      {/* Disc region (bottom-left) */}
      <path
        d="M10,250 L80,130 L150,250 Z"
        fill="var(--color-accent, #d97706)"
        opacity="0.05"
      />
      {/* Sphere region (bottom-right) */}
      <path
        d="M290,250 L220,130 L150,250 Z"
        fill="#22c55e"
        opacity="0.05"
      />

      {/* Triangle outline */}
      <polygon
        points="150,10 10,250 290,250"
        fill="none"
        stroke="var(--color-border, #d1cfc9)"
        strokeWidth="1.5"
      />

      {/* Grid lines — midpoint lines for visual reference */}
      <line x1="80" y1="130" x2="220" y2="130" stroke="var(--color-border, #d1cfc9)" strokeWidth="0.5" strokeDasharray="3,3" opacity="0.6" />
      <line x1="80" y1="130" x2="150" y2="250" stroke="var(--color-border, #d1cfc9)" strokeWidth="0.5" strokeDasharray="3,3" opacity="0.6" />
      <line x1="220" y1="130" x2="150" y2="250" stroke="var(--color-border, #d1cfc9)" strokeWidth="0.5" strokeDasharray="3,3" opacity="0.6" />

      {/* Vertex labels */}
      <text
        x="150"
        y="6"
        textAnchor="middle"
        dominantBaseline="hanging"
        fontSize="11"
        fontWeight="600"
        fill="var(--color-text-secondary, #6b7280)"
        fontFamily="inherit"
      >
        Rod
      </text>
      <text
        x="4"
        y="254"
        textAnchor="start"
        dominantBaseline="hanging"
        fontSize="11"
        fontWeight="600"
        fill="var(--color-text-secondary, #6b7280)"
        fontFamily="inherit"
      >
        Disc
      </text>
      <text
        x="296"
        y="254"
        textAnchor="end"
        dominantBaseline="hanging"
        fontSize="11"
        fontWeight="600"
        fill="var(--color-text-secondary, #6b7280)"
        fontFamily="inherit"
      >
        Sphere
      </text>

      {/* Molecule point — animated entrance per UI-SPEC Animation Contract */}
      <motion.circle
        cx={point.x}
        cy={point.y}
        r={6}
        fill="var(--color-primary, #c41e3a)"
        stroke="white"
        strokeWidth="2"
        initial={{ opacity: 0, scale: 0 }}
        animate={{ opacity: 1, scale: 1 }}
        transition={{ duration: 0.5, type: 'spring', stiffness: 200, damping: 15 }}
      />
    </svg>
  );
}
