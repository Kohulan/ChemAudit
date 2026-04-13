import { motion } from 'framer-motion';
import { ClayCard } from '../ui/ClayCard';
import { Badge } from '../ui/Badge';
import type { RoundTripLoss } from '../../types/diagnostics';

interface LostInfoCardProps {
  losses: RoundTripLoss[];
}

/**
 * Card listing specific round-trip information losses (stereo, charge, isotope).
 *
 * Per UI-SPEC Round-Trip Pipeline Flow contract (D-08):
 * - Renders as ClayCard with status-error 30% alpha border
 * - Groups losses by type with before/after counts
 * - Omits rows where loss.before === loss.after (no actual loss)
 * - Scale-in animation (scale 0.96->1, opacity 0->1, 0.3s ease-out)
 */
export function LostInfoCard({ losses }: LostInfoCardProps) {
  // Filter to only actual losses
  const actualLosses = losses.filter((l) => l.before !== l.after);

  if (actualLosses.length === 0) return null;

  function getLossBadge(type: string) {
    switch (type) {
      case 'stereo':
        return <Badge variant="warning" size="sm">Stereo</Badge>;
      case 'charge':
        return <Badge variant="warning" size="sm">Charge</Badge>;
      case 'isotope':
        return <Badge variant="warning" size="sm">Isotope</Badge>;
      default:
        return <Badge size="sm">{type}</Badge>;
    }
  }

  return (
    <motion.div
      initial={{ scale: 0.96, opacity: 0 }}
      animate={{ scale: 1, opacity: 1 }}
      transition={{ duration: 0.3, ease: 'easeOut' }}
    >
      <ClayCard
        variant="default"
        size="sm"
        className="border border-[rgba(239,68,68,0.3)] mt-4"
      >
        <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display mb-3">
          Lost Information
        </p>
        <div className="space-y-2">
          {actualLosses.map((loss, idx) => (
            <div
              key={idx}
              className="flex items-center justify-between gap-3 min-h-[44px] py-1"
            >
              <div className="flex items-center gap-2 flex-1 min-w-0">
                {getLossBadge(loss.type)}
                <span className="text-sm text-[var(--color-text-secondary)] truncate">
                  {loss.description}
                </span>
              </div>
              <span className="text-xs font-mono text-[var(--color-text-muted)] shrink-0 whitespace-nowrap">
                Before: {loss.before} | After: {loss.after}
              </span>
            </div>
          ))}
        </div>
      </ClayCard>
    </motion.div>
  );
}
