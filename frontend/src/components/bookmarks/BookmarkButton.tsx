/**
 * Star icon button that toggles bookmark state for a molecule.
 */

import { useState, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Star } from 'lucide-react';
import { bookmarksApi } from '../../services/api';
import { cn } from '../../lib/utils';

interface BookmarkButtonProps {
  smiles: string;
  name?: string;
  source?: string;
  jobId?: string;
  isBookmarked?: boolean;
  onBookmark?: (bookmarked: boolean) => void;
  className?: string;
}

export function BookmarkButton({
  smiles,
  name,
  source,
  jobId,
  isBookmarked: initialBookmarked = false,
  onBookmark,
  className,
}: BookmarkButtonProps) {
  const [bookmarked, setBookmarked] = useState(initialBookmarked);
  const [isLoading, setIsLoading] = useState(false);
  const [showToast, setShowToast] = useState(false);

  const handleClick = useCallback(async () => {
    if (isLoading || !smiles.trim()) return;
    setIsLoading(true);
    try {
      if (!bookmarked) {
        await bookmarksApi.createBookmark({
          smiles: smiles.trim(),
          name: name ?? undefined,
          source: source ?? 'single_validation',
          job_id: jobId ?? undefined,
        });
        setBookmarked(true);
        onBookmark?.(true);
        setShowToast(true);
        setTimeout(() => setShowToast(false), 2000);
      }
      // For unbookmark, the user would use the Bookmarks page to remove
    } catch {
      // Silently fail - could already be bookmarked
    } finally {
      setIsLoading(false);
    }
  }, [bookmarked, smiles, name, source, jobId, onBookmark, isLoading]);

  return (
    <div className="relative inline-flex">
      <motion.button
        onClick={handleClick}
        disabled={isLoading || bookmarked}
        className={cn(
          'relative p-2 rounded-lg transition-all duration-200',
          'border shadow-sm',
          bookmarked
            ? 'bg-amber-500/10 text-amber-500 border-amber-500/30 cursor-default'
            : 'bg-[var(--color-surface-elevated)] text-[var(--color-text-secondary)] border-[var(--color-border)] hover:border-amber-500/50 hover:text-amber-500',
          className
        )}
        whileHover={!bookmarked ? { scale: 1.1 } : undefined}
        whileTap={!bookmarked ? { scale: 0.9 } : undefined}
        title={bookmarked ? 'Bookmarked' : 'Bookmark this molecule'}
      >
        <Star
          className={cn('w-4 h-4', bookmarked && 'fill-amber-500')}
        />
      </motion.button>

      {/* Toast */}
      <AnimatePresence>
        {showToast && (
          <motion.div
            initial={{ opacity: 0, y: 5, scale: 0.9 }}
            animate={{ opacity: 1, y: 0, scale: 1 }}
            exit={{ opacity: 0, scale: 0.9 }}
            className="absolute -top-8 left-1/2 -translate-x-1/2 whitespace-nowrap px-2.5 py-1 rounded-lg bg-amber-500 text-white text-xs font-medium shadow-lg"
          >
            Bookmarked!
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
