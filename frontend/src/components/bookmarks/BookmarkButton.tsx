/**
 * Labelled bookmark button that toggles bookmark state for a molecule.
 * Renders as a ClayButton with a star icon and visible label text.
 */

import { useState, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Star } from 'lucide-react';
import { bookmarksApi } from '../../services/api';
import { ClayButton } from '../ui/ClayButton';

interface BookmarkButtonProps {
  smiles: string;
  name?: string;
  source?: string;
  jobId?: string;
  isBookmarked?: boolean;
  onBookmark?: (bookmarked: boolean) => void;
  onAfterBookmark?: (bookmarkId: number) => void;
  className?: string;
}

export function BookmarkButton({
  smiles,
  name,
  source,
  jobId,
  isBookmarked: initialBookmarked = false,
  onBookmark,
  onAfterBookmark,
  className,
}: BookmarkButtonProps) {
  const [bookmarked, setBookmarked] = useState(initialBookmarked);
  const [isLoading, setIsLoading] = useState(false);
  const [showToast, setShowToast] = useState(false);

  const handleClick = useCallback(async () => {
    if (isLoading || !smiles.trim() || bookmarked) return;
    setIsLoading(true);
    try {
      const bookmark = await bookmarksApi.createBookmark({
        smiles: smiles.trim(),
        name: name ?? undefined,
        source: source ?? 'single_validation',
        job_id: jobId ?? undefined,
      });
      setBookmarked(true);
      onBookmark?.(true);
      onAfterBookmark?.(bookmark.id);
      setShowToast(true);
      setTimeout(() => setShowToast(false), 2000);
    } catch {
      // Silently fail - could already be bookmarked
    } finally {
      setIsLoading(false);
    }
  }, [bookmarked, smiles, name, source, jobId, onBookmark, onAfterBookmark, isLoading]);

  return (
    <div className="relative inline-flex">
      <ClayButton
        variant="stone"
        onClick={handleClick}
        disabled={isLoading || bookmarked}
        loading={isLoading}
        leftIcon={<Star className={bookmarked ? 'w-4 h-4 fill-amber-300' : 'w-4 h-4'} />}
        className={className}
        title={bookmarked ? 'Bookmarked' : 'Bookmark this molecule'}
      >
        {bookmarked ? 'Bookmarked' : 'Bookmark Result'}
      </ClayButton>

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
