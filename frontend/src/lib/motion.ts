import type { Variants } from 'framer-motion';

/**
 * Stagger container for child animations.
 *
 * Used with a parent motion element whose children are motion elements
 * that will enter with a small staggered delay.
 */
export const staggerContainer: Variants = {
  initial: {},
  animate: {
    transition: {
      staggerChildren: 0.1,
      delayChildren: 0.1,
    },
  },
};

/**
 * Hover lift animation for cards.
 *
 * Raises the element on hover and resets on tap. Designed to be spread
 * onto a motion element via {...hoverLift}.
 */
export const hoverLift = {
  whileHover: {
    y: -2,
    transition: { duration: 0.2 },
  },
  whileTap: {
    y: 0,
    transition: { duration: 0.1 },
  },
};
