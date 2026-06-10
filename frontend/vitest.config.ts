import { defineConfig } from 'vitest/config';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  test: {
    globals: true,
    environment: 'jsdom',
    setupFiles: ['./src/tests/setup.ts'],
    include: ['src/**/*.{test,spec}.{js,mjs,cjs,ts,mts,cts,jsx,tsx}'],
    coverage: {
      provider: 'v8',
      reporter: ['text', 'json', 'html'],
      exclude: [
        'node_modules/',
        'src/tests/',
        '**/*.d.ts',
        '**/*.config.*',
      ],
      // Ratchet floor — set just below the current measured totals (lines ~31%,
      // branches ~34%). Only enforced when coverage is collected (`--coverage`),
      // so plain `npm test` runs are unaffected. Raise as page-level tests land.
      thresholds: {
        statements: 30,
        branches: 33,
        functions: 30,
        lines: 30,
      },
    },
  },
  resolve: {
    alias: {
      '@': '/src',
    },
  },
});
