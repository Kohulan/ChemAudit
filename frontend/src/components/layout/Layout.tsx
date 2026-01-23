import { Header } from './Header';

interface LayoutProps {
  children: React.ReactNode;
}

/**
 * Main layout wrapper with chemistry-themed styling
 */
export function Layout({ children }: LayoutProps) {
  return (
    <div className="min-h-screen bg-gradient-chem-light bg-molecular-pattern">
      <Header />
      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {children}
      </main>
      {/* Footer */}
      <footer className="border-t border-chem-dark/5 bg-white/50 backdrop-blur-sm">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
          <div className="flex flex-col sm:flex-row items-center justify-between gap-4">
            <p className="text-sm text-chem-dark/50">
              ChemStructVal - Chemical Structure Validation Suite
            </p>
            <div className="flex items-center gap-6">
              <a
                href="http://localhost:8000/docs"
                target="_blank"
                rel="noopener noreferrer"
                className="text-sm text-chem-dark/50 hover:text-chem-primary transition-colors"
              >
                API Documentation
              </a>
              <a
                href="https://github.com"
                target="_blank"
                rel="noopener noreferrer"
                className="text-sm text-chem-dark/50 hover:text-chem-primary transition-colors"
              >
                GitHub
              </a>
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
}
