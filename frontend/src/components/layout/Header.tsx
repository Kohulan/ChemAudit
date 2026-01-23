/**
 * Application header with chemistry-themed gradient and branding
 */
export function Header() {
  return (
    <header className="bg-gradient-chem shadow-chem sticky top-0 z-50">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex items-center justify-between h-16">
          {/* Logo and branding */}
          <div className="flex items-center gap-3">
            {/* Molecule icon */}
            <div className="w-10 h-10 bg-white/20 rounded-xl flex items-center justify-center backdrop-blur-sm">
              <svg
                className="w-6 h-6 text-white"
                viewBox="0 0 24 24"
                fill="currentColor"
              >
                {/* Benzene-like hexagon with atoms */}
                <circle cx="12" cy="12" r="3" />
                <circle cx="12" cy="4" r="2" />
                <circle cx="18.9" cy="8" r="2" />
                <circle cx="18.9" cy="16" r="2" />
                <circle cx="12" cy="20" r="2" />
                <circle cx="5.1" cy="16" r="2" />
                <circle cx="5.1" cy="8" r="2" />
                {/* Bonds */}
                <line x1="12" y1="9" x2="12" y2="6" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="14.6" y1="10.5" x2="17" y2="9" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="14.6" y1="13.5" x2="17" y2="15" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="12" y1="15" x2="12" y2="18" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="9.4" y1="13.5" x2="7" y2="15" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
                <line x1="9.4" y1="10.5" x2="7" y2="9" stroke="currentColor" strokeWidth="1.5" opacity="0.7" />
              </svg>
            </div>
            <div>
              <h1 className="text-xl font-bold text-white tracking-tight">
                ChemStructVal
              </h1>
              <p className="text-xs text-white/70 -mt-0.5">
                Structure Validation Suite
              </p>
            </div>
          </div>

          {/* Navigation */}
          <nav className="flex items-center gap-2">
            <a
              href="/"
              className="px-4 py-2 text-sm font-medium text-white/90 hover:text-white hover:bg-white/10 rounded-lg transition-all"
            >
              Validate
            </a>
            <a
              href="http://localhost:8000/docs"
              target="_blank"
              rel="noopener noreferrer"
              className="px-4 py-2 text-sm font-medium text-white/90 hover:text-white hover:bg-white/10 rounded-lg transition-all flex items-center gap-2"
            >
              <svg className="w-4 h-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M9 12h6M9 16h6M17 21H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
              API Docs
              <svg className="w-3 h-3 opacity-60" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M18 13v6a2 2 0 01-2 2H5a2 2 0 01-2-2V8a2 2 0 012-2h6M15 3h6v6M10 14L21 3" />
              </svg>
            </a>
          </nav>
        </div>
      </div>
    </header>
  );
}
