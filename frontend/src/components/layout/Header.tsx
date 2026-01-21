export function Header() {
  return (
    <header className="bg-white border-b border-gray-200 px-6 py-4">
      <div className="max-w-7xl mx-auto flex items-center justify-between">
        <div className="flex items-center space-x-3">
          <div className="w-8 h-8 bg-blue-600 rounded-lg flex items-center justify-center">
            <span className="text-white font-bold text-sm">CV</span>
          </div>
          <h1 className="text-xl font-semibold text-gray-900">ChemStructVal</h1>
        </div>
        <nav className="flex items-center space-x-6">
          <a href="/" className="text-gray-600 hover:text-gray-900">Validate</a>
          <a href="/docs" className="text-gray-600 hover:text-gray-900">API Docs</a>
        </nav>
      </div>
    </header>
  );
}
