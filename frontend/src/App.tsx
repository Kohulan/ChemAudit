import { BrowserRouter as Router, Routes, Route } from 'react-router-dom'

function HomePage() {
  return (
    <div className="min-h-screen bg-gray-50 flex items-center justify-center">
      <div className="text-center">
        <h1 className="text-4xl font-bold text-gray-900 mb-4">
          ChemStructVal
        </h1>
        <p className="text-lg text-gray-600">
          Chemical Structure Validation Suite
        </p>
      </div>
    </div>
  )
}

function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<HomePage />} />
      </Routes>
    </Router>
  )
}

export default App
