import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
import { ThemeProvider } from './contexts/ThemeContext';
import { ValidationCacheProvider } from './contexts/ValidationCacheContext';
import { BatchCacheProvider } from './contexts/BatchCacheContext';
import './index.css';

// Hidden console signature for the developer-chemists who open DevTools.
// The citation is the original SMILES paper that defines every input this
// tool accepts. Quoted line is verbatim from the paper.
if (typeof window !== 'undefined' && import.meta.env.MODE !== 'test') {
  // eslint-disable-next-line no-console
  console.log(
    '%cChemAudit%c — Chemical Structure Validation\n' +
      '%cWeininger, D. J. Chem. Inf. Comput. Sci. 1989, 29(2), 97-101\n' +
      '%c"The SMILES language was designed to be as compact as possible."\n\n' +
      '%cBackend%c FastAPI + RDKit (Python)   %cFrontend%c React + RDKit.js (Wasm)',
    'font-weight:600;color:#c41e3a;font-size:13px',
    'color:#5c5650',
    'color:#9c958d;font-style:italic',
    'color:#5c5650',
    'font-weight:600;color:#1a1815',
    'color:#5c5650',
    'font-weight:600;color:#1a1815',
    'color:#5c5650',
  );
}

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <ThemeProvider>
      <ValidationCacheProvider>
        <BatchCacheProvider>
          <App />
        </BatchCacheProvider>
      </ValidationCacheProvider>
    </ThemeProvider>
  </React.StrictMode>
);
