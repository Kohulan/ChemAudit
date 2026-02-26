import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
import { ThemeProvider } from './contexts/ThemeContext';
import { ValidationCacheProvider } from './contexts/ValidationCacheContext';
import { BatchCacheProvider } from './contexts/BatchCacheContext';
import './index.css';

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
