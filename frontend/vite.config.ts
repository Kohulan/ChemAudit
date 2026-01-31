import { defineConfig, loadEnv } from 'vite'
import react from '@vitejs/plugin-react'

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  // Load env file based on `mode`
  const env = loadEnv(mode, process.cwd(), '')

  // Parse allowed hosts from env (comma-separated)
  // Check both loadEnv result and process.env (for Docker)
  const allowedHostsStr = env.VITE_ALLOWED_HOSTS || process.env.VITE_ALLOWED_HOSTS
  const allowedHosts = allowedHostsStr
    ? allowedHostsStr.split(',').map(h => h.trim())
    : true  // Allow all hosts if not specified (safe for dev server)

  return {
    plugins: [react()],
    server: {
      host: true,
      port: 3002,
      allowedHosts,
      proxy: {
        '/api': {
          target: 'http://backend:8000',
          changeOrigin: true,
        },
      },
    },
    define: {
      'process.env.VITE_API_URL': JSON.stringify(
        process.env.VITE_API_URL || '/api/v1'
      ),
    },
  }
})
