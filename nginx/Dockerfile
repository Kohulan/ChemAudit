# Multi-stage build: Frontend + Nginx reverse proxy
# Frontend files are baked into the image, eliminating shared volume issues

# Stage 1: Build the React frontend
FROM node:20-alpine AS frontend-builder

ARG HTTP_PROXY
ARG HTTPS_PROXY
ARG NO_PROXY

WORKDIR /app

COPY frontend/package*.json ./
RUN npm ci --prefer-offline --no-audit

COPY frontend/ .
RUN npm run build

# Stage 2: Nginx with built frontend and reverse proxy config
FROM nginx:1.25-alpine

COPY nginx/nginx.conf /etc/nginx/conf.d/default.conf
COPY --from=frontend-builder /app/dist /usr/share/nginx/html

RUN chown -R nginx:nginx /usr/share/nginx/html && \
    chmod -R 755 /usr/share/nginx/html

EXPOSE 80

HEALTHCHECK --interval=30s --timeout=10s --start-period=10s --retries=3 \
    CMD wget --no-verbose --tries=1 --spider http://localhost/health || exit 1

CMD ["nginx", "-g", "daemon off;"]
