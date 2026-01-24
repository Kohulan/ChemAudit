# Docker Deployment Guide

Deploy ChemStructVal using Docker Compose for development or production environments.

## Prerequisites

- **Docker**: 24.0+ ([install](https://docs.docker.com/get-docker/))
- **Docker Compose**: 2.20+ (included with Docker Desktop)
- **System Requirements**:
  - 4GB RAM minimum (8GB recommended)
  - 10GB disk space
  - Linux, macOS, or Windows with WSL2

## Quick Start (Development)

For local development:

```bash
# Clone repository
git clone https://github.com/yourusername/chemstructval.git
cd chemstructval

# Start all services
docker-compose up -d

# View logs
docker-compose logs -f

# Access application
# Frontend: http://localhost:3000
# Backend API: http://localhost:8000
# API Docs: http://localhost:8000/docs
# Grafana: http://localhost:3001 (admin/chemstructval)
# Prometheus: http://localhost:9090
```

## Service Architecture

```
┌──────────────────────────────────────────────────────┐
│                   ChemStructVal Stack                │
├──────────────────────────────────────────────────────┤
│                                                      │
│  [Frontend :3000] ────────┐                         │
│                           ▼                          │
│                    [Backend :8000]                   │
│                           │                          │
│         ┌─────────────────┼─────────────────┐       │
│         ▼                 ▼                 ▼        │
│  [PostgreSQL :5432]  [Redis :6379]  [Metrics]       │
│                           │                          │
│                           ▼                          │
│                  [Celery Workers x4]                 │
│                                                      │
│  [Prometheus :9090] ◄──── Scrapes metrics           │
│         │                                            │
│         ▼                                            │
│  [Grafana :3001] ──── Dashboards                    │
│                                                      │
└──────────────────────────────────────────────────────┘
```

## Development Configuration

### docker-compose.yml

The development configuration includes:

- **Hot reload** for backend (uvicorn --reload)
- **Hot reload** for frontend (Vite dev server)
- **Volume mounts** for code changes
- **Exposed ports** for direct access
- **Monitoring stack** (Prometheus + Grafana)

### Environment Variables

Development defaults are set in `docker-compose.yml`:

```yaml
DATABASE_URL: postgresql+asyncpg://chemstructval:chemstructval@postgres:5432/chemstructval
REDIS_URL: redis://redis:6379/0
ENABLE_METRICS: true
```

### Starting Services

```bash
# Start all services
docker-compose up -d

# Start specific service
docker-compose up -d backend

# View logs (all services)
docker-compose logs -f

# View logs (specific service)
docker-compose logs -f backend

# Check status
docker-compose ps

# Stop services
docker-compose down

# Stop and remove volumes (clean slate)
docker-compose down -v
```

## Production Deployment

Production deployment uses optimized build configurations.

### Step 1: Configure Environment

Create production environment file:

```bash
cp .env.example .env
```

Edit `.env` with production values:

```bash
# Database
DATABASE_URL=postgresql+asyncpg://user:password@postgres:5432/chemstructval

# Redis
REDIS_URL=redis://redis:6379/0
REDIS_RATE_LIMIT_URL=redis://redis:6379/1

# Security
SECRET_KEY=your-secret-key-here  # Generate with: openssl rand -hex 32

# CORS
CORS_ORIGINS=["https://yourdomain.com"]

# Monitoring
ENABLE_METRICS=true
GRAFANA_PASSWORD=secure-password-here

# Optional: Email for Let's Encrypt
LETSENCRYPT_EMAIL=admin@yourdomain.com

# Optional: Domain for SSL
DOMAIN=yourdomain.com
```

### Step 2: Create Production Compose File

Create `docker-compose.prod.yml`:

```yaml
version: '3.8'

services:
  postgres:
    image: postgres:16-alpine
    container_name: chemstructval-postgres
    environment:
      POSTGRES_USER: chemstructval
      POSTGRES_PASSWORD: ${DB_PASSWORD}
      POSTGRES_DB: chemstructval
    volumes:
      - postgres_data:/var/lib/postgresql/data
    restart: unless-stopped
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U chemstructval"]
      interval: 10s
      timeout: 5s
      retries: 5

  redis:
    image: redis:7-alpine
    container_name: chemstructval-redis
    command: redis-server --maxmemory 512mb --maxmemory-policy allkeys-lru
    volumes:
      - redis_data:/data
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "redis-cli", "ping"]
      interval: 10s
      timeout: 5s
      retries: 5

  backend:
    build:
      context: ./backend
      dockerfile: Dockerfile
    container_name: chemstructval-backend
    command: uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4
    environment:
      - DATABASE_URL=${DATABASE_URL}
      - REDIS_URL=${REDIS_URL}
      - SECRET_KEY=${SECRET_KEY}
      - CORS_ORIGINS=${CORS_ORIGINS}
      - ENABLE_METRICS=true
    depends_on:
      postgres:
        condition: service_healthy
      redis:
        condition: service_healthy
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/api/v1/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  celery-worker:
    build:
      context: ./backend
      dockerfile: Dockerfile
    command: celery -A app.celery_app worker --loglevel=info --concurrency=4
    environment:
      - DATABASE_URL=${DATABASE_URL}
      - REDIS_URL=${REDIS_URL}
    depends_on:
      redis:
        condition: service_healthy
    restart: unless-stopped
    deploy:
      replicas: 2  # Multiple workers for high throughput

  nginx:
    image: nginx:alpine
    container_name: chemstructval-nginx
    volumes:
      - ./nginx/nginx.conf:/etc/nginx/nginx.conf:ro
      - ./frontend/dist:/usr/share/nginx/html:ro
      - ./certbot/conf:/etc/letsencrypt:ro
      - ./certbot/www:/var/www/certbot:ro
    ports:
      - "80:80"
      - "443:443"
    depends_on:
      - backend
    restart: unless-stopped

  prometheus:
    image: prom/prometheus:v2.47.0
    container_name: chemstructval-prometheus
    command:
      - '--config.file=/etc/prometheus/prometheus.yml'
      - '--storage.tsdb.path=/prometheus'
      - '--storage.tsdb.retention.time=30d'
    volumes:
      - ./prometheus.yml:/etc/prometheus/prometheus.yml:ro
      - prometheus_data:/prometheus
    restart: unless-stopped

  grafana:
    image: grafana/grafana:10.2.0
    container_name: chemstructval-grafana
    environment:
      - GF_SECURITY_ADMIN_USER=admin
      - GF_SECURITY_ADMIN_PASSWORD=${GRAFANA_PASSWORD}
      - GF_USERS_ALLOW_SIGN_UP=false
      - GF_SERVER_ROOT_URL=https://${DOMAIN}/grafana
    volumes:
      - ./grafana/provisioning:/etc/grafana/provisioning:ro
      - ./grafana/dashboards:/var/lib/grafana/dashboards:ro
      - grafana_data:/var/lib/grafana
    restart: unless-stopped

volumes:
  postgres_data:
  redis_data:
  prometheus_data:
  grafana_data:
```

### Step 3: Build Frontend

```bash
cd frontend
npm install
npm run build
cd ..
```

This creates `frontend/dist/` with optimized static files.

### Step 4: Build Images

```bash
docker-compose -f docker-compose.prod.yml build
```

### Step 5: Start Production Stack

```bash
docker-compose -f docker-compose.prod.yml up -d
```

### Step 6: Verify Deployment

```bash
# Check all services are running
docker-compose -f docker-compose.prod.yml ps

# Check backend health
curl http://localhost:8000/api/v1/health

# Check frontend (via nginx)
curl http://localhost/

# Check metrics
curl http://localhost:8000/metrics
```

## Nginx Configuration

Create `nginx/nginx.conf`:

```nginx
events {
    worker_connections 1024;
}

http {
    include /etc/nginx/mime.types;
    default_type application/octet-stream;

    # Logging
    access_log /var/log/nginx/access.log;
    error_log /var/log/nginx/error.log;

    # Gzip
    gzip on;
    gzip_types text/plain text/css application/json application/javascript text/xml application/xml application/xml+rss text/javascript;

    # Rate limiting
    limit_req_zone $binary_remote_addr zone=api:10m rate=100r/m;

    upstream backend {
        server backend:8000;
    }

    server {
        listen 80;
        server_name _;

        # Let's Encrypt challenge
        location /.well-known/acme-challenge/ {
            root /var/www/certbot;
        }

        # Redirect to HTTPS (enable after SSL setup)
        # location / {
        #     return 301 https://$host$request_uri;
        # }

        # Frontend
        location / {
            root /usr/share/nginx/html;
            try_files $uri $uri/ /index.html;
        }

        # Backend API
        location /api/ {
            limit_req zone=api burst=20 nodelay;

            proxy_pass http://backend;
            proxy_set_header Host $host;
            proxy_set_header X-Real-IP $remote_addr;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Forwarded-Proto $scheme;

            # WebSocket support
            proxy_http_version 1.1;
            proxy_set_header Upgrade $http_upgrade;
            proxy_set_header Connection "upgrade";

            # Timeouts for batch processing
            proxy_connect_timeout 600s;
            proxy_send_timeout 600s;
            proxy_read_timeout 600s;
        }

        # Metrics (internal only)
        location /metrics {
            deny all;
        }
    }

    # HTTPS server (enable after SSL setup)
    # server {
    #     listen 443 ssl http2;
    #     server_name yourdomain.com;
    #
    #     ssl_certificate /etc/letsencrypt/live/yourdomain.com/fullchain.pem;
    #     ssl_certificate_key /etc/letsencrypt/live/yourdomain.com/privkey.pem;
    #     include /etc/nginx/ssl-params.conf;
    #
    #     # Same locations as HTTP server above
    # }
}
```

## Scaling

### Horizontal Scaling

**Backend instances:**
```bash
docker-compose -f docker-compose.prod.yml up -d --scale backend=3
```

Nginx load balances automatically.

**Celery workers:**
```yaml
# In docker-compose.prod.yml
celery-worker:
  deploy:
    replicas: 4  # Increase for higher batch throughput
```

### Vertical Scaling

**Increase worker concurrency:**
```yaml
celery-worker:
  command: celery -A app.celery_app worker --loglevel=info --concurrency=8
```

**Increase Gunicorn workers:**
```yaml
backend:
  command: uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 8
```

Rule of thumb: `workers = (CPU cores × 2) + 1`

## Backup and Restore

### Database Backup

```bash
# Backup
docker exec chemstructval-postgres pg_dump -U chemstructval chemstructval > backup_$(date +%Y%m%d).sql

# Automated daily backups (cron)
0 2 * * * cd /path/to/chemstructval && docker exec chemstructval-postgres pg_dump -U chemstructval chemstructval > backups/backup_$(date +\%Y\%m\%d).sql
```

### Database Restore

```bash
# Restore from backup
docker exec -i chemstructval-postgres psql -U chemstructval chemstructval < backup_20260124.sql
```

### Volume Backup

```bash
# Backup all volumes
docker run --rm \
  -v chemstructval_postgres_data:/data \
  -v $(pwd)/backups:/backup \
  alpine tar czf /backup/postgres_data.tar.gz /data

# Restore volume
docker run --rm \
  -v chemstructval_postgres_data:/data \
  -v $(pwd)/backups:/backup \
  alpine tar xzf /backup/postgres_data.tar.gz -C /
```

## Monitoring

### Container Health

```bash
# Check container status
docker-compose -f docker-compose.prod.yml ps

# Check resource usage
docker stats

# Check logs for errors
docker-compose -f docker-compose.prod.yml logs --tail=100 backend
```

### Application Health

```bash
# Health endpoint
curl http://localhost:8000/api/v1/health

# Expected response
{
  "status": "healthy",
  "database": "connected",
  "redis": "connected"
}
```

### Metrics

Access Prometheus: `http://your-server:9090`
Access Grafana: `http://your-server:3001`

See [Monitoring Guide](monitoring.md) for dashboard setup.

## Updates and Maintenance

### Updating Application

```bash
# Pull latest code
git pull origin main

# Rebuild images
docker-compose -f docker-compose.prod.yml build

# Restart services (zero-downtime)
docker-compose -f docker-compose.prod.yml up -d --no-deps backend

# Verify
docker-compose -f docker-compose.prod.yml ps
```

### Database Migrations

```bash
# Run Alembic migrations
docker-compose -f docker-compose.prod.yml exec backend alembic upgrade head
```

### Log Rotation

Configure Docker logging:

```yaml
# In docker-compose.prod.yml
services:
  backend:
    logging:
      driver: "json-file"
      options:
        max-size: "10m"
        max-file: "3"
```

## Security Considerations

### Environment Variables

- Never commit `.env` to git
- Use strong passwords (min 16 chars)
- Rotate `SECRET_KEY` periodically
- Use different credentials for each environment

### Network Security

- Use firewall to restrict port access
- Only expose port 80/443 publicly
- Keep other ports (5432, 6379, 9090) internal
- Use VPN for administrative access

### SSL/TLS

See [SSL Setup Guide](ssl-setup.md) for HTTPS configuration.

### Container Security

```bash
# Scan images for vulnerabilities
docker scan chemstructval-backend

# Update base images regularly
docker-compose -f docker-compose.prod.yml pull
docker-compose -f docker-compose.prod.yml up -d
```

## Troubleshooting

### Services won't start

```bash
# Check logs
docker-compose -f docker-compose.prod.yml logs

# Check system resources
docker info
df -h  # Disk space
free -h  # Memory
```

### Database connection failed

```bash
# Check PostgreSQL is running
docker-compose -f docker-compose.prod.yml ps postgres

# Check connection
docker exec -it chemstructval-postgres psql -U chemstructval -d chemstructval

# Verify DATABASE_URL in .env
```

### High memory usage

```bash
# Check container resource usage
docker stats

# Reduce Celery concurrency
# Edit docker-compose.prod.yml: --concurrency=2

# Add memory limits
services:
  backend:
    mem_limit: 2g
```

See [FAQ](../troubleshooting/faq.md) for more issues.

## Production Checklist

Before going live:

- [ ] Set strong passwords in `.env`
- [ ] Generate new `SECRET_KEY`
- [ ] Configure SSL certificates (see [SSL Setup](ssl-setup.md))
- [ ] Set up automated backups
- [ ] Configure monitoring alerts
- [ ] Test backup restoration
- [ ] Set up log aggregation
- [ ] Configure firewall rules
- [ ] Test health endpoints
- [ ] Perform load testing
- [ ] Document deployment procedures
- [ ] Set up status page (optional)

## Resources

- [Docker Documentation](https://docs.docker.com/)
- [Docker Compose Reference](https://docs.docker.com/compose/compose-file/)
- [Nginx Documentation](https://nginx.org/en/docs/)
- [PostgreSQL Docker](https://hub.docker.com/_/postgres)
- [Redis Docker](https://hub.docker.com/_/redis)

## Getting Help

- Check [FAQ](../troubleshooting/faq.md)
- Review container logs
- Check Docker daemon status
- Verify system resources
- Open GitHub issue with logs
