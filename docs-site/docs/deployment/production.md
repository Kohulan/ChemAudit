---
sidebar_position: 2
title: Production
description: Deploy ChemAudit to production with Nginx, SSL, and deployment profiles
---

# Production Deployment

Deploy ChemAudit to production using deployment profiles for configurable batch limits and worker counts.

## Deployment Profiles

ChemAudit includes pre-configured profiles that set batch limits, worker counts, and memory allocation:

| Profile | Max Batch | Max File Size | Celery Workers | Use Case |
|---------|-----------|---------------|----------------|----------|
| **small** | 1,000 | 100 MB | 2 | Development, testing |
| **medium** | 10,000 | 500 MB | 4 | Standard production |
| **large** | 50,000 | 500 MB | 8 | High-throughput labs |
| **xl** | 100,000 | 1 GB | 12 | Enterprise scale |
| **coconut** | 1,000,000 | 1 GB | 16 | Full COCONUT database |

## Quick Deploy

```bash
# Interactive deployment script
./deploy.sh

# Or specify profile directly
./deploy.sh medium
```

The script handles frontend building, environment configuration, and Docker Compose orchestration.

## Manual Production Deploy

### 1. Build Frontend

```bash
cd frontend
npm ci --production
npm run build
cd ..
cp -r frontend/dist frontend-dist
```

### 2. Configure Environment

```bash
# Option A: Automated (recommended) — deploy.sh generates all secrets
cp .env.prod.example .env
./deploy.sh medium  # Generates secure secrets, then deploys

# Option B: Manual
cp .env.prod.example .env
# Generate and set each secret:
openssl rand -hex 64   # For SECRET_KEY
openssl rand -hex 32   # For API_KEY_ADMIN_SECRET, CSRF_SECRET_KEY
openssl rand -hex 24   # For POSTGRES_PASSWORD, REDIS_PASSWORD, GRAFANA_PASSWORD
```

Edit `.env` with generated secrets.

:::danger Startup Validation
When `DEBUG=false`, the application rejects startup if `SECRET_KEY`, `API_KEY_ADMIN_SECRET`, or `CSRF_SECRET_KEY` still contain placeholder values. This prevents accidental deployment with insecure defaults.
:::

### 3. Deploy Services

```bash
# Start production stack
docker-compose -f docker-compose.prod.yml up -d

# Or with monitoring
docker-compose -f docker-compose.prod.yml --profile monitoring up -d

# Run database migrations
docker-compose -f docker-compose.prod.yml exec backend alembic upgrade head
```

**Access:** http://localhost (all services behind Nginx)

## SSL Configuration

### Let's Encrypt

```bash
# Install Certbot
sudo apt-get install certbot -y

# Stop nginx
docker-compose -f docker-compose.prod.yml stop nginx

# Obtain certificate
sudo certbot certonly --standalone \
  -d your-domain.com \
  -d www.your-domain.com

# Copy certificates
mkdir -p nginx/ssl
sudo cp /etc/letsencrypt/live/your-domain.com/fullchain.pem nginx/ssl/cert.pem
sudo cp /etc/letsencrypt/live/your-domain.com/privkey.pem nginx/ssl/key.pem
sudo chown $USER:$USER nginx/ssl/*.pem

# Generate DH parameters
openssl dhparam -out nginx/ssl/dhparam.pem 4096

# Restart nginx
docker-compose -f docker-compose.prod.yml start nginx
```

### Auto-Renewal

```bash
# Add to crontab
sudo crontab -e

# Add this line
0 3 * * * certbot renew --quiet --post-hook "docker-compose -f /path/to/docker-compose.prod.yml restart nginx"
```

## Scaling

### Horizontal Scaling

```bash
# Scale Celery workers
docker-compose -f docker-compose.prod.yml up -d --scale celery-worker=8

# Verify workers
docker-compose -f docker-compose.prod.yml exec backend \
  celery -A app.celery_app inspect active
```

### Vertical Scaling

Use larger deployment profile:

```bash
./deploy.sh xl  # 12 workers, 100K batch limit
```

## Maintenance

### Database Backup

```bash
# Create backup
docker-compose -f docker-compose.prod.yml exec postgres \
  pg_dump -U chemaudit chemaudit > backup_$(date +%Y%m%d).sql

# Restore backup
cat backup_20260204.sql | docker-compose -f docker-compose.prod.yml exec -T postgres \
  psql -U chemaudit chemaudit
```

### Application Updates

```bash
# Pull latest changes
git pull origin main

# Rebuild frontend
cd frontend && npm ci && npm run build && cd ..
cp -r frontend/dist frontend-dist

# Rebuild and restart
docker-compose -f docker-compose.prod.yml build
docker-compose -f docker-compose.prod.yml up -d

# Run migrations
docker-compose -f docker-compose.prod.yml exec backend alembic upgrade head
```

## Security Checklist

**Secrets & Authentication:**

- [ ] `SECRET_KEY`, `API_KEY_ADMIN_SECRET`, `CSRF_SECRET_KEY` generated (enforced at startup)
- [ ] `POSTGRES_PASSWORD` set to strong value
- [ ] `REDIS_PASSWORD` set to strong value (Redis auth enforced)
- [ ] `GRAFANA_PASSWORD` set to strong value
- [ ] `DEBUG=false` in production

**Network & SSL:**

- [ ] SSL certificate configured (see nginx/nginx-ssl.conf template)
- [ ] Firewall configured (ports 80, 443 only)
- [ ] Database not publicly accessible (bound to internal network)
- [ ] Redis authenticated and not publicly accessible

**Application Security (enabled by default):**

- [ ] CSRF protection active for browser sessions
- [ ] API rate limiting enabled
- [ ] Batch job ownership enforced (session-based)
- [ ] /metrics restricted to internal IPs
- [ ] OpenAPI docs disabled (`DEBUG=false`)
- [ ] Nginx security headers (CSP, Permissions-Policy, server_tokens off)
- [ ] Containers run as non-root with cap_drop ALL

**Operations:**

- [ ] Regular database backups configured
- [ ] Backup restoration tested
- [ ] Monitoring configured (Prometheus + Grafana)

## Next Steps

- **[Monitoring](/docs/deployment/monitoring)** - Set up Prometheus and Grafana
- **[Docker](/docs/deployment/docker)** - Development setup
- **[Troubleshooting](/docs/troubleshooting)** - Common issues
