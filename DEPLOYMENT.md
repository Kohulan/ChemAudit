# ChemStructVal Production Deployment Guide

## Overview

This guide covers production deployment of ChemStructVal using Docker Compose with Nginx reverse proxy, SSL/TLS encryption, and monitoring infrastructure.

## Architecture

```
Internet → Nginx (Port 80/443)
           ├─→ Frontend (React SPA served by Nginx)
           ├─→ Backend API (Gunicorn + Uvicorn workers)
           │   ├─→ PostgreSQL (Database)
           │   └─→ Redis (Cache + Message Broker)
           └─→ Celery Workers (Batch Processing)
```

## Prerequisites

- Docker Engine 24.0+
- Docker Compose 2.20+
- 4GB RAM minimum (8GB recommended)
- 20GB disk space
- Domain name (for production SSL)
- Ports 80 and 443 available

## Quick Start (Development-like Production)

For testing production build locally without SSL:

```bash
# 1. Build frontend
cd frontend
npm install
npm run build
cd ..

# 2. Copy frontend build to nginx directory
cp -r frontend/dist frontend-dist

# 3. Set environment variables
cp .env.prod.example .env
# Edit .env and change passwords

# 4. Start services
docker-compose -f docker-compose.prod.yml up -d

# 5. View logs
docker-compose -f docker-compose.prod.yml logs -f

# 6. Access application
# Open http://localhost
```

## Production Deployment with SSL

### Step 1: Server Setup

```bash
# Update system
sudo apt-get update && sudo apt-get upgrade -y

# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER

# Install Docker Compose
sudo curl -L "https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# Configure firewall
sudo ufw allow 22/tcp   # SSH
sudo ufw allow 80/tcp   # HTTP
sudo ufw allow 443/tcp  # HTTPS
sudo ufw enable
```

### Step 2: SSL Certificate Setup (Let's Encrypt)

```bash
# Install Certbot
sudo apt-get install certbot -y

# Obtain certificate (standalone mode - requires port 80 available)
sudo certbot certonly --standalone -d your-domain.com -d www.your-domain.com

# Certificates will be in:
# /etc/letsencrypt/live/your-domain.com/fullchain.pem
# /etc/letsencrypt/live/your-domain.com/privkey.pem

# Create nginx SSL directory
mkdir -p nginx/ssl

# Copy certificates
sudo cp /etc/letsencrypt/live/your-domain.com/fullchain.pem nginx/ssl/cert.pem
sudo cp /etc/letsencrypt/live/your-domain.com/privkey.pem nginx/ssl/key.pem

# Generate Diffie-Hellman parameters (takes several minutes)
sudo openssl dhparam -out nginx/ssl/dhparam.pem 4096

# Set certificate renewal cron job
sudo crontab -e
# Add: 0 3 * * * certbot renew --quiet && cp /etc/letsencrypt/live/your-domain.com/*.pem /path/to/nginx/ssl/
```

### Step 3: Configure Nginx for SSL

Edit `nginx/nginx.conf`:

```nginx
# Uncomment the HTTPS server block
server {
    listen 443 ssl http2;
    server_name your-domain.com;

    ssl_certificate /etc/nginx/ssl/cert.pem;
    ssl_certificate_key /etc/nginx/ssl/key.pem;

    include /etc/nginx/conf.d/ssl-params.conf;
    include /etc/nginx/conf.d/locations.conf;
}

# Enable HTTP to HTTPS redirect
server {
    listen 80;
    server_name your-domain.com;

    location = /health {
        access_log off;
        return 200 "healthy\n";
        add_header Content-Type text/plain;
    }

    location / {
        return 301 https://$host$request_uri;
    }
}
```

Update `docker-compose.prod.yml` to mount SSL directory:

```yaml
nginx:
  volumes:
    - ./nginx/ssl:/etc/nginx/ssl:ro  # Add this line
```

Edit `nginx/ssl-params.conf` and uncomment:

```nginx
ssl_dhparam /etc/nginx/ssl/dhparam.pem;
add_header Strict-Transport-Security "max-age=31536000; includeSubDomains; preload" always;
```

### Step 4: Environment Configuration

```bash
# Copy production environment template
cp .env.prod.example .env

# Edit with secure passwords
nano .env
```

**IMPORTANT**: Change these values in `.env`:
- `POSTGRES_PASSWORD` - Strong random password
- `GRAFANA_PASSWORD` - Strong random password
- Update `DATABASE_URL` with the new password

### Step 5: Build and Deploy

```bash
# Clone repository
git clone https://github.com/yourusername/chemstructval.git
cd chemstructval

# Build frontend
cd frontend
npm install
npm run build
cd ..

# Copy to nginx directory
cp -r frontend/dist frontend-dist

# Start services (without monitoring)
docker-compose -f docker-compose.prod.yml up -d

# Or with monitoring
docker-compose -f docker-compose.prod.yml --profile monitoring up -d

# Check logs
docker-compose -f docker-compose.prod.yml logs -f

# Verify all services are healthy
docker-compose -f docker-compose.prod.yml ps
```

### Step 6: Database Initialization

```bash
# Run database migrations
docker-compose -f docker-compose.prod.yml exec backend alembic upgrade head

# Verify database connection
docker-compose -f docker-compose.prod.yml exec backend python -c "from app.database import engine; print('Database connected')"
```

### Step 7: Verify Deployment

Test each component:

```bash
# Health check
curl https://your-domain.com/api/v1/health

# Test validation API
curl -X POST https://your-domain.com/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{"structure": "CCO", "format": "smiles"}'

# Check SSL rating
# Visit: https://www.ssllabs.com/ssltest/analyze.html?d=your-domain.com
```

## Monitoring

Access monitoring services (if started with `--profile monitoring`):

- **Application**: https://your-domain.com
- **Grafana**: http://your-domain.com:3001 (configure reverse proxy for HTTPS)
- **Prometheus**: http://your-domain.com:9090 (restrict access in production)

### Restrict Monitoring Access

Edit `nginx/locations.conf` to restrict Prometheus metrics:

```nginx
location /metrics {
    allow 10.0.0.0/8;      # Internal network
    allow YOUR_IP_HERE;    # Your IP
    deny all;

    proxy_pass http://backend;
}
```

## Maintenance

### Viewing Logs

```bash
# All services
docker-compose -f docker-compose.prod.yml logs -f

# Specific service
docker-compose -f docker-compose.prod.yml logs -f backend

# Last 100 lines
docker-compose -f docker-compose.prod.yml logs --tail=100 backend
```

### Updating Application

```bash
# Pull latest code
git pull

# Rebuild frontend
cd frontend
npm install
npm run build
cd ..
cp -r frontend/dist frontend-dist

# Rebuild and restart services
docker-compose -f docker-compose.prod.yml build
docker-compose -f docker-compose.prod.yml up -d

# Run migrations if needed
docker-compose -f docker-compose.prod.yml exec backend alembic upgrade head
```

### Database Backup

```bash
# Backup database
docker-compose -f docker-compose.prod.yml exec postgres pg_dump -U chemstructval chemstructval > backup_$(date +%Y%m%d).sql

# Restore database
cat backup_20260124.sql | docker-compose -f docker-compose.prod.yml exec -T postgres psql -U chemstructval chemstructval
```

### Scaling Celery Workers

```bash
# Scale to 8 workers
docker-compose -f docker-compose.prod.yml up -d --scale celery-worker=8

# Verify workers
docker-compose -f docker-compose.prod.yml exec backend celery -A app.celery_app inspect active
```

## Performance Tuning

### Gunicorn Workers

Adjust in `backend/Dockerfile.prod`:

```dockerfile
CMD ["gunicorn", "app.main:app", \
    "--workers", "8", \  # 2-4 × CPU cores
    ...
]
```

### PostgreSQL

Edit `docker-compose.prod.yml` to add performance settings:

```yaml
postgres:
  command:
    - postgres
    - -c
    - max_connections=200
    - -c
    - shared_buffers=256MB
    - -c
    - effective_cache_size=1GB
```

### Redis Memory

Adjust in `docker-compose.prod.yml`:

```yaml
redis:
  command: redis-server --maxmemory 1gb --maxmemory-policy allkeys-lru
```

## Troubleshooting

### Service Won't Start

```bash
# Check service status
docker-compose -f docker-compose.prod.yml ps

# View detailed logs
docker-compose -f docker-compose.prod.yml logs backend

# Check health status
docker inspect chemstructval-backend-prod | grep -A 10 Health
```

### Database Connection Issues

```bash
# Verify PostgreSQL is running
docker-compose -f docker-compose.prod.yml exec postgres pg_isready

# Check connection from backend
docker-compose -f docker-compose.prod.yml exec backend python -c "from app.database import engine; print(engine.url)"

# View PostgreSQL logs
docker-compose -f docker-compose.prod.yml logs postgres
```

### Nginx Configuration Errors

```bash
# Test nginx configuration
docker-compose -f docker-compose.prod.yml exec nginx nginx -t

# Reload nginx
docker-compose -f docker-compose.prod.yml exec nginx nginx -s reload
```

### SSL Certificate Issues

```bash
# Verify certificate files exist
ls -la nginx/ssl/

# Check certificate expiration
openssl x509 -in nginx/ssl/cert.pem -noout -dates

# Test SSL configuration
openssl s_client -connect your-domain.com:443 -servername your-domain.com
```

## Security Checklist

- [ ] Changed all default passwords
- [ ] SSL/TLS certificates configured
- [ ] HSTS enabled in nginx
- [ ] Firewall configured (only 80, 443 exposed)
- [ ] Database not exposed to public internet
- [ ] Redis not exposed to public internet
- [ ] Prometheus metrics access restricted
- [ ] Grafana access secured
- [ ] Regular backups configured
- [ ] Log monitoring set up
- [ ] Security headers configured in nginx
- [ ] Rate limiting enabled
- [ ] Non-root users in Docker containers

## High Availability Considerations

For production at scale:

1. **Load Balancer**: Use AWS ALB, Google Cloud Load Balancer, or HAProxy
2. **Database**: Use managed PostgreSQL (AWS RDS, Google Cloud SQL)
3. **Redis**: Use managed Redis (AWS ElastiCache, Redis Cloud)
4. **Container Orchestration**: Consider Kubernetes for multi-node deployments
5. **CDN**: Use CloudFlare or AWS CloudFront for static assets
6. **Monitoring**: Use DataDog, New Relic, or Grafana Cloud
7. **Logging**: Centralized logging with ELK stack or Grafana Loki

## Cost Optimization

- Use reserved instances for long-term deployments
- Enable auto-scaling for Celery workers based on queue depth
- Use spot instances for batch processing workers
- Implement cache warming to reduce database load
- Optimize Docker image sizes with multi-stage builds
- Use CDN for static assets to reduce bandwidth

## Support

For issues and questions:
- GitHub Issues: https://github.com/yourusername/chemstructval/issues
- Documentation: https://github.com/yourusername/chemstructval/docs

## License

[Your License Here]
