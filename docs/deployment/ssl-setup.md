# SSL Certificate Setup with Let's Encrypt

Secure your ChemStructVal deployment with free SSL certificates from Let's Encrypt.

## Overview

This guide walks through:
- Setting up Let's Encrypt with Certbot
- Configuring Nginx for HTTPS
- Enabling HTTP to HTTPS redirect
- Automating certificate renewal
- Security best practices

## Prerequisites

- ChemStructVal deployed with Docker Compose
- Domain name pointing to your server's IP address
- Ports 80 and 443 open and accessible
- Root or sudo access to server

### Verify DNS

Before proceeding, verify your domain points to your server:

```bash
# Check DNS resolution
nslookup yourdomain.com

# Should return your server's IP address
```

## Setup Process

### Step 1: Create Certbot Directories

```bash
cd /path/to/chemstructval
mkdir -p certbot/conf certbot/www
```

These directories will hold:
- `certbot/conf`: SSL certificates and Let's Encrypt configuration
- `certbot/www`: Files for HTTP-01 challenge verification

### Step 2: Configure Nginx for HTTP Challenge

Ensure your `nginx/nginx.conf` includes the ACME challenge location:

```nginx
server {
    listen 80;
    server_name yourdomain.com;

    # Let's Encrypt challenge
    location /.well-known/acme-challenge/ {
        root /var/www/certbot;
    }

    # Temporary: serve app over HTTP
    location / {
        root /usr/share/nginx/html;
        try_files $uri $uri/ /index.html;
    }

    location /api/ {
        proxy_pass http://backend:8000;
        # ... other proxy settings
    }
}
```

### Step 3: Start Nginx

```bash
docker-compose -f docker-compose.prod.yml up -d nginx
```

Verify Nginx is accessible:

```bash
curl http://yourdomain.com
```

### Step 4: Request Certificate

Use Certbot in a Docker container to request the certificate:

```bash
docker run --rm \
  -v $(pwd)/certbot/conf:/etc/letsencrypt \
  -v $(pwd)/certbot/www:/var/www/certbot \
  certbot/certbot certonly \
  --webroot \
  --webroot-path=/var/www/certbot \
  --email your@email.com \
  --agree-tos \
  --no-eff-email \
  -d yourdomain.com
```

**Parameters explained:**
- `--webroot`: Use webroot plugin (doesn't need to stop Nginx)
- `--webroot-path`: Where to place challenge files
- `--email`: For renewal and security notices
- `--agree-tos`: Accept Let's Encrypt terms of service
- `--no-eff-email`: Don't share email with EFF
- `-d`: Domain name (can specify multiple with additional `-d` flags)

**For multiple domains:**
```bash
docker run --rm \
  -v $(pwd)/certbot/conf:/etc/letsencrypt \
  -v $(pwd)/certbot/www:/var/www/certbot \
  certbot/certbot certonly \
  --webroot \
  --webroot-path=/var/www/certbot \
  --email your@email.com \
  --agree-tos \
  --no-eff-email \
  -d yourdomain.com \
  -d www.yourdomain.com \
  -d api.yourdomain.com
```

**Expected output:**
```
Successfully received certificate.
Certificate is saved at: /etc/letsencrypt/live/yourdomain.com/fullchain.pem
Key is saved at:         /etc/letsencrypt/live/yourdomain.com/privkey.pem
```

### Step 5: Verify Certificate Files

```bash
ls -la certbot/conf/live/yourdomain.com/

# Should show:
# fullchain.pem (certificate + intermediate)
# privkey.pem (private key)
# cert.pem (certificate only)
# chain.pem (intermediate certificate)
```

### Step 6: Configure SSL Parameters

Create `nginx/ssl-params.conf` with security headers:

```nginx
# SSL protocols and ciphers
ssl_protocols TLSv1.2 TLSv1.3;
ssl_prefer_server_ciphers on;
ssl_ciphers ECDHE-ECDSA-AES128-GCM-SHA256:ECDHE-RSA-AES128-GCM-SHA256:ECDHE-ECDSA-AES256-GCM-SHA384:ECDHE-RSA-AES256-GCM-SHA384;

# SSL session settings
ssl_session_cache shared:SSL:10m;
ssl_session_timeout 10m;
ssl_session_tickets off;

# OCSP Stapling
ssl_stapling on;
ssl_stapling_verify on;
ssl_trusted_certificate /etc/letsencrypt/live/yourdomain.com/chain.pem;
resolver 8.8.8.8 8.8.4.4 valid=300s;
resolver_timeout 5s;

# Security headers
add_header Strict-Transport-Security "max-age=31536000; includeSubDomains" always;
add_header X-Frame-Options "SAMEORIGIN" always;
add_header X-Content-Type-Options "nosniff" always;
add_header X-XSS-Protection "1; mode=block" always;
add_header Referrer-Policy "no-referrer-when-downgrade" always;
```

### Step 7: Update Nginx Configuration for HTTPS

Edit `nginx/nginx.conf` to add HTTPS server block:

```nginx
events {
    worker_connections 1024;
}

http {
    include /etc/nginx/mime.types;
    default_type application/octet-stream;

    # ... other settings ...

    upstream backend {
        server backend:8000;
    }

    # HTTP server - redirect to HTTPS
    server {
        listen 80;
        server_name yourdomain.com;

        # Let's Encrypt challenge
        location /.well-known/acme-challenge/ {
            root /var/www/certbot;
        }

        # Redirect all other traffic to HTTPS
        location / {
            return 301 https://$host$request_uri;
        }
    }

    # HTTPS server
    server {
        listen 443 ssl http2;
        server_name yourdomain.com;

        # SSL certificates
        ssl_certificate /etc/letsencrypt/live/yourdomain.com/fullchain.pem;
        ssl_certificate_key /etc/letsencrypt/live/yourdomain.com/privkey.pem;

        # SSL parameters
        include /etc/nginx/ssl-params.conf;

        # Gzip compression
        gzip on;
        gzip_types text/plain text/css application/json application/javascript text/xml application/xml application/xml+rss text/javascript;

        # Rate limiting
        limit_req_zone $binary_remote_addr zone=api:10m rate=100r/m;

        # Frontend
        location / {
            root /usr/share/nginx/html;
            try_files $uri $uri/ /index.html;

            # Cache static assets
            location ~* \.(js|css|png|jpg|jpeg|gif|ico|svg|woff|woff2|ttf|eot)$ {
                expires 1y;
                add_header Cache-Control "public, immutable";
            }
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
}
```

### Step 8: Update Docker Compose

Ensure `docker-compose.prod.yml` mounts SSL certificates:

```yaml
nginx:
  image: nginx:alpine
  container_name: chemstructval-nginx
  volumes:
    - ./nginx/nginx.conf:/etc/nginx/nginx.conf:ro
    - ./nginx/ssl-params.conf:/etc/nginx/ssl-params.conf:ro
    - ./frontend/dist:/usr/share/nginx/html:ro
    - ./certbot/conf:/etc/letsencrypt:ro
    - ./certbot/www:/var/www/certbot:ro
  ports:
    - "80:80"
    - "443:443"
  depends_on:
    - backend
  restart: unless-stopped
```

### Step 9: Restart Nginx

```bash
docker-compose -f docker-compose.prod.yml restart nginx
```

### Step 10: Verify HTTPS

```bash
# Test HTTPS connection
curl -I https://yourdomain.com

# Expected: HTTP/2 200 OK

# Test HTTP redirect
curl -I http://yourdomain.com

# Expected: 301 Moved Permanently
# Location: https://yourdomain.com/
```

## Certificate Renewal

Let's Encrypt certificates expire after **90 days**. Set up automatic renewal.

### Manual Renewal

```bash
docker run --rm \
  -v $(pwd)/certbot/conf:/etc/letsencrypt \
  -v $(pwd)/certbot/www:/var/www/certbot \
  certbot/certbot renew
```

Reload Nginx after renewal:

```bash
docker-compose -f docker-compose.prod.yml restart nginx
```

### Automated Renewal with Cron

Create renewal script `scripts/renew-certs.sh`:

```bash
#!/bin/bash
cd /path/to/chemstructval

# Renew certificates
docker run --rm \
  -v $(pwd)/certbot/conf:/etc/letsencrypt \
  -v $(pwd)/certbot/www:/var/www/certbot \
  certbot/certbot renew --quiet

# Reload Nginx if renewal occurred
docker-compose -f docker-compose.prod.yml restart nginx

# Log renewal
echo "Certificate renewal checked: $(date)" >> logs/cert-renewal.log
```

Make executable:

```bash
chmod +x scripts/renew-certs.sh
```

Add to crontab:

```bash
crontab -e

# Add this line (runs daily at 2 AM)
0 2 * * * /path/to/chemstructval/scripts/renew-certs.sh
```

### Test Renewal (Dry Run)

```bash
docker run --rm \
  -v $(pwd)/certbot/conf:/etc/letsencrypt \
  -v $(pwd)/certbot/www:/var/www/certbot \
  certbot/certbot renew --dry-run
```

Expected output:
```
Congratulations, all simulated renewals succeeded:
  /etc/letsencrypt/live/yourdomain.com/fullchain.pem (success)
```

## Monitoring Certificate Expiry

### Check Expiration Date

```bash
# Check certificate expiry
openssl x509 -in certbot/conf/live/yourdomain.com/cert.pem -noout -enddate

# Example output:
# notAfter=Apr 24 10:30:00 2026 GMT
```

### Set Up Expiry Alerts

Create monitoring script `scripts/check-cert-expiry.sh`:

```bash
#!/bin/bash
DOMAIN="yourdomain.com"
CERT_FILE="/path/to/chemstructval/certbot/conf/live/$DOMAIN/cert.pem"
DAYS_WARN=14

# Get expiry date
EXPIRY_DATE=$(openssl x509 -in "$CERT_FILE" -noout -enddate | cut -d= -f2)
EXPIRY_EPOCH=$(date -d "$EXPIRY_DATE" +%s)
NOW_EPOCH=$(date +%s)
DAYS_LEFT=$(( ($EXPIRY_EPOCH - $NOW_EPOCH) / 86400 ))

if [ $DAYS_LEFT -le $DAYS_WARN ]; then
    echo "WARNING: SSL certificate for $DOMAIN expires in $DAYS_LEFT days!"
    # Send alert (email, Slack, etc.)
fi
```

Add to crontab (runs daily):

```bash
0 9 * * * /path/to/chemstructval/scripts/check-cert-expiry.sh
```

## Wildcard Certificates

For subdomains, use wildcard certificate:

```bash
docker run --rm \
  -v $(pwd)/certbot/conf:/etc/letsencrypt \
  certbot/certbot certonly \
  --manual \
  --preferred-challenges dns \
  --email your@email.com \
  --agree-tos \
  -d "*.yourdomain.com" \
  -d yourdomain.com
```

**Note:** DNS-01 challenge requires manual DNS TXT record creation or API automation.

## Security Best Practices

### SSL/TLS Configuration

✅ **Do:**
- Use TLS 1.2 and 1.3 only
- Enable HSTS with long max-age
- Enable OCSP stapling
- Use strong cipher suites
- Disable SSL compression

❌ **Don't:**
- Use SSLv3 or TLS 1.0/1.1 (deprecated)
- Allow weak ciphers
- Skip HSTS header
- Use self-signed certificates in production

### Certificate Management

✅ **Do:**
- Set up automated renewal
- Monitor expiration dates
- Test renewal process regularly
- Keep backups of certificates
- Use strong permissions on private keys

```bash
# Secure private key permissions
chmod 600 certbot/conf/live/yourdomain.com/privkey.pem
```

❌ **Don't:**
- Commit certificates to git
- Share private keys
- Skip renewal testing
- Ignore expiry warnings

### Nginx Security

✅ **Do:**
- Keep Nginx updated
- Enable security headers
- Set appropriate timeouts
- Configure rate limiting
- Use HTTP/2

❌ **Don't:**
- Expose internal ports
- Disable security headers
- Allow unlimited connections
- Skip access logging

## Testing SSL Configuration

### OpenSSL Test

```bash
# Test SSL connection
openssl s_client -connect yourdomain.com:443 -servername yourdomain.com

# Check certificate chain
openssl s_client -connect yourdomain.com:443 -showcerts

# Test specific TLS version
openssl s_client -connect yourdomain.com:443 -tls1_2
openssl s_client -connect yourdomain.com:443 -tls1_3
```

### Online Testing

**SSL Labs Test:**
Visit: https://www.ssllabs.com/ssltest/analyze.html?d=yourdomain.com

Target grade: **A** or **A+**

**SecurityHeaders.com:**
Visit: https://securityheaders.com/?q=yourdomain.com

Target grade: **A** or better

## Troubleshooting

### Certificate Request Failed

**Error:** Challenge validation failed

**Solutions:**
1. Verify domain DNS points to server
2. Check port 80 is accessible
3. Ensure `/.well-known/acme-challenge/` is reachable
4. Check Nginx logs: `docker-compose logs nginx`

**Test challenge directory:**
```bash
# Create test file
echo "test" > certbot/www/test.txt

# Access via HTTP
curl http://yourdomain.com/.well-known/acme-challenge/test.txt
```

### Certificate Not Trusted

**Error:** "Your connection is not private" or "NET::ERR_CERT_AUTHORITY_INVALID"

**Solutions:**
1. Verify `fullchain.pem` is used (not `cert.pem`)
2. Check intermediate certificate is included
3. Clear browser cache
4. Verify certificate matches domain

### Renewal Failed

**Error:** Unable to renew certificate

**Solutions:**
1. Check disk space: `df -h`
2. Verify certbot directory permissions
3. Test with `--dry-run`
4. Check certbot logs: `certbot/conf/letsencrypt.log`
5. Manually renew with verbose output: `certbot renew --verbose`

### Mixed Content Warnings

**Error:** Browser blocks HTTP content on HTTPS page

**Solutions:**
1. Update frontend API URL to use HTTPS
2. Check for hardcoded HTTP URLs in HTML/JS
3. Use protocol-relative URLs: `//yourdomain.com/api/`
4. Enable Content-Security-Policy header

## Advanced Configurations

### HSTS Preloading

For maximum security, submit to HSTS preload list:

```nginx
add_header Strict-Transport-Security "max-age=63072000; includeSubDomains; preload" always;
```

Submit at: https://hstspreload.org/

### Certificate Pinning

For high-security applications:

```nginx
# Get public key hash
openssl x509 -in certbot/conf/live/yourdomain.com/cert.pem -pubkey -noout | \
  openssl pkey -pubin -outform der | \
  openssl dgst -sha256 -binary | \
  base64

# Add to Nginx
add_header Public-Key-Pins 'pin-sha256="base64hash"; max-age=5184000' always;
```

### Custom DH Parameters

```bash
# Generate strong DH parameters (takes several minutes)
openssl dhparam -out nginx/dhparam.pem 4096

# Add to ssl-params.conf
ssl_dhparam /etc/nginx/dhparam.pem;
```

## Resources

- [Let's Encrypt Documentation](https://letsencrypt.org/docs/)
- [Certbot Documentation](https://certbot.eff.org/docs/)
- [Mozilla SSL Configuration Generator](https://ssl-config.mozilla.org/)
- [SSL Labs](https://www.ssllabs.com/)
- [Security Headers](https://securityheaders.com/)

## Getting Help

- Check [FAQ](../troubleshooting/faq.md)
- Review Certbot logs
- Test with `--dry-run`
- Verify DNS configuration
- Check Nginx logs
- Open GitHub issue with details
