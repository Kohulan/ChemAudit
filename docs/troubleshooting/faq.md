# Troubleshooting FAQ

Common issues and solutions for ChemStructVal.

## Table of Contents

- [Installation Issues](#installation-issues)
- [Validation Issues](#validation-issues)
- [Batch Processing Issues](#batch-processing-issues)
- [Performance Issues](#performance-issues)
- [API Issues](#api-issues)
- [Monitoring Issues](#monitoring-issues)
- [Frontend Issues](#frontend-issues)
- [Database Issues](#database-issues)

---

## Installation Issues

### Docker fails to start services

**Symptom:** `docker-compose up` fails with connection errors or services crash immediately.

**Solutions:**

1. **Ensure Docker daemon is running:**
   ```bash
   docker info
   # Should show system info, not "Cannot connect to Docker daemon"
   ```

2. **Check available memory:**
   ```bash
   docker stats
   # ChemStructVal needs minimum 4GB RAM
   ```

3. **Remove existing containers and volumes:**
   ```bash
   docker-compose down -v
   docker-compose build --no-cache
   docker-compose up -d
   ```

4. **Check logs for specific errors:**
   ```bash
   docker-compose logs backend
   docker-compose logs postgres
   ```

5. **Verify Docker Compose version:**
   ```bash
   docker-compose version
   # Should be 2.20.0 or higher
   ```

### Port already in use

**Symptom:** `Error: bind: address already in use` when starting services.

**Common ports:**
- 3000: Frontend
- 8000: Backend
- 5432: PostgreSQL
- 6379: Redis
- 9090: Prometheus
- 3001: Grafana

**Solutions:**

1. **Find process using port:**
   ```bash
   # macOS/Linux
   lsof -i :8000

   # Windows
   netstat -ano | findstr :8000
   ```

2. **Kill the process:**
   ```bash
   # macOS/Linux
   kill -9 <PID>

   # Windows
   taskkill /PID <PID> /F
   ```

3. **Change port in docker-compose.yml:**
   ```yaml
   backend:
     ports:
       - "8001:8000"  # Changed from 8000:8000
   ```

### RDKit.js fails to load

**Symptom:** "Failed to load RDKit.js" error in browser console.

**Solutions:**

1. **Clear browser cache:**
   - Chrome: Ctrl+Shift+Delete → Clear browsing data
   - Firefox: Ctrl+Shift+Delete → Clear Data

2. **Check browser console for WASM errors:**
   ```
   F12 → Console tab
   Look for: "wasm streaming compile failed"
   ```

3. **Verify RDKit WASM file exists:**
   ```bash
   ls frontend/node_modules/@rdkit/rdkit/dist/RDKit_minimal.wasm
   # Should exist and be ~15MB
   ```

4. **Check memory:**
   - RDKit WASM needs ~100MB
   - Close other browser tabs
   - Check Chrome Task Manager (Shift+Esc)

5. **Try different browser:**
   - WASM support required
   - Chrome 89+, Firefox 78+, Safari 14+

### npm install fails

**Symptom:** Errors during `npm install` in frontend.

**Solutions:**

1. **Clear npm cache:**
   ```bash
   npm cache clean --force
   rm -rf node_modules package-lock.json
   npm install
   ```

2. **Use correct Node version:**
   ```bash
   node --version
   # Should be 18.x or 20.x

   # Use nvm to switch versions
   nvm install 20
   nvm use 20
   ```

3. **Check disk space:**
   ```bash
   df -h
   # Need at least 1GB free
   ```

### Python dependencies fail to install

**Symptom:** Errors during `poetry install` or `pip install`.

**Solutions:**

1. **Use correct Python version:**
   ```bash
   python --version
   # Should be 3.11 or higher
   ```

2. **Install system dependencies for RDKit:**
   ```bash
   # Ubuntu/Debian
   sudo apt-get install -y libboost-all-dev libcairo2-dev

   # macOS
   brew install boost cairo
   ```

3. **Clear Poetry cache:**
   ```bash
   poetry cache clear pypi --all
   poetry install
   ```

---

## Validation Issues

### "Invalid SMILES" for valid molecule

**Symptom:** Known valid SMILES is rejected as invalid.

**Possible causes:**

1. **Encoding issues:**
   ```bash
   # Check file encoding
   file -I input.csv
   # Should show: charset=utf-8

   # Convert to UTF-8
   iconv -f ISO-8859-1 -t UTF-8 input.csv > output.csv
   ```

2. **Leading/trailing whitespace:**
   ```python
   # Trim whitespace
   smiles = smiles.strip()
   ```

3. **Invisible characters:**
   ```bash
   # Show all characters
   cat -A input.csv | head
   # Look for ^M (carriage return) or other oddities
   ```

4. **RDKit version compatibility:**
   ```bash
   # Check RDKit version
   docker exec chemstructval-backend python -c "import rdkit; print(rdkit.__version__)"
   # Should be 2024.3.0+
   ```

**Debug steps:**

```python
from rdkit import Chem

# Test parsing
smiles = "YOUR_SMILES"
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    print(f"Failed to parse: {smiles}")
else:
    print(f"Canonical: {Chem.MolToSmiles(mol)}")
```

### Stereochemistry warnings on intended racemates

**Symptom:** "Undefined stereocenter" warning for molecules where stereochemistry is intentionally unspecified.

**Explanation:** This is **expected behavior** for racemates.

**Interpretation:**
- WARNING severity = informational, not an error
- Undefined stereo is intentional for racemates
- Different from ERROR which indicates structural problems

**Actions:**

1. **For racemates:** Accept the warning, document as intentional
2. **For single enantiomers:** Add stereochemistry with `@` or `@@`
3. **For ML datasets:** Consider enumerating stereoisomers

**Example:**
```
Input: CC(O)N
Warning: "Undefined stereocenter at atom 1"

Options:
- Keep as racemate: CC(O)N (both R and S)
- Specify S: C[C@H](O)N
- Specify R: C[C@@H](O)N
```

### Standardization removes desired stereochemistry

**Symptom:** E/Z double bond stereochemistry is lost after standardization.

**Cause:** Tautomer canonicalization can remove E/Z stereo.

**Solution:** Tautomer canonicalization is **OFF by default** to protect stereochemistry.

If you enabled it and seeing issues:

```python
# Check standardization options
standardization_options = {
    "remove_salts": True,
    "neutralize": True,
    "canonicalize_tautomers": False,  # Should be False
}
```

**Workaround if tautomers needed:**
1. Run standardization without tautomer canonicalization
2. Manually review tautomer changes
3. Preserve stereo where critical

### PAINS alert on known safe compound

**Symptom:** PAINS warning on FDA-approved drug or well-known compound.

**Explanation:** PAINS are **warnings, not disqualifications**.

**Context:**
- 87 FDA-approved drugs contain PAINS patterns
- Context matters more than pattern presence
- Different assays have different sensitivities

**Actions:**

1. **Review the specific pattern:**
   ```json
   {
     "alert": "PAINS_B",
     "pattern": "catechol",
     "affected_atoms": [3, 4]
   }
   ```

2. **Consider your use case:**
   - High-throughput screening: May want to filter
   - Lead optimization: Review in context
   - Known drugs: Usually acceptable

3. **Document decision:**
   - Note why PAINS alert is acceptable
   - Reference literature if available
   - Include in compound metadata

**Example:** L-DOPA contains catechol (PAINS_B) but is FDA-approved for Parkinson's disease.

### Multiple fragments detected

**Symptom:** "Connectivity warning: 2 fragments detected"

**Causes:**
- Salt form: `CCO.Cl` (ethanol + chloride)
- Counter ions: `[Na+].CC(=O)[O-]` (sodium acetate)
- Mixtures: `CCO.CCCC` (ethanol + butane)

**When it's acceptable:**
- Intentional salts (HCl, Na+, etc.)
- Ionic compounds
- Coordination complexes

**When to fix:**
- Unintended mixtures
- Data entry errors
- ML training data (usually want neutral forms)

**Solutions:**

1. **Keep salt form (database storage):**
   - Document as salt
   - Acceptable for most use cases

2. **Strip salts (ML workflows):**
   ```python
   # Use standardization with salt stripping
   result = client.standardize(smiles, options={"remove_salts": True})
   ```

3. **Manual review:**
   - Check if both fragments are intentional
   - Verify data source

---

## Batch Processing Issues

### Batch stuck at 0%

**Symptom:** Progress bar doesn't move, stays at 0% indefinitely.

**Solutions:**

1. **Check Celery worker is running:**
   ```bash
   docker-compose logs celery-worker
   # Should show: "celery@worker ready"
   ```

2. **Check Redis connection:**
   ```bash
   docker exec chemstructval-redis redis-cli ping
   # Should return: PONG
   ```

3. **Verify batch was uploaded:**
   ```bash
   docker-compose logs backend | grep "Batch job"
   # Should show: "Batch job created: <job_id>"
   ```

4. **Check Celery task queue:**
   ```bash
   docker exec chemstructval-redis redis-cli llen celery
   # Should show task count
   ```

5. **Restart Celery worker:**
   ```bash
   docker-compose restart celery-worker
   docker-compose logs -f celery-worker
   ```

### Individual molecules fail in batch

**Symptom:** Some molecules show "failed" status in results.

**Explanation:** This is **expected and normal** for invalid structures.

**Behavior:**
- Individual failures don't stop the batch
- Failed molecules are marked with status="failed"
- Error message explains the specific issue
- Other molecules continue processing

**Common failure reasons:**

1. **Parse errors:**
   ```
   Molecule 42: "Invalid SMILES syntax"
   Error: "Unmatched parenthesis"
   ```

2. **Encoding issues:**
   ```
   Molecule 108: "Cannot parse input"
   Error: "Invalid UTF-8 sequence"
   ```

3. **Empty rows:**
   ```
   Molecule 256: "Missing SMILES"
   Error: "Empty input"
   ```

**Actions:**

1. **Export results to see error messages:**
   ```bash
   # Download CSV with error column
   # Filter for status="failed"
   ```

2. **Review first 10 failures:**
   - Look for patterns
   - Check common error types
   - Fix source data if systematic

3. **Re-process failed molecules:**
   - Extract failed rows
   - Fix issues
   - Upload as new batch

### WebSocket disconnects

**Symptom:** Progress updates stop, but batch continues processing.

**Causes:**
- Long-running batches (>10 minutes)
- Network issues
- Nginx WebSocket timeout
- Browser tab inactive

**Impact:**
- Progress bar freezes
- Batch continues on server
- Results are not lost

**Solutions:**

1. **Refresh page to reconnect:**
   ```
   F5 or Ctrl+R
   Progress will resume from current point
   ```

2. **Check Nginx WebSocket timeout:**
   ```nginx
   # In nginx.conf
   location /api/ {
       proxy_read_timeout 86400s;  # 24 hours
   }
   ```

3. **Monitor in backend logs:**
   ```bash
   docker-compose logs -f celery-worker
   # Shows actual progress
   ```

4. **Check job status via API:**
   ```bash
   curl http://localhost:8000/api/v1/batch/<job_id>/status
   ```

### File upload fails

**Symptom:** "Upload failed" error when submitting batch.

**Solutions:**

1. **Check file size limits:**
   ```
   CSV: Maximum 50MB
   SDF: Maximum 100MB
   ```

2. **Verify file format:**
   ```bash
   # Check CSV structure
   head input.csv

   # Check for BOM (Byte Order Mark)
   xxd input.csv | head
   # Should not start with: efbbbf (UTF-8 BOM)
   ```

3. **Check Nginx upload limits:**
   ```nginx
   # In nginx.conf
   client_max_body_size 100M;
   ```

4. **Check backend logs:**
   ```bash
   docker-compose logs backend | grep -i upload
   ```

---

## Performance Issues

### Validation takes too long

**Symptom:** Single molecule validation takes >10 seconds.

**Expected performance:**
- Simple molecules: <100ms
- Complex molecules: 100-500ms
- Very large molecules (>100 atoms): 1-3s

**If slower than expected:**

1. **Check molecule size:**
   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles(smiles)
   print(f"Atoms: {mol.GetNumAtoms()}")
   print(f"Bonds: {mol.GetNumBonds()}")
   # >200 atoms will be slow
   ```

2. **Verify Redis is running (for caching):**
   ```bash
   docker-compose ps redis
   curl http://localhost:8000/metrics | grep cache_hits
   ```

3. **Check Prometheus metrics for bottlenecks:**
   ```
   Visit: http://localhost:3001
   Dashboard: ChemStructVal Overview
   Panel: Validation Duration
   ```

4. **Profile with detailed logging:**
   ```bash
   # Enable debug logging
   docker-compose exec backend python -c "
   import logging
   logging.basicConfig(level=logging.DEBUG)
   "
   ```

5. **Check system resources:**
   ```bash
   docker stats
   # Look for high CPU or memory usage
   ```

### High memory usage

**Symptom:** Container killed with OOM (Out Of Memory) error.

**Solutions:**

1. **Check container resource usage:**
   ```bash
   docker stats --no-stream
   ```

2. **Increase container memory limits:**
   ```yaml
   # In docker-compose.prod.yml
   backend:
     mem_limit: 4g
     mem_reservation: 2g
   ```

3. **Reduce Celery concurrency:**
   ```yaml
   celery-worker:
     command: celery -A app.celery_app worker --concurrency=2
   ```

4. **Check for memory leaks:**
   ```bash
   # Monitor over time
   watch -n 5 docker stats --no-stream
   # Memory should stabilize, not grow indefinitely
   ```

5. **Clear Redis cache:**
   ```bash
   docker exec chemstructval-redis redis-cli FLUSHDB
   ```

### Batch processing is slow

**Symptom:** Batch of 1000 molecules takes >2 minutes.

**Expected performance:**
- Target: >100 molecules/second
- 1,000 molecules: <10 seconds
- 10,000 molecules: <2 minutes

**Optimization steps:**

1. **Check Celery worker count:**
   ```yaml
   # Increase workers
   docker-compose up -d --scale celery-worker=4
   ```

2. **Increase worker concurrency:**
   ```yaml
   celery-worker:
     command: celery -A app.celery_app worker --concurrency=8
   ```

3. **Check Redis performance:**
   ```bash
   docker exec chemstructval-redis redis-cli INFO stats
   # Look for: instantaneous_ops_per_sec
   ```

4. **Enable caching:**
   ```
   Verify REDIS_URL is set in environment
   Check cache hit rate in Grafana
   ```

5. **Profile batch processing:**
   ```bash
   docker-compose logs celery-worker | grep -i "duration"
   ```

---

## API Issues

### Rate limit exceeded (429)

**Symptom:** HTTP 429 "Too Many Requests" response.

**Rate limits:**
- **Anonymous**: 10 requests/minute
- **With API key**: 300 requests/minute

**Solutions:**

1. **Create API key for higher limits:**
   ```bash
   curl -X POST http://localhost:8000/api/v1/api-keys \
     -H "Content-Type: application/json" \
     -d '{"name": "My Application"}'
   ```

   Response:
   ```json
   {
     "key": "csv_abc123...",
     "name": "My Application"
   }
   ```

2. **Use API key in requests:**
   ```bash
   curl http://localhost:8000/api/v1/validate \
     -H "X-API-Key: csv_abc123..." \
     -H "Content-Type: application/json" \
     -d '{"smiles": "CCO"}'
   ```

3. **Implement retry logic:**
   ```python
   import time
   from requests.exceptions import HTTPError

   def validate_with_retry(smiles, max_retries=3):
       for attempt in range(max_retries):
           try:
               response = requests.post(url, json={"smiles": smiles})
               response.raise_for_status()
               return response.json()
           except HTTPError as e:
               if e.response.status_code == 429:
                   wait_time = 2 ** attempt  # Exponential backoff
                   time.sleep(wait_time)
               else:
                   raise
       raise Exception("Max retries exceeded")
   ```

4. **Use Python client (has built-in retry):**
   ```python
   from chemstructval import ChemStructValClient

   client = ChemStructValClient(
       base_url="http://localhost:8000",
       api_key="csv_abc123..."
   )

   # Automatically handles rate limiting with retries
   result = client.validate("CCO")
   ```

### CORS errors

**Symptom:** Browser shows "blocked by CORS policy" error.

**Causes:**
- Frontend and backend on different domains
- CORS_ORIGINS not configured correctly

**Solutions:**

1. **Check CORS_ORIGINS setting:**
   ```bash
   # In .env or docker-compose.yml
   CORS_ORIGINS=["http://localhost:3000", "https://yourdomain.com"]
   ```

2. **For development (local testing):**
   ```yaml
   # docker-compose.yml
   backend:
     environment:
       - CORS_ORIGINS=["http://localhost:3000", "http://127.0.0.1:3000"]
   ```

3. **For production:**
   ```yaml
   # docker-compose.prod.yml
   backend:
     environment:
       - CORS_ORIGINS=["https://yourdomain.com"]
   ```

4. **Verify CORS headers:**
   ```bash
   curl -I -X OPTIONS http://localhost:8000/api/v1/validate \
     -H "Origin: http://localhost:3000" \
     -H "Access-Control-Request-Method: POST"

   # Should include:
   # Access-Control-Allow-Origin: http://localhost:3000
   ```

5. **Restart backend after changes:**
   ```bash
   docker-compose restart backend
   ```

### Authentication fails

**Symptom:** API returns 401 Unauthorized even with valid API key.

**Solutions:**

1. **Verify API key format:**
   ```
   Key should start with: csv_
   Example: csv_1234567890abcdef...
   ```

2. **Check header name:**
   ```bash
   # Correct
   curl -H "X-API-Key: csv_abc123..."

   # Wrong
   curl -H "Authorization: Bearer csv_abc123..."
   ```

3. **Verify key exists in database:**
   ```bash
   docker exec -it chemstructval-postgres psql -U chemstructval -d chemstructval

   SELECT id, name, created_at FROM api_keys;
   ```

4. **Check for whitespace:**
   ```python
   api_key = api_key.strip()  # Remove leading/trailing spaces
   ```

---

## Monitoring Issues

### Grafana shows "No Data"

**Symptom:** Dashboard panels are empty.

**Solutions:**

1. **Check Prometheus is scraping:**
   ```
   Visit: http://localhost:9090/targets
   Status should be: UP for chemstructval-backend
   ```

2. **Verify backend exposes /metrics:**
   ```bash
   curl http://localhost:8000/metrics
   # Should return Prometheus metrics in text format
   ```

3. **Check ENABLE_METRICS environment variable:**
   ```yaml
   # docker-compose.yml
   backend:
     environment:
       - ENABLE_METRICS=true
   ```

4. **Generate some traffic:**
   ```bash
   # Make requests to generate metrics
   for i in {1..10}; do
     curl http://localhost:8000/api/v1/validate \
       -H "Content-Type: application/json" \
       -d '{"smiles": "CCO"}'
   done
   ```

5. **Check Grafana data source:**
   ```
   Grafana → Settings → Data Sources → Prometheus
   Click "Test" → Should show "Data source is working"
   ```

6. **Check time range:**
   ```
   Top right in Grafana: Change to "Last 1 hour"
   If no recent data, try "Last 24 hours"
   ```

### Prometheus target down

**Symptom:** Prometheus shows backend target as DOWN.

**Solutions:**

1. **Check backend health:**
   ```bash
   curl http://localhost:8000/api/v1/health
   # Should return: {"status": "healthy"}
   ```

2. **Verify network connectivity:**
   ```bash
   # From Prometheus container
   docker exec chemstructval-prometheus wget -O- http://backend:8000/metrics
   ```

3. **Check backend logs:**
   ```bash
   docker-compose logs backend | grep -i error
   ```

4. **Verify prometheus.yml configuration:**
   ```bash
   docker exec chemstructval-prometheus cat /etc/prometheus/prometheus.yml
   # Should include:
   # - targets: ['backend:8000']
   ```

5. **Restart Prometheus:**
   ```bash
   docker-compose restart prometheus
   ```

### High Prometheus memory usage

**Symptom:** Prometheus container uses excessive memory (>4GB).

**Solutions:**

1. **Check TSDB size:**
   ```bash
   docker exec chemstructval-prometheus du -sh /prometheus
   ```

2. **Reduce retention time:**
   ```yaml
   # docker-compose.prod.yml
   prometheus:
     command:
       - '--storage.tsdb.retention.time=15d'  # Reduce from 30d
   ```

3. **Check metric cardinality:**
   ```bash
   curl http://localhost:9090/api/v1/status/tsdb
   # Look for high cardinality metrics
   ```

4. **Clear old data:**
   ```bash
   docker-compose down
   docker volume rm chemstructval_prometheus_data
   docker-compose up -d
   ```

---

## Frontend Issues

### Molecule visualization not showing

**Symptom:** Blank space where molecule should render.

**Solutions:**

1. **Check RDKit.js loaded:**
   ```javascript
   // Browser console
   window.RDKit
   // Should be defined, not undefined
   ```

2. **Check for JavaScript errors:**
   ```
   F12 → Console
   Look for errors mentioning RDKit or molecule
   ```

3. **Verify SMILES is valid:**
   ```javascript
   // Browser console
   const mol = window.RDKit.get_mol("CCO");
   console.log(mol);  // Should not be null
   mol.delete();  // Clean up
   ```

4. **Check SVG generation:**
   ```javascript
   // Browser console
   const mol = window.RDKit.get_mol("CCO");
   const svg = mol.get_svg();
   console.log(svg);  // Should contain <svg>...</svg>
   mol.delete();
   ```

### Frontend build fails

**Symptom:** `npm run build` fails with errors.

**Solutions:**

1. **Clear build cache:**
   ```bash
   rm -rf frontend/dist frontend/.vite
   npm run build
   ```

2. **Check for TypeScript errors:**
   ```bash
   npm run type-check
   ```

3. **Verify dependencies:**
   ```bash
   npm ci  # Clean install from lock file
   npm run build
   ```

4. **Check Node version:**
   ```bash
   node --version
   # Should be 18.x or 20.x
   ```

### Vite dev server hot reload not working

**Symptom:** Changes to code don't reflect in browser.

**Solutions:**

1. **Check Vite is watching files:**
   ```bash
   # Look for "page reload" in terminal
   npm run dev
   ```

2. **Manually reload browser:**
   ```
   Ctrl+R or F5
   ```

3. **Check file permissions:**
   ```bash
   ls -la frontend/src/
   # Files should be readable
   ```

4. **Restart Vite:**
   ```bash
   # Kill process (Ctrl+C)
   npm run dev
   ```

---

## Database Issues

### Database connection failed

**Symptom:** Backend logs show "could not connect to PostgreSQL".

**Solutions:**

1. **Check PostgreSQL is running:**
   ```bash
   docker-compose ps postgres
   # State should be: Up
   ```

2. **Verify connection string:**
   ```bash
   # In docker-compose.yml or .env
   DATABASE_URL=postgresql+asyncpg://chemstructval:chemstructval@postgres:5432/chemstructval
   ```

3. **Test connection from backend container:**
   ```bash
   docker exec -it chemstructval-backend python -c "
   from sqlalchemy import create_engine
   engine = create_engine('postgresql+psycopg2://chemstructval:chemstructval@postgres:5432/chemstructval')
   with engine.connect() as conn:
       print('Connected successfully')
   "
   ```

4. **Check PostgreSQL logs:**
   ```bash
   docker-compose logs postgres
   ```

5. **Reset database:**
   ```bash
   docker-compose down
   docker volume rm chemstructval_postgres_data
   docker-compose up -d
   ```

### Migrations fail

**Symptom:** Alembic migration errors.

**Solutions:**

1. **Check current revision:**
   ```bash
   docker-compose exec backend alembic current
   ```

2. **Run migrations:**
   ```bash
   docker-compose exec backend alembic upgrade head
   ```

3. **Reset migrations (development only):**
   ```bash
   # WARNING: This deletes all data
   docker-compose down -v
   docker-compose up -d postgres
   sleep 5
   docker-compose exec backend alembic upgrade head
   ```

### Database disk space full

**Symptom:** "No space left on device" errors.

**Solutions:**

1. **Check disk usage:**
   ```bash
   docker system df
   df -h
   ```

2. **Clean up Docker:**
   ```bash
   docker system prune -a
   docker volume prune
   ```

3. **Vacuum database:**
   ```bash
   docker exec chemstructval-postgres vacuumdb -U chemstructval -d chemstructval --full
   ```

---

## Getting Help

If your issue isn't listed here:

1. **Check logs:**
   ```bash
   docker-compose logs backend
   docker-compose logs celery-worker
   docker-compose logs frontend
   ```

2. **Search existing issues:**
   - GitHub Issues: https://github.com/yourusername/chemstructval/issues

3. **Open new issue with:**
   - Steps to reproduce
   - Full error messages
   - Environment details:
     ```bash
     docker --version
     docker-compose --version
     uname -a  # OS info
     ```
   - Relevant logs (attach as file)

4. **Community support:**
   - Discord: (if available)
   - Stack Overflow: Tag with `chemstructval`

5. **Commercial support:**
   - Email: support@chemstructval.com (if applicable)
