---
sidebar_position: 6
title: Rate Limits
description: API rate limits, headers, and best practices
---

# Rate Limits

ChemAudit implements rate limiting to ensure fair usage and prevent abuse. Rate limits are set per endpoint.

## Rate Limit Tiers

### Per-Endpoint Limits

Requests without an API key are counted per client IP address:

| Endpoint | Limit |
|----------|-------|
| `POST /validate` | 10 req/min |
| `POST /validate/async` | 30 req/min |
| `GET /checks` | 10 req/min |
| `POST /alerts` | 10 req/min |
| `POST /alerts/quick-check` | 10 req/min |
| `GET /alerts/catalogs` | 10 req/min |
| `POST /score` | 10 req/min |
| `POST /score/compare` | 10 req/min |
| `POST /standardize` | 10 req/min |
| `POST /batch/upload` | 3 req/min |
| `GET /batch/{job_id}` | 60 req/min |
| `GET /batch/{job_id}/status` | 10 req/min |
| `GET /batch/{job_id}/stats` | 10 req/min |
| `DELETE /batch/{job_id}` | 10 req/min |
| `POST /batch/detect-columns` | 10 req/min |
| `POST /qsar-ready/single` | 30 req/min |
| `POST /qsar-ready/batch/upload` | 3 req/min |
| `GET /qsar-ready/batch/*/status` | 60 req/min |
| `GET /qsar-ready/batch/*/results` | 30 req/min |
| `GET /qsar-ready/batch/*/download/*` | 10 req/min |
| `POST /structure-filter/filter` | 20 req/min |
| `POST /structure-filter/score` | 30 req/min |
| `POST /structure-filter/reinvent-score` | 30 req/min |
| `POST /structure-filter/batch/upload` | 3 req/min |
| `GET /structure-filter/batch/*/status` | 60 req/min |
| `GET /structure-filter/batch/*/results` | 30 req/min |
| `GET /structure-filter/batch/*/download/*` | 10 req/min |
| `POST /dataset/upload` | 3 req/min |
| `GET /dataset/*/status` | 60 req/min |
| `GET /dataset/*/results` | 30 req/min |
| `POST /dataset/*/diff` | 10 req/min |
| `GET /dataset/*/download/*` | 10 req/min |
| `POST /integrations/*` | 30 req/min |
| `POST /integrations/compare` | 10 req/min |
| `POST /resolve` | 30 req/min |

### With an API Key

Requests with an `X-API-Key` header are counted in a separate bucket for that key instead of the client IP address. The per-endpoint limits above still apply; an API key does not raise them.

## Rate Limit Headers

Only `429` responses carry rate limit headers. Successful responses do not include them.

| Header | Description |
|--------|-------------|
| `Retry-After` | Seconds to wait before retrying |
| `X-RateLimit-Limit` | The limit that was exceeded, e.g. `10 per 1 minute` |

## When Rate Limited

If you exceed the rate limit, you'll receive a 429 response:

```json
{
  "error": "rate_limit_exceeded",
  "message": "Rate limit exceeded: 10 per 1 minute",
  "retry_after": 60,
  "violations": 1,
  "ban_threshold": 100
}
```

**Headers:**

```
Retry-After: 60
X-RateLimit-Limit: 10 per 1 minute
```

## IP Banning

Repeated rate limit violations (100+ within the tracking window) result in a temporary IP ban:

- **Ban duration**: 60 minutes (default)
- **Response**: 403 Forbidden
- **Message**: "Your IP has been temporarily banned due to repeated rate limit violations"

:::warning Avoid Bans
Implement proper rate limiting and error handling in your code to avoid IP bans.
:::

## Best Practices

### 1. Respect Retry-After

The server does not report remaining requests. When you get a 429, wait for the number of seconds in `Retry-After`:

```python
import requests
import time

response = requests.post(url, json=data)
if response.status_code == 429:
    wait_seconds = int(response.headers.get('Retry-After', 60))
    print(f"Rate limited. Waiting {wait_seconds} seconds...")
    time.sleep(wait_seconds)
```

### 2. Retry When Rate Limited

Retry after the delay the server asks for:

```python
import time

def api_call_with_backoff(url, data, max_retries=3):
    for attempt in range(max_retries):
        response = requests.post(url, json=data)

        if response.status_code == 200:
            return response.json()
        elif response.status_code == 429:
            # Rate limited
            wait_time = int(response.headers.get('Retry-After', 60))

            print(f"Rate limited. Waiting {wait_time} seconds...")
            time.sleep(wait_time)
        else:
            response.raise_for_status()

    raise Exception("Max retries exceeded")
```

### 3. Batch Requests

Instead of many single validations, use batch processing:

```python
# Bad: Many API calls
for smiles in smiles_list:
    response = requests.post("/validate", json={"molecule": smiles})
    # Likely to hit rate limit

# Good: Single batch upload
with open('molecules.csv', 'rb') as f:
    response = requests.post("/batch/upload", files={"file": f})
job_id = response.json()['job_id']
# Monitor via WebSocket
```

### 4. Use API Keys

An API key gives your requests their own rate limit bucket, separate from other clients on the same IP address. The limits are the same:

```python
import os

api_key = os.environ['CHEMAUDIT_API_KEY']

response = requests.post(
    url,
    headers={"X-API-Key": api_key},
    json=data
)
```

### 5. Distribute Requests

Spread requests evenly across time:

```python
import time

def rate_limited_validate(molecules, max_per_minute=10):
    delay = 60.0 / max_per_minute  # Seconds between requests

    results = []
    for molecule in molecules:
        response = requests.post("/validate", json={"molecule": molecule})
        results.append(response.json())

        time.sleep(delay)  # Ensure we don't exceed limit

    return results
```

### 6. Cache Results

Avoid redundant requests by caching:

```python
cache = {}

def validate_with_cache(smiles):
    if smiles in cache:
        return cache[smiles]

    response = requests.post("/validate", json={"molecule": smiles})
    result = response.json()

    cache[smiles] = result
    return result
```

:::info Server-Side Caching
ChemAudit caches validation results by InChIKey for 1 hour, so redundant requests are already fast.
:::

## Rate Limit Strategies

### Strategy 1: Token Bucket

Maintain a local token bucket to ensure you never exceed limits:

```python
import time
from collections import deque

class RateLimiter:
    def __init__(self, max_requests, time_window):
        self.max_requests = max_requests
        self.time_window = time_window
        self.requests = deque()

    def wait_if_needed(self):
        now = time.time()

        # Remove old requests outside time window
        while self.requests and self.requests[0] < now - self.time_window:
            self.requests.popleft()

        # If at limit, wait
        if len(self.requests) >= self.max_requests:
            sleep_time = self.time_window - (now - self.requests[0])
            if sleep_time > 0:
                time.sleep(sleep_time)

        self.requests.append(time.time())

# Usage
limiter = RateLimiter(max_requests=10, time_window=60)

for molecule in molecules:
    limiter.wait_if_needed()
    response = requests.post("/validate", json={"molecule": molecule})
```

### Strategy 2: Adaptive Rate Limiting

Adjust request rate based on server responses:

```python
class AdaptiveRateLimiter:
    def __init__(self):
        self.delay = 0.1  # Start with 100ms delay

    def on_success(self, response):
        # Request went through, speed up slightly
        self.delay = max(0.05, self.delay * 0.9)

    def on_rate_limit(self, response):
        # Hit limit: wait as long as the server asks, then slow down
        time.sleep(int(response.headers.get('Retry-After', 60)))
        self.delay = min(30.0, self.delay * 2.0)

    def wait(self):
        time.sleep(self.delay)
```

## Testing Rate Limits

Test your rate limiting implementation:

```python
import time
import requests

# Test hitting the limit
for i in range(15):  # More than limit
    start = time.time()
    response = requests.post("/validate", json={"molecule": "CCO"})

    print(f"Request {i+1}: Status {response.status_code}, "
          f"Retry-After: {response.headers.get('Retry-After')}")

    if response.status_code == 429:
        print("Rate limited as expected")
        break

    time.sleep(0.1)  # Small delay between requests
```

## Next Steps

- **[Authentication](/docs/api/authentication)** - Get API keys
- **[Error Handling](/docs/api/error-handling)** - Handle 429 responses
- **[Batch Processing](/docs/user-guide/batch-processing)** - Process many molecules efficiently
