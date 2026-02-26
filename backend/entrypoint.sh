#!/bin/bash
set -e

echo "Running database migrations..."
python /app/run_migrations.py
echo "Starting application..."
exec "$@"
