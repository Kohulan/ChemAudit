FROM python:3.11-slim

# Install system dependencies required for RDKit and WeasyPrint
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libsm6 \
    libgomp1 \
    libexpat1 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf-2.0-0 \
    libffi-dev \
    shared-mime-info \
    curl \
    g++ \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Download OPSIN JAR for IUPAC name conversion
RUN mkdir -p /app/data && curl -L -o /app/data/opsin.jar \
    https://github.com/dan2097/opsin/releases/download/2.8.0/opsin-cli-2.8.0-jar-with-dependencies.jar

# Copy pyproject.toml and install dependencies
COPY pyproject.toml ./
RUN pip install --no-cache-dir -e .

# Copy application code
COPY . .

# Run as non-root user for security (matches Dockerfile.prod pattern)
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# Expose port
EXPOSE 8000

# Entrypoint runs Alembic migrations before starting the app
ENTRYPOINT ["bash", "entrypoint.sh"]

# Run uvicorn
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
