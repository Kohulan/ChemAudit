"""
ChemStructVal API - Chemical Structure Validation Suite
"""
from contextlib import asynccontextmanager
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware

from app.core.config import settings
from app.core.exceptions import (
    ChemStructValException,
    chemstructval_exception_handler,
    generic_exception_handler,
)
from app.api.routes import alerts, batch, health, scoring, standardization, validation
from app.websockets import manager


@asynccontextmanager
async def lifespan(app: FastAPI):
    """
    Application lifespan manager.

    Handles startup and shutdown events.
    """
    # Startup
    print(f"Starting {settings.APP_NAME} v{settings.APP_VERSION}")

    # Initialize WebSocket manager Redis connection
    await manager.init_redis()
    print("WebSocket manager initialized with Redis")

    yield

    # Shutdown
    await manager.close_redis()
    print("Shutting down...")


app = FastAPI(
    title="ChemStructVal API",
    description="Chemical Structure Validation and Standardization API",
    version=settings.APP_VERSION,
    lifespan=lifespan,
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Register exception handlers
app.add_exception_handler(ChemStructValException, chemstructval_exception_handler)
app.add_exception_handler(Exception, generic_exception_handler)

# Include routers
app.include_router(health.router, prefix="/api/v1", tags=["health"])
app.include_router(validation.router, prefix="/api/v1", tags=["validation"])
app.include_router(alerts.router, prefix="/api/v1", tags=["alerts"])
app.include_router(standardization.router, prefix="/api/v1", tags=["standardization"])
app.include_router(scoring.router, prefix="/api/v1", tags=["scoring"])
app.include_router(batch.router, prefix="/api/v1", tags=["batch"])


@app.websocket("/ws/batch/{job_id}")
async def batch_progress_websocket(websocket: WebSocket, job_id: str):
    """
    WebSocket endpoint for real-time batch progress updates.

    Connect to /ws/batch/{job_id} after uploading a file to receive
    progress updates in real-time.

    Message format:
    {
        "job_id": "...",
        "status": "processing|complete|failed|cancelled",
        "progress": 0-100,
        "processed": int,
        "total": int,
        "eta_seconds": int|null
    }
    """
    await manager.connect(job_id, websocket)

    # Send initial status
    await manager.send_initial_status(job_id, websocket)

    try:
        # Keep connection alive, waiting for close
        while True:
            try:
                # Wait for client messages (pings, close)
                data = await websocket.receive_text()
                # Client can send "ping" to keep alive
                if data == "ping":
                    await websocket.send_text("pong")
            except WebSocketDisconnect:
                break
    finally:
        manager.disconnect(job_id, websocket)


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "status": "running",
    }
