"""
WebSocket Connection Manager

Manages WebSocket connections with Redis pub/sub for horizontal scaling.
Forwards batch progress updates to connected clients.
"""
import asyncio
import json
from typing import Dict, List, Optional

import redis.asyncio as redis
from fastapi import WebSocket, WebSocketDisconnect

from app.core.config import settings


class ConnectionManager:
    """
    Manages WebSocket connections with Redis pub/sub integration.

    Features:
    - Multiple clients can connect to same job_id
    - Redis pub/sub enables horizontal scaling across multiple backend instances
    - Automatic cleanup on disconnect
    """

    def __init__(self):
        self.active_connections: Dict[str, List[WebSocket]] = {}
        self._redis: Optional[redis.Redis] = None
        self._redis_url = settings.REDIS_URL
        self._subscriber_tasks: Dict[str, asyncio.Task] = {}

    async def init_redis(self, redis_url: Optional[str] = None) -> None:
        """Initialize Redis connection for pub/sub."""
        url = redis_url or self._redis_url
        self._redis = redis.from_url(url)

    async def close_redis(self) -> None:
        """Close Redis connection."""
        if self._redis:
            await self._redis.close()
            self._redis = None

    async def connect(self, job_id: str, websocket: WebSocket) -> None:
        """
        Accept WebSocket connection and register for job updates.

        Args:
            job_id: Job identifier to subscribe to
            websocket: WebSocket connection to register
        """
        await websocket.accept()

        if job_id not in self.active_connections:
            self.active_connections[job_id] = []
        self.active_connections[job_id].append(websocket)

        # Start subscriber task for this job if not already running
        if job_id not in self._subscriber_tasks:
            task = asyncio.create_task(self._subscribe_to_job(job_id))
            self._subscriber_tasks[job_id] = task

    def disconnect(self, job_id: str, websocket: WebSocket) -> None:
        """
        Remove WebSocket connection from job subscription.

        Args:
            job_id: Job identifier
            websocket: WebSocket connection to remove
        """
        if job_id in self.active_connections:
            try:
                self.active_connections[job_id].remove(websocket)
            except ValueError:
                pass

            # Clean up if no more connections for this job
            if not self.active_connections[job_id]:
                del self.active_connections[job_id]

                # Cancel subscriber task
                if job_id in self._subscriber_tasks:
                    self._subscriber_tasks[job_id].cancel()
                    del self._subscriber_tasks[job_id]

    async def _subscribe_to_job(self, job_id: str) -> None:
        """
        Subscribe to Redis channel for job updates and forward to WebSockets.

        Args:
            job_id: Job identifier to subscribe to
        """
        if not self._redis:
            return

        pubsub = self._redis.pubsub()
        channel = f"batch:progress:{job_id}"

        try:
            await pubsub.subscribe(channel)

            while job_id in self.active_connections:
                try:
                    message = await asyncio.wait_for(
                        pubsub.get_message(ignore_subscribe_messages=True),
                        timeout=1.0,
                    )

                    if message and message["type"] == "message":
                        data = json.loads(message["data"])
                        await self._broadcast_to_job(job_id, data)

                        # If job is complete/failed/cancelled, stop subscribing
                        status = data.get("status", "")
                        if status in ("complete", "failed", "cancelled"):
                            break

                except asyncio.TimeoutError:
                    # No message, check if still have connections
                    if job_id not in self.active_connections:
                        break
                except asyncio.CancelledError:
                    break

        except Exception:
            pass
        finally:
            try:
                await pubsub.unsubscribe(channel)
                await pubsub.close()
            except Exception:
                pass

    async def _broadcast_to_job(self, job_id: str, data: dict) -> None:
        """
        Send message to all WebSocket connections for a job.

        Args:
            job_id: Job identifier
            data: Data to broadcast
        """
        if job_id not in self.active_connections:
            return

        dead_connections = []

        for websocket in self.active_connections[job_id]:
            try:
                await websocket.send_json(data)
            except Exception:
                dead_connections.append(websocket)

        # Clean up dead connections
        for ws in dead_connections:
            self.disconnect(job_id, ws)

    async def send_initial_status(
        self, job_id: str, websocket: WebSocket
    ) -> None:
        """
        Send current job status when client first connects.

        Args:
            job_id: Job identifier
            websocket: WebSocket to send status to
        """
        if not self._redis:
            return

        try:
            # Get current progress from Redis
            data = await self._redis.get(f"batch:job:{job_id}")
            if data:
                progress = json.loads(data)
                await websocket.send_json(progress)
        except Exception:
            pass


# Singleton instance
manager = ConnectionManager()
