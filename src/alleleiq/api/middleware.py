"""ASGI middleware: request context/logging, rate limiting, and request size limiting.

Known limitation: :class:`RateLimitMiddleware` is an in-memory,
per-process fixed-window-ish limiter (a sliding window over a deque per
client IP). It is adequate for a single-process deployment but does
**not** share state across multiple workers or replicas — a real
multi-worker production deployment needs a shared store (e.g. Redis).
This is a known, documented scope limitation for v0.1.0-alpha, not an
oversight.
"""

from __future__ import annotations

import time
from collections import defaultdict, deque
from uuid import uuid4

from starlette.middleware.base import BaseHTTPMiddleware, RequestResponseEndpoint
from starlette.requests import Request
from starlette.responses import JSONResponse, Response
from starlette.types import ASGIApp

from alleleiq.logging import get_logger

logger = get_logger(__name__)

_SECONDS_PER_HOUR = 3600


class RequestContextMiddleware(BaseHTTPMiddleware):
    """Attach a request ID to every response and log request/response metadata.

    Logged fields are limited to request metadata (method, path, status
    code, duration, request ID) — never request or response bodies,
    which may contain genomic data.
    """

    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        request_id = str(uuid4())
        started_at = time.monotonic()

        response = await call_next(request)

        duration_ms = round((time.monotonic() - started_at) * 1000, 1)
        response.headers["X-Request-ID"] = request_id
        logger.info(
            "http request",
            extra={
                "request_id": request_id,
                "method": request.method,
                "path": request.url.path,
                "status_code": response.status_code,
                "duration_ms": duration_ms,
            },
        )
        return response


class _SlidingWindowRateLimiter:
    """Per-client sliding-window request counter."""

    def __init__(self, limit_per_hour: int) -> None:
        self._limit = limit_per_hour
        self._requests: dict[str, deque[float]] = defaultdict(deque)

    def is_allowed(self, client_id: str) -> bool:
        now = time.monotonic()
        window_start = now - _SECONDS_PER_HOUR
        timestamps = self._requests[client_id]
        while timestamps and timestamps[0] < window_start:
            timestamps.popleft()
        if len(timestamps) >= self._limit:
            return False
        timestamps.append(now)
        return True


class RateLimitMiddleware(BaseHTTPMiddleware):
    """Reject requests beyond a configured per-client-IP hourly limit with HTTP 429."""

    def __init__(self, app: ASGIApp, *, limit_per_hour: int) -> None:
        super().__init__(app)
        self._limiter = _SlidingWindowRateLimiter(limit_per_hour)

    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        client_id = request.client.host if request.client else "unknown"
        if not self._limiter.is_allowed(client_id):
            return JSONResponse(
                status_code=429,
                content={"detail": "rate limit exceeded; try again later"},
            )
        return await call_next(request)


class RequestSizeLimitMiddleware(BaseHTTPMiddleware):
    """Reject requests whose declared Content-Length exceeds a configured maximum.

    Checked against the ``Content-Length`` header before the body is
    read, so an oversized upload is rejected without buffering it into
    memory first. This is a complementary, earlier check to the
    read-time size check in the variants router (which remains
    necessary for chunked transfer-encoded requests with no
    ``Content-Length`` header).
    """

    def __init__(self, app: ASGIApp, *, max_bytes: int) -> None:
        super().__init__(app)
        self._max_bytes = max_bytes

    async def dispatch(self, request: Request, call_next: RequestResponseEndpoint) -> Response:
        content_length = request.headers.get("content-length")
        if content_length is not None:
            try:
                declared_size = int(content_length)
            except ValueError:
                declared_size = None
            if declared_size is not None and declared_size > self._max_bytes:
                return JSONResponse(
                    status_code=413,
                    content={"detail": f"request body exceeds the {self._max_bytes}-byte limit"},
                )
        return await call_next(request)
