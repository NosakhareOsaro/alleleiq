"""Abstract base class for external evidence-source HTTP clients.

Every client in :mod:`alleleiq.clients` (gnomAD, ClinVar, VEP) shares the
same operational concerns — retrying transient network failures, being a
polite citizen of a shared public API, and logging what was requested
without logging genomic content — so those concerns live here once
rather than being reimplemented per client.
"""

from __future__ import annotations

import asyncio
import time
from abc import ABC, abstractmethod
from types import TracebackType
from typing import Any, Self

import httpx
from tenacity import (
    AsyncRetrying,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential_jitter,
)

from alleleiq.logging import get_logger

logger = get_logger(__name__)

_RETRYABLE_EXCEPTIONS = (httpx.TransportError, httpx.HTTPStatusError)


class BaseClient(ABC):
    """Shared HTTP plumbing for evidence-source API clients.

    Subclasses are responsible for building request payloads, parsing
    responses, and translating failures into their own
    :class:`~alleleiq.exceptions.DatabaseQueryError` subclass; this base
    class only handles connection lifecycle, retry, rate limiting, and
    request-level logging.
    """

    def __init__(
        self,
        base_url: str,
        timeout_seconds: float,
        max_retries: int,
        min_request_interval_seconds: float = 0.0,
    ) -> None:
        """Initialise shared client state.

        Args:
            base_url: The API's base URL.
            timeout_seconds: Per-request timeout.
            max_retries: Number of retries after an initial failed attempt
                (so ``max_retries=3`` means up to 4 total attempts).
            min_request_interval_seconds: Minimum time to wait between the
                start of consecutive requests, enforced regardless of
                concurrency, as a courtesy to shared public endpoints.
        """
        self._base_url = base_url
        self._timeout_seconds = timeout_seconds
        self._max_retries = max_retries
        self._min_request_interval_seconds = min_request_interval_seconds
        self._last_request_started_at: float | None = None
        self._throttle_lock = asyncio.Lock()
        self._client: httpx.AsyncClient | None = None

    async def __aenter__(self) -> Self:
        self._client = httpx.AsyncClient(base_url=self._base_url, timeout=self._timeout_seconds)
        return self

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        if self._client is not None:
            await self._client.aclose()
            self._client = None

    @property
    def client(self) -> httpx.AsyncClient:
        """The underlying httpx client. Only valid inside an ``async with`` block."""
        if self._client is None:
            raise RuntimeError(
                f"{type(self).__name__} must be used as an async context manager "
                "before making requests, e.g. `async with SomeClient(...) as c: ...`"
            )
        return self._client

    async def _throttle(self) -> None:
        """Enforce the minimum interval between consecutive request starts."""
        if self._min_request_interval_seconds <= 0:
            return
        async with self._throttle_lock:
            now = time.monotonic()
            if self._last_request_started_at is not None:
                elapsed = now - self._last_request_started_at
                remaining = self._min_request_interval_seconds - elapsed
                if remaining > 0:
                    await asyncio.sleep(remaining)
            self._last_request_started_at = time.monotonic()

    async def _request(
        self,
        method: str,
        url: str,
        *,
        raise_for_status: bool = True,
        **kwargs: Any,
    ) -> httpx.Response:
        """Issue an HTTP request with throttling, retry, and structured logging.

        Args:
            method: HTTP method, e.g. ``"GET"`` or ``"POST"``.
            url: Path or full URL, resolved against the client's base URL.
            raise_for_status: Whether to raise ``httpx.HTTPStatusError`` on
                4xx/5xx responses (and thus make them eligible for retry).
            **kwargs: Forwarded to :meth:`httpx.AsyncClient.request`.

        Returns:
            The successful HTTP response.

        Raises:
            httpx.TransportError: On a network failure that persisted
                through all retry attempts.
            httpx.HTTPStatusError: On a persistent 4xx/5xx response, if
                ``raise_for_status`` is True.
        """
        await self._throttle()
        attempt_number = 0

        async for attempt in AsyncRetrying(
            stop=stop_after_attempt(self._max_retries + 1),
            wait=wait_exponential_jitter(initial=0.5, max=10),
            retry=retry_if_exception_type(_RETRYABLE_EXCEPTIONS),
            reraise=True,
        ):
            with attempt:
                attempt_number += 1
                started_at = time.monotonic()
                logger.debug(
                    "issuing external API request",
                    extra={
                        "client": type(self).__name__,
                        "method": method,
                        "url": url,
                        "attempt": attempt_number,
                    },
                )
                response = await self.client.request(method, url, **kwargs)
                duration_ms = round((time.monotonic() - started_at) * 1000, 1)
                if raise_for_status:
                    response.raise_for_status()
                logger.debug(
                    "external API request completed",
                    extra={
                        "client": type(self).__name__,
                        "method": method,
                        "url": url,
                        "status_code": response.status_code,
                        "duration_ms": duration_ms,
                        "attempt": attempt_number,
                    },
                )
                return response

        # Unreachable: AsyncRetrying with reraise=True either returns from the
        # loop body above or raises. This satisfies mypy's control-flow check.
        raise AssertionError("unreachable")  # pragma: no cover

    @abstractmethod
    async def health_check(self) -> bool:
        """Check whether the upstream service is reachable and responding.

        Returns:
            True if the service responded successfully, False otherwise.
            Implementations should not raise for a failed health check.
        """
