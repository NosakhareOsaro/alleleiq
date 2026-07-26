"""A minimal process-wide TTL cache, used by GnomADClient to avoid redundant queries.

Deliberately process-wide (a module-level singleton, not per-instance
state): :mod:`alleleiq.pipeline` constructs a fresh ``GnomADClient`` per
variant analysed, so a per-instance cache would never help across the
multiple variants in a single multi-variant VCF — the case this cache
exists to speed up.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

_MISSING = object()


@dataclass
class TTLCache:
    """A simple time-based cache with hit/miss statistics."""

    _store: dict[str, tuple[float, Any]] = field(default_factory=dict)
    hits: int = 0
    misses: int = 0

    def get(self, key: str) -> Any:
        """Return the cached value for ``key``, or the sentinel ``_MISSING``."""
        entry = self._store.get(key)
        if entry is None:
            self.misses += 1
            return _MISSING

        expires_at, value = entry
        if time.monotonic() > expires_at:
            del self._store[key]
            self.misses += 1
            return _MISSING

        self.hits += 1
        return value

    def set(self, key: str, value: Any, ttl_seconds: float) -> None:
        """Cache ``value`` under ``key`` for ``ttl_seconds``."""
        self._store[key] = (time.monotonic() + ttl_seconds, value)

    def clear(self) -> None:
        """Remove all entries and reset statistics. Primarily for test isolation."""
        self._store.clear()
        self.hits = 0
        self.misses = 0

    @property
    def size(self) -> int:
        """Number of entries currently cached (including possibly-expired ones)."""
        return len(self._store)


def is_cache_miss(value: Any) -> bool:
    """True if ``value`` is the sentinel returned by :meth:`TTLCache.get` on a miss."""
    return value is _MISSING
