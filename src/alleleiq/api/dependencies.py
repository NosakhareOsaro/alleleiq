"""Shared FastAPI dependencies.

Includes an authentication stub. AlleleIQ v0.2.0 has no
authentication implemented — this dependency is a documented
placeholder for the access-control layer a real clinical deployment
would require (e.g. NHS Care Identity Service / OAuth2 integration),
made explicit here rather than left as a silent gap. It is wired into
the routers now so adding real authentication later is a one-function
change, not an API-wide refactor.
"""

from __future__ import annotations

from alleleiq.config import Settings, get_settings


def get_settings_dependency() -> Settings:
    """FastAPI dependency wrapping :func:`alleleiq.config.get_settings`.

    ``get_settings`` is itself process-wide cached (``lru_cache``), so
    this dependency needs no caching of its own — it exists purely so
    routers depend on an injectable function rather than importing the
    module-level getter directly, which keeps them testable via
    FastAPI's dependency-override mechanism.
    """
    return get_settings()


async def require_authenticated_caller() -> None:
    """Authentication stub.

    Raises nothing and returns nothing in v0.2.0: every caller is
    currently treated as authenticated. This is a placeholder, not a
    security control — do not deploy this release against real patient
    data behind a network boundary that does not enforce access control
    at another layer (e.g. an authenticating reverse proxy).
    """
    return None
