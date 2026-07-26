"""FastAPI application entry point: app factory, lifespan, and router wiring."""

from __future__ import annotations

from collections.abc import AsyncIterator
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from alleleiq.__version__ import VERSION
from alleleiq.api.middleware import (
    RateLimitMiddleware,
    RequestContextMiddleware,
    RequestSizeLimitMiddleware,
)
from alleleiq.api.routers import health, variants
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.config import get_settings
from alleleiq.logging import configure_logging, get_logger

logger = get_logger(__name__)

_API_DESCRIPTION = """
AlleleIQ aggregates evidence from population frequency databases (gnomAD v4),
in silico prediction tools (REVEL, SpliceAI via Ensembl VEP), and curated
variant databases (ClinVar); applies ACMG/AMP criteria as codified by ACGS
2024; and produces structured, ancestry-stratified evidence reports.

**PRELIMINARY EVIDENCE SUMMARY — FOR CLINICAL SCIENTIST REVIEW ONLY.**
This API does not classify variants. Every report requires review by a
qualified Clinical Genetic Scientist following ACGS 2024 guidelines with
full consideration of clinical phenotype, family history, and all available
evidence.

Rate limiting: requests are limited per client IP (see
`ALLELEIQ_API_RATE_LIMIT_PER_HOUR`, default 100/hour). A rate-limited
request receives HTTP 429.
"""


@asynccontextmanager
async def lifespan(_app: FastAPI) -> AsyncIterator[None]:
    """Configure logging on startup; log gnomAD cache statistics on shutdown."""
    settings = get_settings()
    configure_logging(settings.log_level)
    logger.info("AlleleIQ API starting", extra={"version": VERSION})
    yield
    logger.info("gnomAD client cache statistics", extra=GnomADClient.cache_stats())
    logger.info("AlleleIQ API shutting down")


def create_app() -> FastAPI:
    """Build the FastAPI application: metadata, middleware, and routers."""
    settings = get_settings()
    app = FastAPI(
        title="AlleleIQ API",
        description=_API_DESCRIPTION,
        version=VERSION,
        lifespan=lifespan,
    )

    # Middleware order: CORS must be outermost so even error/rate-limited responses
    # carry CORS headers, hence it is added last (Starlette applies the
    # most-recently-added middleware first / outermost).
    app.add_middleware(RequestContextMiddleware)
    app.add_middleware(RateLimitMiddleware, limit_per_hour=settings.api_rate_limit_per_hour)
    app.add_middleware(RequestSizeLimitMiddleware, max_bytes=settings.api_max_upload_bytes)

    is_permissive_dev_cors = settings.api_cors_allow_origins == ["*"]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.api_cors_allow_origins,
        # Wildcard origins and credentialed CORS are mutually exclusive per the
        # Fetch spec (browsers reject the combination); ALLELEIQ_API_CORS_ALLOW_ORIGINS=*
        # is meant for permissive local development only, so credentials are disabled
        # in that case rather than silently producing a header browsers ignore.
        allow_credentials=not is_permissive_dev_cors,
        allow_methods=["GET", "POST"],
        allow_headers=["*"],
    )

    app.include_router(health.router, prefix="/api/v1")
    app.include_router(variants.router, prefix="/api/v1")

    return app


app = create_app()
