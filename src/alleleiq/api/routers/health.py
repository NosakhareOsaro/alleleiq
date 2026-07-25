"""Health and readiness endpoints."""

from __future__ import annotations

from fastapi import APIRouter, Depends, Response, status

from alleleiq.__version__ import VERSION
from alleleiq.api.dependencies import get_settings_dependency
from alleleiq.api.schemas import HealthResponse, ReadinessResponse
from alleleiq.clients.clinvar import ClinVarClient
from alleleiq.clients.gnomad import GnomADClient
from alleleiq.clients.vep import VEPClient
from alleleiq.config import Settings

router = APIRouter(prefix="/health", tags=["health"])


@router.get("", response_model=HealthResponse, summary="Liveness check")
async def health() -> HealthResponse:
    """Return process liveness. Always 200 if the server can respond at all."""
    return HealthResponse(status="healthy", version=VERSION)


@router.get(
    "/ready",
    response_model=ReadinessResponse,
    summary="Readiness check: upstream evidence-source connectivity",
)
async def ready(
    response: Response, settings: Settings = Depends(get_settings_dependency)
) -> ReadinessResponse:
    """Check connectivity to gnomAD, VEP, and ClinVar.

    Returns HTTP 503 (with the same body) if any upstream service is
    unreachable, so this can be used directly as a load-balancer or
    orchestrator readiness probe.
    """
    async with (
        GnomADClient(settings) as gnomad_client,
        VEPClient(settings) as vep_client,
        ClinVarClient(settings) as clinvar_client,
    ):
        gnomad_ok = await gnomad_client.health_check()
        vep_ok = await vep_client.health_check()
        clinvar_ok = await clinvar_client.health_check()

    all_ready = gnomad_ok and vep_ok and clinvar_ok
    if not all_ready:
        response.status_code = status.HTTP_503_SERVICE_UNAVAILABLE

    return ReadinessResponse(
        status="ready" if all_ready else "degraded",
        gnomad=gnomad_ok,
        vep=vep_ok,
        clinvar=clinvar_ok,
    )
