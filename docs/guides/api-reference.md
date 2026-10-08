# API Reference

The full interactive OpenAPI documentation is always available at
`/docs` (Swagger UI) or `/redoc` on a running instance. This page is a
static companion covering the same endpoints.

Base path: `/api/v1`

!!! warning "No authentication in v0.2.0"
    `alleleiq.api.dependencies.require_authenticated_caller` is a
    documented stub — every caller is currently treated as
    authenticated. Do not expose this release directly to an untrusted
    network without an authenticating layer in front of it (e.g. a
    reverse proxy).

## `POST /api/v1/variants/analyse`

Analyse a single-variant VCF upload.

**Request** — `multipart/form-data`:

| Field | Type | Default | Description |
|---|---|---|---|
| `vcf_file` | file | required | A VCF containing exactly one variant. |
| `genome_build` | string | `GRCh38` | Genome build; GRCh38 only. |
| `include_clinvar` | bool | `true` | Whether to query ClinVar. PS1, PM5, and PM1's variant-density signal depend on it; disabling it leaves those criteria unassessed and can change the preliminary classification. |

**Response** — `200 OK`: an `EvidenceReport` (see
[Data Models](../architecture/data-models.md)).

**Errors**:

- `422` — malformed VCF, unsupported genome build, or a VCF containing
  anything other than exactly one variant (use the batch endpoint
  instead).
- `502` — an upstream evidence source (gnomAD, VEP, or ClinVar) could
  not be reached.
- `413` — upload exceeds `ALLELEIQ_API_MAX_UPLOAD_BYTES` (default 50MB).
- `429` — rate limit exceeded (see below).

```bash
curl -X POST http://localhost:8000/api/v1/variants/analyse \
  -F "vcf_file=@variant.vcf" \
  -F "genome_build=GRCh38"
```

## `POST /api/v1/variants/analyse/batch`

Analyse every variant in a multi-variant VCF upload. Same request shape
as above.

**Response** — `200 OK`:

```json
{
  "total_variants": 5,
  "successful": 4,
  "reports": ["<EvidenceReport>", "..."],
  "errors": ["17-43094692-G-A: gnomAD request failed: ..."]
}
```

A per-variant evidence-source failure is recorded in `errors` rather
than aborting the whole batch — one unreachable upstream query cannot
silently drop the entire batch's results.

## `GET /api/v1/health`

Liveness check. Always `200` if the process can respond at all.

```json
{"status": "healthy", "version": "0.2.0"}
```

## `GET /api/v1/health/ready`

Readiness check: live connectivity to gnomAD, VEP, ClinVar, and UniProt.

**Response** — `200` if all reachable, **`503`** (same body) if any
upstream service is unreachable:

```json
{"status": "ready", "gnomad": true, "vep": true, "clinvar": true, "uniprot": true}
```

Suitable for direct use as a load-balancer or orchestrator readiness
probe (this is what the Docker image's `HEALTHCHECK` calls internally).

## Rate limiting

100 requests/hour per client IP by default
(`ALLELEIQ_API_RATE_LIMIT_PER_HOUR`). Exceeding it returns `429` with
`{"detail": "rate limit exceeded; try again later"}`. The limiter is
in-memory and per-process — a multi-worker deployment needs a shared
store (documented limitation in `alleleiq.api.middleware`).

## Request headers

Every response carries an `X-Request-ID` header (a fresh UUID per
request), logged alongside method, path, status code, and duration —
useful for correlating a client-side error with server-side logs.

## CORS

Configured via `ALLELEIQ_API_CORS_ALLOW_ORIGINS` (comma-separated).
Setting it to `*` for permissive local development automatically
disables credentialed CORS, since browsers reject the combination of a
wildcard origin and credentials.
