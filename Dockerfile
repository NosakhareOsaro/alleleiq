# syntax=docker/dockerfile:1

# ─────────────────────────────────────────────────────────────
# Stage 1: builder — install dependencies with uv into a venv
# ─────────────────────────────────────────────────────────────
FROM python:3.12-slim AS builder

# Build tools + htslib dependencies, in case cyvcf2/pysam need to compile
# from source on this platform (no manylinux wheel available). Not carried
# into the runtime image.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

WORKDIR /app
COPY pyproject.toml uv.lock README.md ./
COPY src ./src

RUN uv sync --frozen --no-dev --no-editable

# ─────────────────────────────────────────────────────────────
# Stage 2: runtime — minimal image, non-root user
# ─────────────────────────────────────────────────────────────
FROM python:3.12-slim AS runtime

# Shared-library runtime equivalents of the builder's -dev packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    libcurl4 \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid 1000 alleleiq \
    && useradd --uid 1000 --gid alleleiq --shell /usr/sbin/nologin --create-home alleleiq

WORKDIR /app
COPY --from=builder /app/.venv /app/.venv
COPY src ./src
COPY pyproject.toml ./

ENV PATH="/app/.venv/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN mkdir -p /app/reports && chown -R alleleiq:alleleiq /app

USER alleleiq

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import httpx; httpx.get('http://localhost:8000/api/v1/health', timeout=3).raise_for_status()" || exit 1

EXPOSE 8000

# Exec form: uvicorn runs as PID 1 and receives SIGTERM directly for a
# clean shutdown, rather than being wrapped by a shell that swallows it.
ENTRYPOINT ["uvicorn", "alleleiq.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
