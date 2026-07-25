#!/usr/bin/env bash
# Bootstrap a local AlleleIQ development environment.
set -euo pipefail

if ! command -v uv >/dev/null 2>&1; then
    echo "uv is required: https://docs.astral.sh/uv/getting-started/installation/" >&2
    exit 1
fi

uv python pin 3.12
uv sync --all-extras
uv run pre-commit install

if [ ! -f .env ]; then
    cp .env.example .env
    echo "Created .env from .env.example — review and adjust before running the API or CLI."
fi

echo "AlleleIQ development environment ready. Run 'make test' to verify."
