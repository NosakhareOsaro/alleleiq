.PHONY: install test lint type-check format docs docker-build clean run-api run-cli

UV := uv

install:
	$(UV) sync --all-extras
	$(UV) run pre-commit install

test:
	$(UV) run pytest

lint:
	$(UV) run ruff check src tests
	$(UV) run black --check src tests

format:
	$(UV) run ruff check --fix src tests
	$(UV) run black src tests

type-check:
	$(UV) run mypy --config-file=mypy.ini src

docs:
	$(UV) run mkdocs build --strict

docs-serve:
	$(UV) run mkdocs serve

docker-build:
	docker build -t alleleiq:latest .

run-api:
	$(UV) run uvicorn alleleiq.api.main:app --reload --port 8000

run-cli:
	$(UV) run alleleiq --help

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} +
	rm -rf .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage coverage.xml dist build site
