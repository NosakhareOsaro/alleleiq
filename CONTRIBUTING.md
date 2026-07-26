# Contributing to AlleleIQ

Thank you for considering a contribution. AlleleIQ is clinical decision
support software — contributions that touch evidence-gathering or
scoring logic are held to a higher bar than typical open-source
tooling, because a subtle bug there could affect a patient's
interpretation. This document explains what that means in practice.

## Development environment setup

```bash
git clone https://github.com/NosakhareOsaro/alleleiq.git
cd alleleiq
./scripts/setup_dev.sh
```

This pins Python 3.12 (needed for reliable `cyvcf2`/`pysam` wheel
availability), installs all dependencies via `uv sync --all-extras`,
installs pre-commit hooks, and copies `.env.example` to `.env`.

Verify your setup:

```bash
make test        # pytest
make lint         # ruff + black --check
make type-check   # mypy --strict
```

All three must pass before opening a PR — the CI pipeline
(`.github/workflows/ci.yml`) runs the same commands and will reject a
PR that doesn't.

## Commit message format

This project uses [Conventional Commits](https://www.conventionalcommits.org/),
enforced by `python-semantic-release` for automated versioning:

```
<type>: <description>

[optional body]
```

Types: `feat`, `fix`, `docs`, `test`, `refactor`, `chore`, `perf`, `style`, `ci`.

- `feat` and `fix` drive automatic version bumps (minor and patch
  respectively) — see `[tool.semantic_release]` in `pyproject.toml`.
- Write the description in the imperative mood ("add", not "added" or
  "adds").
- Reference the clinical or architectural *why* in the body when it
  isn't obvious from the diff — see recent commit history for examples
  of the level of detail expected for scoring-logic changes
  specifically.

## Pull request process

1. Open an issue first for anything beyond a small fix — especially
   changes to ACMG/AMP criteria thresholds, the Bayesian point weights,
   or coverage assessment logic — so the clinical rationale can be
   discussed before code is written.
2. One logical change per PR. A new criterion evaluator, a bug fix, and
   a docs update are three PRs, not one.
3. PRs must pass CI (lint, type-check, test with coverage gate,
   security scan, docs build) before review.
4. At least one maintainer review is required before merge. For changes
   to `alleleiq.scoring.*`, review must include verification against
   `scripts/validate_known_variants.py` (see below).

## Testing requirements

- **Every new or modified ACMG/AMP criterion evaluator requires unit
  tests covering both its "met" and "not met" paths, plus every
  `requires_clinical_input` branch.** See
  `tests/unit/test_pvs1_criterion.py` for the level of decision-tree
  coverage expected — it exercises every branch of PVS1's logic, not
  just the happy path.
- New external API integrations (a new client, or a new endpoint on an
  existing one) require unit tests with mocked responses
  (`httpx.MockTransport`) covering success, not-found, malformed
  response, and network-failure cases — see any file in
  `tests/unit/test_*_client.py` for the pattern.
- Run `uv run python scripts/validate_known_variants.py` after any
  change to `alleleiq.scoring.*` or `alleleiq.clients.*`. A change that
  introduces a **dangerous** outcome (an opposite-direction
  classification versus a ClinVar expert-panel call) must not be
  merged without an explicit, reviewed explanation of why the change is
  correct despite that result.
- Coverage must not regress below the CI gate (`--cov-fail-under=85`).

## Documentation requirements

- Every public class, method, and function needs a Google-style
  docstring (see `.claude`/existing modules for the house style: lead
  with what it does, explain non-obvious *why* in the body, document
  `Args`/`Returns`/`Raises`).
- A new or materially changed criterion, data source, or scoring rule
  needs a corresponding update to the relevant page under `docs/`
  (`docs/architecture/acmg-framework.md` for criteria changes,
  `docs/architecture/data-models.md` for model changes, etc.).
- Run `make docs` (`mkdocs build --strict`) before submitting — a
  broken internal link or missing nav entry fails CI.

## Clinical accuracy review process

Changes to `alleleiq.scoring.criteria.*`, `alleleiq.scoring.bayesian`,
or `alleleiq.scoring.coverage` should, where practical, be reviewed by
someone with clinical genomics domain knowledge in addition to the
usual code review — the correctness bar here is "does this reflect
ACGS 2024 guidance accurately," not just "does this pass tests." If you
are not confident the change is clinically accurate, say so explicitly
in the PR description and request a domain-knowledgeable reviewer
rather than merging on code review alone.

## Code style

- `black` (line length 100) and `ruff` are enforced by pre-commit and
  CI. Run `make format` to apply both.
- `mypy --strict` must pass with zero errors — no `# type: ignore`
  without a comment explaining why the suppression is safe (see
  `alleleiq/models/report.py` for an example of a justified one).
- No `print()` statements — use `alleleiq.logging.get_logger`.
- No magic numbers in clinical logic — every threshold should be a
  named constant or, if it's something a deployment might reasonably
  need to tune, a `Settings` field with a docstring explaining its
  clinical or operational rationale.
