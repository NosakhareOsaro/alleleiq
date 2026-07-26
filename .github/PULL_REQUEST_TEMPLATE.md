## Summary

<!-- What does this PR change, and why? -->

## Type of change

- [ ] `feat` — new feature
- [ ] `fix` — bug fix
- [ ] `docs` — documentation only
- [ ] `test` — tests only
- [ ] `refactor` — no behaviour change
- [ ] `chore` / `ci` — tooling
- [ ] `perf` — performance improvement

## Clinical/scoring impact

- [ ] This PR touches `alleleiq.scoring.*`, `alleleiq.clients.*`, or
      criterion thresholds in `alleleiq.config`.
- [ ] If checked above: I ran `uv run python scripts/validate_known_variants.py`
      after this change and confirmed no new **dangerous** (opposite-
      direction) outcomes versus the ClinVar gold-standard set.

## Testing

- [ ] `make test` passes locally.
- [ ] `make lint` and `make type-check` pass locally.
- [ ] New/changed criterion evaluators have unit tests covering both
      the "met" and "not met" paths, plus `requires_clinical_input`
      branches (see CONTRIBUTING.md).
- [ ] New external API integrations have unit tests with mocked
      responses covering success, not-found, malformed response, and
      network-failure cases.

## Documentation

- [ ] Relevant `docs/` pages updated (if this changes criteria, data
      models, or the API surface).
- [ ] `make docs` (`mkdocs build --strict`) passes locally.

## Checklist

- [ ] Commit messages follow [Conventional Commits](https://www.conventionalcommits.org/).
- [ ] No new `mypy --strict` errors or `# type: ignore` without an
      explanatory comment.
- [ ] No magic numbers introduced in clinical logic (named constants or
      `Settings` fields only).
