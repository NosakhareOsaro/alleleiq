# Security Policy

## Genomic data sensitivity

AlleleIQ processes genomic variant data that may be linked to real
patients in a clinical deployment. Treat any security issue that could
expose variant data, patient identifiers, or database query patterns
that reveal what variants were analysed as **high severity**, even if
it would normally be classified lower in a non-clinical context.

Specific areas of concern for this codebase:

- **Logging**: `alleleiq.logging` is designed so that INFO-level and
  above logs never carry allele frequencies, classifications, or
  criterion outcomes (see `docs/clinical/safety-design.md` §6). A
  regression that leaks this into INFO-level logs is a security issue,
  not just a design deviation.
- **Report files**: JSON/HTML reports written to disk (`--output-dir`,
  or the API's response body) contain the full evidence picture for a
  variant. Deployments are responsible for access-controlling the
  filesystem/network path these are written to or transmitted over —
  AlleleIQ itself does not encrypt report output at rest.
- **API authentication**: `alleleiq.api.dependencies.require_authenticated_caller`
  is a documented no-op in v0.1.0-alpha. Do not expose the API directly
  to an untrusted network without an authenticating layer in front of
  it. This is a known, intentional gap (see `docs/guides/api-reference.md`),
  not something to report as a new vulnerability — but deployment
  guidance questions about it are welcome via the contact below.
- **XML parsing**: ClinVar responses are parsed with `defusedxml`
  specifically to prevent XXE-class attacks from a compromised or
  spoofed upstream response — do not reintroduce `xml.etree.ElementTree`
  for parsing external XML.

## Supported versions

AlleleIQ is in active alpha development (`0.x`). Security fixes are
applied to the latest released version only; there is no long-term
support branch at this stage.

## Reporting a vulnerability

Please report security issues privately rather than opening a public
GitHub issue: **osaronosakhare@gmail.com**.

Include:

- A description of the issue and its potential impact (with particular
  attention to whether it could expose genomic/patient data, per
  above).
- Steps to reproduce, or a minimal proof of concept.
- The AlleleIQ version affected (`uv run alleleiq version`).

We aim to acknowledge reports within 5 business days and to provide a
remediation timeline once the issue is confirmed. Please allow a
reasonable disclosure window before publishing details publicly.

## Dependency vulnerabilities

Dependency vulnerabilities are scanned on every CI run via `pip-audit`
(see `.github/workflows/ci.yml`). If you discover a vulnerable
dependency that CI hasn't caught, please report it the same way as
above rather than opening a public issue, in case it's exploitable
against a deployed instance before a fix is released.
