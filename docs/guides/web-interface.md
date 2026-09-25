# Web Interface

A minimal, server-rendered upload-and-review interface, mounted at the
application root alongside the versioned JSON REST API
(`/api/v1/...`, see [API Reference](api-reference.md)). It is a
presentation layer over the same `alleleiq.pipeline.build_evidence_report`
entry point the CLI and JSON API use — it cannot compute a different
result than either of them.

!!! warning "No authentication in v0.2.0"
    Like the JSON API, the web interface has no authentication layer.
    Do not expose it directly to an untrusted network without an
    authenticating reverse proxy in front of it.

!!! note "No screenshots yet"
    This page describes the interface in words rather than images.
    Capturing real screenshots requires a browser, which the
    environment this interface was built in did not have available.
    If you're a maintainer with a browser handy: run `make run-api`,
    visit `http://localhost:8000/`, and add screenshots of the upload
    form and a completed result page here — see
    [CONTRIBUTING.md](https://github.com/NosakhareOsaro/alleleiq/blob/main/CONTRIBUTING.md).
    Please don't add placeholder or stock images in the meantime; an
    honest textual description is preferable to an image that isn't
    actually of this interface.

## The flow

1. **`GET /`** — an upload form, plus a live connectivity check against
   gnomAD, VEP, and ClinVar (the same three checks
   `/api/v1/health/ready` performs). If any are unreachable, a warning
   explains that analysis will still run but affected criteria will be
   reported as `requires_clinical_input` rather than silently skipped.

   A permanent banner reading **"CLINICAL DECISION SUPPORT ONLY — NOT
   A CLINICAL CLASSIFICATION TOOL"** is part of the shared page layout
   every page in this interface extends — there is no page that can
   render without it, the same non-optional guarantee the JSON/HTML
   reports enforce at the model level (see
   [Safety Design](../clinical/safety-design.md)).

2. **Upload** — a drag-and-drop zone, backed by a real
   `<form method="post" action="/analyse" enctype="multipart/form-data">`.
   The `hx-*` attributes on it are progressive enhancement: with
   JavaScript disabled, the same form still submits and the browser
   navigates to a full result page. With HTMX loaded, the identical
   submission swaps the result into the page in place instead.

3. **Progress** (JavaScript-enhanced only) — while a multi-variant VCF
   is being analysed, a Server-Sent Events feed at
   `GET /analyse/{job_id}/progress` streams the *actual* pipeline
   stage each variant is currently at (`querying gnomAD` →
   `querying VEP` → `querying ClinVar` → `querying REVEL` →
   `querying ClinGen` → `scoring`, per variant, then a final
   `complete` event) — sourced directly from
   `alleleiq.pipeline.build_evidence_report`'s optional `on_progress`
   hook, not a simulated or time-based progress bar. This is
   convenience-only: the final result is identical whether or not
   anything was subscribed to this feed.

4. **Result** — one card per analysed variant: its preliminary
   classification, the clinical disclaimer, the same rendered HTML
   evidence report embedded inline via an iframe (the identical file
   the download link points at — nothing is rendered twice), and
   download links for both the JSON and HTML report files. Any
   variant that failed to analyse is listed separately with the
   reason, rather than silently dropped from the batch.

5. **Download** — `GET /reports/{job_id}/{filename}` serves a
   previously generated report file. `job_id` must be a syntactically
   valid UUID and `filename` must be one of the two files that
   analysis run actually produced; anything else — a malformed
   `job_id`, a path-traversal attempt, an unrecognised extension —
   returns `404` rather than touching the filesystem outside that
   job's own report directory.

## Errors

A VCF that fails to parse, declares an unsupported genome build, or
exceeds the configured upload size limit renders a clinical-grade
error message — specific and actionable, but never a raw exception
string or stack trace — instead of a generic failure page.

## Configuration

`ALLELEIQ_WEB_REPORTS_DIR` (default `./reports/web`) controls where
the web interface writes generated reports, one subdirectory per
analysis job. This is independent of the CLI's `--output-dir`, which
defaults to `./reports`.

## Known limitations

- The Server-Sent Events progress feed uses an in-memory,
  per-process queue keyed by a client-generated job ID — the same
  kind of single-process scope limitation already documented for the
  API's rate limiter (see `alleleiq.api.middleware`). A multi-worker
  deployment would need a shared pub/sub layer (e.g. Redis) for a
  progress subscription to reliably reach a worker other than the one
  handling the corresponding `POST /analyse`.
- The interface has no authentication, matching the JSON API's
  current scope (see `SECURITY.md`).
