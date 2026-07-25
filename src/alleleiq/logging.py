"""Structured JSON logging for AlleleIQ.

All application logs are emitted as single-line JSON objects to
``stderr``. ``stdout`` is reserved exclusively for report output and CLI
results, so logs never contaminate machine-readable output that a
downstream LIMS or shell pipeline might parse.

Clinical safety note
---------------------
At ``INFO`` level and above, log records must contain only process
metadata: variant identifiers (chrom/pos/ref/alt or HGVS notation),
pipeline stage names, counts, and timings. They must **never** contain
allele frequencies, ACMG criterion outcomes, preliminary classifications,
or anything else that could reveal a patient's genomic interpretation
result. That category of detail is permitted only at ``DEBUG``, which is
expected to be enabled solely in controlled development environments and
never in a deployment that processes real patient data.
"""

from __future__ import annotations

import json
import logging
import sys
from datetime import UTC, datetime
from typing import Any

_RESERVED_LOG_RECORD_ATTRS = frozenset(logging.LogRecord("", 0, "", 0, "", (), None).__dict__)


class JSONFormatter(logging.Formatter):
    """Render log records as single-line JSON objects.

    Any keyword arguments passed via ``extra={...}`` on the logging call
    are nested under a ``context`` key, keeping the structured payload
    separate from the fixed set of standard fields.
    """

    def format(self, record: logging.LogRecord) -> str:
        """Serialise a log record to a JSON string.

        Args:
            record: The log record to format.

        Returns:
            A single line of JSON containing timestamp, level, module,
            function, line number, message, and any extra context.
        """
        payload: dict[str, Any] = {
            "timestamp": datetime.fromtimestamp(record.created, tz=UTC).isoformat(),
            "level": record.levelname,
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
            "message": record.getMessage(),
        }

        context = {
            key: value
            for key, value in record.__dict__.items()
            if key not in _RESERVED_LOG_RECORD_ATTRS
        }
        if context:
            payload["context"] = context

        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)

        return json.dumps(payload, default=str, sort_keys=True)


def configure_logging(level: str = "INFO") -> None:
    """Configure the root ``alleleiq`` logger for structured JSON output.

    Idempotent: safe to call multiple times (e.g. once at CLI startup and
    once in API lifespan) without producing duplicate handlers.

    Args:
        level: One of ``DEBUG``, ``INFO``, ``WARNING``, ``ERROR``.
    """
    root = logging.getLogger("alleleiq")
    root.setLevel(level.upper())
    root.propagate = False

    for handler in list(root.handlers):
        root.removeHandler(handler)

    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setFormatter(JSONFormatter())
    root.addHandler(handler)


def get_logger(name: str) -> logging.Logger:
    """Return a namespaced logger under the ``alleleiq`` hierarchy.

    Args:
        name: Typically ``__name__`` of the calling module.

    Returns:
        A standard library :class:`logging.Logger` configured to emit
        structured JSON via the handler installed by
        :func:`configure_logging`.
    """
    if not name.startswith("alleleiq"):
        name = f"alleleiq.{name}"
    return logging.getLogger(name)
