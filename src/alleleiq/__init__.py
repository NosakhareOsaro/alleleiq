"""AlleleIQ: population-inclusive evidence aggregation for clinical genomic
variant interpretation.

AlleleIQ aggregates evidence from population frequency databases, in
silico prediction tools, and curated variant databases; applies ACMG/AMP
and ACGS 2024 classification criteria; and produces structured,
ancestry-stratified evidence reports for review by a qualified Clinical
Genetic Scientist. It does not itself classify variants.
"""

from alleleiq.__version__ import VERSION
from alleleiq.exceptions import (
    AlleleIQError,
    ClinVarError,
    ConfigurationError,
    DatabaseQueryError,
    GnomADError,
    InsufficientEvidenceError,
    ReportGenerationError,
    ScoringError,
    UniProtError,
    VariantParseError,
    VEPError,
)

__version__ = VERSION

__all__ = [
    "AlleleIQError",
    "ClinVarError",
    "ConfigurationError",
    "DatabaseQueryError",
    "GnomADError",
    "InsufficientEvidenceError",
    "ReportGenerationError",
    "ScoringError",
    "UniProtError",
    "VEPError",
    "VariantParseError",
    "__version__",
]
