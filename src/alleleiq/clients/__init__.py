"""Async HTTP clients for external clinical genomics evidence sources."""

from alleleiq.clients.base import BaseClient
from alleleiq.clients.clinvar import ClinVarClient, ClinVarSubmission
from alleleiq.clients.gnomad import GeneConstraint, GnomADClient
from alleleiq.clients.revel import REVELClient, REVELScore
from alleleiq.clients.uniprot import ProteinRegion, UniProtClient
from alleleiq.clients.vep import VEPAnnotation, VEPClient

__all__ = [
    "BaseClient",
    "ClinVarClient",
    "ClinVarSubmission",
    "GeneConstraint",
    "GnomADClient",
    "ProteinRegion",
    "REVELClient",
    "REVELScore",
    "UniProtClient",
    "VEPAnnotation",
    "VEPClient",
]
