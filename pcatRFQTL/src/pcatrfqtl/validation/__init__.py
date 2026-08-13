"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/__init__.py

Description:
    Validation utilities and dataset validators for the
    PCa-tRFQTL research pipeline.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.validation.genomic import (
    GenomicCoordinate,
    is_valid_genomic_coordinate,
    parse_genomic_coordinate,
)
from pcatrfqtl.validation.validators.trfqtl import (
    TRFQTLValidator,
)

from pcatrfqtl.validation.anomaly import (
    Anomaly,
    AnomalyFactory,
)

__all__ = [
    "GenomicCoordinate",
    "TRFQTLValidator",
    "is_valid_genomic_coordinate",
    "parse_genomic_coordinate",
    "Anomaly",
    "AnomalyFactory",
]