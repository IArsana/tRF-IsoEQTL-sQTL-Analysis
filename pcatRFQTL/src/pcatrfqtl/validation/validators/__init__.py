"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/validators/__init__.py

Description:
    Validation classes for project datasets.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.validation.validators.trfqtl import TRFQTLValidator

__all__ = [
    "TRFQTLValidator",
]