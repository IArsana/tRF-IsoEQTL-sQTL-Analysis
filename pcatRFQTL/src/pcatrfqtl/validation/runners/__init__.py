"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/runners/__init__.py

Description:
    Validation runners for project datasets.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.validation.runners.trfqtl_workbook import (
    TRFQTLWorkbookValidator,
)

__all__ = [
    "TRFQTLWorkbookValidator",
]