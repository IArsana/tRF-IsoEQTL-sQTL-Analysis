"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/__init__.py

Description:
    Execution runners for M3 scientific analyses.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.analysis.m3.runners.extract_gwas import (
    M31GWASInputs,
    M31GWASRunner,
)
from pcatrfqtl.analysis.m3.runners.extract_trfqtl import (
    M32TRFQTLInputs,
    M32TRFQTLRunner,
)
from pcatrfqtl.analysis.m3.runners.direct_overlap import (
    M33OverlapInputs,
    M33OverlapRunner,
)
from pcatrfqtl.analysis.m3.runners.build_candidates import (
    M34CandidateInputs,
    M34CandidateRunner,
)
from pcatrfqtl.analysis.m3.runners.summarize_statistics import (
    M35StatisticsInputs,
    M35StatisticsRunner,
)
from pcatrfqtl.analysis.m3.runners.build_figures import (
    M36FigureInputs,
    M36FigureRunner,
)

__all__ = [
    "M31GWASInputs",
    "M31GWASRunner",
    "M32TRFQTLInputs",
    "M32TRFQTLRunner",
    "M31GWASInputs",
    "M31GWASRunner",
    "M32TRFQTLInputs",
    "M32TRFQTLRunner",
    "M33OverlapInputs",
    "M33OverlapRunner",
    "M31GWASInputs",
    "M31GWASRunner",
    "M32TRFQTLInputs",
    "M32TRFQTLRunner",
    "M33OverlapInputs",
    "M33OverlapRunner",
    "M34CandidateInputs",
    "M34CandidateRunner",
    "M31GWASInputs",
    "M31GWASRunner",
    "M32TRFQTLInputs",
    "M32TRFQTLRunner",
    "M33OverlapInputs",
    "M33OverlapRunner",
    "M34CandidateInputs",
    "M34CandidateRunner",
    "M35StatisticsInputs",
    "M35StatisticsRunner",
    "M36FigureInputs",
    "M36FigureRunner",
]