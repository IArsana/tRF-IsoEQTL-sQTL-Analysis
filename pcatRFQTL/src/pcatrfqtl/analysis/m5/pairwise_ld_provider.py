"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/pairwise_ld_provider.py

Description:
    Provider abstraction for M5.3C.3G pairwise LD queries.

    The interface retrieves LD for one canonical rsID pair and one reference
    population.

    Provider-specific implementations should return explicit unavailable or
    failed states rather than silently treating missing data as r² = 0.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True)
class PairwiseLDRecord:
    """Pairwise LD result for one variant pair and population."""

    variant_a: str

    variant_b: str

    population: str

    r2: float | None

    d_prime: float | None

    status: str

    provider: str

    reason: str | None = None


class PairwiseLDProvider(Protocol):
    """Protocol implemented by pairwise LD providers."""

    def query_pair(
        self,
        *,
        variant_a: str,
        variant_b: str,
        population: str,
    ) -> PairwiseLDRecord:
        """Query pairwise LD for one rsID pair."""
        ...