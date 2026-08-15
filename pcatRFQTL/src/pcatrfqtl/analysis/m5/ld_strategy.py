"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/ld_strategy.py

Description:
    Defines the M5.1 linkage disequilibrium reference strategy.

    M5.1 is a methodological configuration stage only. It records the
    assumptions and safeguards required before performing locus-aware
    variant integration.

    No LD coefficients are calculated in this module.

    Core principles:

        - Cancer-tRFQTL and Moradi coordinates are treated as hg19 /
          GRCh37 unless otherwise documented by locked harmonization.
        - rsID remains the preferred cross-dataset variant identity.
        - genomic windows are search regions, not LD blocks.
        - source-reported LD and independently calculated LD must remain
          distinguishable.
        - population ancestry of the LD reference panel must be recorded.
        - LD thresholds must be explicit and configurable.
        - no causal interpretation is made from LD alone.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum
from typing import Any


class LDReferenceStatus(StrEnum):
    """Availability state for an external LD reference."""

    REQUIRED = "REQUIRED"
    AVAILABLE = "AVAILABLE"
    UNRESOLVED = "UNRESOLVED"


class LDMetric(StrEnum):
    """Supported LD statistics."""

    R2 = "R2"
    D_PRIME = "D_PRIME"


@dataclass(frozen=True)
class LDReferenceStrategy:
    """
    Analysis policy for M5 locus-aware integration.

    The default screening window is intentionally independent from the
    actual LD threshold.

    Example:
        A variant inside ±500 kb is not considered linked until an
        external/reference LD calculation supports that relationship.
    """

    source_build: str = "hg19"
    reference_assembly: str = "GRCh37"

    preferred_variant_identifier: str = "rsID"

    screening_window_bp: int = 500_000

    primary_ld_metric: LDMetric = LDMetric.R2

    primary_r2_threshold: float = 0.8

    secondary_r2_threshold: float = 0.5

    reference_panel_status: LDReferenceStatus = (
        LDReferenceStatus.REQUIRED
    )

    reference_panel_name: str | None = None

    reference_panel_population: str | None = None

    population_matching_required: bool = True

    multi_population_sensitivity_allowed: bool = True

    source_reported_ld_is_external_evidence: bool = True

    source_reported_ld_recomputed: bool = False

    coordinate_matching_is_ld_evidence: bool = False

    physical_proximity_is_ld_evidence: bool = False

    liftover_performed: bool = False

    ld_calculation_performed: bool = False

    colocalization_performed: bool = False

    causal_inference_performed: bool = False

    def validate(self) -> None:
        """Validate strategy parameters."""

        if self.screening_window_bp <= 0:
            raise ValueError(
                "screening_window_bp must be greater than zero."
            )

        if not (
            0.0
            <= self.primary_r2_threshold
            <= 1.0
        ):
            raise ValueError(
                "primary_r2_threshold must lie between 0 and 1."
            )

        if not (
            0.0
            <= self.secondary_r2_threshold
            <= 1.0
        ):
            raise ValueError(
                "secondary_r2_threshold must lie between 0 and 1."
            )

        if (
            self.secondary_r2_threshold
            > self.primary_r2_threshold
        ):
            raise ValueError(
                "secondary_r2_threshold cannot exceed "
                "primary_r2_threshold."
            )

        if (
            self.reference_panel_status
            == LDReferenceStatus.AVAILABLE
            and self.reference_panel_name is None
        ):
            raise ValueError(
                "An available LD reference panel requires "
                "reference_panel_name."
            )

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-safe strategy representation."""

        self.validate()

        data = asdict(
            self
        )

        data[
            "primary_ld_metric"
        ] = self.primary_ld_metric.value

        data[
            "reference_panel_status"
        ] = self.reference_panel_status.value

        return data


def default_m5_ld_strategy(
    *,
    screening_window_bp: int = 500_000,
    primary_r2_threshold: float = 0.8,
    secondary_r2_threshold: float = 0.5,
) -> LDReferenceStrategy:
    """
    Return conservative M5 LD strategy.

    A specific external reference panel is deliberately not selected at
    M5.1 unless population/reference provenance has been explicitly
    supplied.
    """

    strategy = LDReferenceStrategy(
        screening_window_bp=(
            screening_window_bp
        ),
        primary_r2_threshold=(
            primary_r2_threshold
        ),
        secondary_r2_threshold=(
            secondary_r2_threshold
        ),
    )

    strategy.validate()

    return strategy