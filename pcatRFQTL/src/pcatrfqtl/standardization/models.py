"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/models.py

Description:
    Shared models for source-independent data standardization.

    These models provide stable representations of standardization
    status, provenance, and transformed values across GWAS Catalog,
    Cancer-tRFQTL, and Moradi-derived datasets.

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
from enum import Enum
from typing import Any


class StandardizationStatus(str, Enum):
    """Outcome of a single standardization operation."""

    STANDARDIZED = "STANDARDIZED"
    PRESERVED = "PRESERVED"
    MISSING = "MISSING"
    UNSUPPORTED = "UNSUPPORTED"


@dataclass(frozen=True)
class StandardizationResult:
    """
    Represent the result of one standardization operation.

    Attributes
    ----------
    raw_value:
        Original source value.

    standardized_value:
        Derived standardized representation.

    status:
        Standardization outcome.

    source:
        Description of how the standardized value was derived.

    usable:
        Whether the standardized representation is suitable for
        downstream operations requiring that field.

    note:
        Optional provenance or source-data note.
    """

    raw_value: Any

    standardized_value: Any

    status: StandardizationStatus

    source: str

    usable: bool

    note: str | None = None

    def to_dict(
        self,
    ) -> dict[str, Any]:
        """Return a JSON-serializable dictionary."""

        result = asdict(
            self
        )

        result[
            "status"
        ] = self.status.value

        return result