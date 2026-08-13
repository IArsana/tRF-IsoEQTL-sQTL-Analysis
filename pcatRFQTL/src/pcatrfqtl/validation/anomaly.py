"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/anomaly.py

Description:
    Provides standardized anomaly representations for validation
    results throughout the pcatRFQTL research pipeline.

    Anomaly identifiers are generated using a consistent sequential
    format (ANOM-000001, ANOM-000002, ...).

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
from typing import Any


@dataclass(slots=True)
class Anomaly:
    """
    Represent a single data-quality anomaly.
    """

    anomaly_id: str
    type: str
    category: str
    severity: str
    dataset: str
    file: str | None
    sheet: str | None
    row: int | None
    column: str | None
    value: Any
    message: str
    action: str

    def to_dict(self) -> dict[str, Any]:
        """Convert the anomaly into a serializable dictionary."""
        return asdict(self)


class AnomalyFactory:
    """
    Generate standardized anomaly identifiers.
    """

    def __init__(self, start: int = 1) -> None:
        self._counter = start

    def create(
        self,
        *,
        type: str,
        category: str,
        severity: str,
        dataset: str,
        file: str | None = None,
        sheet: str | None = None,
        row: int | None = None,
        column: str | None = None,
        value: Any = None,
        message: str,
        action: str,
    ) -> Anomaly:
        """
        Create a new anomaly with a unique identifier.
        """

        anomaly_id = f"ANOM-{self._counter:06d}"
        self._counter += 1

        return Anomaly(
            anomaly_id=anomaly_id,
            type=type,
            category=category,
            severity=severity,
            dataset=dataset,
            file=file,
            sheet=sheet,
            row=row,
            column=column,
            value=value,
            message=message,
            action=action,
        )