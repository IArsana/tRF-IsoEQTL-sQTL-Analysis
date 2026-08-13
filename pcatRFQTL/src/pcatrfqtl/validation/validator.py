"""
Base dataset validator.

This module provides the base validation interface used by
dataset-specific validators throughout the PCa-tRFQTL pipeline.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from pcatrfqtl.validation.result import (
    ValidationReport,
)


class BaseValidator(ABC):
    """
    Abstract base class for dataset validators.
    """

    def __init__(
        self,
        dataset_id: str,
        dataset_config: dict[str, Any],
    ) -> None:

        self.dataset_id = dataset_id
        self.dataset_config = dataset_config

    @abstractmethod
    def validate(self) -> ValidationReport:
        """
        Execute dataset-specific validation rules.

        Returns
        -------
        ValidationReport
            Validation results for the dataset.
        """

        raise NotImplementedError

    @staticmethod
    def resolve_path(
        project_root: Path,
        relative_path: str,
    ) -> Path:
        """
        Resolve a project-relative dataset path.
        """

        return project_root / relative_path