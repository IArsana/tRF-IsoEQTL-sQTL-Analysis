"""
Data models for dataset inspection.

This module defines structured representations of files and datasets
discovered during the M1.1 dataset inspection stage.

Author:
    I Putu Indra Arsana
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class FileInspection:
    """
    Store inspection results for a single data file.
    """

    dataset_id: str
    file_id: str
    path: str
    exists: bool
    format: str | None = None
    size_bytes: int | None = None
    size_mb: float | None = None
    readable: bool = False
    row_count: int | None = None
    column_count: int | None = None
    columns: list[str] = field(default_factory=list)
    sheets: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """
        Convert the inspection result into a serializable dictionary.
        """

        return asdict(self)