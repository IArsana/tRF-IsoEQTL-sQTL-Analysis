"""
Dataset registry loader.

This module loads dataset metadata from the central
configs/datasets.yaml configuration file.

Author:
    I Putu Indra Arsana
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


class DatasetRegistry:
    """
    Load and expose dataset metadata from datasets.yaml.
    """

    def __init__(self, config_path: Path) -> None:
        self.config_path = config_path.resolve()

        if not self.config_path.exists():
            raise FileNotFoundError(
                f"Dataset configuration not found: {self.config_path}"
            )

        self.project_root = self.config_path.parent.parent

        with self.config_path.open("r", encoding="utf-8") as handle:
            self.config: dict[str, Any] = yaml.safe_load(handle) or {}

        self.datasets: dict[str, Any] = self.config.get("datasets", {})

    def get_dataset(self, dataset_id: str) -> dict[str, Any]:
        """
        Return metadata for a specific dataset.
        """

        if dataset_id not in self.datasets:
            raise KeyError(
                f"Dataset '{dataset_id}' is not defined in datasets.yaml."
            )

        return self.datasets[dataset_id]

    def iter_datasets(self):
        """
        Iterate through registered datasets.
        """

        yield from self.datasets.items()

    def resolve_path(self, relative_path: str) -> Path:
        """
        Resolve a project-relative dataset path.
        """

        return self.project_root / relative_path