"""
Dataset inspection orchestrator.

This module coordinates file-level inspection across all datasets
registered in configs/datasets.yaml.

Author:
    I Putu Indra Arsana
"""

from __future__ import annotations

import json
from pathlib import Path

from .file_inspector import (
    detect_format,
    inspect_delimited_file,
    inspect_excel_file,
)
from .models import FileInspection
from .yaml_registry import DatasetRegistry


class DatasetInspector:
    """
    Inspect all files registered in the dataset registry.
    """

    def __init__(self, registry: DatasetRegistry) -> None:
        self.registry = registry

    def inspect_file(
        self,
        dataset_id: str,
        file_id: str,
        file_metadata: dict,
    ) -> FileInspection:

        relative_path = file_metadata["path"]
        path = self.registry.resolve_path(relative_path)

        result = FileInspection(
            dataset_id=dataset_id,
            file_id=file_id,
            path=relative_path,
            exists=path.exists(),
        )

        if not path.exists():
            result.errors.append(
                f"File does not exist: {path}"
            )
            return result

        result.size_bytes = path.stat().st_size
        result.size_mb = round(
            result.size_bytes / (1024 * 1024),
            3,
        )

        detected_format = detect_format(path)
        result.format = detected_format

        if detected_format is None:
            result.warnings.append(
                f"Unsupported file extension: {path.suffix}"
            )
            return result

        try:
            if detected_format == "csv":
                rows, columns = inspect_delimited_file(
                    path,
                    delimiter=",",
                )

                result.row_count = rows
                result.column_count = len(columns)
                result.columns = columns

            elif detected_format == "tsv":
                rows, columns = inspect_delimited_file(
                    path,
                    delimiter="\t",
                )

                result.row_count = rows
                result.column_count = len(columns)
                result.columns = columns

            elif detected_format == "xlsx":
                sheets, sheet_metadata = inspect_excel_file(path)

                result.sheets = sheets

                result.warnings.extend(
                    self._validate_excel_metadata(
                        dataset_id,
                        file_id,
                        file_metadata,
                        sheets,
                    )
                )

                if len(sheets) == 1:
                    metadata = sheet_metadata[sheets[0]]

                    result.row_count = metadata["rows"]
                    result.column_count = metadata["columns"]
                    result.columns = metadata["column_names"]

            result.readable = True

        except Exception as exc:
            result.errors.append(
                f"Failed to inspect file: {exc}"
            )

        return result

    def _validate_excel_metadata(
        self,
        dataset_id: str,
        file_id: str,
        file_metadata: dict,
        actual_sheets: list[str],
    ) -> list[str]:

        warnings: list[str] = []

        expected_sheets = file_metadata.get("expected_sheets")

        if expected_sheets is None:
            return warnings

        missing = set(expected_sheets) - set(actual_sheets)

        if missing:
            warnings.append(
                f"Missing expected sheets: {sorted(missing)}"
            )

        return warnings

    def inspect_all(self) -> list[FileInspection]:
        """
        Inspect all registered dataset files.
        """

        results: list[FileInspection] = []

        for dataset_id, dataset in self.registry.iter_datasets():

            files = dataset.get("files", {})

            for file_id, file_metadata in files.items():

                results.append(
                    self.inspect_file(
                        dataset_id=dataset_id,
                        file_id=file_id,
                        file_metadata=file_metadata,
                    )
                )

        return results

    @staticmethod
    def save_json(
        results: list[FileInspection],
        output_path: Path,
    ) -> None:
        """
        Save inspection results as JSON.
        """

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        payload = [
            result.to_dict()
            for result in results
        ]

        with output_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                payload,
                handle,
                indent=2,
                ensure_ascii=False,
            )