"""
Low-level file inspection utilities.

This module inspects raw data files without modifying them. It supports
CSV, TSV, XLSX, and XLS-style tabular sources where supported by the
installed Python libraries.

Author:
    I Putu Indra Arsana
"""

from __future__ import annotations

import csv
from pathlib import Path

from openpyxl import load_workbook

SUPPORTED_FORMATS = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".xlsx": "xlsx",
}


def detect_format(path: Path) -> str | None:
    """
    Detect the supported file format from the file extension.
    """

    return SUPPORTED_FORMATS.get(path.suffix.lower())


def inspect_delimited_file(
    path: Path,
    delimiter: str,
) -> tuple[int, list[str]]:
    """
    Inspect a CSV/TSV file.

    Returns:
        Tuple containing row count and column names.
    """

    row_count = 0

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.reader(handle, delimiter=delimiter)

        try:
            header = next(reader)
        except StopIteration:
            return 0, []

        for _ in reader:
            row_count += 1

    return row_count, [column.strip() for column in header]


def inspect_excel_file(
    path: Path,
) -> tuple[list[str], dict[str, dict[str, object]]]:
    """
    Inspect an Excel workbook without loading cell data into memory.

    Returns:
        A tuple containing:
            - sheet names
            - sheet metadata
    """

    workbook = load_workbook(
        filename=path,
        read_only=True,
        data_only=False,
    )

    sheets = workbook.sheetnames
    sheet_metadata: dict[str, dict[str, object]] = {}

    for sheet_name in sheets:
        worksheet = workbook[sheet_name]

        rows = worksheet.iter_rows(values_only=True)

        try:
            header_row = next(rows)
        except StopIteration:
            header_row = ()

        columns = [
            str(value).strip()
            for value in header_row
            if value is not None
        ]

        row_count = max(worksheet.max_row - 1, 0)

        sheet_metadata[sheet_name] = {
            "rows": row_count,
            "columns": len(columns),
            "column_names": columns,
            "max_column": worksheet.max_column,
            "max_row": worksheet.max_row,
        }

    workbook.close()

    return sheets, sheet_metadata