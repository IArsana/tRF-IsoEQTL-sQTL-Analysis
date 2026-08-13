"""
Generic validation rules.

This module contains reusable validation rules that can be applied
to multiple dataset types in the PCa-tRFQTL pipeline.

Dataset-specific rules should be implemented in dedicated validator
modules rather than added to this file.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

from pcatrfqtl.validation.result import (
    ValidationCheck,
    ValidationSeverity,
    ValidationStatus,
)


def validate_file_exists(
    path: Path,
) -> ValidationCheck:
    """
    Validate that a file exists.
    """

    if path.exists():
        return ValidationCheck(
            name="file_exists",
            status=ValidationStatus.PASS,
            severity=ValidationSeverity.INFO,
            message=f"File exists: {path}",
        )

    return ValidationCheck(
        name="file_exists",
        status=ValidationStatus.FAIL,
        severity=ValidationSeverity.ERROR,
        message=f"File does not exist: {path}",
    )


def validate_required_columns(
    actual_columns: Iterable[str],
    required_columns: Iterable[str],
) -> ValidationCheck:
    """
    Validate that all required columns are present.
    """

    actual = {
        str(column).strip()
        for column in actual_columns
    }

    required = {
        str(column).strip()
        for column in required_columns
    }

    missing = sorted(required - actual)

    if not missing:
        return ValidationCheck(
            name="required_columns",
            status=ValidationStatus.PASS,
            severity=ValidationSeverity.INFO,
            message="All required columns are present.",
            details={
                "required_columns": sorted(required),
            },
        )

    return ValidationCheck(
        name="required_columns",
        status=ValidationStatus.FAIL,
        severity=ValidationSeverity.ERROR,
        message=(
            "Required columns are missing: "
            f"{missing}"
        ),
        details={
            "missing_columns": missing,
            "required_columns": sorted(required),
            "available_columns": sorted(actual),
        },
    )


def validate_non_empty_columns(
    columns: Iterable[str],
) -> ValidationCheck:
    """
    Validate that column names are not empty.
    """

    empty_columns = [
        column
        for column in columns
        if not str(column).strip()
    ]

    if not empty_columns:
        return ValidationCheck(
            name="non_empty_columns",
            status=ValidationStatus.PASS,
            severity=ValidationSeverity.INFO,
            message="All column names are non-empty.",
        )

    return ValidationCheck(
        name="non_empty_columns",
        status=ValidationStatus.FAIL,
        severity=ValidationSeverity.ERROR,
        message="One or more column names are empty.",
        details={
            "empty_columns": empty_columns,
        },
    )


def validate_numeric_range(
    values: Iterable[Any],
    minimum: float,
    maximum: float,
    name: str,
) -> ValidationCheck:
    """
    Validate that numeric values fall within an allowed range.

    Non-numeric values are considered invalid.
    """

    invalid_values: list[Any] = []
    total_values = 0

    for value in values:
        total_values += 1

        if value is None:
            continue

        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            invalid_values.append(value)
            continue

        if not minimum <= numeric_value <= maximum:
            invalid_values.append(value)

    if not invalid_values:
        return ValidationCheck(
            name=name,
            status=ValidationStatus.PASS,
            severity=ValidationSeverity.INFO,
            message=(
                f"All values in '{name}' are within "
                f"[{minimum}, {maximum}]."
            ),
            details={
                "total_values": total_values,
                "minimum": minimum,
                "maximum": maximum,
            },
        )

    return ValidationCheck(
        name=name,
        status=ValidationStatus.FAIL,
        severity=ValidationSeverity.ERROR,
        message=(
            f"Invalid values detected in '{name}'."
        ),
        details={
            "invalid_count": len(invalid_values),
            "total_values": total_values,
            "minimum": minimum,
            "maximum": maximum,
        },
    )


def validate_no_duplicates(
    values: Iterable[Any],
    name: str,
) -> ValidationCheck:
    """
    Validate that a collection contains no duplicate values.
    """

    values_list = list(values)
    unique_values = set(values_list)

    duplicate_count = (
        len(values_list) - len(unique_values)
    )

    if duplicate_count == 0:
        return ValidationCheck(
            name=name,
            status=ValidationStatus.PASS,
            severity=ValidationSeverity.INFO,
            message=f"No duplicate values detected in '{name}'.",
            details={
                "total_values": len(values_list),
            },
        )

    return ValidationCheck(
        name=name,
        status=ValidationStatus.WARNING,
        severity=ValidationSeverity.WARNING,
        message=(
            f"Duplicate values detected in '{name}'."
        ),
        details={
            "total_values": len(values_list),
            "unique_values": len(unique_values),
            "duplicate_count": duplicate_count,
        },
    )