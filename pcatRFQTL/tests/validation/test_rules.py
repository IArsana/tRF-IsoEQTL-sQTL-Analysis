"""
Tests for generic validation rules.

Author:
    I Putu Indra Arsana
"""

from pathlib import Path

from pcatrfqtl.validation.result import ValidationStatus
from pcatrfqtl.validation.rules import (
    validate_file_exists,
    validate_no_duplicates,
    validate_numeric_range,
    validate_required_columns,
)


def test_validate_file_exists(tmp_path: Path) -> None:
    """
    Test file existence validation.
    """

    test_file = tmp_path / "test.txt"
    test_file.write_text("test", encoding="utf-8")

    result = validate_file_exists(test_file)

    assert result.status == ValidationStatus.PASS


def test_validate_required_columns() -> None:
    """
    Test required column validation.
    """

    result = validate_required_columns(
        actual_columns=[
            "SNP",
            "P-VALUE",
            "GENE",
        ],
        required_columns=[
            "SNP",
            "P-VALUE",
        ],
    )

    assert result.status == ValidationStatus.PASS


def test_validate_required_columns_missing() -> None:
    """
    Test missing required columns.
    """

    result = validate_required_columns(
        actual_columns=[
            "SNP",
        ],
        required_columns=[
            "SNP",
            "P-VALUE",
        ],
    )

    assert result.status == ValidationStatus.FAIL


def test_validate_numeric_range() -> None:
    """
    Test numeric range validation.
    """

    result = validate_numeric_range(
        values=[
            0.001,
            0.05,
            1.0,
        ],
        minimum=0.0,
        maximum=1.0,
        name="p_value",
    )

    assert result.status == ValidationStatus.PASS


def test_validate_numeric_range_invalid() -> None:
    """
    Test invalid numeric range values.
    """

    result = validate_numeric_range(
        values=[
            0.001,
            2.0,
        ],
        minimum=0.0,
        maximum=1.0,
        name="p_value",
    )

    assert result.status == ValidationStatus.FAIL


def test_validate_no_duplicates() -> None:
    """
    Test duplicate detection.
    """

    result = validate_no_duplicates(
        values=[
            "rs1",
            "rs2",
            "rs3",
        ],
        name="snp_id",
    )

    assert result.status == ValidationStatus.PASS