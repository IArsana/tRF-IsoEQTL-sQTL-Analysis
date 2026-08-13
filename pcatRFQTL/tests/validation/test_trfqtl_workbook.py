"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_trfqtl_workbook.py

Description:
    Unit tests for the Cancer-tRFQTL workbook validation runner.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pathlib import Path

import pandas as pd

from pcatrfqtl.validation.runners.trfqtl_workbook import (
    TRFQTLWorkbookValidator,
)


def test_clean_dataframe() -> None:
    """Test dataframe cleaning."""

    dataframe = pd.DataFrame(
        [
            ["A", 1],
            [None, None],
            ["B", 2],
        ],
        columns=["type", "value"],
    )

    cleaned = (
        TRFQTLWorkbookValidator
        ._clean_dataframe(dataframe)
    )

    assert len(cleaned) == 2
    assert list(cleaned.columns) == [
        "type",
        "value",
    ]


def test_missing_value_statistics() -> None:
    """Test missing-value statistics."""

    dataframe = pd.DataFrame(
        {
            "A": [1, None, 3],
            "B": ["x", "y", None],
        }
    )

    result = (
        TRFQTLWorkbookValidator
        ._missing_value_statistics(dataframe)
    )

    assert result["A"] == 1
    assert result["B"] == 1


def test_duplicate_statistics() -> None:
    """Test duplicate-row statistics."""

    dataframe = pd.DataFrame(
        {
            "SNP": [
                "rs1",
                "rs2",
                "rs1",
            ],
        }
    )

    result = (
        TRFQTLWorkbookValidator
        ._duplicate_statistics(dataframe)
    )

    assert result["duplicate_rows"] == 2
    assert result["duplicate_groups"] == 1


def test_unique_statistics() -> None:
    """Test unique-value statistics."""

    dataframe = pd.DataFrame(
        {
            "Cancer": [
                "PRAD",
                "PRAD",
                "BRCA",
            ],
            "SNP": [
                "rs1",
                "rs2",
                "rs1",
            ],
        }
    )

    result = (
        TRFQTLWorkbookValidator
        ._unique_statistics(
            dataframe,
            ["Cancer", "SNP"],
        )
    )

    assert result["Cancer"] == 2
    assert result["SNP"] == 2


def test_missing_workbook() -> None:
    """Test missing workbook handling."""

    validator = TRFQTLWorkbookValidator(
        Path("does-not-exist.xlsx")
    )

    try:
        validator.validate_all()
    except FileNotFoundError:
        return

    raise AssertionError(
        "Expected FileNotFoundError"
    )