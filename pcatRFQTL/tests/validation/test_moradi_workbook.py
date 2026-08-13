"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_moradi_workbook.py

Description:
    Unit tests for Moradi workbook validation helper methods.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import pandas as pd

from pcatrfqtl.validation.runners.moradi_workbook import (
    MoradiWorkbookValidator,
)


from pcatrfqtl.validation.validators.moradi import (
    MoradiQTLValidator,
)


def test_clean_dataframe() -> None:
    """Completely empty rows should be removed."""

    dataframe = pd.DataFrame(
        [
            ["rs1", 1],
            [None, None],
            ["rs2", 2],
        ],
        columns=["SNP", "value"],
    )

    cleaned = (
        MoradiWorkbookValidator
        ._clean_dataframe(dataframe)
    )

    assert len(cleaned) == 2


def test_missing_value_statistics() -> None:
    """Missing values should be counted per column."""

    dataframe = pd.DataFrame(
        {
            "A": [1, None, 3],
            "B": ["x", "y", None],
        }
    )

    result = (
        MoradiWorkbookValidator
        ._missing_value_statistics(dataframe)
    )

    assert result == {
        "A": 1,
        "B": 1,
    }


def test_unique_statistics() -> None:
    """Unique values should be counted correctly."""

    dataframe = pd.DataFrame(
        {
            "SNP": [
                "rs1",
                "rs2",
                "rs1",
            ],
            "Cancer": [
                "prostate_cancer",
                "prostate_cancer",
                "breast_cancer",
            ],
        }
    )

    result = (
        MoradiWorkbookValidator
        ._unique_statistics(
            dataframe,
            [
                "SNP",
                "Cancer",
            ],
        )
    )

    assert result == {
        "SNP": 2,
        "Cancer": 2,
    }


def test_duplicate_statistics() -> None:
    """Exact duplicate rows should be counted."""

    dataframe = pd.DataFrame(
        {
            "A": [1, 1, 2],
            "B": ["x", "x", "y"],
        }
    )

    result = (
        MoradiWorkbookValidator
        ._duplicate_statistics(dataframe)
    )

    assert result["duplicate_rows"] == 2
    assert result["duplicate_excess_rows"] == 1


def test_required_columns() -> None:
    """Missing required columns should be reported."""

    dataframe = pd.DataFrame(
        {
            "SNP": ["rs1"],
            "SNP_pos": ["1:100"],
        }
    )

    result = (
        MoradiWorkbookValidator
        ._validate_required_columns(
            dataframe,
            [
                "SNP",
                "SNP_pos",
                "beta",
            ],
        )
    )

    assert result == ["beta"]

def test_status_from_anomalies_pass() -> None:
    """No anomalies should produce PASS."""

    result = (
        MoradiWorkbookValidator
        ._status_from_anomalies([])
    )

    assert result == "PASS"


def test_status_from_anomalies_warning() -> None:
    """Warning-only anomalies should produce WARNING."""

    anomalies = [
        {
            "severity": "WARNING",
        }
    ]

    result = (
        MoradiWorkbookValidator
        ._status_from_anomalies(
            anomalies
        )
    )

    assert result == "WARNING"


def test_status_from_anomalies_fail() -> None:
    """Any error anomaly should produce FAIL."""

    anomalies = [
        {
            "severity": "WARNING",
        },
        {
            "severity": "ERROR",
        },
    ]

    result = (
        MoradiWorkbookValidator
        ._status_from_anomalies(
            anomalies
        )
    )

    assert result == "FAIL"


def test_known_malformed_intron_event() -> None:
    """Scientific-notation-like intron IDs should be recognized."""

    assert (
        MoradiWorkbookValidator
        ._is_known_malformed_intron_event(
            "INT1e+05"
        )
    )

    assert not (
        MoradiWorkbookValidator
        ._is_known_malformed_intron_event(
            "INT132623"
        )
    )

def test_s3_trans_exon_uses_t_as_snp_column() -> None:
    """S3 trans-exon source schema should preserve the T column."""

    schema = (
        MoradiWorkbookValidator
        .S3_SHEETS[
            "Trans-exon-sQTl"
        ]
    )

    assert schema["snp_column"] == "T"


def test_s4_uses_composite_gwas_tag_ids() -> None:
    """S4 tag SNP identifiers may contain composite dbSNP IDs."""

    validator = MoradiQTLValidator()

    result = (
        validator
        .validate_composite_snp_ids(
            [
                "rs376652476:rs792460",
                "rs544592782:rs3081977",
            ],
            "tag_SNP_pos",
        )
    )

    assert result == []