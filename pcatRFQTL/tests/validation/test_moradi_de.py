"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_moradi_de.py

Description:
    Unit tests for Moradi differential-expression validation utilities
    and dataset-level runner behavior.

    Tests cover:
        - exon feature identifiers
        - intron feature identifiers
        - transcript identifiers
        - baseMean validation
        - generic numeric validation
        - probability validation
        - missing-value behavior
        - required schema handling
        - duplicate statistics
        - S5-S10 dataset configuration
        - single-file runner validation
        - full runner validation

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pcatrfqtl.validation.runners.moradi_de import (
    MoradiDERunner,
)
from pcatrfqtl.validation.validators.moradi_de import (
    MoradiDEValidator,
)


def create_validator() -> MoradiDEValidator:
    """Create a Moradi differential-expression field validator."""

    return MoradiDEValidator()


# ===========================================================================
# Missing-value handling
# ===========================================================================


def test_is_missing() -> None:
    """Common missing representations should be detected."""

    validator = create_validator()

    assert validator._is_missing(None)
    assert validator._is_missing(np.nan)
    assert validator._is_missing(float("nan"))
    assert validator._is_missing("")
    assert validator._is_missing(" ")
    assert validator._is_missing("NA")
    assert validator._is_missing("nan")
    assert validator._is_missing("None")
    assert validator._is_missing("null")

    assert not validator._is_missing("EX123")
    assert not validator._is_missing(0)
    assert not validator._is_missing(1.0)


# ===========================================================================
# Feature validation
# ===========================================================================


def test_valid_exon_features() -> None:
    """EX-prefixed exon identifiers should pass."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "EX1",
            "EX123",
            "ex999",
        ],
        "exon",
    )

    assert result == []


def test_invalid_exon_features() -> None:
    """Non-EX identifiers should fail exon validation."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "INT1",
            "ENST000003",
            "EXABC",
            "1",
        ],
        "exon",
    )

    assert len(result) == 4


def test_valid_intron_features() -> None:
    """INT-prefixed intron identifiers should pass."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "INT1",
            "INT123",
            "int999",
        ],
        "intron",
    )

    assert result == []


def test_invalid_intron_features() -> None:
    """Non-INT identifiers should fail intron validation."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "EX1",
            "ENST000003",
            "INTABC",
            "1",
        ],
        "intron",
    )

    assert len(result) == 4


def test_valid_transcript_features() -> None:
    """Canonical ENST identifiers with optional version should pass."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "ENST00000318325",
            "ENST00000318325.6",
            "enst00000412345.1",
        ],
        "transcript",
    )

    assert result == []


def test_invalid_transcript_features() -> None:
    """Non-ENST or malformed transcript identifiers should fail."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "EX1",
            "INT1",
            "ENSG000001",
            "ENSTABC",
        ],
        "transcript",
    )

    assert len(result) == 4


def test_missing_feature_invalid_by_default() -> None:
    """Missing feature IDs should fail by default."""

    validator = create_validator()

    result = validator.validate_features(
        [
            None,
            "",
            np.nan,
        ],
        "exon",
    )

    assert len(result) == 3


def test_missing_feature_allowed_when_configured() -> None:
    """Missing feature IDs should be accepted when explicitly allowed."""

    validator = create_validator()

    result = validator.validate_features(
        [
            None,
            "",
            np.nan,
            "EX1",
        ],
        "exon",
        allow_missing=True,
    )

    assert result == []


def test_unknown_feature_type_raises() -> None:
    """Unsupported feature types should raise ValueError."""

    validator = create_validator()

    with pytest.raises(ValueError):
        validator.validate_features(
            ["EX1"],
            "unknown",
        )


def test_feature_error_contains_rule_id() -> None:
    """Feature validation messages should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_features(
        [
            "INVALID",
        ],
        "exon",
        "feature_id",
    )

    assert len(result) == 1
    assert "MORADI-DE-FEATURE-001" in result[0]
    assert "feature_id[0]" in result[0]


# ===========================================================================
# baseMean validation
# ===========================================================================


def test_valid_base_mean() -> None:
    """Non-negative finite baseMean values should pass."""

    validator = create_validator()

    result = validator.validate_base_mean(
        [
            0,
            1,
            10.5,
            "100",
            "0.001",
        ]
    )

    assert result == []


def test_invalid_base_mean() -> None:
    """Negative, non-numeric, or non-finite baseMean values should fail."""

    validator = create_validator()

    result = validator.validate_base_mean(
        [
            -1,
            "invalid",
            math.inf,
            -math.inf,
        ]
    )

    assert len(result) == 4


def test_missing_base_mean_invalid_by_default() -> None:
    """Missing baseMean should fail in strict mode."""

    validator = create_validator()

    result = validator.validate_base_mean(
        [
            None,
            "",
            np.nan,
        ]
    )

    assert len(result) == 3


def test_base_mean_error_contains_rule_id() -> None:
    """baseMean validation errors should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_base_mean(
        [
            -1,
        ]
    )

    assert len(result) == 1
    assert "MORADI-DE-BASEMEAN-001" in result[0]
    assert "baseMean[0]" in result[0]


# ===========================================================================
# Generic numeric validation
# ===========================================================================


def test_valid_numeric_values() -> None:
    """Finite numeric DE statistics should pass."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            -10,
            -0.5,
            0,
            1,
            3.14,
            "5.0",
            "1e-4",
        ],
        "stat",
    )

    assert result == []


def test_invalid_numeric_values() -> None:
    """Non-numeric and non-finite values should fail."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            "invalid",
            math.inf,
            -math.inf,
        ],
        "stat",
    )

    assert len(result) == 3


def test_missing_numeric_invalid_by_default() -> None:
    """Missing numeric values should fail by default."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            None,
            "",
            np.nan,
        ],
        "stat",
    )

    assert len(result) == 3


def test_missing_numeric_allowed() -> None:
    """Missing numeric values should pass when explicitly allowed."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            None,
            "",
            np.nan,
            1.0,
        ],
        "stat",
        allow_missing=True,
    )

    assert result == []


def test_numeric_error_contains_rule_id() -> None:
    """Numeric validation errors should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            "invalid",
        ],
        "stat",
    )

    assert len(result) == 1
    assert "MORADI-DE-NUMERIC-001" in result[0]
    assert "stat[0]" in result[0]


# ===========================================================================
# Probability validation
# ===========================================================================


def test_valid_probability_values() -> None:
    """Probability-like values in 0 <= x <= 1 should pass."""

    validator = create_validator()

    result = validator.validate_probability(
        [
            0,
            1e-300,
            0.05,
            0.5,
            1,
            "0.01",
        ],
        "pvalue",
    )

    assert result == []


def test_invalid_probability_values() -> None:
    """Out-of-range and malformed probability values should fail."""

    validator = create_validator()

    result = validator.validate_probability(
        [
            -0.01,
            1.01,
            "invalid",
            math.inf,
            -math.inf,
        ],
        "pvalue",
    )

    assert len(result) == 5


def test_missing_probability_invalid_by_default() -> None:
    """Missing p-values should fail in strict Moradi DE validation."""

    validator = create_validator()

    result = validator.validate_probability(
        [
            None,
            "",
            np.nan,
        ],
        "pvalue",
    )

    assert len(result) == 3


def test_probability_error_contains_rule_id() -> None:
    """Probability errors should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_probability(
        [
            2.0,
        ],
        "padj",
    )

    assert len(result) == 1
    assert "MORADI-DE-PROBABILITY-001" in result[0]
    assert "padj[0]" in result[0]


# ===========================================================================
# Runner configuration
# ===========================================================================


def test_runner_contains_all_s5_to_s10_datasets() -> None:
    """Runner should define all six Moradi DE supplementary datasets."""

    expected = {
        "S5",
        "S6",
        "S7",
        "S8",
        "S9",
        "S10",
    }

    assert set(
        MoradiDERunner.DATASETS
    ) == expected


def test_runner_feature_type_mapping() -> None:
    """S5-S10 should use the correct biological feature validator."""

    datasets = (
        MoradiDERunner.DATASETS
    )

    assert (
        datasets["S5"][
            "feature_type"
        ]
        == "exon"
    )

    assert (
        datasets["S6"][
            "feature_type"
        ]
        == "exon"
    )

    assert (
        datasets["S7"][
            "feature_type"
        ]
        == "intron"
    )

    assert (
        datasets["S8"][
            "feature_type"
        ]
        == "intron"
    )

    assert (
        datasets["S9"][
            "feature_type"
        ]
        == "transcript"
    )

    assert (
        datasets["S10"][
            "feature_type"
        ]
        == "transcript"
    )


def test_runner_analysis_scope_mapping() -> None:
    """cis/trans metadata should match each supplementary dataset."""

    datasets = (
        MoradiDERunner.DATASETS
    )

    assert (
        datasets["S5"][
            "analysis_scope"
        ]
        == "cis"
    )

    assert (
        datasets["S6"][
            "analysis_scope"
        ]
        == "trans"
    )

    assert (
        datasets["S7"][
            "analysis_scope"
        ]
        == "cis"
    )

    assert (
        datasets["S8"][
            "analysis_scope"
        ]
        == "trans"
    )

    assert (
        datasets["S9"][
            "analysis_scope"
        ]
        == "cis"
    )

    assert (
        datasets["S10"][
            "analysis_scope"
        ]
        == "trans"
    )


# ===========================================================================
# Runner helpers
# ===========================================================================


def test_missing_required_columns() -> None:
    """Schema helper should identify missing required columns."""

    dataframe = pd.DataFrame(
        {
            "Unnamed: 0": [
                "EX1"
            ],
            "baseMean": [
                10.0
            ],
        }
    )

    missing = (
        MoradiDERunner
        ._missing_required_columns(
            dataframe
        )
    )

    assert (
        "log2FoldChange"
        in missing
    )

    assert (
        "pvalue"
        in missing
    )

    assert (
        "padj"
        in missing
    )


def test_duplicate_statistics() -> None:
    """Runner should distinguish duplicate rows from excess duplicates."""

    dataframe = pd.DataFrame(
        {
            "a": [
                1,
                1,
                2,
            ],
            "b": [
                "x",
                "x",
                "y",
            ],
        }
    )

    result = (
        MoradiDERunner
        ._duplicate_statistics(
            dataframe
        )
    )

    assert (
        result[
            "duplicate_rows"
        ]
        == 2
    )

    assert (
        result[
            "duplicate_excess_rows"
        ]
        == 1
    )


# ===========================================================================
# Temporary source-file helpers
# ===========================================================================


def write_de_file(
    directory: Path,
    filename: str,
    feature_ids: list[str],
) -> None:
    """Write a small valid Moradi-like DE CSV fixture."""

    dataframe = pd.DataFrame(
        {
            "Unnamed: 0":
                feature_ids,
            "baseMean":
                [
                    10.0
                    + index
                    for index
                    in range(
                        len(
                            feature_ids
                        )
                    )
                ],
            "log2FoldChange":
                [
                    1.0
                    for _
                    in feature_ids
                ],
            "lfcSE":
                [
                    0.2
                    for _
                    in feature_ids
                ],
            "stat":
                [
                    5.0
                    for _
                    in feature_ids
                ],
            "pvalue":
                [
                    0.001
                    for _
                    in feature_ids
                ],
            "padj":
                [
                    0.01
                    for _
                    in feature_ids
                ],
        }
    )

    dataframe.to_csv(
        directory
        / filename,
        index=False,
    )


def create_complete_test_dataset(
    directory: Path,
) -> None:
    """Create valid temporary S5-S10 source files."""

    feature_map = {
        "S5": [
            "EX1",
            "EX2",
        ],
        "S6": [
            "EX3",
            "EX4",
        ],
        "S7": [
            "INT1",
            "INT2",
        ],
        "S8": [
            "INT3",
            "INT4",
        ],
        "S9": [
            "ENST00000318325.6",
            "ENST00000412345.1",
        ],
        "S10": [
            "ENST00000511111.2",
            "ENST00000622222.1",
        ],
    }

    for (
        dataset_name,
        specification,
    ) in (
        MoradiDERunner
        .DATASETS
        .items()
    ):
        write_de_file(
            directory,
            specification[
                "filename"
            ],
            feature_map[
                dataset_name
            ],
        )


# ===========================================================================
# Runner integration tests
# ===========================================================================


def test_validate_complete_temporary_dataset(
    tmp_path: Path,
) -> None:
    """A valid S5-S10 fixture should produce PASS."""

    create_complete_test_dataset(
        tmp_path
    )

    runner = MoradiDERunner(
        tmp_path
    )

    report = runner.validate()

    assert (
        report[
            "status"
        ]
        == "PASS"
    )

    assert (
        report[
            "files_validated"
        ]
        == 6
    )

    assert (
        report[
            "total_rows"
        ]
        == 12
    )

    assert (
        report[
            "validation_summary"
        ][
            "severity_counts"
        ][
            "ERROR"
        ]
        == 0
    )

    for result in (
        report[
            "files"
        ]
        .values()
    ):
        assert (
            result[
                "status"
            ]
            == "PASS"
        )


def test_missing_source_file_causes_fail(
    tmp_path: Path,
) -> None:
    """A missing required S5-S10 source file should fail validation."""

    create_complete_test_dataset(
        tmp_path
    )

    missing_path = (
        tmp_path
        / MoradiDERunner.DATASETS[
            "S10"
        ][
            "filename"
        ]
    )

    missing_path.unlink()

    runner = MoradiDERunner(
        tmp_path
    )

    report = runner.validate()

    assert (
        report[
            "status"
        ]
        == "FAIL"
    )

    assert (
        report[
            "files"
        ][
            "S10"
        ][
            "status"
        ]
        == "FAIL"
    )


def test_invalid_feature_causes_fail(
    tmp_path: Path,
) -> None:
    """A malformed biological feature ID should fail its source file."""

    create_complete_test_dataset(
        tmp_path
    )

    specification = (
        MoradiDERunner
        .DATASETS[
            "S5"
        ]
    )

    file_path = (
        tmp_path
        / specification[
            "filename"
        ]
    )

    dataframe = pd.read_csv(
        file_path
    )

    dataframe.loc[
        0,
        "Unnamed: 0",
    ] = "INT999"

    dataframe.to_csv(
        file_path,
        index=False,
    )

    runner = MoradiDERunner(
        tmp_path
    )

    report = runner.validate()

    assert (
        report[
            "status"
        ]
        == "FAIL"
    )

    assert (
        report[
            "files"
        ][
            "S5"
        ][
            "validation_error_counts"
        ][
            "features"
        ]
        == 1
    )


def test_duplicate_rows_produce_warning(
    tmp_path: Path,
) -> None:
    """Exact duplicate DE records should produce WARNING, not ERROR."""

    create_complete_test_dataset(
        tmp_path
    )

    specification = (
        MoradiDERunner
        .DATASETS[
            "S5"
        ]
    )

    file_path = (
        tmp_path
        / specification[
            "filename"
        ]
    )

    dataframe = pd.read_csv(
        file_path
    )

    dataframe = pd.concat(
        [
            dataframe,
            dataframe.iloc[
                [0]
            ],
        ],
        ignore_index=True,
    )

    dataframe.to_csv(
        file_path,
        index=False,
    )

    runner = MoradiDERunner(
        tmp_path
    )

    report = runner.validate()

    s5 = (
        report[
            "files"
        ][
            "S5"
        ]
    )

    assert (
        s5[
            "status"
        ]
        == "WARNING"
    )

    assert (
        s5[
            "duplicates"
        ][
            "duplicate_excess_rows"
        ]
        == 1
    )

    assert (
        report[
            "status"
        ]
        == "WARNING"
    )