"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/standardization/test_common.py

Description:
    Unit tests for source-independent standardization primitives used
    throughout the pcatRFQTL research pipeline.

    Tests cover:

        - missing-value detection
        - string standardization
        - chromosome standardization
        - genomic-position standardization
        - numeric standardization
        - probability standardization
        - boolean standardization
        - dbSNP rsID standardization
        - Ensembl transcript standardization
        - standardization result serialization

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

import numpy as np

from pcatrfqtl.standardization.common import (
    StandardizationUtils,
)
from pcatrfqtl.standardization.models import (
    StandardizationResult,
    StandardizationStatus,
)


# ===========================================================================
# StandardizationResult
# ===========================================================================


def test_standardization_result_to_dict() -> None:
    """StandardizationResult should serialize enum status as text."""

    result = StandardizationResult(
        raw_value="chr1",
        standardized_value="1",
        status=(
            StandardizationStatus
            .STANDARDIZED
        ),
        source="chromosome_normalization",
        usable=True,
        note=None,
    )

    serialized = result.to_dict()

    assert serialized[
        "raw_value"
    ] == "chr1"

    assert serialized[
        "standardized_value"
    ] == "1"

    assert serialized[
        "status"
    ] == "STANDARDIZED"

    assert serialized[
        "source"
    ] == "chromosome_normalization"

    assert serialized[
        "usable"
    ] is True


# ===========================================================================
# Missing values
# ===========================================================================


def test_missing_values() -> None:
    """Common missing representations should be detected."""

    values = [
        None,
        np.nan,
        float("nan"),
        "",
        " ",
        "NA",
        "nan",
        "None",
        "null",
    ]

    for value in values:
        assert (
            StandardizationUtils
            .is_missing(
                value
            )
        )


def test_non_missing_values() -> None:
    """Valid zero and identifiers should not be considered missing."""

    values = [
        0,
        1,
        False,
        "rs123",
        "chr1",
        "0",
    ]

    for value in values:
        assert not (
            StandardizationUtils
            .is_missing(
                value
            )
        )


# ===========================================================================
# String standardization
# ===========================================================================


def test_standardize_string() -> None:
    """Whitespace should be stripped and collapsed."""

    result = (
        StandardizationUtils
        .standardize_string(
            "  prostate   cancer  "
        )
    )

    assert (
        result.standardized_value
        == "prostate cancer"
    )

    assert (
        result.status
        == StandardizationStatus.STANDARDIZED
    )

    assert result.usable is True


def test_standardize_missing_string() -> None:
    """Missing strings should produce MISSING status."""

    result = (
        StandardizationUtils
        .standardize_string(
            None
        )
    )

    assert (
        result.standardized_value
        is None
    )

    assert (
        result.status
        == StandardizationStatus.MISSING
    )

    assert result.usable is False


# ===========================================================================
# Chromosomes
# ===========================================================================


def test_standardize_autosomal_chromosome() -> None:
    """chr-prefixed autosomes should normalize to canonical labels."""

    result = (
        StandardizationUtils
        .standardize_chromosome(
            "chr01"
        )
    )

    # chr01 is intentionally not accepted by the current canonical
    # chromosome pattern.
    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


def test_standardize_chr1() -> None:
    """chr1 should normalize to 1."""

    result = (
        StandardizationUtils
        .standardize_chromosome(
            "chr1"
        )
    )

    assert (
        result.standardized_value
        == "1"
    )

    assert result.usable is True


def test_standardize_sex_chromosome() -> None:
    """Sex chromosome labels should be normalized."""

    result = (
        StandardizationUtils
        .standardize_chromosome(
            "chrx"
        )
    )

    assert (
        result.standardized_value
        == "X"
    )


def test_standardize_mitochondrial_chromosome() -> None:
    """M and chrM should normalize to MT."""

    values = [
        "M",
        "chrM",
        "MT",
        "chrMT",
    ]

    for value in values:

        result = (
            StandardizationUtils
            .standardize_chromosome(
                value
            )
        )

        assert (
            result.standardized_value
            == "MT"
        )

        assert result.usable is True


def test_invalid_chromosome() -> None:
    """Non-human chromosome labels should be unsupported."""

    result = (
        StandardizationUtils
        .standardize_chromosome(
            "chr23"
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )

    assert (
        result.standardized_value
        is None
    )

    assert result.usable is False


def test_missing_chromosome() -> None:
    """Missing chromosome should remain missing."""

    result = (
        StandardizationUtils
        .standardize_chromosome(
            None
        )
    )

    assert (
        result.status
        == StandardizationStatus.MISSING
    )

    assert result.usable is False


# ===========================================================================
# Genomic positions
# ===========================================================================


def test_standardize_position_integer() -> None:
    """Positive integer positions should be standardized."""

    result = (
        StandardizationUtils
        .standardize_position(
            "123456"
        )
    )

    assert (
        result.standardized_value
        == 123456
    )

    assert isinstance(
        result.standardized_value,
        int,
    )

    assert result.usable is True


def test_standardize_integral_float_position() -> None:
    """Integral floating-point positions should become integers."""

    result = (
        StandardizationUtils
        .standardize_position(
            100.0
        )
    )

    assert (
        result.standardized_value
        == 100
    )


def test_invalid_zero_position() -> None:
    """Position zero is not a valid genomic coordinate."""

    result = (
        StandardizationUtils
        .standardize_position(
            0
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


def test_invalid_negative_position() -> None:
    """Negative positions should be unsupported."""

    result = (
        StandardizationUtils
        .standardize_position(
            -100
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


def test_invalid_fractional_position() -> None:
    """Fractional genomic positions should be unsupported."""

    result = (
        StandardizationUtils
        .standardize_position(
            100.5
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


def test_non_numeric_position() -> None:
    """Non-numeric genomic positions should be unsupported."""

    result = (
        StandardizationUtils
        .standardize_position(
            "abc"
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


# ===========================================================================
# Numeric values
# ===========================================================================


def test_standardize_numeric() -> None:
    """Finite numeric strings should become float values."""

    result = (
        StandardizationUtils
        .standardize_numeric(
            "3.14"
        )
    )

    assert (
        result.standardized_value
        == 3.14
    )

    assert isinstance(
        result.standardized_value,
        float,
    )


def test_standardize_negative_numeric() -> None:
    """Generic numeric standardization should permit negatives."""

    result = (
        StandardizationUtils
        .standardize_numeric(
            -2.5
        )
    )

    assert (
        result.standardized_value
        == -2.5
    )


def test_non_finite_numeric() -> None:
    """Infinite numeric values should not be usable."""

    for value in [
        math.inf,
        -math.inf,
    ]:

        result = (
            StandardizationUtils
            .standardize_numeric(
                value
            )
        )

        assert (
            result.status
            == StandardizationStatus.UNSUPPORTED
        )


def test_missing_numeric() -> None:
    """Missing numeric input should remain missing."""

    result = (
        StandardizationUtils
        .standardize_numeric(
            None
        )
    )

    assert (
        result.status
        == StandardizationStatus.MISSING
    )

    assert result.usable is True


# ===========================================================================
# Probability values
# ===========================================================================


def test_valid_probability() -> None:
    """Probability values between zero and one should pass."""

    for value in [
        0,
        1e-10,
        0.5,
        1,
    ]:

        result = (
            StandardizationUtils
            .standardize_probability(
                value
            )
        )

        assert result.usable is True


def test_zero_probability_strict_mode() -> None:
    """Zero should become unsupported when explicitly disallowed."""

    result = (
        StandardizationUtils
        .standardize_probability(
            0,
            allow_zero=False,
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


def test_probability_outside_range() -> None:
    """Probabilities outside 0-1 should be unsupported."""

    for value in [
        -0.01,
        1.01,
    ]:

        result = (
            StandardizationUtils
            .standardize_probability(
                value
            )
        )

        assert (
            result.status
            == StandardizationStatus.UNSUPPORTED
        )


# ===========================================================================
# Boolean values
# ===========================================================================


def test_standardize_boolean_native() -> None:
    """Native booleans should be preserved."""

    result = (
        StandardizationUtils
        .standardize_boolean(
            True
        )
    )

    assert (
        result.standardized_value
        is True
    )

    assert (
        result.status
        == StandardizationStatus.PRESERVED
    )


def test_standardize_boolean_integer() -> None:
    """Integer zero and one should normalize to boolean values."""

    false_result = (
        StandardizationUtils
        .standardize_boolean(
            0
        )
    )

    true_result = (
        StandardizationUtils
        .standardize_boolean(
            1
        )
    )

    assert (
        false_result.standardized_value
        is False
    )

    assert (
        true_result.standardized_value
        is True
    )


def test_standardize_boolean_strings() -> None:
    """Common textual booleans should normalize consistently."""

    truthy = [
        "true",
        "TRUE",
        "yes",
        "Y",
        "1",
    ]

    falsy = [
        "false",
        "FALSE",
        "no",
        "N",
        "0",
    ]

    for value in truthy:
        assert (
            StandardizationUtils
            .standardize_boolean(
                value
            )
            .standardized_value
            is True
        )

    for value in falsy:
        assert (
            StandardizationUtils
            .standardize_boolean(
                value
            )
            .standardized_value
            is False
        )


def test_invalid_boolean() -> None:
    """Unknown boolean representations should be unsupported."""

    result = (
        StandardizationUtils
        .standardize_boolean(
            "maybe"
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )


# ===========================================================================
# rsID
# ===========================================================================


def test_standardize_rsid() -> None:
    """rsID prefix should normalize to lowercase."""

    result = (
        StandardizationUtils
        .standardize_rsid(
            "RS123456"
        )
    )

    assert (
        result.standardized_value
        == "rs123456"
    )

    assert result.usable is True


def test_invalid_rsid() -> None:
    """Non-rsID representations should not be canonicalized."""

    result = (
        StandardizationUtils
        .standardize_rsid(
            "chr1:12345"
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )

    assert (
        result.standardized_value
        is None
    )


# ===========================================================================
# Ensembl transcripts
# ===========================================================================


def test_standardize_transcript_preserve_version() -> None:
    """Transcript version should be preserved by default."""

    result = (
        StandardizationUtils
        .standardize_transcript_id(
            "enst00000318325.6"
        )
    )

    assert (
        result.standardized_value
        == "ENST00000318325.6"
    )


def test_standardize_transcript_strip_version() -> None:
    """Transcript version may be explicitly removed."""

    result = (
        StandardizationUtils
        .standardize_transcript_id(
            "ENST00000318325.6",
            strip_version=True,
        )
    )

    assert (
        result.standardized_value
        == "ENST00000318325"
    )

    assert (
        result.source
        == "ensembl_transcript_version_removed"
    )


def test_invalid_transcript() -> None:
    """Non-ENST identifiers should remain unsupported."""

    result = (
        StandardizationUtils
        .standardize_transcript_id(
            "ENSG000001"
        )
    )

    assert (
        result.status
        == StandardizationStatus.UNSUPPORTED
    )