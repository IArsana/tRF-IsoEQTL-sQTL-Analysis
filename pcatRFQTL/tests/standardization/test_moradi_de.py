"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/standardization/test_moradi_de.py

Description:
    Unit tests for Moradi differential-expression standardization.

    Tests cover:

        - exon feature identifiers
        - intron feature identifiers
        - Ensembl transcript identifiers
        - known INT1e+05 source anomaly
        - DESeq2-like numeric fields
        - probability fields
        - source table metadata
        - cis/trans scope metadata
        - record usability
        - DataFrame cardinality preservation

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

from pcatrfqtl.standardization.moradi_de import (
    MoradiDEStandardizer,
)


# ===========================================================================
# Feature identifiers
# ===========================================================================


def test_standardize_exon_feature() -> None:
    """Canonical exon IDs should remain usable."""

    result = (
        MoradiDEStandardizer
        .standardize_feature(
            "EX123",
            "exon",
        )
    )

    assert (
        result[
            "feature_id"
        ]
        == "EX123"
    )

    assert (
        result[
            "feature_usable"
        ]
        is True
    )

    assert (
        result[
            "source_feature_anomaly"
        ]
        is False
    )


def test_standardize_intron_feature() -> None:
    """Canonical intron IDs should remain usable."""

    result = (
        MoradiDEStandardizer
        .standardize_feature(
            "INT123",
            "intron",
        )
    )

    assert (
        result[
            "feature_id"
        ]
        == "INT123"
    )

    assert (
        result[
            "feature_usable"
        ]
        is True
    )


def test_standardize_transcript_feature() -> None:
    """Ensembl transcript IDs should normalize to uppercase."""

    result = (
        MoradiDEStandardizer
        .standardize_feature(
            "enst00000318325",
            "transcript",
        )
    )

    assert (
        result[
            "feature_id"
        ]
        == "ENST00000318325"
    )

    assert (
        result[
            "feature_usable"
        ]
        is True
    )


def test_transcript_version_is_preserved() -> None:
    """Transcript version suffix should be retained."""

    result = (
        MoradiDEStandardizer
        .standardize_feature(
            "ENST00000318325.6",
            "transcript",
        )
    )

    assert (
        result[
            "feature_id"
        ]
        == "ENST00000318325.6"
    )


def test_scientific_notation_like_intron_is_source_anomaly() -> None:
    """
    INT1e+05 must not be automatically converted to INT100000.
    """

    result = (
        MoradiDEStandardizer
        .standardize_feature(
            "INT1e+05",
            "intron",
        )
    )

    assert (
        result[
            "feature_id"
        ]
        is None
    )

    assert (
        result[
            "feature_usable"
        ]
        is False
    )

    assert (
        result[
            "source_feature_anomaly"
        ]
        is True
    )

    assert (
        result[
            "feature_source"
        ]
        == "scientific_notation_like_source_identifier"
    )


def test_noncanonical_feature_is_not_source_anomaly() -> None:
    """Generic malformed IDs should remain distinguishable."""

    result = (
        MoradiDEStandardizer
        .standardize_feature(
            "BAD123",
            "intron",
        )
    )

    assert (
        result[
            "feature_usable"
        ]
        is False
    )

    assert (
        result[
            "source_feature_anomaly"
        ]
        is False
    )

    assert (
        result[
            "feature_source"
        ]
        == "non_canonical_feature_id"
    )


# ===========================================================================
# Numeric fields
# ===========================================================================


def test_standardize_numeric() -> None:
    """Finite numeric values should become floats."""

    result = (
        MoradiDEStandardizer
        .standardize_numeric(
            "2.5"
        )
    )

    assert result == 2.5


def test_negative_log2_fold_change_is_valid() -> None:
    """Negative fold changes are biologically valid."""

    result = (
        MoradiDEStandardizer
        .standardize_numeric(
            -3.2
        )
    )

    assert result == -3.2


def test_nonnegative_base_mean() -> None:
    """baseMean must be non-negative."""

    valid = (
        MoradiDEStandardizer
        .standardize_nonnegative_numeric(
            10.0
        )
    )

    invalid = (
        MoradiDEStandardizer
        .standardize_nonnegative_numeric(
            -1.0
        )
    )

    assert valid == 10.0
    assert invalid is None


def test_probability() -> None:
    """P-values and adjusted P-values should remain within 0-1."""

    result = (
        MoradiDEStandardizer
        .standardize_probability(
            0.001
        )
    )

    assert result == 0.001


def test_invalid_probability() -> None:
    """Out-of-range probability should not be standardized."""

    result = (
        MoradiDEStandardizer
        .standardize_probability(
            1.5
        )
    )

    assert result is None


# ===========================================================================
# Complete records
# ===========================================================================


def test_standardize_s5_exon_record() -> None:
    """A valid S5 exon record should standardize fully."""

    record = {
        "Unnamed: 0":
            "EX123",

        "baseMean":
            100.5,

        "log2FoldChange":
            -2.1,

        "lfcSE":
            0.3,

        "stat":
            -7.0,

        "pvalue":
            1e-8,

        "padj":
            1e-6,
    }

    result = (
        MoradiDEStandardizer
        .standardize_record(
            record,
            source_table="S5",
            feature_type="exon",
            analysis_scope="cis",
        )
    )

    assert (
        result[
            "source_table"
        ]
        == "S5"
    )

    assert (
        result[
            "analysis_scope"
        ]
        == "cis"
    )

    assert (
        result[
            "data_type"
        ]
        == "differential_expression"
    )

    assert (
        result[
            "feature_id"
        ]
        == "EX123"
    )

    assert (
        result[
            "base_mean"
        ]
        == 100.5
    )

    assert (
        result[
            "log2_fold_change"
        ]
        == -2.1
    )

    assert (
        result[
            "p_value"
        ]
        == 1e-8
    )

    assert (
        result[
            "adjusted_p_value"
        ]
        == 1e-6
    )

    assert (
        result[
            "de_statistics_usable"
        ]
        is True
    )

    assert (
        result[
            "de_record_usable"
        ]
        is True
    )

    assert (
        result[
            "standardization_status"
        ]
        == "STANDARDIZED"
    )


def test_standardize_s9_transcript_record() -> None:
    """A valid isoform-expression record should remain DE data."""

    record = {
        "Unnamed: 0":
            "ENST00000318325",

        "baseMean":
            250.0,

        "log2FoldChange":
            1.5,

        "lfcSE":
            0.2,

        "stat":
            7.5,

        "pvalue":
            1e-10,

        "padj":
            1e-8,
    }

    result = (
        MoradiDEStandardizer
        .standardize_record(
            record,
            source_table="S9",
            feature_type="transcript",
            analysis_scope="cis",
        )
    )

    assert (
        result[
            "feature_type"
        ]
        == "transcript"
    )

    assert (
        result[
            "feature_id"
        ]
        == "ENST00000318325"
    )

    assert (
        result[
            "data_type"
        ]
        == "differential_expression"
    )

    assert (
        result[
            "de_record_usable"
        ]
        is True
    )


# ===========================================================================
# Known S7 anomaly
# ===========================================================================


def test_s7_known_source_anomaly_is_partial() -> None:
    """
    The known INT1e+05 record must preserve raw data and remain excluded
    from canonical feature-ID analyses.
    """

    record = {
        "Unnamed: 0":
            "INT1e+05",

        "baseMean":
            25.0,

        "log2FoldChange":
            -1.2,

        "lfcSE":
            0.2,

        "stat":
            -6.0,

        "pvalue":
            1e-5,

        "padj":
            0.001,
    }

    result = (
        MoradiDEStandardizer
        .standardize_record(
            record,
            source_table="S7",
            feature_type="intron",
            analysis_scope="cis",
        )
    )

    assert (
        result[
            "feature_raw"
        ]
        == "INT1e+05"
    )

    assert (
        result[
            "feature_id"
        ]
        is None
    )

    assert (
        result[
            "source_feature_anomaly"
        ]
        is True
    )

    assert (
        result[
            "de_statistics_usable"
        ]
        is True
    )

    assert (
        result[
            "de_record_usable"
        ]
        is False
    )

    assert (
        result[
            "standardization_status"
        ]
        == "PARTIAL"
    )


# ===========================================================================
# Missing adjusted P-value
# ===========================================================================


def test_missing_adjusted_p_value_can_still_standardize_record() -> None:
    """
    padj may be unavailable while core DE statistics remain usable.
    """

    record = {
        "Unnamed: 0":
            "EX123",

        "baseMean":
            100.0,

        "log2FoldChange":
            1.0,

        "lfcSE":
            0.25,

        "stat":
            4.0,

        "pvalue":
            0.001,

        "padj":
            None,
    }

    result = (
        MoradiDEStandardizer
        .standardize_record(
            record,
            source_table="S5",
            feature_type="exon",
            analysis_scope="cis",
        )
    )

    assert (
        result[
            "adjusted_p_value"
        ]
        is None
    )

    assert (
        result[
            "adjusted_p_value_available"
        ]
        is False
    )

    assert (
        result[
            "de_record_usable"
        ]
        is True
    )


# ===========================================================================
# Partial record
# ===========================================================================


def test_invalid_statistics_create_partial_record() -> None:
    """Valid feature with invalid DE statistics should remain partial."""

    record = {
        "Unnamed: 0":
            "INT123",

        "baseMean":
            -1,

        "log2FoldChange":
            1.0,

        "lfcSE":
            0.2,

        "stat":
            5.0,

        "pvalue":
            0.001,

        "padj":
            0.01,
    }

    result = (
        MoradiDEStandardizer
        .standardize_record(
            record,
            source_table="S7",
            feature_type="intron",
            analysis_scope="cis",
        )
    )

    assert (
        result[
            "feature_usable"
        ]
        is True
    )

    assert (
        result[
            "de_statistics_usable"
        ]
        is False
    )

    assert (
        result[
            "de_record_usable"
        ]
        is False
    )

    assert (
        result[
            "standardization_status"
        ]
        == "PARTIAL"
    )


# ===========================================================================
# DataFrame
# ===========================================================================


def test_standardize_dataframe_preserves_cardinality() -> None:
    """DataFrame standardization must preserve source row count."""

    dataframe = pd.DataFrame(
        [
            {
                "Unnamed: 0":
                    "INT1",
                "baseMean":
                    10,
                "log2FoldChange":
                    1.0,
                "lfcSE":
                    0.1,
                "stat":
                    10.0,
                "pvalue":
                    0.001,
                "padj":
                    0.01,
            },
            {
                "Unnamed: 0":
                    "INT2",
                "baseMean":
                    20,
                "log2FoldChange":
                    -1.0,
                "lfcSE":
                    0.2,
                "stat":
                    -5.0,
                "pvalue":
                    0.002,
                "padj":
                    0.02,
            },
        ]
    )

    result = (
        MoradiDEStandardizer
        .standardize_dataframe(
            dataframe,
            source_table="S7",
            feature_type="intron",
            analysis_scope="cis",
        )
    )

    assert len(
        result
    ) == 2

    assert (
        "feature_id"
        in result.columns
    )

    assert (
        "log2_fold_change"
        in result.columns
    )

    assert (
        "de_record_usable"
        in result.columns
    )

    assert (
        "standardization_status"
        in result.columns
    )


def test_standardize_dataframe_preserves_source_table_metadata() -> None:
    """Every standardized row should retain dataset provenance."""

    dataframe = pd.DataFrame(
        [
            {
                "Unnamed: 0":
                    "ENST000001",
                "baseMean":
                    10,
                "log2FoldChange":
                    1,
                "lfcSE":
                    0.2,
                "stat":
                    5,
                "pvalue":
                    0.01,
                "padj":
                    0.05,
            }
        ]
    )

    result = (
        MoradiDEStandardizer
        .standardize_dataframe(
            dataframe,
            source_table="S10",
            feature_type="transcript",
            analysis_scope="trans",
        )
    )

    assert (
        result.iloc[
            0
        ][
            "source_table"
        ]
        == "S10"
    )

    assert (
        result.iloc[
            0
        ][
            "analysis_scope"
        ]
        == "trans"
    )