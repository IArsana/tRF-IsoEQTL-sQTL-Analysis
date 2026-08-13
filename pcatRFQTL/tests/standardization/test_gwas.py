"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/standardization/test_gwas.py

Description:
    Unit tests for GWAS Catalog source-specific standardization.

    Tests cover:

        - GWAS variant classification
        - canonical rsID derivation
        - multi-rsID handling
        - SNP interaction handling
        - chromosome standardization
        - genomic-position standardization
        - P-value handling
        - P-value reconstruction from PVALUE_MLOG
        - risk-allele extraction
        - row-level usability flags
        - DataFrame standardization

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
import pytest

from pcatrfqtl.standardization.gwas import (
    GWASStandardizer,
)


# ===========================================================================
# Variant classification
# ===========================================================================


def test_classify_single_rsid() -> None:
    """Single rsID should retain GWAS validation classification."""

    assert (
        GWASStandardizer
        .classify_variant(
            "rs12345"
        )
        == "single_rsid"
    )


def test_classify_interaction() -> None:
    """SNP interaction descriptors should remain identifiable."""

    assert (
        GWASStandardizer
        .classify_variant(
            "rs123 x rs456"
        )
        == "snp_interaction"
    )


def test_classify_coordinate_variant() -> None:
    """Coordinate variants should be recognized."""

    assert (
        GWASStandardizer
        .classify_variant(
            "chr12:9098995"
        )
        == "coordinate_variant"
    )


# ===========================================================================
# Canonical rsID
# ===========================================================================


def test_single_rsid_becomes_canonical() -> None:
    """Single rsID should produce canonical_rsid."""

    result = (
        GWASStandardizer
        .standardize_canonical_rsid(
            "RS12345"
        )
    )

    assert (
        result[
            "canonical_rsid"
        ]
        == "rs12345"
    )

    assert (
        result[
            "canonical_rsid_usable"
        ]
        is True
    )


def test_multi_rsid_not_collapsed() -> None:
    """Multiple variants must not be collapsed to one rsID."""

    result = (
        GWASStandardizer
        .standardize_canonical_rsid(
            "rs123; rs456"
        )
    )

    assert (
        result[
            "canonical_rsid"
        ]
        is None
    )

    assert (
        result[
            "canonical_rsid_usable"
        ]
        is False
    )


# ===========================================================================
# Multi-rsID parsing
# ===========================================================================


def test_split_multi_rsid() -> None:
    """Multi-rsID fields should produce standardized rsID lists."""

    result = (
        GWASStandardizer
        .split_multi_rsid(
            "RS123; rs456,rs789"
        )
    )

    assert result == [
        "rs123",
        "rs456",
        "rs789",
    ]


# ===========================================================================
# Chromosomes
# ===========================================================================


def test_standardize_single_chromosome() -> None:
    """Single chromosome should become canonical."""

    result = (
        GWASStandardizer
        .standardize_chromosome(
            "chr6"
        )
    )

    assert (
        result[
            "chromosome"
        ]
        == "6"
    )

    assert (
        result[
            "interaction_chromosomes"
        ]
        is None
    )

    assert (
        result[
            "chromosome_usable"
        ]
        is True
    )


def test_standardize_interaction_chromosomes() -> None:
    """Interaction chromosome fields should remain multi-coordinate."""

    result = (
        GWASStandardizer
        .standardize_chromosome(
            "6 x 6"
        )
    )

    assert (
        result[
            "chromosome"
        ]
        is None
    )

    assert (
        result[
            "interaction_chromosomes"
        ]
        == [
            "6",
            "6",
        ]
    )

    assert (
        result[
            "chromosome_usable"
        ]
        is False
    )


# ===========================================================================
# Positions
# ===========================================================================


def test_standardize_single_position() -> None:
    """Single position should become integer."""

    result = (
        GWASStandardizer
        .standardize_position(
            "123456"
        )
    )

    assert (
        result[
            "position"
        ]
        == 123456
    )

    assert (
        result[
            "position_usable"
        ]
        is True
    )


def test_standardize_interaction_positions() -> None:
    """Interaction coordinates should be retained as a list."""

    result = (
        GWASStandardizer
        .standardize_position(
            "31268092 x 31463174"
        )
    )

    assert (
        result[
            "position"
        ]
        is None
    )

    assert (
        result[
            "interaction_positions"
        ]
        == [
            31268092,
            31463174,
        ]
    )

    assert (
        result[
            "position_usable"
        ]
        is False
    )


# ===========================================================================
# P-values
# ===========================================================================


def test_reported_positive_p_value() -> None:
    """Positive source P-value should be used directly."""

    result = (
        GWASStandardizer
        .standardize_p_value(
            5e-8,
            7.30103,
        )
    )

    assert (
        result[
            "p_value_standardized"
        ]
        == pytest.approx(
            5e-8
        )
    )

    assert (
        result[
            "p_value_source"
        ]
        == "reported_p"
    )

    assert (
        result[
            "p_value_reconstructed"
        ]
        is False
    )


def test_zero_p_value_reconstructed_from_mlog() -> None:
    """Zero P-values should use mlog when representable."""

    result = (
        GWASStandardizer
        .standardize_p_value(
            0,
            10,
        )
    )

    assert (
        result[
            "p_value_standardized"
        ]
        == pytest.approx(
            1e-10
        )
    )

    assert (
        result[
            "p_value_source"
        ]
        == "reconstructed_from_mlog"
    )

    assert (
        result[
            "p_value_reconstructed"
        ]
        is True
    )

    assert (
        result[
            "p_value_usable"
        ]
        is True
    )


def test_zero_p_value_without_mlog() -> None:
    """Zero without supporting mlog should remain provenance-limited."""

    result = (
        GWASStandardizer
        .standardize_p_value(
            0,
            None,
        )
    )

    assert (
        result[
            "p_value_standardized"
        ]
        == 0.0
    )

    assert (
        result[
            "p_value_source"
        ]
        == "reported_zero"
    )

    assert (
        result[
            "p_value_usable"
        ]
        is False
    )


def test_missing_p_value() -> None:
    """Missing P-value should remain missing."""

    result = (
        GWASStandardizer
        .standardize_p_value(
            None,
            None,
        )
    )

    assert (
        result[
            "p_value_standardized"
        ]
        is None
    )

    assert (
        result[
            "p_value_source"
        ]
        == "missing"
    )

    assert (
        result[
            "p_value_usable"
        ]
        is False
    )


# ===========================================================================
# Risk allele
# ===========================================================================


def test_standardize_risk_allele() -> None:
    """Risk allele should be extracted from GWAS Catalog source field."""

    result = (
        GWASStandardizer
        .standardize_risk_allele(
            "rs12345-A"
        )
    )

    assert (
        result[
            "risk_allele"
        ]
        == "A"
    )

    assert (
        result[
            "risk_allele_usable"
        ]
        is True
    )


def test_unknown_risk_allele() -> None:
    """Unknown risk allele should remain unavailable."""

    result = (
        GWASStandardizer
        .standardize_risk_allele(
            "rs12345-?"
        )
    )

    assert (
        result[
            "risk_allele"
        ]
        is None
    )

    assert (
        result[
            "risk_allele_usable"
        ]
        is False
    )


# ===========================================================================
# Full record
# ===========================================================================


def test_standardize_single_variant_record() -> None:
    """A conventional single-variant record should standardize cleanly."""

    record = {
        "STUDY ACCESSION":
            "GCST000001",
        "PUBMEDID":
            "12345678",
        "DISEASE/TRAIT":
            "Prostate cancer",
        "MAPPED_TRAIT":
            "prostate carcinoma",
        "MAPPED_TRAIT_URI":
            "example",
        "SNPS":
            "RS12345",
        "CHR_ID":
            "chr6",
        "CHR_POS":
            "123456",
        "STRONGEST SNP-RISK ALLELE":
            "rs12345-A",
        "P-VALUE":
            5e-8,
        "PVALUE_MLOG":
            7.30103,
        "MAPPED_GENE":
            "GENE1",
        "REPORTED GENE(S)":
            "GENE1",
    }

    result = (
        GWASStandardizer
        .standardize_record(
            record
        )
    )

    assert (
        result[
            "canonical_rsid"
        ]
        == "rs12345"
    )

    assert (
        result[
            "chromosome"
        ]
        == "6"
    )

    assert (
        result[
            "position"
        ]
        == 123456
    )

    assert (
        result[
            "coordinate_usable"
        ]
        is True
    )

    assert (
        result[
            "variant_usable"
        ]
        is True
    )

    assert (
        result[
            "standardization_status"
        ]
        == "STANDARDIZED"
    )


def test_interaction_record_not_coordinate_usable() -> None:
    """SNP interactions should remain outside single-coordinate analysis."""

    record = {
        "SNPS":
            "rs123 x rs456",
        "CHR_ID":
            "6 x 6",
        "CHR_POS":
            "100 x 200",
        "P-VALUE":
            1e-8,
        "PVALUE_MLOG":
            8,
    }

    result = (
        GWASStandardizer
        .standardize_record(
            record
        )
    )

    assert (
        result[
            "variant_class"
        ]
        == "snp_interaction"
    )

    assert (
        result[
            "coordinate_usable"
        ]
        is False
    )

    assert (
        result[
            "variant_usable"
        ]
        is False
    )


# ===========================================================================
# DataFrame
# ===========================================================================


def test_standardize_dataframe() -> None:
    """DataFrame standardization should preserve row cardinality."""

    dataframe = pd.DataFrame(
        [
            {
                "SNPS":
                    "rs123",
                "CHR_ID":
                    "1",
                "CHR_POS":
                    100,
                "P-VALUE":
                    1e-8,
                "PVALUE_MLOG":
                    8,
            },
            {
                "SNPS":
                    "chr2:200",
                "CHR_ID":
                    "2",
                "CHR_POS":
                    200,
                "P-VALUE":
                    0.05,
                "PVALUE_MLOG":
                    1.30103,
            },
        ]
    )

    result = (
        GWASStandardizer
        .standardize_dataframe(
            dataframe
        )
    )

    assert len(
        result
    ) == 2

    assert (
        "source_variant"
        in result.columns
    )

    assert (
        "variant_class"
        in result.columns
    )

    assert (
        "coordinate_usable"
        in result.columns
    )