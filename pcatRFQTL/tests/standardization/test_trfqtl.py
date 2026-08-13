"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/standardization/test_trfqtl.py

Description:
    Unit tests for Cancer-tRFQTL source-specific standardization.

    Tests cover:

        - rsID standardization
        - hg19 coordinate parsing
        - allele parsing
        - cancer-type standardization
        - tRF identifier preservation
        - numeric and probability fields
        - S2 standardization
        - S10 standardization
        - known non-rsID S10 source behavior
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

from pcatrfqtl.standardization.trfqtl import (
    TRFQTLStandardizer,
)


# ===========================================================================
# rsID
# ===========================================================================


def test_standardize_rsid() -> None:
    """Canonical dbSNP identifiers should normalize consistently."""

    result = (
        TRFQTLStandardizer
        .standardize_rsid(
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
            "rsid_usable"
        ]
        is True
    )


def test_non_rsid_preserved_as_unusable() -> None:
    """Known non-rsID values must not be converted to artificial rsIDs."""

    result = (
        TRFQTLStandardizer
        .standardize_rsid(
            "6"
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
            "rsid_usable"
        ]
        is False
    )


# ===========================================================================
# Coordinates
# ===========================================================================


def test_standardize_hg19_coordinate() -> None:
    """hg19 coordinate should split into chromosome and position."""

    result = (
        TRFQTLStandardizer
        .standardize_coordinate(
            "6:28958399"
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
            "position"
        ]
        == 28958399
    )

    assert (
        result[
            "genome_build"
        ]
        == "hg19"
    )

    assert (
        result[
            "coordinate_usable"
        ]
        is True
    )


def test_standardize_chr_prefixed_coordinate() -> None:
    """chr-prefixed coordinates should normalize identically."""

    result = (
        TRFQTLStandardizer
        .standardize_coordinate(
            "chr6:28958399"
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
            "position"
        ]
        == 28958399
    )


def test_invalid_coordinate() -> None:
    """Malformed coordinates should remain unusable."""

    result = (
        TRFQTLStandardizer
        .standardize_coordinate(
            "invalid"
        )
    )

    assert (
        result[
            "coordinate_usable"
        ]
        is False
    )

    assert (
        result[
            "chromosome"
        ]
        is None
    )


# ===========================================================================
# Alleles
# ===========================================================================


def test_standardize_allele_pair() -> None:
    """Reference/alternate allele pairs should be parsed."""

    result = (
        TRFQTLStandardizer
        .standardize_alleles(
            "A/T"
        )
    )

    assert (
        result[
            "reference_allele"
        ]
        == "A"
    )

    assert (
        result[
            "alternate_allele"
        ]
        == "T"
    )

    assert (
        result[
            "allele_pair_usable"
        ]
        is True
    )


def test_standardize_sequence_allele_pair() -> None:
    """Multi-base allele pairs should remain valid."""

    result = (
        TRFQTLStandardizer
        .standardize_alleles(
            "AT/GC"
        )
    )

    assert (
        result[
            "reference_allele"
        ]
        == "AT"
    )

    assert (
        result[
            "alternate_allele"
        ]
        == "GC"
    )


def test_single_allele_preserved() -> None:
    """Single alleles should be preserved without inventing a pair."""

    result = (
        TRFQTLStandardizer
        .standardize_alleles(
            "A"
        )
    )

    assert (
        result[
            "allele_value"
        ]
        == "A"
    )

    assert (
        result[
            "reference_allele"
        ]
        is None
    )

    assert (
        result[
            "allele_pair_usable"
        ]
        is False
    )


# ===========================================================================
# Cancer / tRF
# ===========================================================================


def test_standardize_cancer_type() -> None:
    """Cancer-type strings should have normalized whitespace."""

    result = (
        TRFQTLStandardizer
        .standardize_cancer_type(
            "  prostate   cancer "
        )
    )

    assert (
        result[
            "cancer_type"
        ]
        == "prostate cancer"
    )


def test_standardize_trf() -> None:
    """Source tRF identifiers should be preserved."""

    result = (
        TRFQTLStandardizer
        .standardize_trf(
            "tRF-Example-1"
        )
    )

    assert (
        result[
            "trf_id"
        ]
        == "tRF-Example-1"
    )

    assert (
        result[
            "trf_usable"
        ]
        is True
    )


# ===========================================================================
# S2
# ===========================================================================


def test_standardize_s2_record() -> None:
    """A complete S2 record should become fully standardized."""

    record = {
        "Cancer type":
            "PRAD",

        "SNP ID":
            "RS12345",

        "SNP position (hg19)":
            "6:28958399",

        "Alleles":
            "A/T",

        "tRF":
            "tRF-test",

        "statistic":
            4.2,

        "P-value":
            1e-6,

        "FDR":
            0.01,

        "GWAS tagSnp":
            "rs999",

        "LD (r2)":
            0.8,
    }

    result = (
        TRFQTLStandardizer
        .standardize_s2_record(
            record
        )
    )

    assert (
        result[
            "source_table"
        ]
        == "S2"
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
        == 28958399
    )

    assert (
        result[
            "gwas_tag_rsid"
        ]
        == "rs999"
    )

    assert (
        result[
            "ld_r2"
        ]
        == 0.8
    )

    assert (
        result[
            "variant_usable"
        ]
        is True
    )

    assert (
        result[
            "qtl_usable"
        ]
        is True
    )

    assert (
        result[
            "standardization_status"
        ]
        == "STANDARDIZED"
    )


def test_s2_without_rsid_but_coordinate_still_variant_usable() -> None:
    """Coordinate information may keep an otherwise non-rsID variant usable."""

    record = {
        "Cancer type":
            "PRAD",

        "SNP ID":
            "invalid",

        "SNP position (hg19)":
            "6:28958399",

        "Alleles":
            "A/T",

        "tRF":
            "tRF-test",

        "P-value":
            0.01,
    }

    result = (
        TRFQTLStandardizer
        .standardize_s2_record(
            record
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


# ===========================================================================
# S10
# ===========================================================================


def test_standardize_s10_record() -> None:
    """A valid S10 GWAS-effect record should standardize fully."""

    record = {
        "Cancer type":
            "PRAD",

        "SNP ID":
            "rs12345",

        "Position (hg19)":
            "6:28958399",

        "A1 (effect allele)":
            "a",

        "Effect":
            0.25,

        "GWAS P-value":
            5e-8,
    }

    result = (
        TRFQTLStandardizer
        .standardize_s10_record(
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
            "effect_allele"
        ]
        == "A"
    )

    assert (
        result[
            "effect"
        ]
        == 0.25
    )

    assert (
        result[
            "effect_usable"
        ]
        is True
    )

    assert (
        result[
            "standardization_status"
        ]
        == "STANDARDIZED"
    )


def test_known_s10_non_rsid_source_value() -> None:
    """
    Known S10 source value '6' must remain non-rsID without automatic
    correction.
    """

    record = {
        "Cancer type":
            "PRAD",

        "SNP ID":
            "6",

        "Position (hg19)":
            "6:28958399",

        "A1 (effect allele)":
            "A",

        "Effect":
            0.25,

        "GWAS P-value":
            5e-8,
    }

    result = (
        TRFQTLStandardizer
        .standardize_s10_record(
            record
        )
    )

    assert (
        result[
            "snp_id_raw"
        ]
        == "6"
    )

    assert (
        result[
            "canonical_rsid"
        ]
        is None
    )

    assert (
        result[
            "rsid_usable"
        ]
        is False
    )

    # The coordinate remains independently usable.
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


# ===========================================================================
# DataFrames
# ===========================================================================


def test_standardize_s2_dataframe() -> None:
    """S2 DataFrame standardization should preserve cardinality."""

    dataframe = pd.DataFrame(
        [
            {
                "Cancer type":
                    "PRAD",
                "SNP ID":
                    "rs1",
                "SNP position (hg19)":
                    "1:100",
                "Alleles":
                    "A/G",
                "tRF":
                    "tRF-1",
                "statistic":
                    1.0,
                "P-value":
                    0.01,
                "FDR":
                    0.05,
                "GWAS tagSnp":
                    "rs2",
                "LD (r2)":
                    0.8,
            },
            {
                "Cancer type":
                    "BRCA",
                "SNP ID":
                    "rs3",
                "SNP position (hg19)":
                    "2:200",
                "Alleles":
                    "C/T",
                "tRF":
                    "tRF-2",
                "statistic":
                    -1.0,
                "P-value":
                    0.02,
                "FDR":
                    0.06,
                "GWAS tagSnp":
                    "rs4",
                "LD (r2)":
                    0.7,
            },
        ]
    )

    result = (
        TRFQTLStandardizer
        .standardize_s2_dataframe(
            dataframe
        )
    )

    assert len(
        result
    ) == 2

    assert (
        "canonical_rsid"
        in result.columns
    )

    assert (
        "genome_build"
        in result.columns
    )

    assert (
        "qtl_usable"
        in result.columns
    )


def test_standardize_s10_dataframe() -> None:
    """S10 DataFrame standardization should preserve cardinality."""

    dataframe = pd.DataFrame(
        [
            {
                "Cancer type":
                    "PRAD",
                "SNP ID":
                    "rs1",
                "Position (hg19)":
                    "1:100",
                "A1 (effect allele)":
                    "A",
                "Effect":
                    0.2,
                "GWAS P-value":
                    1e-8,
            },
            {
                "Cancer type":
                    "PRAD",
                "SNP ID":
                    "6",
                "Position (hg19)":
                    "6:28958399",
                "A1 (effect allele)":
                    "G",
                "Effect":
                    -0.1,
                "GWAS P-value":
                    0.01,
            },
        ]
    )

    result = (
        TRFQTLStandardizer
        .standardize_s10_dataframe(
            dataframe
        )
    )

    assert len(
        result
    ) == 2

    assert (
        "canonical_rsid"
        in result.columns
    )

    assert (
        "effect_usable"
        in result.columns
    )