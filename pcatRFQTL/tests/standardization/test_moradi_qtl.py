"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/standardization/test_moradi_qtl.py

Description:
    Unit tests for Moradi QTL source-specific standardization.

    Tests cover:

        - canonical rsID standardization
        - composite rsID standardization
        - genomic-coordinate parsing
        - allele standardization
        - exon identifiers
        - intron identifiers
        - transcript identifiers
        - known INT1e+05 source anomaly
        - S1 standardization
        - S2/S4 semantic column inversion
        - S3 source column T handling
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

from pcatrfqtl.standardization.moradi_qtl import (
    MoradiQTLStandardizer,
)


# ===========================================================================
# rsID
# ===========================================================================


def test_standardize_rsid() -> None:
    """Single rsIDs should normalize to canonical lowercase prefix."""

    result = (
        MoradiQTLStandardizer
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


def test_standardize_composite_rsid() -> None:
    """Colon-separated rsIDs should remain explicitly composite."""

    result = (
        MoradiQTLStandardizer
        .standardize_composite_rsid(
            "RS123:rs456"
        )
    )

    assert result[
        "rsids"
    ] == [
        "rs123",
        "rs456",
    ]

    assert (
        result[
            "primary_rsid"
        ]
        is None
    )

    assert (
        result[
            "composite_rsid"
        ]
        is True
    )

    assert (
        result[
            "rsid_usable"
        ]
        is True
    )


def test_single_composite_field_rsid() -> None:
    """A single rsID in a composite-capable field remains single."""

    result = (
        MoradiQTLStandardizer
        .standardize_composite_rsid(
            "rs123"
        )
    )

    assert result[
        "rsids"
    ] == [
        "rs123"
    ]

    assert (
        result[
            "primary_rsid"
        ]
        == "rs123"
    )

    assert (
        result[
            "composite_rsid"
        ]
        is False
    )


# ===========================================================================
# Coordinate
# ===========================================================================


def test_standardize_coordinate() -> None:
    """Chromosome-position strings should be decomposed."""

    result = (
        MoradiQTLStandardizer
        .standardize_coordinate(
            "11:69561656"
        )
    )

    assert (
        result[
            "chromosome"
        ]
        == "11"
    )

    assert (
        result[
            "position"
        ]
        == 69561656
    )

    assert (
        result[
            "coordinate_usable"
        ]
        is True
    )

    assert (
        result[
            "genome_build"
        ]
        == "to_be_verified"
    )


def test_chr_prefixed_coordinate() -> None:
    """chr-prefixed coordinates should also standardize."""

    result = (
        MoradiQTLStandardizer
        .standardize_coordinate(
            "chr11:69561656"
        )
    )

    assert (
        result[
            "chromosome"
        ]
        == "11"
    )


# ===========================================================================
# Features
# ===========================================================================


def test_standardize_exon_feature() -> None:
    """Canonical exon identifiers should remain usable."""

    result = (
        MoradiQTLStandardizer
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


def test_standardize_intron_feature() -> None:
    """Canonical intron identifiers should remain usable."""

    result = (
        MoradiQTLStandardizer
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


def test_scientific_notation_like_intron_not_corrected() -> None:
    """
    Known source anomaly INT1e+05 must not be silently converted to
    INT100000.
    """

    result = (
        MoradiQTLStandardizer
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
            "feature_source"
        ]
        == "non_canonical_feature_id"
    )


def test_standardize_transcript_feature() -> None:
    """Ensembl transcript identifiers should remain canonical."""

    result = (
        MoradiQTLStandardizer
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


# ===========================================================================
# Alleles
# ===========================================================================


def test_standardize_biallelic_field() -> None:
    """Comma-separated alleles should become a list."""

    result = (
        MoradiQTLStandardizer
        .standardize_alleles(
            "A,G"
        )
    )

    assert result[
        "alleles"
    ] == [
        "A",
        "G",
    ]

    assert (
        result[
            "allele_count"
        ]
        == 2
    )

    assert (
        result[
            "alleles_usable"
        ]
        is True
    )


def test_standardize_multiallelic_field() -> None:
    """Multi-allelic source fields should remain supported."""

    result = (
        MoradiQTLStandardizer
        .standardize_alleles(
            "A,C,G"
        )
    )

    assert result[
        "alleles"
    ] == [
        "A",
        "C",
        "G",
    ]

    assert (
        result[
            "allele_source"
        ]
        == "multi_allelic"
    )


# ===========================================================================
# S1
# ===========================================================================


def test_standardize_s1_intron_record() -> None:
    """A complete S1 intron-QTL record should standardize."""

    record = {
        "SNP":
            "rs12786544",

        "R":
            "A",

        "A":
            "G",

        "SNP_pos":
            "11:69561656",

        "splicing_event":
            "INT123",

        "Stat":
            3.2,

        "p_value":
            1e-5,

        "FDR":
            0.01,

        "beta":
            0.25,

        "MAF":
            0.2,

        "AvgCall":
            0.99,

        "Rsq":
            0.95,
    }

    result = (
        MoradiQTLStandardizer
        .standardize_s1_record(
            record,
            feature_type="intron",
        )
    )

    assert (
        result[
            "canonical_rsid"
        ]
        == "rs12786544"
    )

    assert (
        result[
            "chromosome"
        ]
        == "11"
    )

    assert (
        result[
            "position"
        ]
        == 69561656
    )

    assert (
        result[
            "feature_id"
        ]
        == "INT123"
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


def test_s1_known_source_feature_anomaly() -> None:
    """INT1e+05 should result in partial standardization."""

    record = {
        "SNP":
            "rs12786544",

        "SNP_pos":
            "11:69561656",

        "splicing_event":
            "INT1e+05",

        "p_value":
            1e-5,
    }

    result = (
        MoradiQTLStandardizer
        .standardize_s1_record(
            record,
            feature_type="intron",
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
            "qtl_usable"
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
# S2 / S4 semantic inversion
# ===========================================================================


def test_standardize_s2_ld_record() -> None:
    """S2 semantic inversion should be handled explicitly."""

    record = {
        "sQTL_SNP":
            "11:69561656",

        "sQTL_SNP-pos":
            "rs123",

        "splicing_event":
            "INT123",

        "tag_SNP":
            "11:69570000",

        "tag_SNP_pos":
            "rs456:rs789",

        "LD":
            0.8,

        "gwas_cancer":
            "Prostate cancer",
    }

    result = (
        MoradiQTLStandardizer
        .standardize_ld_record(
            record,
            source_table="S2",
            feature_type="intron",
            qtl_scope="cis",
            feature_column="splicing_event",
        )
    )

    assert (
        result[
            "qtl_chromosome"
        ]
        == "11"
    )

    assert (
        result[
            "qtl_position"
        ]
        == 69561656
    )

    assert (
        result[
            "qtl_rsid"
        ]
        == "rs123"
    )

    assert result[
        "tag_rsids"
    ] == [
        "rs456",
        "rs789",
    ]

    assert (
        result[
            "tag_is_composite_rsid"
        ]
        is True
    )

    assert (
        result[
            "qtl_usable"
        ]
        is True
    )


def test_standardize_s4_trans_iso_source_quirk() -> None:
    """
    S4 trans-iso uses source splicing_event despite its sheet-level
    iso-eQTL label, so the source feature representation is respected.
    """

    record = {
        "sQTL_SNP":
            "1:100",

        "sQTL_SNP-pos":
            "rs1",

        "splicing_event":
            "INT55",

        "tag_SNP":
            "1:200",

        "tag_SNP_pos":
            "rs2",

        "LD":
            0.7,

        "gwas_cancer":
            "PRAD",
    }

    result = (
        MoradiQTLStandardizer
        .standardize_ld_record(
            record,
            source_table="S4",
            feature_type="intron",
            qtl_scope="trans",
            feature_column="splicing_event",
        )
    )

    assert (
        result[
            "source_table"
        ]
        == "S4"
    )

    assert (
        result[
            "feature_id"
        ]
        == "INT55"
    )

    assert (
        result[
            "qtl_scope"
        ]
        == "trans"
    )


# ===========================================================================
# S3
# ===========================================================================


def test_standardize_s3_standard_variant_column() -> None:
    """Ordinary S3 datasets should support SNP as variant column."""

    record = {
        "SNP":
            "rs100",

        "SNP_pos":
            "1:100",

        "splicing_event":
            "INT123",

        "Stat":
            2.0,

        "p_value":
            0.001,

        "FDR":
            0.01,

        "beta":
            0.2,

        "MAF":
            0.1,

        "AvgCall":
            0.98,

        "Rsq":
            0.9,
    }

    result = (
        MoradiQTLStandardizer
        .standardize_s3_record(
            record,
            feature_type="intron",
            feature_column="splicing_event",
        )
    )

    assert (
        result[
            "canonical_rsid"
        ]
        == "rs100"
    )

    assert (
        result[
            "source_variant_column"
        ]
        == "SNP"
    )

    assert (
        result[
            "qtl_usable"
        ]
        is True
    )


def test_standardize_s3_trans_exon_t_column() -> None:
    """
    S3 trans-exon source column T must be explicitly supported without
    renaming the raw source schema.
    """

    record = {
        "T":
            "rs200",

        "SNP_pos":
            "2:200",

        "splicing_event":
            "EX55",

        "p_value":
            0.001,
    }

    result = (
        MoradiQTLStandardizer
        .standardize_s3_record(
            record,
            feature_type="exon",
            feature_column="splicing_event",
            variant_column="T",
        )
    )

    assert (
        result[
            "snp_id_raw"
        ]
        == "rs200"
    )

    assert (
        result[
            "canonical_rsid"
        ]
        == "rs200"
    )

    assert (
        result[
            "source_variant_column"
        ]
        == "T"
    )

    assert (
        result[
            "feature_id"
        ]
        == "EX55"
    )


# ===========================================================================
# DataFrame
# ===========================================================================


def test_standardize_s1_dataframe_preserves_cardinality() -> None:
    """S1 DataFrame standardization must preserve row count."""

    dataframe = pd.DataFrame(
        [
            {
                "SNP":
                    "rs1",
                "SNP_pos":
                    "1:100",
                "splicing_event":
                    "INT1",
                "p_value":
                    0.01,
            },
            {
                "SNP":
                    "rs2",
                "SNP_pos":
                    "1:200",
                "splicing_event":
                    "INT2",
                "p_value":
                    0.02,
            },
        ]
    )

    result = (
        MoradiQTLStandardizer
        .standardize_s1_dataframe(
            dataframe,
            feature_type="intron",
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
        "feature_id"
        in result.columns
    )

    assert (
        "qtl_usable"
        in result.columns
    )


def test_standardize_ld_dataframe_preserves_cardinality() -> None:
    """S2/S4 standardization must preserve source row count."""

    dataframe = pd.DataFrame(
        [
            {
                "sQTL_SNP":
                    "1:100",
                "sQTL_SNP-pos":
                    "rs1",
                "splicing_event":
                    "INT1",
                "tag_SNP":
                    "1:200",
                "tag_SNP_pos":
                    "rs2",
                "LD":
                    0.8,
                "gwas_cancer":
                    "PRAD",
            }
        ]
    )

    result = (
        MoradiQTLStandardizer
        .standardize_ld_dataframe(
            dataframe,
            source_table="S2",
            feature_type="intron",
            qtl_scope="cis",
            feature_column="splicing_event",
        )
    )

    assert len(
        result
    ) == 1

    assert (
        result.iloc[
            0
        ][
            "qtl_rsid"
        ]
        == "rs1"
    )


def test_standardize_s3_dataframe_with_t_column() -> None:
    """S3 DataFrame API must support the trans-exon T source column."""

    dataframe = pd.DataFrame(
        [
            {
                "T":
                    "rs1",
                "SNP_pos":
                    "1:100",
                "splicing_event":
                    "EX1",
                "p_value":
                    0.01,
            }
        ]
    )

    result = (
        MoradiQTLStandardizer
        .standardize_s3_dataframe(
            dataframe,
            feature_type="exon",
            feature_column="splicing_event",
            variant_column="T",
        )
    )

    assert len(
        result
    ) == 1

    assert (
        result.iloc[
            0
        ][
            "canonical_rsid"
        ]
        == "rs1"
    )