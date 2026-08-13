"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_moradi.py

Description:
    Unit tests for Moradi QTL field-level validation utilities.

    Tests cover SNP identifiers, genomic coordinates, SNV/indel/
    multi-allelic representations, splicing-event identifiers,
    Ensembl transcript identifiers, probability fields, LD values,
    MAF, genotype call rate, imputation Rsq, and generic numerical
    statistics.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pcatrfqtl.validation.validators.moradi import (
    MoradiQTLValidator,
)


def create_validator() -> MoradiQTLValidator:
    """Create a Moradi QTL validator."""

    return MoradiQTLValidator()


# ===========================================================================
# SNP identifiers
# ===========================================================================


def test_valid_snp_ids() -> None:
    """Valid dbSNP identifiers should pass."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "rs2153655",
            "rs4612711",
            "RS11186696",
        ]
    )

    assert result == []


def test_invalid_snp_ids() -> None:
    """Invalid SNP identifiers should be detected."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "rs2153655",
            "2153655",
            "invalid",
            "",
        ]
    )

    assert len(result) == 3


# ===========================================================================
# Genomic coordinates
# ===========================================================================


def test_valid_coordinates() -> None:
    """Valid genomic coordinates should pass."""

    validator = create_validator()

    result = validator.validate_coordinates(
        [
            "10:93563056",
            "1:100",
            "22:999999",
            "X:12345",
            "Y:42",
            "MT:100",
            "chr10:93563056",
        ]
    )

    assert result == []


def test_invalid_coordinates() -> None:
    """Malformed genomic coordinates should be detected."""

    validator = create_validator()

    result = validator.validate_coordinates(
        [
            "10:93563056",
            "10",
            "chr23:100",
            "1:0",
            "invalid",
        ]
    )

    assert len(result) == 4


# ===========================================================================
# Alleles
# ===========================================================================


def test_valid_alleles() -> None:
    """
    Moradi alleles may include SNVs, indels, and multi-allelic values.
    """

    validator = create_validator()

    result = validator.validate_alleles(
        [
            "A",
            "C",
            "G",
            "T",
            "CT",
            "GGT",
            "A,T",
            "C,T",
        ],
        "A",
    )

    assert result == []


def test_invalid_alleles() -> None:
    """Invalid allele representations should be detected."""

    validator = create_validator()

    result = validator.validate_alleles(
        [
            "A",
            "X",
            "A/T",
            "A,,T",
            "",
        ],
        "A",
    )

    assert len(result) == 4


# ===========================================================================
# Splicing-event identifiers
# ===========================================================================


def test_valid_splicing_events() -> None:
    """Valid intron and exon event identifiers should pass."""

    validator = create_validator()

    result = validator.validate_splicing_events(
        [
            "INT132623",
            "INT78964",
            "EX342488",
            "EX420915",
        ]
    )

    assert result == []


def test_invalid_splicing_events() -> None:
    """Malformed splicing-event identifiers should fail."""

    validator = create_validator()

    result = validator.validate_splicing_events(
        [
            "INT132623",
            "ENST00000318325.6",
            "EXABC",
            "",
        ]
    )

    assert len(result) == 3


# ===========================================================================
# Transcript identifiers
# ===========================================================================


def test_valid_transcripts() -> None:
    """Valid Ensembl transcript identifiers should pass."""

    validator = create_validator()

    result = validator.validate_transcripts(
        [
            "ENST00000318325.6",
            "ENST00000370575.5",
            "ENST00000257915.10",
            "ENST00000603017",
        ]
    )

    assert result == []


def test_invalid_transcripts() -> None:
    """Malformed transcript identifiers should be detected."""

    validator = create_validator()

    result = validator.validate_transcripts(
        [
            "ENST00000318325.6",
            "ENSG00000123456",
            "INT12345",
            "",
        ]
    )

    assert len(result) == 3


# ===========================================================================
# P-values / FDR
# ===========================================================================


def test_valid_p_values() -> None:
    """Valid P-values should pass."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0,
            2.130980e-55,
            0.05,
            1,
        ]
    )

    assert result == []


def test_invalid_p_values() -> None:
    """Out-of-range and non-numeric P-values should fail."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            -0.1,
            1.1,
            "invalid",
        ]
    )

    assert len(result) == 3


def test_valid_fdr() -> None:
    """Valid FDR values should pass."""

    validator = create_validator()

    result = validator.validate_fdr(
        [
            0,
            8.643118e-49,
            0.05,
            1,
        ]
    )

    assert result == []


def test_invalid_fdr() -> None:
    """Invalid FDR values should fail."""

    validator = create_validator()

    result = validator.validate_fdr(
        [
            -0.01,
            1.01,
        ]
    )

    assert len(result) == 2


# ===========================================================================
# LD
# ===========================================================================


def test_valid_ld() -> None:
    """Valid LD r² values should pass."""

    validator = create_validator()

    result = validator.validate_ld(
        [
            0,
            0.540508,
            0.696408,
            0.977277,
            1,
        ]
    )

    assert result == []


def test_invalid_ld() -> None:
    """Invalid LD values should fail."""

    validator = create_validator()

    result = validator.validate_ld(
        [
            -0.1,
            1.1,
        ]
    )

    assert len(result) == 2


# ===========================================================================
# MAF
# ===========================================================================


def test_valid_maf() -> None:
    """Valid MAF values should pass."""

    validator = create_validator()

    result = validator.validate_maf(
        [
            0,
            0.07348,
            0.48173,
            1,
        ]
    )

    assert result == []


def test_invalid_maf() -> None:
    """Invalid MAF values should fail."""

    validator = create_validator()

    result = validator.validate_maf(
        [
            -0.01,
            1.01,
        ]
    )

    assert len(result) == 2


# ===========================================================================
# Call rate
# ===========================================================================


def test_valid_call_rate() -> None:
    """Valid genotype call-rate values should pass."""

    validator = create_validator()

    result = validator.validate_call_rate(
        [
            0,
            0.99570,
            0.99998,
            1,
        ]
    )

    assert result == []


def test_invalid_call_rate() -> None:
    """Invalid call-rate values should fail."""

    validator = create_validator()

    result = validator.validate_call_rate(
        [
            -0.1,
            1.2,
        ]
    )

    assert len(result) == 2


# ===========================================================================
# Rsq
# ===========================================================================


def test_valid_rsq() -> None:
    """Valid imputation Rsq values should pass."""

    validator = create_validator()

    result = validator.validate_rsq(
        [
            0,
            0.82606,
            0.99198,
            1,
        ]
    )

    assert result == []


def test_invalid_rsq() -> None:
    """Invalid Rsq values should fail."""

    validator = create_validator()

    result = validator.validate_rsq(
        [
            -0.1,
            1.01,
        ]
    )

    assert len(result) == 2


# ===========================================================================
# Generic numerical fields
# ===========================================================================


def test_valid_numeric_fields() -> None:
    """Statistical and effect values may be positive or negative."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            -18.660559,
            0,
            410.594200,
            -519.560942,
        ],
        "beta",
    )

    assert result == []


def test_invalid_numeric_fields() -> None:
    """Non-numeric values should be detected."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            1.0,
            "invalid",
        ],
        "beta",
    )

    assert len(result) == 1

def test_valid_composite_snp_ids() -> None:
    """Single and colon-separated dbSNP identifiers should pass."""

    validator = create_validator()

    result = validator.validate_composite_snp_ids(
        [
            "rs7103835",
            "rs559911988:rs11348515",
            "rs545285266:rs568161688:rs150973389",
        ],
        "tag_SNP_pos",
    )

    assert result == []


def test_invalid_composite_snp_ids() -> None:
    """Malformed composite dbSNP identifiers should fail."""

    validator = create_validator()

    result = validator.validate_composite_snp_ids(
        [
            "rs12345",
            "rs12345:",
            ":rs12345",
            "rs12345:123456",
            "rs12345::rs67890",
            "invalid",
        ],
        "tag_SNP_pos",
    )

    assert len(result) == 5