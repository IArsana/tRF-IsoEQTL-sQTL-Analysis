"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_trfqtl.py

Description:
    Unit tests for Cancer-tRFQTL validation.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.validation.validators.trfqtl import (
    TRFQTLValidator,
)


def create_validator() -> TRFQTLValidator:
    """Create a tRFQTL validator."""

    return TRFQTLValidator()


def test_trfqtl_snp_ids() -> None:
    """Test valid SNP identifiers."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "rs123",
            "rs456789",
            "RS123",
        ]
    )

    assert result == []


def test_trfqtl_invalid_snp_ids() -> None:
    """Test invalid SNP identifiers."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "rs123",
            "123456",
            "chr1:12345",
        ]
    )

    assert len(result) == 2


def test_trfqtl_coordinates() -> None:
    """Test valid genomic coordinates."""

    validator = create_validator()

    result = validator.validate_coordinates(
        [
            "chr1:147308207",
            "15:25826064",
            "chrX:123456",
        ]
    )

    assert result == []


def test_trfqtl_invalid_coordinates() -> None:
    """Test invalid genomic coordinates."""

    validator = create_validator()

    result = validator.validate_coordinates(
        [
            "chr1",
            "chr1:0",
            "abc:123",
        ]
    )

    assert len(result) == 3


def test_trfqtl_alleles() -> None:
    """Test valid allele representations."""

    validator = create_validator()

    result = validator.validate_alleles(
        [
            "A/T",
            "G/C",
            "C",
            "A",
            "AT/GC",
        ]
    )

    assert result == []


def test_trfqtl_invalid_alleles() -> None:
    """Test invalid allele representations."""

    validator = create_validator()

    result = validator.validate_alleles(
        [
            "A/T",
            "ATGC",
            "X/Y",
        ]
    )

    assert len(result) == 2


def test_trfqtl_ids() -> None:
    """Test valid tRF identifiers."""

    validator = create_validator()

    result = validator.validate_trf_ids(
        [
            "tRF-30-34HWH3RXSINH",
            "tRF-23-V47PU9XW0N",
        ]
    )

    assert result == []


def test_trfqtl_invalid_ids() -> None:
    """Test invalid tRF identifiers."""

    validator = create_validator()

    result = validator.validate_trf_ids(
        [
            "tRF-30-34HWH3RXSINH",
            "miR-21",
        ]
    )

    assert len(result) == 1


def test_trfqtl_p_values() -> None:
    """Test valid P-values."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0.05,
            1e-10,
            1.0,
        ]
    )

    assert result == []


def test_trfqtl_invalid_p_values() -> None:
    """Test invalid P-values."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0.05,
            0,
            1.5,
        ]
    )

    assert len(result) == 2


def test_trfqtl_fdr() -> None:
    """Test valid FDR values."""

    validator = create_validator()

    result = validator.validate_fdr(
        [
            0,
            0.001,
            1.0,
        ]
    )

    assert result == []


def test_trfqtl_ld_r2() -> None:
    """Test valid LD r2 values."""

    validator = create_validator()

    result = validator.validate_ld_r2(
        [
            0,
            0.5,
            1,
        ]
    )

    assert result == []


def test_trfqtl_invalid_ld_r2() -> None:
    """Test invalid LD r2 values."""

    validator = create_validator()

    result = validator.validate_ld_r2(
        [
            -0.1,
            0.5,
            1.1,
        ]
    )

    assert len(result) == 2