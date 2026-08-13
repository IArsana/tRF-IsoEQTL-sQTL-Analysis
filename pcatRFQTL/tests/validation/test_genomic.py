"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_genomic.py

Description:
    Unit tests for genomic coordinate parsing and validation.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.validation.genomic import (
    GenomicCoordinate,
    is_valid_genomic_coordinate,
    parse_genomic_coordinate,
)


def test_parse_chr_coordinate() -> None:
    """Test parsing coordinates with a chr prefix."""

    result = parse_genomic_coordinate(
        "chr1:147308207"
    )

    assert result == GenomicCoordinate(
        chromosome="1",
        position=147308207,
    )


def test_parse_coordinate_without_chr() -> None:
    """Test parsing coordinates without a chr prefix."""

    result = parse_genomic_coordinate(
        "15:25826064"
    )

    assert result == GenomicCoordinate(
        chromosome="15",
        position=25826064,
    )


def test_parse_sex_chromosome() -> None:
    """Test parsing sex chromosome coordinates."""

    result = parse_genomic_coordinate(
        "chrX:123456"
    )

    assert result == GenomicCoordinate(
        chromosome="X",
        position=123456,
    )


def test_invalid_coordinate() -> None:
    """Test rejection of invalid genomic coordinates."""

    assert not is_valid_genomic_coordinate(
        "chr1"
    )

    assert not is_valid_genomic_coordinate(
        "chr1:0"
    )

    assert not is_valid_genomic_coordinate(
        "abc:123"
    )


def test_non_string_coordinate() -> None:
    """Test rejection of non-string coordinates."""

    assert not is_valid_genomic_coordinate(
        123456  # type: ignore[arg-type]
    )