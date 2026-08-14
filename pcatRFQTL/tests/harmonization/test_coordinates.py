"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/harmonization/test_coordinates.py

Description:
    Unit tests for M3.2 canonical genomic coordinate harmonization.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pcatrfqtl.harmonization.coordinates import (
    CoordinateHarmonizer,
)


def test_normalize_autosomal_chromosome() -> None:
    """Autosomal chromosome labels should be canonicalized."""

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "chr6"
        )
        == "6"
    )

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "06"
        )
        == "6"
    )


def test_normalize_sex_chromosome() -> None:
    """Sex chromosomes should be canonicalized."""

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "chrX"
        )
        == "X"
    )

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "y"
        )
        == "Y"
    )


def test_normalize_mitochondrial_chromosome() -> None:
    """Mitochondrial aliases should normalize to MT."""

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "chrM"
        )
        == "MT"
    )

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "M"
        )
        == "MT"
    )


def test_invalid_chromosome_returns_none() -> None:
    """Unsupported chromosome labels should not be corrected."""

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "chr23"
        )
        is None
    )

    assert (
        CoordinateHarmonizer
        .normalize_chromosome(
            "GL000220.1"
        )
        is None
    )


def test_normalize_position() -> None:
    """Positive positions should normalize to integers."""

    assert (
        CoordinateHarmonizer
        .normalize_position(
            "28958399"
        )
        == 28958399
    )


def test_invalid_position() -> None:
    """Invalid positions should return None."""

    assert (
        CoordinateHarmonizer
        .normalize_position(
            0
        )
        is None
    )

    assert (
        CoordinateHarmonizer
        .normalize_position(
            -1
        )
        is None
    )

    assert (
        CoordinateHarmonizer
        .normalize_position(
            "abc"
        )
        is None
    )


def test_parse_coordinate() -> None:
    """Chromosome-position string should be decomposed."""

    result = (
        CoordinateHarmonizer
        .parse_coordinate(
            "chr6:28958399"
        )
    )

    assert result == (
        "6",
        28958399,
    )


def test_parse_moradi_coordinate() -> None:
    """Moradi hg19 coordinate format should parse directly."""

    result = (
        CoordinateHarmonizer
        .parse_coordinate(
            "11:69561656"
        )
    )

    assert result == (
        "11",
        69561656,
    )


def test_invalid_coordinate_returns_none() -> None:
    """Malformed coordinates must remain unusable."""

    result = (
        CoordinateHarmonizer
        .parse_coordinate(
            "chr6"
        )
    )

    assert result == (
        None,
        None,
    )


def test_coordinate_key() -> None:
    """Canonical coordinate key should omit build information."""

    result = (
        CoordinateHarmonizer
        .make_coordinate_key(
            "6",
            28958399,
        )
    )

    assert (
        result
        == "6:28958399"
    )


def test_build_aware_coordinate_key() -> None:
    """Build-aware key should protect against assembly ambiguity."""

    result = (
        CoordinateHarmonizer
        .make_build_aware_key(
            genome_build="hg19",
            chromosome="6",
            position=28958399,
        )
    )

    assert (
        result
        == "hg19:6:28958399"
    )


def test_harmonize_cancer_trfqtl_coordinate() -> None:
    """Cancer-tRFQTL coordinates should be usable hg19 coordinates."""

    result = (
        CoordinateHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            coordinate="6:28958399",
        )
    )

    assert (
        result.chromosome
        == "6"
    )

    assert (
        result.position
        == 28958399
    )

    assert (
        result.coordinate_key
        == "6:28958399"
    )

    assert (
        result.build_aware_key
        == "hg19:6:28958399"
    )

    assert (
        result.source_genome_build
        == "hg19"
    )

    assert (
        result.source_assembly
        == "GRCh37"
    )

    assert (
        result.coordinate_usable
        is True
    )

    assert (
        result.coordinate_join_allowed
        is True
    )

    assert (
        result.liftover_required
        is True
    )

    assert (
        result.liftover_performed
        is False
    )


def test_harmonize_moradi_coordinate() -> None:
    """Moradi coordinates should produce hg19 build-aware keys."""

    result = (
        CoordinateHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            coordinate="11:69561656",
        )
    )

    assert (
        result.build_aware_key
        == "hg19:11:69561656"
    )

    assert (
        result.coordinate_join_allowed
        is True
    )


def test_gwas_coordinate_join_is_blocked() -> None:
    """GWAS coordinates must remain blocked while build is unresolved."""

    result = (
        CoordinateHarmonizer
        .harmonize(
            dataset="gwas_catalog",
            chromosome="6",
            position=28958399,
        )
    )

    assert (
        result.coordinate_usable
        is True
    )

    assert (
        result.source_genome_build
        == "source_dependent"
    )

    assert (
        result.coordinate_join_allowed
        is False
    )

    assert (
        result.liftover_performed
        is False
    )


def test_trfqtl_and_moradi_are_matchable() -> None:
    """Verified hg19 sources should permit coordinate comparison."""

    result = (
        CoordinateHarmonizer
        .are_coordinates_matchable(
            left_dataset=(
                "cancer_trfqtl"
            ),
            left_chromosome="6",
            left_position=28958399,
            right_dataset=(
                "moradi_qtl"
            ),
            right_chromosome="6",
            right_position=28958399,
        )
    )

    assert result is True


def test_trfqtl_and_gwas_are_not_matchable() -> None:
    """Unresolved GWAS build should block coordinate comparison."""

    result = (
        CoordinateHarmonizer
        .are_coordinates_matchable(
            left_dataset=(
                "cancer_trfqtl"
            ),
            left_chromosome="6",
            left_position=28958399,
            right_dataset=(
                "gwas_catalog"
            ),
            right_chromosome="6",
            right_position=28958399,
        )
    )

    assert result is False


def test_equal_hg19_coordinates() -> None:
    """Identical coordinates on compatible builds should match."""

    result = (
        CoordinateHarmonizer
        .coordinates_equal(
            left_dataset=(
                "cancer_trfqtl"
            ),
            left_chromosome="6",
            left_position=28958399,
            right_dataset=(
                "moradi_qtl"
            ),
            right_chromosome="6",
            right_position=28958399,
        )
    )

    assert result is True


def test_different_hg19_coordinates_are_not_equal() -> None:
    """Different coordinates must not match."""

    result = (
        CoordinateHarmonizer
        .coordinates_equal(
            left_dataset=(
                "cancer_trfqtl"
            ),
            left_chromosome="6",
            left_position=28958399,
            right_dataset=(
                "moradi_qtl"
            ),
            right_chromosome="6",
            right_position=28958400,
        )
    )

    assert result is False


def test_gwas_coordinate_cannot_be_declared_equal() -> None:
    """Same numeric coordinate is insufficient with unresolved build."""

    result = (
        CoordinateHarmonizer
        .coordinates_equal(
            left_dataset=(
                "cancer_trfqtl"
            ),
            left_chromosome="6",
            left_position=28958399,
            right_dataset=(
                "gwas_catalog"
            ),
            right_chromosome="6",
            right_position=28958399,
        )
    )

    assert result is False


def test_invalid_coordinate_is_not_joinable() -> None:
    """Invalid genomic positions should disable joins."""

    result = (
        CoordinateHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            chromosome="6",
            position=0,
        )
    )

    assert (
        result.coordinate_usable
        is False
    )

    assert (
        result.coordinate_key
        is None
    )

    assert (
        result.build_aware_key
        is None
    )

    assert (
        result.coordinate_join_allowed
        is False
    )


def test_coordinate_to_dict() -> None:
    """Coordinate model should be serializable."""

    result = (
        CoordinateHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            coordinate="1:100",
        )
    )

    payload = result.to_dict()

    assert (
        payload[
            "coordinate_key"
        ]
        == "1:100"
    )

    assert (
        payload[
            "build_aware_key"
        ]
        == "hg19:1:100"
    )

def test_coordinate_key_alone_is_not_cross_build_identity() -> None:
    """
    Numeric coordinate identity must not override genome-build
    compatibility.
    """

    trfqtl = (
        CoordinateHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            coordinate="6:28958399",
        )
    )

    gwas = (
        CoordinateHarmonizer
        .harmonize(
            dataset="gwas_catalog",
            coordinate="6:28958399",
        )
    )

    assert (
        trfqtl.coordinate_key
        == gwas.coordinate_key
    )

    assert (
        trfqtl.build_aware_key
        != gwas.build_aware_key
    )

    assert (
        CoordinateHarmonizer
        .coordinates_equal(
            left_dataset="cancer_trfqtl",
            left_chromosome="6",
            left_position=28958399,
            right_dataset="gwas_catalog",
            right_chromosome="6",
            right_position=28958399,
        )
        is False
    )