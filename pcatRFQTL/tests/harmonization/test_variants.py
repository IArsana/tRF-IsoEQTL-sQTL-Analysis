"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/harmonization/test_variants.py

Description:
    Unit tests for M3.3 variant identity harmonization.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pcatrfqtl.harmonization.variants import (
    VariantHarmonizer,
    VariantIdentityStatus,
    VariantMatchMethod,
)


def test_normalize_rsid() -> None:
    """Single rsIDs should normalize to lowercase rs prefix."""

    assert (
        VariantHarmonizer
        .normalize_rsid(
            "RS123"
        )
        == "rs123"
    )

    assert (
        VariantHarmonizer
        .normalize_rsid(
            "rs999"
        )
        == "rs999"
    )


def test_invalid_rsid_returns_none() -> None:
    """Unsupported variant descriptors must not become rsIDs."""

    assert (
        VariantHarmonizer
        .normalize_rsid(
            "6"
        )
        is None
    )

    assert (
        VariantHarmonizer
        .normalize_rsid(
            "rs1:rs2"
        )
        is None
    )

    assert (
        VariantHarmonizer
        .normalize_rsid(
            "rs1 x rs2"
        )
        is None
    )


def test_harmonize_variant_with_rsid_and_coordinate() -> None:
    """Complete hg19 variant identity should be fully resolved."""

    result = (
        VariantHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            source_variant_id="rs12786544",
            rsid="rs12786544",
            coordinate="11:69561656",
        )
    )

    assert (
        result.canonical_rsid
        == "rs12786544"
    )

    assert (
        result.rsid_usable
        is True
    )

    assert (
        result.coordinate.build_aware_key
        == "hg19:11:69561656"
    )

    assert (
        result.identity_key
        == "rsid:rs12786544"
    )

    assert (
        result.identity_method
        is VariantMatchMethod.RSID_AND_COORDINATE
    )

    assert (
        result.status
        is VariantIdentityStatus.RESOLVED
    )


def test_harmonize_variant_with_rsid_only() -> None:
    """A valid rsID remains usable without genomic coordinates."""

    result = (
        VariantHarmonizer
        .harmonize(
            dataset="gwas_catalog",
            source_variant_id="rs123",
            rsid="rs123",
        )
    )

    assert (
        result.canonical_rsid
        == "rs123"
    )

    assert (
        result.identity_key
        == "rsid:rs123"
    )

    assert (
        result.identity_method
        is VariantMatchMethod.RSID
    )

    assert (
        result.status
        is VariantIdentityStatus.RESOLVED
    )


def test_coordinate_fallback_for_non_rsid_trfqtl_variant() -> None:
    """
    Cancer-tRFQTL source value '6' should remain non-rsID but retain
    usable hg19 coordinate identity.
    """

    result = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            source_variant_id="6",
            rsid="6",
            coordinate="6:28958399",
        )
    )

    assert (
        result.source_variant_id
        == "6"
    )

    assert (
        result.canonical_rsid
        is None
    )

    assert (
        result.rsid_usable
        is False
    )

    assert (
        result.coordinate.build_aware_key
        == "hg19:6:28958399"
    )

    assert (
        result.identity_key
        == "coord:hg19:6:28958399"
    )

    assert (
        result.identity_method
        is VariantMatchMethod.COORDINATE
    )

    assert (
        result.status
        is VariantIdentityStatus.PARTIAL
    )


def test_unresolved_variant() -> None:
    """Variant without usable rsID or coordinate should remain unresolved."""

    result = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            source_variant_id="unknown",
            rsid="unknown",
        )
    )

    assert (
        result.identity_key
        is None
    )

    assert (
        result.identity_method
        is VariantMatchMethod.NONE
    )

    assert (
        result.status
        is VariantIdentityStatus.UNRESOLVED
    )


def test_same_rsid_matches_across_datasets() -> None:
    """Canonical rsID should support cross-source matching."""

    left = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            rsid="rs12345",
            coordinate="1:100",
        )
    )

    right = (
        VariantHarmonizer
        .harmonize(
            dataset="gwas_catalog",
            rsid="RS12345",
            chromosome="1",
            position=999999,
        )
    )

    assert (
        VariantHarmonizer
        .variants_equal(
            left,
            right,
        )
        is True
    )

    assert (
        VariantHarmonizer
        .match_method(
            left,
            right,
        )
        is VariantMatchMethod.RSID
    )


def test_same_verified_hg19_coordinate_matches_without_rsid() -> None:
    """Coordinate fallback should work between verified hg19 sources."""

    left = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            source_variant_id="6",
            rsid="6",
            coordinate="6:28958399",
        )
    )

    right = (
        VariantHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            coordinate="6:28958399",
        )
    )

    assert (
        VariantHarmonizer
        .variants_equal(
            left,
            right,
        )
        is True
    )

    assert (
        VariantHarmonizer
        .match_method(
            left,
            right,
        )
        is VariantMatchMethod.COORDINATE
    )


def test_gwas_coordinate_fallback_is_blocked() -> None:
    """Unresolved GWAS build must block coordinate-only matching."""

    left = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            coordinate="6:28958399",
        )
    )

    right = (
        VariantHarmonizer
        .harmonize(
            dataset="gwas_catalog",
            coordinate="6:28958399",
        )
    )

    assert (
        VariantHarmonizer
        .variants_equal(
            left,
            right,
        )
        is False
    )

    assert (
        VariantHarmonizer
        .match_method(
            left,
            right,
        )
        is VariantMatchMethod.NONE
    )


def test_rsid_match_is_independent_of_coordinate_build() -> None:
    """
    rsID matching may still identify the same dbSNP variant while
    coordinate comparison remains unavailable.
    """

    left = (
        VariantHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            rsid="rs555",
            coordinate="1:1000",
        )
    )

    right = (
        VariantHarmonizer
        .harmonize(
            dataset="gwas_catalog",
            rsid="rs555",
            coordinate="1:1000",
        )
    )

    assert (
        VariantHarmonizer
        .rsids_equal(
            left,
            right,
        )
        is True
    )

    assert (
        VariantHarmonizer
        .coordinates_equal(
            left,
            right,
        )
        is False
    )

    assert (
        VariantHarmonizer
        .variants_equal(
            left,
            right,
        )
        is True
    )


def test_rsid_and_coordinate_match() -> None:
    """Two verified records may agree on both identity paths."""

    left = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            rsid="rs100",
            coordinate="2:200",
        )
    )

    right = (
        VariantHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            rsid="rs100",
            coordinate="2:200",
        )
    )

    assert (
        VariantHarmonizer
        .match_method(
            left,
            right,
        )
        is VariantMatchMethod.RSID_AND_COORDINATE
    )


def test_same_rsid_different_verified_coordinate_still_matches_by_rsid() -> None:
    """
    M3.3 records identity evidence but does not yet adjudicate rsID /
    coordinate discordance.
    """

    left = (
        VariantHarmonizer
        .harmonize(
            dataset="cancer_trfqtl",
            rsid="rs100",
            coordinate="2:200",
        )
    )

    right = (
        VariantHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            rsid="rs100",
            coordinate="2:201",
        )
    )

    assert (
        VariantHarmonizer
        .rsids_equal(
            left,
            right,
        )
        is True
    )

    assert (
        VariantHarmonizer
        .coordinates_equal(
            left,
            right,
        )
        is False
    )

    assert (
        VariantHarmonizer
        .variants_equal(
            left,
            right,
        )
        is True
    )

    assert (
        VariantHarmonizer
        .match_method(
            left,
            right,
        )
        is VariantMatchMethod.RSID
    )


def test_variant_to_dict() -> None:
    """Variant model should be JSON-serializable."""

    result = (
        VariantHarmonizer
        .harmonize(
            dataset="moradi_qtl",
            rsid="rs10",
            coordinate="1:10",
        )
    )

    payload = (
        result.to_dict()
    )

    assert (
        payload[
            "canonical_rsid"
        ]
        == "rs10"
    )

    assert (
        payload[
            "identity_method"
        ]
        == "RSID_AND_COORDINATE"
    )

    assert (
        payload[
            "status"
        ]
        == "RESOLVED"
    )

    assert (
        payload[
            "coordinate"
        ][
            "build_aware_key"
        ]
        == "hg19:1:10"
    )