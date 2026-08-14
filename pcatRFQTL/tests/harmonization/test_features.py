"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/harmonization/test_features.py

Description:
    Unit tests for M3.4 biological feature identity harmonization.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pcatrfqtl.harmonization.features import (
    FeatureHarmonizer,
    FeatureIdentifierSystem,
    FeatureIdentityStatus,
    FeatureType,
)


def test_classify_transcript() -> None:
    """ENST identifiers should classify as transcript isoforms."""

    assert (
        FeatureHarmonizer.classify(
            "ENST00000264639.9"
        )
        is FeatureType.TRANSCRIPT_ISOFORM
    )


def test_classify_exon() -> None:
    """EX identifiers should classify as exon events."""

    assert (
        FeatureHarmonizer.classify(
            "EX12345"
        )
        is FeatureType.EXON_EVENT
    )


def test_classify_intron() -> None:
    """INT identifiers should classify as intron events."""

    assert (
        FeatureHarmonizer.classify(
            "INT12345"
        )
        is FeatureType.INTRON_EVENT
    )


def test_scientific_intron_still_classifies_as_intron() -> None:
    """
    Scientific-notation-like source anomaly still belongs to the
    intron identifier class without becoming canonical.
    """

    assert (
        FeatureHarmonizer.classify(
            "INT1e+05"
        )
        is FeatureType.INTRON_EVENT
    )


def test_harmonize_transcript_with_version() -> None:
    """Transcript version should be preserved."""

    result = (
        FeatureHarmonizer
        .harmonize_transcript(
            "ENST00000264639.9"
        )
    )

    assert (
        result.canonical_feature_id
        == "ENST00000264639.9"
    )

    assert (
        result.base_identifier
        == "ENST00000264639"
    )

    assert (
        result.version
        == "9"
    )

    assert (
        result.identity_key
        == "transcript:ENST00000264639.9"
    )

    assert result.usable is True

    assert (
        result.status
        is FeatureIdentityStatus.RESOLVED
    )


def test_harmonize_transcript_without_version() -> None:
    """Unversioned ENST identifiers remain valid."""

    result = (
        FeatureHarmonizer
        .harmonize_transcript(
            "ENST00000264639"
        )
    )

    assert (
        result.base_identifier
        == "ENST00000264639"
    )

    assert result.version is None
    assert result.usable is True


def test_harmonize_exon_event() -> None:
    """Canonical exon events should remain usable."""

    result = (
        FeatureHarmonizer
        .harmonize_exon_event(
            "EX123"
        )
    )

    assert (
        result.identity_key
        == "exon:EX123"
    )

    assert (
        result.identifier_system
        is FeatureIdentifierSystem.MORADI_EXON_EVENT
    )

    assert result.usable is True


def test_harmonize_intron_event() -> None:
    """Canonical intron events should remain usable."""

    result = (
        FeatureHarmonizer
        .harmonize_intron_event(
            "INT456"
        )
    )

    assert (
        result.identity_key
        == "intron:INT456"
    )

    assert result.usable is True


def test_malformed_intron_is_preserved() -> None:
    """Known INT1e+05 source anomaly must never be auto-corrected."""

    result = (
        FeatureHarmonizer
        .harmonize_intron_event(
            "INT1e+05"
        )
    )

    assert (
        result.source_feature_id
        == "INT1e+05"
    )

    assert (
        result.canonical_feature_id
        is None
    )

    assert result.usable is False

    assert (
        result.source_anomaly
        is True
    )

    assert (
        result.status
        is FeatureIdentityStatus.PARTIAL
    )


def test_trf_identifier_is_preserved() -> None:
    """tRF identifiers should remain in their source namespace."""

    result = (
        FeatureHarmonizer
        .harmonize(
            "tRF-example",
            feature_type=FeatureType.TRF,
        )
    )

    assert (
        result.canonical_feature_id
        == "tRF-example"
    )

    assert (
        result.identifier_system
        is FeatureIdentifierSystem.SOURCE_TRF
    )

    assert (
        result.identity_key
        == "trf:source:tRF-example"
    )

    assert result.usable is True

    assert (
        result.status
        is FeatureIdentityStatus.PRESERVED
    )


def test_ensembl_gene_identifier() -> None:
    """Ensembl gene IDs should be canonicalized."""

    result = (
        FeatureHarmonizer
        .harmonize_gene(
            "ENSG00000141510.18"
        )
    )

    assert (
        result.canonical_feature_id
        == "ENSG00000141510.18"
    )

    assert (
        result.base_identifier
        == "ENSG00000141510"
    )

    assert (
        result.version
        == "18"
    )

    assert (
        result.identifier_system
        is FeatureIdentifierSystem.ENSEMBL_GENE
    )


def test_gene_symbol_is_preserved() -> None:
    """Gene symbols must not be silently mapped to Ensembl."""

    result = (
        FeatureHarmonizer
        .harmonize_gene(
            "TP53"
        )
    )

    assert (
        result.canonical_feature_id
        == "TP53"
    )

    assert (
        result.identity_key
        == "gene:symbol:TP53"
    )

    assert (
        result.identifier_system
        is FeatureIdentifierSystem.GENE_SYMBOL
    )

    assert (
        result.status
        is FeatureIdentityStatus.PRESERVED
    )


def test_unknown_identifier_is_preserved_but_not_usable() -> None:
    """Unknown feature classes should not participate in joins."""

    result = (
        FeatureHarmonizer
        .harmonize(
            "UNKNOWN_FEATURE_123"
        )
    )

    assert (
        result.feature_type
        is FeatureType.UNKNOWN
    )

    assert (
        result.source_feature_id
        == "UNKNOWN_FEATURE_123"
    )

    assert result.usable is False


def test_same_transcript_features_match() -> None:
    """Identical transcript identities should compare equal."""

    left = (
        FeatureHarmonizer
        .harmonize(
            "ENST00000264639.9"
        )
    )

    right = (
        FeatureHarmonizer
        .harmonize(
            "ENST00000264639.9"
        )
    )

    assert (
        FeatureHarmonizer
        .features_equal(
            left,
            right,
        )
        is True
    )


def test_different_transcript_versions_do_not_match() -> None:
    """
    Transcript versions remain distinct during conservative M3.4
    matching.
    """

    left = (
        FeatureHarmonizer
        .harmonize(
            "ENST00000264639.8"
        )
    )

    right = (
        FeatureHarmonizer
        .harmonize(
            "ENST00000264639.9"
        )
    )

    assert (
        FeatureHarmonizer
        .features_equal(
            left,
            right,
        )
        is False
    )


def test_feature_classes_do_not_cross_match() -> None:
    """Different biological feature classes must remain distinct."""

    exon = (
        FeatureHarmonizer
        .harmonize(
            "EX123"
        )
    )

    intron = (
        FeatureHarmonizer
        .harmonize(
            "INT123"
        )
    )

    assert (
        FeatureHarmonizer
        .features_equal(
            exon,
            intron,
        )
        is False
    )


def test_feature_to_dict() -> None:
    """Feature models should be JSON serializable."""

    result = (
        FeatureHarmonizer
        .harmonize(
            "ENST00000264639.9"
        )
    )

    payload = result.to_dict()

    assert (
        payload["feature_type"]
        == "TRANSCRIPT_ISOFORM"
    )

    assert (
        payload["identifier_system"]
        == "ENSEMBL_TRANSCRIPT"
    )

    assert (
        payload["status"]
        == "RESOLVED"
    )