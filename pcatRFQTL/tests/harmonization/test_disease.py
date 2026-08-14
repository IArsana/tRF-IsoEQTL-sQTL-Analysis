"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/harmonization/test_disease.py

Description:
    Unit tests for M3.5 disease-context harmonization.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pcatrfqtl.harmonization.disease import (
    DiseaseContextType,
    DiseaseHarmonizer,
    DiseaseIdentifierSystem,
    DiseaseIdentityStatus,
)


def test_normalize_prostate_cancer() -> None:
    """Disease labels should normalize for matching."""

    assert (
        DiseaseHarmonizer
        .normalize_for_matching(
            "Prostate Cancer"
        )
        == "prostate cancer"
    )


def test_normalize_underscore_label() -> None:
    """Underscore source labels should normalize to spaces."""

    assert (
        DiseaseHarmonizer
        .normalize_for_matching(
            "prostate_cancer"
        )
        == "prostate cancer"
    )


def test_prad_is_primary_disease() -> None:
    """TCGA PRAD code should resolve to prostate cancer."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "PRAD"
        )
    )

    assert (
        result.canonical_disease_id
        == "prostate_cancer"
    )

    assert (
        result.canonical_name
        == "Prostate cancer"
    )

    assert (
        result.tcga_code
        == "PRAD"
    )

    assert (
        result.is_primary_disease
        is True
    )

    assert (
        result.context_type
        is DiseaseContextType.PRIMARY_DISEASE
    )

    assert (
        result.status
        is DiseaseIdentityStatus.RESOLVED
    )


def test_prca_is_primary_disease() -> None:
    """Moradi PrCa terminology should resolve to prostate cancer."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "PrCa"
        )
    )

    assert (
        result.is_primary_disease
        is True
    )

    assert (
        result.identity_key
        == "disease:prostate_cancer"
    )


def test_prostate_cancer_label_is_primary() -> None:
    """GWAS prostate cancer label should resolve canonically."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "Prostate cancer"
        )
    )

    assert (
        result.is_primary_disease
        is True
    )

    assert (
        result.canonical_disease_id
        == "prostate_cancer"
    )


def test_prostate_carcinoma_is_primary() -> None:
    """Alternative prostate carcinoma terminology should resolve."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "prostate carcinoma"
        )
    )

    assert (
        result.is_primary_disease
        is True
    )


def test_prostate_cancer_underscore_is_primary() -> None:
    """Machine-readable prostate_cancer label should resolve."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "prostate_cancer"
        )
    )

    assert (
        result.is_primary_disease
        is True
    )


def test_non_primary_cancer_is_preserved() -> None:
    """Other cancer labels should not be converted to prostate cancer."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "BRCA"
        )
    )

    assert (
        result.source_value
        == "BRCA"
    )

    assert (
        result.canonical_disease_id
        is None
    )

    assert (
        result.is_primary_disease
        is False
    )

    assert (
        result.context_type
        is DiseaseContextType.OTHER_DISEASE
    )

    assert (
        result.identifier_system
        is DiseaseIdentifierSystem.SOURCE_LABEL
    )

    assert (
        result.status
        is DiseaseIdentityStatus.PRESERVED
    )

    assert (
        result.usable
        is False
    )


def test_unknown_disease_is_unresolved() -> None:
    """Missing disease context should remain unresolved."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            None
        )
    )

    assert (
        result.canonical_disease_id
        is None
    )

    assert (
        result.context_type
        is DiseaseContextType.UNKNOWN
    )

    assert (
        result.status
        is DiseaseIdentityStatus.UNRESOLVED
    )

    assert (
        result.usable
        is False
    )


def test_empty_disease_is_unresolved() -> None:
    """Empty labels should remain unresolved."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "   "
        )
    )

    assert (
        result.status
        is DiseaseIdentityStatus.UNRESOLVED
    )


def test_source_value_is_preserved() -> None:
    """Original source disease representation should remain available."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "  Prostate Cancer  "
        )
    )

    assert (
        result.source_value
        == "Prostate Cancer"
    )

    assert (
        result.normalized_source_value
        == "prostate cancer"
    )


def test_same_primary_disease_matches_across_sources() -> None:
    """PRAD and Prostate cancer should have the same identity."""

    trfqtl = (
        DiseaseHarmonizer
        .harmonize(
            "PRAD"
        )
    )

    gwas = (
        DiseaseHarmonizer
        .harmonize(
            "Prostate cancer"
        )
    )

    assert (
        DiseaseHarmonizer
        .diseases_equal(
            trfqtl,
            gwas,
        )
        is True
    )


def test_moradi_and_trfqtl_primary_context_match() -> None:
    """PrCa and PRAD should resolve to the same disease."""

    moradi = (
        DiseaseHarmonizer
        .harmonize(
            "PrCa"
        )
    )

    trfqtl = (
        DiseaseHarmonizer
        .harmonize(
            "PRAD"
        )
    )

    assert (
        DiseaseHarmonizer
        .diseases_equal(
            moradi,
            trfqtl,
        )
        is True
    )


def test_non_primary_disease_does_not_match_primary() -> None:
    """Other cancers must not match prostate cancer."""

    prostate = (
        DiseaseHarmonizer
        .harmonize(
            "PRAD"
        )
    )

    breast = (
        DiseaseHarmonizer
        .harmonize(
            "BRCA"
        )
    )

    assert (
        DiseaseHarmonizer
        .diseases_equal(
            prostate,
            breast,
        )
        is False
    )


def test_both_primary_disease() -> None:
    """Different prostate cancer aliases should share primary context."""

    left = (
        DiseaseHarmonizer
        .harmonize(
            "PRAD"
        )
    )

    right = (
        DiseaseHarmonizer
        .harmonize(
            "prostate carcinoma"
        )
    )

    assert (
        DiseaseHarmonizer
        .both_primary_disease(
            left,
            right,
        )
        is True
    )


def test_disease_to_dict() -> None:
    """Disease model should be JSON serializable."""

    result = (
        DiseaseHarmonizer
        .harmonize(
            "PRAD"
        )
    )

    payload = result.to_dict()

    assert (
        payload[
            "canonical_disease_id"
        ]
        == "prostate_cancer"
    )

    assert (
        payload[
            "identifier_system"
        ]
        == "CANONICAL_INTERNAL"
    )

    assert (
        payload[
            "context_type"
        ]
        == "PRIMARY_DISEASE"
    )

    assert (
        payload[
            "status"
        ]
        == "RESOLVED"
    )