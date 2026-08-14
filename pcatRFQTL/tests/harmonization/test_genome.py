"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/harmonization/test_genome.py

Description:
    Unit tests for M3.1 genome-build verification.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import pytest

from pcatrfqtl.harmonization.genome import (
    GenomeBuildRegistry,
)
from pcatrfqtl.harmonization.models import (
    GenomeAssembly,
    GenomeBuildStatus,
)


def test_cancer_trfqtl_is_verified_hg19() -> None:
    """Cancer-tRFQTL should be registered as hg19 / GRCh37."""

    result = GenomeBuildRegistry.get(
        "cancer_trfqtl"
    )

    assert (
        result.genome_build
        == "hg19"
    )

    assert (
        result.assembly
        is GenomeAssembly.GRCH37
    )

    assert (
        result.status
        is GenomeBuildStatus.VERIFIED
    )

    assert (
        result.coordinate_join_allowed
        is True
    )


def test_moradi_qtl_is_verified_hg19() -> None:
    """Moradi QTL coordinates should be registered as hg19."""

    result = GenomeBuildRegistry.get(
        "moradi_qtl"
    )

    assert (
        result.genome_build
        == "hg19"
    )

    assert (
        result.assembly
        is GenomeAssembly.GRCH37
    )

    assert (
        result.status
        is GenomeBuildStatus.VERIFIED
    )


def test_moradi_de_not_coordinate_dataset() -> None:
    """Moradi DE outputs should not participate in variant joins."""

    result = GenomeBuildRegistry.get(
        "moradi_de"
    )

    assert (
        result.coordinate_join_allowed
        is False
    )

    assert (
        result.liftover_required
        is False
    )


def test_gwas_catalog_remains_source_dependent() -> None:
    """GWAS Catalog should not receive an assumed genome build."""

    result = GenomeBuildRegistry.get(
        "gwas_catalog"
    )

    assert (
        result.genome_build
        == "source_dependent"
    )

    assert (
        result.assembly
        is GenomeAssembly.SOURCE_DEPENDENT
    )

    assert (
        result.status
        is GenomeBuildStatus.SOURCE_DEPENDENT
    )

    assert (
        result.coordinate_join_allowed
        is False
    )


def test_trfqtl_and_moradi_are_coordinate_compatible() -> None:
    """Verified hg19 sources should permit direct coordinate matching."""

    assert (
        GenomeBuildRegistry
        .are_coordinate_compatible(
            "cancer_trfqtl",
            "moradi_qtl",
        )
        is True
    )


def test_gwas_and_moradi_not_yet_coordinate_compatible() -> None:
    """Source-dependent GWAS coordinates must not be joined yet."""

    assert (
        GenomeBuildRegistry
        .are_coordinate_compatible(
            "gwas_catalog",
            "moradi_qtl",
        )
        is False
    )


def test_hg19_sources_require_liftover_for_hg38_target() -> None:
    """Verified hg19 datasets require future conversion to hg38."""

    assert (
        GenomeBuildRegistry
        .requires_liftover_to_target(
            "cancer_trfqtl"
        )
        is True
    )

    assert (
        GenomeBuildRegistry
        .requires_liftover_to_target(
            "moradi_qtl"
        )
        is True
    )


def test_unknown_dataset_raises_key_error() -> None:
    """Unregistered datasets should fail explicitly."""

    with pytest.raises(
        KeyError
    ):
        GenomeBuildRegistry.get(
            "unknown_dataset"
        )


def test_summary_is_serializable() -> None:
    """Registry summary should contain plain serializable metadata."""

    result = (
        GenomeBuildRegistry
        .summary()
    )

    assert (
        result[
            "moradi_qtl"
        ][
            "genome_build"
        ]
        == "hg19"
    )

    assert (
        result[
            "moradi_qtl"
        ][
            "assembly"
        ]
        == "GRCh37"
    )

    assert (
        result[
            "moradi_qtl"
        ][
            "status"
        ]
        == "VERIFIED"
    )