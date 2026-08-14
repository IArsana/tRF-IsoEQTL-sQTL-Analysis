"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/test_config.py

Description:
    Unit tests for project dataset configuration loading.

    configs/datasets.yaml is the single source of truth for dataset
    metadata, raw data paths, integration settings, and genomic
    harmonization metadata.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pcatrfqtl.config.settings import load_datasets


def test_datasets_loads() -> None:
    """Dataset configuration should load the expected root structure."""

    config = load_datasets()

    assert isinstance(
        config,
        dict,
    )

    assert (
        "project"
        in config
    )

    assert (
        "directories"
        in config
    )

    assert (
        "datasets"
        in config
    )

    assert (
        "integration"
        in config
    )


def test_expected_datasets_exist() -> None:
    """Core source datasets should be registered."""

    config = load_datasets()

    datasets = config[
        "datasets"
    ]

    assert (
        "gwas_catalog"
        in datasets
    )

    assert (
        "cancer_trfqtl"
        in datasets
    )

    assert (
        "moradi"
        in datasets
    )


def test_gwas_catalog_configuration() -> None:
    """GWAS Catalog should have the expected source registration."""

    config = load_datasets()

    gwas = config[
        "datasets"
    ][
        "gwas_catalog"
    ]

    assert (
        gwas[
            "id"
        ]
        == "gwas_catalog_associations"
    )

    assert (
        gwas[
            "source"
        ][
            "database"
        ]
        == "GWAS Catalog"
    )


def test_primary_disease_is_prostate_cancer() -> None:
    """Integration configuration should target prostate cancer."""

    config = load_datasets()

    primary_disease = (
        config[
            "integration"
        ][
            "primary_disease"
        ]
    )

    assert (
        primary_disease[
            "name"
        ]
        == "Prostate cancer"
    )

    assert (
        primary_disease[
            "tcga_code"
        ]
        == "PRAD"
    )


def test_genome_build_metadata() -> None:
    """Genome-build metadata should remain explicit and structured."""

    config = load_datasets()

    reference_build = (
        config[
            "integration"
        ][
            "genomic_harmonization"
        ][
            "reference_build"
        ]
    )

    trfqtl = (
        reference_build[
            "trfqtl"
        ]
    )

    moradi = (
        reference_build[
            "moradi"
        ]
    )

    gwas = (
        reference_build[
            "gwas_catalog"
        ]
    )

    assert (
        trfqtl[
            "build"
        ]
        == "hg19"
    )

    assert (
        trfqtl[
            "assembly"
        ]
        == "GRCh37"
    )

    assert (
        trfqtl[
            "status"
        ]
        == "verified"
    )

    assert (
        moradi[
            "build"
        ]
        == "hg19"
    )

    assert (
        moradi[
            "assembly"
        ]
        == "GRCh37"
    )

    assert (
        moradi[
            "status"
        ]
        == "verified"
    )

    assert (
        gwas[
            "build"
        ]
        == "source_dependent"
    )

    assert (
        gwas[
            "assembly"
        ]
        == "source_dependent"
    )

    assert (
        gwas[
            "status"
        ]
        == "source_dependent"
    )

def test_harmonization_target_build() -> None:
    """Canonical harmonization target should be hg38."""

    config = load_datasets()

    harmonization = (
        config[
            "integration"
        ][
            "genomic_harmonization"
        ]
    )

    assert (
        harmonization[
            "target_build"
        ][
            "build"
        ]
        == "hg38"
    )

    assert (
        harmonization[
            "target_build"
        ][
            "assembly"
        ]
        == "GRCh38"
    )

    assert (
        harmonization[
            "liftover"
        ][
            "enabled"
        ]
        is False
    )