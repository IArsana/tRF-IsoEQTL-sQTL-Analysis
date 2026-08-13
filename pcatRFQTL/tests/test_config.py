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
    """Known and unresolved genome builds should remain explicit."""

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

    assert (
        reference_build[
            "trfqtl"
        ]
        == "hg19"
    )

    assert (
        reference_build[
            "moradi"
        ]
        == "to_be_verified"
    )