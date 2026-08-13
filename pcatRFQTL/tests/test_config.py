from pcatrfqtl.config.settings import (
    load_config,
    load_datasets,
    load_thresholds,
)


def test_config_loads():
    config = load_config()

    assert "project" in config
    assert "paths" in config


def test_datasets_loads():
    datasets = load_datasets()

    assert "gwas" in datasets
    assert "trfqtl" in datasets


def test_thresholds_loads():
    thresholds = load_thresholds()

    assert "gwas" in thresholds
    assert "qtl" in thresholds
