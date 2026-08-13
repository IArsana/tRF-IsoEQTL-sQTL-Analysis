"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/config/settings.py

Description:
    Provides configuration loading, validation, and project-level
    settings management for the PCa-tRFQTL research pipeline.

    This module is responsible for loading YAML-based configuration
    files, resolving project paths, and providing centralized access
    to dataset configuration, analysis thresholds, and general
    pipeline settings.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CONFIG_DIR = PROJECT_ROOT / "configs"


def load_yaml(filename: str) -> dict:
    path = CONFIG_DIR / filename

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    with path.open("r", encoding="utf-8") as file:
        return yaml.safe_load(file)


def load_config() -> dict:
    return load_yaml("config.yaml")


def load_datasets() -> dict:
    return load_yaml("datasets.yaml")


def load_thresholds() -> dict:
    return load_yaml("thresholds.yaml")
