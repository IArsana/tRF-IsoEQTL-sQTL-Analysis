"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_standardization.py

Description:
    Full-dataset smoke test for M2 data standardization.

    The script standardizes all validated source datasets and writes
    standardized Parquet outputs together with a JSON execution summary.

    Source datasets:

        - GWAS Catalog
        - Cancer-tRFQTL
        - Moradi QTL S1-S4
        - Moradi DE S5-S10

    This script is intended for local pipeline verification.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import json
from pathlib import Path

from pcatrfqtl.standardization.runners.standardize import (
    StandardizationInputs,
    StandardizationRunner,
)


PROJECT_ROOT = (
    Path(
        __file__
    )
    .resolve()
    .parents[1]
)


def main() -> None:
    """Execute complete M2 standardization."""

    inputs = StandardizationInputs(
        gwas_catalog=(
            PROJECT_ROOT
            / "data"
            / "raw"
            / "gwas"
            / "gwas-catalog-download-associations-alt-full.tsv"
        ),

        cancer_trfqtl_workbook=(
            PROJECT_ROOT
            / "data"
            / "raw"
            / "trfqtl"
            / "can-25-1282_supplementary_data_suppst1-11.xlsx"
        ),

        moradi_directory=(
            PROJECT_ROOT
            / "data"
            / "raw"
            / "moradi"
        ),
    )

    runner = StandardizationRunner(
        inputs=inputs,

        output_directory=(
            PROJECT_ROOT
            / "data"
            / "interim"
            / "standardized"
        ),

        gwas_chunk_size=100_000,
    )

    report = runner.run_all()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()