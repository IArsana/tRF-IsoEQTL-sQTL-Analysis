"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3b_whole_catalog_proxy_audit.py

Description:
    Smoke test for M5.3B whole-GWAS-Catalog LD proxy audit.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pcatrfqtl.analysis.m5.runners.audit_whole_catalog_proxies import (
    M53BWholeCatalogProxyAuditRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def parse_arguments() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description=(
            "Run M5.3B whole-Catalog proxy audit."
        )
    )

    parser.add_argument(
        "--population",
        default="EAS",
        choices=[
            "EAS",
            "EUR",
            "SAS",
        ],
    )

    parser.add_argument(
        "--gwas-directory",
        required=True,
        help=(
            "Directory containing all standardized GWAS "
            "Catalog Parquet parts."
        ),
    )

    return parser.parse_args()


def main() -> None:

    args = parse_arguments()

    population = (
        args.population
        .strip()
        .upper()
    )

    runner = (
        M53BWholeCatalogProxyAuditRunner(
            normalized_ld_path=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m5"
                / "ld"
                / population
                / "normalized_ld_evidence.parquet"
            ),
            gwas_directory=(
                Path(
                    args.gwas_directory
                )
            ),
            output_directory=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m5"
            ),
            qc_directory=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m5"
                / "qc"
            ),
            population=population,
            secondary_r2_threshold=0.5,
        )
    )

    report = runner.run()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()