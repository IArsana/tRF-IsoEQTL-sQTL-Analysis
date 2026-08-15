"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_trfqtl_all_harmonization.py

Description:
    Full-dataset smoke test for complete Cancer-tRFQTL supplementary
    harmonization.

    The script backfills harmonization for S1, S3-S9, and S11 while
    auditing the locked S2 and S10 harmonized tables.

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

from pcatrfqtl.harmonization.runners.harmonize_trfqtl_all import (
    CancerTRFQTLAllHarmonizationRunner,
    TRFQTLAllHarmonizationInputs,
)


PROJECT_ROOT = (
    Path(
        __file__
    )
    .resolve()
    .parents[
        1
    ]
)


def main() -> None:
    """Run complete Cancer-tRFQTL harmonization smoke test."""

    standardized_directory = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "standardized"
        / "cancer_trfqtl"
    )

    harmonized_directory = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "harmonized"
        / "cancer_trfqtl"
    )

    qc_directory = (
        PROJECT_ROOT
        / "data"
        / "interim"
        / "harmonized"
        / "qc"
    )

    runner = (
        CancerTRFQTLAllHarmonizationRunner(
            inputs=TRFQTLAllHarmonizationInputs(
                standardized_directory=(
                    standardized_directory
                )
            ),
            output_directory=(
                harmonized_directory
            ),
            qc_directory=(
                qc_directory
            ),
        )
    )

    report = (
        runner.run()
    )

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()