"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_trfqtl_all_standardization.py

Description:
    Full-workbook smoke test for complete Cancer-tRFQTL supplementary
    standardization.

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

from pcatrfqtl.standardization.runners.standardize_trfqtl_all import (
    CancerTRFQTLAllStandardizationRunner,
    TRFQTLAllInputs,
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
    """Run complete Cancer-tRFQTL standardization smoke test."""

    runner = (
        CancerTRFQTLAllStandardizationRunner(
            inputs=TRFQTLAllInputs(
                workbook=(
                    PROJECT_ROOT
                    / "data"
                    / "raw"
                    / "trfqtl"
                    / "can-25-1282_supplementary_data_suppst1-11.xlsx"
                )
            ),
            output_directory=(
                PROJECT_ROOT
                / "data"
                / "interim"
                / "standardized"
                / "cancer_trfqtl"
            ),
            qc_directory=(
                PROJECT_ROOT
                / "data"
                / "interim"
                / "standardized"
                / "qc"
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