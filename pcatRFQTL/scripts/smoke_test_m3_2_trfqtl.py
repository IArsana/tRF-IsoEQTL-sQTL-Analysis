"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m3_2_trfqtl.py

Description:
    Full-dataset smoke test for M3.2 prostate-cancer tRF-QTL
    extraction.

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

from pcatrfqtl.analysis.m3.runners.extract_trfqtl import (
    M32TRFQTLInputs,
    M32TRFQTLRunner,
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
    """Execute M3.2 against harmonized Cancer-tRFQTL S2."""

    runner = (
        M32TRFQTLRunner(
            inputs=M32TRFQTLInputs(
                harmonized_trfqtl_s2=(
                    PROJECT_ROOT
                    / "data"
                    / "interim"
                    / "harmonized"
                    / "cancer_trfqtl"
                    / "s2.parquet"
                )
            ),
            output_directory=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
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