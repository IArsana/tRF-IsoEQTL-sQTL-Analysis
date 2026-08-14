"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_harmonization.py

Description:
    Full-dataset smoke test for M3 data harmonization.

    Reads M2 standardized Parquet outputs and generates the complete M3
    harmonized dataset tree together with harmonization QC metadata.

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

from pcatrfqtl.harmonization.runners.harmonize import (
    HarmonizationInputs,
    HarmonizationRunner,
)


PROJECT_ROOT = (
    Path(
        __file__
    )
    .resolve()
    .parents[1]
)


def main() -> None:
    """Execute complete M3 harmonization."""

    runner = HarmonizationRunner(
        inputs=HarmonizationInputs(
            standardized_directory=(
                PROJECT_ROOT
                / "data"
                / "interim"
                / "standardized"
            )
        ),

        output_directory=(
            PROJECT_ROOT
            / "data"
            / "interim"
            / "harmonized"
        ),
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