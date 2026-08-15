"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m3_1_gwas.py

Description:
    Full-dataset smoke test for M3.1 prostate-cancer GWAS extraction.

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

from pcatrfqtl.analysis.m3.runners.extract_gwas import (
    M31GWASInputs,
    M31GWASRunner,
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
    """Execute M3.1 against the complete harmonized GWAS dataset."""

    runner = M31GWASRunner(
        inputs=M31GWASInputs(
            harmonized_gwas_directory=(
                PROJECT_ROOT
                / "data"
                / "interim"
                / "harmonized"
                / "gwas_catalog"
            )
        ),
        output_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m3"
        ),
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