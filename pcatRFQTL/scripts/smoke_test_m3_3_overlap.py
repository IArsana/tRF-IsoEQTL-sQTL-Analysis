"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m3_3_overlap.py

Description:
    Full-dataset smoke test for M3.3 direct GWAS × tRF-QTL overlap.

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

from pcatrfqtl.analysis.m3.runners.direct_overlap import (
    M33OverlapInputs,
    M33OverlapRunner,
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
    """Execute M3.3 against M3.1 and M3.2 outputs."""

    runner = M33OverlapRunner(
        inputs=M33OverlapInputs(
            prostate_gwas=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
                / "prostate_gwas.parquet"
            ),
            prostate_trfqtl=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
                / "prostate_trfqtl.parquet"
            ),
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