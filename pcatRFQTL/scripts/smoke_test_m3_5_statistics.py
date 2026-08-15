"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m3_5_statistics.py

Description:
    Full-dataset smoke test for M3.5 descriptive statistics.

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

from pcatrfqtl.analysis.m3.runners.summarize_statistics import (
    M35StatisticsInputs,
    M35StatisticsRunner,
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

    m3_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m3"
    )

    runner = M35StatisticsRunner(
        inputs=M35StatisticsInputs(
            prostate_gwas=(
                m3_directory
                / "prostate_gwas.parquet"
            ),
            prostate_trfqtl=(
                m3_directory
                / "prostate_trfqtl.parquet"
            ),
            candidate_snptrf=(
                m3_directory
                / "candidate_snptrf.parquet"
            ),
        ),
        output_directory=(
            m3_directory
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