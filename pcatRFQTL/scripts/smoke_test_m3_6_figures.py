"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m3_6_figures.py

Description:
    Full-dataset smoke test for M3.6 figure generation.

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

from pcatrfqtl.analysis.m3.runners.build_figures import (
    M36FigureInputs,
    M36FigureRunner,
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

    runner = M36FigureRunner(
        inputs=M36FigureInputs(
            statistics_summary=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
                / "qc"
                / "m3_5_statistics_summary.json"
            )
        ),
        figure_directory=(
            PROJECT_ROOT
            / "results"
            / "figures"
            / "m3"
        ),
        qc_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m3"
            / "qc"
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