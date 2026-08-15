"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m4_1_moradi_index.py

Description:
    Smoke test for M4.1 Unified Moradi QTL Index construction.

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

from pcatrfqtl.analysis.m4.runners.build_moradi_index import (
    M41MoradiIndexInputs,
    M41MoradiIndexRunner,
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
    """Run M4.1 full-data smoke test."""

    runner = M41MoradiIndexRunner(
        inputs=M41MoradiIndexInputs(
            harmonized_moradi_directory=(
                PROJECT_ROOT
                / "data"
                / "interim"
                / "harmonized"
                / "moradi_qtl"
            )
        ),
        output_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m4"
            / "moradi_qtl_index"
        ),
        qc_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m4"
            / "qc"
        ),
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