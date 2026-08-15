"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m4_2_trfqtl_moradi_match.py

Description:
    Full-data smoke test for M4.2 exact variant integration between
    PRAD tRF-QTL associations and the Unified Moradi QTL Index.

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

from pcatrfqtl.analysis.m4.runners.match_trfqtl_moradi import (
    M42TRFQTLMoradiInputs,
    M42TRFQTLMoradiRunner,
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
    """Run the M4.2 full-data smoke test."""

    runner = M42TRFQTLMoradiRunner(
        inputs=M42TRFQTLMoradiInputs(
            prostate_trfqtl=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
                / "prostate_trfqtl.parquet"
            ),
            moradi_index_directory=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m4"
                / "moradi_qtl_index"
            ),
        ),
        output_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m4"
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