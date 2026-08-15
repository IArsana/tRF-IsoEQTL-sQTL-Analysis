"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_1_m5_2_candidate_loci.py

Description:
    Full-data smoke test for:

        M5.1 - LD Reference Strategy
        M5.2 - Candidate Locus Construction

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

from pcatrfqtl.analysis.m5.runners.build_candidate_loci import (
    M51M52CandidateLocusRunner,
    M51M52Inputs,
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
    """Run M5.1 + M5.2 smoke test."""

    runner = M51M52CandidateLocusRunner(
        inputs=M51M52Inputs(
            bridge_readiness=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m4"
                / "m4_regulatory_bridge_readiness.parquet"
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
            / "m5"
        ),
        qc_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m5"
            / "qc"
        ),
        screening_window_bp=500_000,
        primary_r2_threshold=0.8,
        secondary_r2_threshold=0.5,
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