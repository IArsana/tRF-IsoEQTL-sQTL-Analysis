"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m3_4_candidates.py

Description:
    Full-dataset smoke test for M3.4 integrated SNP-tRF candidate
    generation.

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

from pcatrfqtl.analysis.m3.runners.build_candidates import (
    M34CandidateInputs,
    M34CandidateRunner,
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
    """Execute M3.4 against M3.1 and M3.2 outputs."""

    m3_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m3"
    )

    runner = M34CandidateRunner(
        inputs=M34CandidateInputs(
            prostate_gwas=(
                m3_directory
                / "prostate_gwas.parquet"
            ),
            prostate_trfqtl=(
                m3_directory
                / "prostate_trfqtl.parquet"
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