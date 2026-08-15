"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_2_prioritize_gwas_sumstats.py

Description:
    Full-data smoke test for M5.3C.2 prostate cancer GWAS
    summary-statistics prioritization.

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

from pcatrfqtl.analysis.m5.runners.prioritize_gwas_sumstats import (
    M53CProstateGWASPrioritizationRunner,
)


# ============================================================================
# Project
# ============================================================================

PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Run M5.3C.2 prioritization smoke test."""

    discovery_path = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "sumstats"
        / "prostate_gwas_sumstats_studies.parquet"
    )

    output_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "sumstats"
    )

    qc_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "qc"
    )

    runner = (
        M53CProstateGWASPrioritizationRunner(
            discovery_path=(
                discovery_path
            ),
            output_directory=(
                output_directory
            ),
            qc_directory=(
                qc_directory
            ),
        )
    )

    report = runner.run()

    print(
        "\n"
        + "=" * 72
    )

    print(
        "M5.3C.2 PRIORITIZATION SMOKE TEST"
    )

    print(
        "=" * 72
    )

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print(
        "=" * 72
    )


if __name__ == "__main__":
    main()