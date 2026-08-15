"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3_m5_4_ld_bridges.py

Description:
    Full-data smoke test for:

        M5.3 - LD-aware GWAS bridge
        M5.4 - LD-aware Moradi regulatory bridge

    Default exploratory population:
        EAS

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pcatrfqtl.analysis.m5.runners.build_ld_bridges import (
    M53M54Inputs,
    M53M54LDBridgeRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def parse_arguments() -> argparse.Namespace:
    """Parse smoke-test arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Run M5.3 + M5.4 LD-aware bridge analysis."
        )
    )

    parser.add_argument(
        "--population",
        default="EAS",
        choices=[
            "EAS",
            "EUR",
            "SAS",
        ],
        help=(
            "Population-specific LD reference to use. "
            "Default: EAS."
        ),
    )

    return parser.parse_args()


def main() -> None:
    """Run M5.3 + M5.4 smoke test."""

    args = parse_arguments()

    runner = M53M54LDBridgeRunner(
        inputs=M53M54Inputs(
            candidate_loci=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m5"
                / "candidate_loci.parquet"
            ),
            ld_root_directory=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m5"
                / "ld"
            ),
            prostate_gwas=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
                / "prostate_gwas.parquet"
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
            / "m5"
        ),
        qc_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m5"
            / "qc"
        ),
        population=args.population,
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