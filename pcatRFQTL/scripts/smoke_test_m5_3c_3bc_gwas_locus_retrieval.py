"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3bc_gwas_locus_retrieval.py

Description:
    Full-data smoke test for M5.3C.3B-C harmonised lead coordinate
    resolution and locus-level GWAS summary-statistics retrieval.

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

from pcatrfqtl.analysis.m5.runners.retrieve_gwas_loci import (
    M53C3LocusRetrievalRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def main() -> None:
    """Run M5.3C.3B-C smoke test."""

    runner = M53C3LocusRetrievalRunner(
        config_path=(
            PROJECT_ROOT
            / "configs"
            / "m5_gwas_locus_studies.yaml"
        ),
        remote_manifest_path=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m5"
            / "sumstats"
            / "gwas_locus_remote_manifest.parquet"
        ),
        output_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m5"
            / "sumstats"
        ),
        qc_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m5"
            / "qc"
        ),
    )

    report = runner.run()

    print()
    print(
        "=" * 72
    )

    print(
        "M5.3C.3B-C GWAS LOCUS RETRIEVAL"
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