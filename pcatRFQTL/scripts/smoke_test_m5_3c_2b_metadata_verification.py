"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_2b_metadata_verification.py

Description:
    Full-data smoke test for M5.3C.2B official GWAS Catalog metadata
    verification.

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

from pcatrfqtl.analysis.m5.runners.verify_gwas_metadata import (
    M53CMetadataVerificationRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def main() -> None:
    """Run M5.3C.2B."""

    runner = M53CMetadataVerificationRunner(
        prioritized_path=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m5"
            / "sumstats"
            / "prostate_gwas_sumstats_prioritized.parquet"
        ),
        studies_metadata_path=(
            PROJECT_ROOT
            / "data"
            / "raw"
            / "gwas"
            / "metadata"
            / "gwas_catalog_all_studies.tsv"
        ),
        ancestry_metadata_path=(
            PROJECT_ROOT
            / "data"
            / "raw"
            / "gwas"
            / "metadata"
            / "gwas_catalog_all_ancestry.tsv"
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
        "M5.3C.2B METADATA VERIFICATION"
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