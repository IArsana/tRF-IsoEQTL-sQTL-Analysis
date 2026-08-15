"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3a_remote_sumstats_discovery.py

Description:
    Full-data smoke test for M5.3C.3A remote GWAS summary-statistics
    resource discovery.

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

from pcatrfqtl.analysis.m5.runners.discover_locus_sumstats import (
    M53C3RemoteDiscoveryRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def main() -> None:
    """Run M5.3C.3A smoke test."""

    runner = M53C3RemoteDiscoveryRunner(
        config_path=(
            PROJECT_ROOT
            / "configs"
            / "m5_gwas_locus_studies.yaml"
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
        "M5.3C.3A REMOTE SUMSTATS DISCOVERY"
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