"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3f_disease_ld_bridge.py

Description:
    Smoke test for M5.3C.3F disease-signal ↔ LD bridge integration.

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

from pcatrfqtl.analysis.m5.runners.integrate_disease_ld_bridges import (
    M53C3FDiseaseLDBridgeRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


def main() -> None:
    """Run M5.3C.3F smoke test."""

    runner = M53C3FDiseaseLDBridgeRunner(
        recurrence_path=(
            ROOT
            / "data/processed/m5/sumstats/"
            "gwas_cross_study_variant_recurrence.parquet"
        ),
        m5_directory=(
            ROOT
            / "data/processed/m5"
        ),
        qc_directory=(
            ROOT
            / "data/processed/m5/qc"
        ),
    )

    report = runner.run()

    print()
    print(
        "="
        * 72
    )

    print(
        "M5.3C.3F DISEASE-SIGNAL ↔ LD BRIDGE INTEGRATION"
    )

    print(
        "="
        * 72
    )

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print(
        "="
        * 72
    )


if __name__ == "__main__":
    main()