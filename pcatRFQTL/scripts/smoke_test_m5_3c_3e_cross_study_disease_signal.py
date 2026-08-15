"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3e_cross_study_disease_signal.py

Description:
    Smoke test for M5.3C.3E cross-study GWAS disease-signal synthesis.

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

from pcatrfqtl.analysis.m5.runners.synthesize_gwas_disease_signal import (
    M53C3ECrossStudyDiseaseSignalRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


def main() -> None:
    """Run M5.3C.3E smoke test."""

    runner = M53C3ECrossStudyDiseaseSignalRunner(
        standardized_path=(
            ROOT
            / "data/processed/m5/sumstats/standardized/"
            "gwas_locus_standardized.parquet"
        ),
        disease_signal_path=(
            ROOT
            / "data/processed/m5/sumstats/"
            "gwas_locus_disease_signal_summary.parquet"
        ),
        output_directory=(
            ROOT
            / "data/processed/m5/sumstats"
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
        "M5.3C.3E CROSS-STUDY DISEASE-SIGNAL SYNTHESIS"
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