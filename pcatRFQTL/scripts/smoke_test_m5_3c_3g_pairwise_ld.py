"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3g_pairwise_ld.py

Description:
    Full smoke/execution script for M5.3C.3G disease-proxy ↔
    regulatory-proxy pairwise LD.

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

from pcatrfqtl.analysis.m5.runners.query_disease_regulatory_pairwise_ld import (
    M53C3GPairwiseLDRunner,
)


ROOT = Path(
    __file__
).resolve().parents[1]


def main() -> None:
    """Run full M5.3C.3G."""

    runner = M53C3GPairwiseLDRunner(
        scaffold_path=(
            ROOT
            / "data/processed/m5/pairwise_ld/"
            "disease_regulatory_pairwise_candidates.parquet"
        ),
        cache_path=(
            ROOT
            / "data/interim/m5/pairwise_ld/"
            "pairwise_ld_cache.sqlite"
        ),
        output_directory=(
            ROOT
            / "data/processed/m5/pairwise_ld"
        ),
        qc_directory=(
            ROOT
            / "data/processed/m5/qc"
        ),
        env_path=(
            ROOT
            / ".env"
        ),
    )

    report = runner.run()

    print()
    print("=" * 72)
    print("M5.3C.3G DISEASE ↔ REGULATORY PAIRWISE LD")
    print("=" * 72)

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print("=" * 72)


if __name__ == "__main__":
    main()