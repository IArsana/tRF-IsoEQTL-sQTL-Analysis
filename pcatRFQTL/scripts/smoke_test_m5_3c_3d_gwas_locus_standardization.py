"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3d_gwas_locus_standardization.py

Description:
    Smoke test for M5.3C.3D GWAS locus standardization and disease-signal QC.

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

from pcatrfqtl.analysis.m5.runners.standardize_gwas_loci import (
    M53C3DGWASLocusStandardizationRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


def main() -> None:
    """Run M5.3C.3D smoke test."""

    runner = (
        M53C3DGWASLocusStandardizationRunner(
            retrieval_manifest_path=(
                ROOT
                / "data/processed/m5/sumstats/"
                "gwas_locus_retrieval_manifest.parquet"
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
    )

    report = runner.run()

    print()
    print(
        "="
        * 72
    )
    print(
        "M5.3C.3D GWAS LOCUS STANDARDIZATION & DISEASE-SIGNAL QC"
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