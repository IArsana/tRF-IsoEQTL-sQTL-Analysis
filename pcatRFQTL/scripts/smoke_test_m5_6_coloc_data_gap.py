"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_6_coloc_data_gap.py

Description:
    Smoke test for M5.6 colocalization data-gap resolution.

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

from pcatrfqtl.analysis.m5.runners.resolve_coloc_data_gap import (
    M56ColocDataGapRunner,
)


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m5_coloc_data_gap.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "coloc_data_gap"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "qc"
)


def main() -> None:
    """Run M5.6."""

    print()
    print("=" * 72)
    print("M5.6 COLOCALIZATION DATA-GAP RESOLUTION")
    print("=" * 72)

    print(
        f"Project root:       {ROOT}"
    )

    print(
        f"Config:             {CONFIG_PATH}"
    )

    print(
        f"Output directory:   {OUTPUT_DIRECTORY}"
    )

    print(
        f"QC directory:       {QC_DIRECTORY}"
    )

    print("=" * 72)

    runner = M56ColocDataGapRunner(
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print()
    print("=" * 72)
    print("M5.6 RESULT")
    print("=" * 72)

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print()
    print("=" * 72)
    print("M5.6 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()