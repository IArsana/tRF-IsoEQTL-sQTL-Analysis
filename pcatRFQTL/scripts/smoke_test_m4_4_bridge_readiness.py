"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m4_4_bridge_readiness.py

Description:
    Smoke test for M4.4 Regulatory Bridge Readiness.

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

from pcatrfqtl.analysis.m4.runners.assess_bridge_readiness import (
    M44BridgeReadinessInputs,
    M44BridgeReadinessRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def main() -> None:

    runner = M44BridgeReadinessRunner(
        inputs=M44BridgeReadinessInputs(
            evidence_summary=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m4"
                / "m4_candidate_evidence_summary.parquet"
            )
        ),
        output_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m4"
        ),
        qc_directory=(
            PROJECT_ROOT
            / "data"
            / "processed"
            / "m4"
            / "qc"
        ),
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