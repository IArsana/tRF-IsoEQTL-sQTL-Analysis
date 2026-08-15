"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m4_3_evidence_summary.py

Description:
    Smoke test for M4.3 Multi-omic Exact-Evidence Summary.

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

from pcatrfqtl.analysis.m4.runners.summarize_m4_evidence import (
    M43EvidenceSummaryInputs,
    M43EvidenceSummaryRunner,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def main() -> None:

    candidate_scaffold = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m3"
        / "candidate_snptrf.parquet"
    )

    runner = M43EvidenceSummaryRunner(
        inputs=M43EvidenceSummaryInputs(
            prostate_trfqtl=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m3"
                / "prostate_trfqtl.parquet"
            ),
            moradi_matches=(
                PROJECT_ROOT
                / "data"
                / "processed"
                / "m4"
                / "trfqtl_moradi_matches.parquet"
            ),
            m3_candidate_scaffold=(
                candidate_scaffold
                if candidate_scaffold.exists()
                else None
            ),
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