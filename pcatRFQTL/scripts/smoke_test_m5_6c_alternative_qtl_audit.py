"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_6c_alternative_qtl_audit.py

Description:
    Smoke test for M5.6C alternative PRAD QTL resource audit.

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

from pcatrfqtl.analysis.m5.runners.audit_alternative_qtl_resources import (
    M56CAlternativeQTLAuditRunner,
)


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m5_alternative_qtl_audit.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "alternative_qtl_audit"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "qc"
)


def main() -> None:
    """Run M5.6C."""

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            f"M5.6C config not found: {CONFIG_PATH}"
        )

    print()
    print("=" * 72)
    print("M5.6C ALTERNATIVE PRAD QTL RESOURCE AUDIT")
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

    runner = M56CAlternativeQTLAuditRunner(
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print()
    print("=" * 72)
    print("M5.6C RESULT SUMMARY")
    print("=" * 72)

    summary = report[
        "summary"
    ]

    print(
        f"Alternative resources assessed:      "
        f"{summary['alternative_resources_assessed']}"
    )

    print(
        f"ABF-ready resources:                 "
        f"{summary['abf_ready_resources']}"
    )

    print(
        f"SuSiE-ready resources:               "
        f"{summary['susie_ready_resources']}"
    )

    print(
        f"Resources requiring inspection:      "
        f"{summary['resources_requiring_dataset_inspection']}"
    )

    print(
        f"Candidate leads assessed:            "
        f"{summary['candidate_leads_assessed']}"
    )

    print(
        f"Formal coloc-ready candidates:       "
        f"{summary['formal_coloc_ready_candidates']}"
    )

    print(
        f"Blocked/pending candidates:          "
        f"{summary['formal_coloc_blocked_or_pending_candidates']}"
    )

    print()
    print(
        "Resource classifications:"
    )

    for (
        classification,
        count,
    ) in report[
        "resource_classification_counts"
    ].items():

        print(
            f"  {classification}: {count}"
        )

    print()
    print(
        "Candidate resolutions:"
    )

    for (
        status,
        count,
    ) in report[
        "candidate_resolution_counts"
    ].items():

        print(
            f"  {status}: {count}"
        )

    print("=" * 72)

    print()
    print(
        "Full QC report:"
    )
    print()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print()
    print("=" * 72)
    print("M5.6C EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()