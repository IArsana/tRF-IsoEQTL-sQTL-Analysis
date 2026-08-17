"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_3c_controlled_access_readiness.py

Description:
    Smoke test for M7.3C Controlled-Access Quantification Readiness.

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

from pcatrfqtl.analysis.m7.runners.assess_controlled_access_readiness import (
    M73CControlledAccessReadinessRunner,
)


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m7_controlled_access_readiness.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "trf_controlled_access_readiness"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


def main() -> None:

    print()
    print("=" * 72)
    print("M7.3C CONTROLLED-ACCESS QUANTIFICATION READINESS")
    print("=" * 72)

    runner = M73CControlledAccessReadinessRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    summary = report[
        "summary"
    ]

    print()
    print("Controlled BAM files assessed:            "
          f"{summary['controlled_bam_files_assessed']}")

    print("Controlled BAM files resource-ready:      "
          f"{summary['controlled_bam_files_resource_ready']}")

    print("Candidate sequences ready:                "
          f"{summary['candidate_sequences_ready']}")

    print()
    print("Token path configured:                    "
          f"{summary['token_path_configured']}")

    print("Token file exists:                        "
          f"{summary['token_file_exists']}")

    print("Token file readable:                      "
          f"{summary['token_file_readable']}")

    print()
    print("Access configurations ready:              "
          f"{summary['candidate_access_configurations_ready']}")

    print("Candidates ready for controlled download: "
          f"{summary['candidates_ready_for_controlled_download']}")

    print()
    print("Remote authorization tested:              "
          f"{summary['remote_authorization_tested']}")

    print("BAM download performed:                   "
          f"{summary['bam_download_performed']}")

    print("Read counting performed:                  "
          f"{summary['read_counting_performed']}")

    print("tRF quantification performed:             "
          f"{summary['trf_quantification_performed']}")

    print()
    print("Overall readiness:")
    print(
        f"  {summary['overall_readiness_status']}"
    )

    print()
    print("Candidate readiness:")

    for row in report[
        "candidate_controlled_access_readiness"
    ]:

        print()
        print(
            f"  {row['lead_rsid']} | {row['trf_id']}"
        )

        print(
            "    sequence ready:             "
            f"{row['sequence_ready_for_quantification']}"
        )

        print(
            "    controlled BAMs:            "
            f"{row['controlled_aligned_bam_count']}"
        )

        print(
            "    access configuration ready: "
            f"{row['access_configuration_ready']}"
        )

        print(
            "    download ready:             "
            f"{row['ready_for_controlled_download']}"
        )

        print(
            "    status:                     "
            f"{row['controlled_access_readiness_status']}"
        )

    print()
    print("Next stage:")
    print(
        f"  {summary['next_stage']}"
    )

    print()
    print("Full QC report:")
    print()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
    )

    print()
    print("=" * 72)
    print("M7.3C EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()