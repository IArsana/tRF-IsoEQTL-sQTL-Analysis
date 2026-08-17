"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_1_validation_readiness.py

Description:
    Smoke test for M7.1 TCGA-PRAD Clinical Validation Readiness Audit.

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
from typing import Any

from pcatrfqtl.analysis.m7.runners.assess_validation_readiness import (
    M71ValidationReadinessRunner,
)


# ============================================================================
# Paths
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[
    1
]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m7_validation_readiness.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "validation_readiness"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


# ============================================================================
# Reporting helpers
# ============================================================================


def print_resource_readiness(
    report: dict[str, Any],
) -> None:
    """Print resource-level validation routes."""

    print()
    print("Validation resource readiness:")
    print()

    for row in report[
        "validation_resource_readiness"
    ]:

        print(
            f"  {row['resource_id']}"
        )

        print(
            "    Domain:                      "
            f"{row['domain']}"
        )

        print(
            "    Public resource available:   "
            f"{row['public_resource_available']}"
        )

        print(
            "    Controlled access required:  "
            f"{row['controlled_access_required']}"
        )

        print(
            "    Direct candidate measure:     "
            f"{row['direct_candidate_measurement']}"
        )

        print(
            "    Capability status:            "
            f"{row['capability_status']}"
        )

        print(
            "    Next action:                  "
            f"{row['next_action']}"
        )

        print()


def print_candidate_readiness(
    report: dict[str, Any],
) -> None:
    """Print candidate-specific M7 readiness."""

    print()
    print("Candidate validation readiness:")
    print()

    for row in report[
        "candidate_validation_readiness"
    ]:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    tRF IDs:                       "
            f"{row['trf_ids']}"
        )

        print(
            "    tRF identity available:        "
            f"{row['trf_identity_available']}"
        )

        print(
            "    tRF sequence resolved:          "
            f"{row['trf_sequence_resolved']}"
        )

        print(
            "    tRF validation status:          "
            f"{row['trf_validation_status']}"
        )

        print(
            "    Direct tRF quantification:      "
            f"{row['direct_trf_quantification_ready']}"
        )

        print(
            "    Variant validation status:      "
            f"{row['variant_validation_status']}"
        )

        print(
            "    Direct variant validation:      "
            f"{row['direct_variant_validation_ready']}"
        )

        print(
            "    Regulatory feature:             "
            f"{row['regulatory_feature_id']}"
        )

        print(
            "    Regulatory-event status:        "
            f"{row['regulatory_event_validation_status']}"
        )

        print(
            "    Clinical metadata route ready:  "
            f"{row['clinical_metadata_route_ready']}"
        )

        print()


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M7.1 smoke test."""

    print()
    print("=" * 72)
    print("M7.1 TCGA-PRAD CLINICAL VALIDATION READINESS AUDIT")
    print("=" * 72)

    runner = M71ValidationReadinessRunner(
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
    print("=" * 72)
    print("M7.1 RESULT SUMMARY")
    print("=" * 72)

    print(
        "Target project:                              "
        f"{summary['target_project']}"
    )

    print(
        "Candidates assessed:                         "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Candidates with tRF identity:                "
        f"{summary['candidates_with_trf_identity']}"
    )

    print(
        "Candidates with resolved tRF sequence:       "
        f"{summary['candidates_with_resolved_trf_sequence']}"
    )

    print(
        "Candidates ready for tRF quantification:     "
        f"{summary['candidates_ready_for_direct_trf_quantification']}"
    )

    print(
        "Candidates ready for variant validation:     "
        f"{summary['candidates_ready_for_direct_variant_validation']}"
    )

    print(
        "Candidates with regulatory feature:          "
        f"{summary['candidates_with_regulatory_feature']}"
    )

    print(
        "Regulatory events ready for validation:      "
        f"{summary['regulatory_events_ready_for_direct_validation']}"
    )

    print()

    print(
        "Clinical metadata route ready:               "
        f"{summary['clinical_metadata_route_ready']}"
    )

    print(
        "Controlled-access domains:                   "
        f"{summary['controlled_access_domains']}"
    )

    print()

    print(
        "Overall readiness status:"
    )

    print(
        f"  {summary['overall_readiness_status']}"
    )

    print_resource_readiness(
        report
    )

    print_candidate_readiness(
        report
    )

    print()
    print("Next routes:")
    print()

    print(
        "  Clinical:"
        f"   {summary['clinical_next_stage']}"
    )

    print(
        "  tRF:"
        f"        {summary['trf_next_stage']}"
    )

    print(
        "  Variant:"
        f"    {summary['variant_next_stage']}"
    )

    print(
        "  Regulatory:"
        f" {summary['regulatory_event_next_stage']}"
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
    print("M7.1 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()