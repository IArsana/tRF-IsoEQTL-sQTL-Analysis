"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_2b_endpoint_resolution.py

Description:
    Smoke test for M7.2B Clinical Endpoint Readiness Resolution.

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

from pcatrfqtl.analysis.m7.runners.resolve_clinical_endpoints import (
    M72BEndpointResolutionRunner,
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
    / "m7_endpoint_resolution.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "endpoint_resolution"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


# ============================================================================
# Reporting
# ============================================================================


def print_time_to_event(
    report: dict[str, Any],
) -> None:
    """Print time-to-event resolution."""

    print()
    print("Time-to-event endpoint resolution:")
    print()

    for row in report[
        "time_to_event_resolution"
    ]:

        print(
            f"  {row['endpoint_id']}"
        )

        print(
            "    Usable cases:                "
            f"{row['usable_cases']}"
        )

        print(
            "    Usable fraction:             "
            f"{row['usable_fraction']:.3f}"
        )

        print(
            "    Events:                      "
            f"{row['analyzable_events']}"
        )

        print(
            "    Censored:                    "
            f"{row['analyzable_censored_cases']}"
        )

        print(
            "    Sufficient events:           "
            f"{row['has_sufficient_events']}"
        )

        print(
            "    Sufficient completeness:     "
            f"{row['has_sufficient_completeness']}"
        )

        print(
            "    Has censoring:               "
            f"{row['has_analyzable_censoring']}"
        )

        print(
            "    Final ready:                 "
            f"{row['final_time_to_event_ready']}"
        )

        print(
            "    Resolution:                  "
            f"{row['resolution_status']}"
        )

        print(
            "    Downstream use:              "
            f"{row['downstream_use']}"
        )

        print()


def print_clinicopathologic(
    report: dict[str, Any],
) -> None:
    """Print clinicopathologic endpoint readiness."""

    print()
    print("Clinicopathologic endpoint readiness:")
    print()

    for row in report[
        "clinicopathologic_endpoint_readiness"
    ]:

        print(
            f"  {row['endpoint_id']}"
        )

        print(
            "    Source column:               "
            f"{row['source_column']}"
        )

        print(
            "    Non-missing:                 "
            f"{row['non_missing_cases']}/"
            f"{row['total_cases']}"
        )

        print(
            "    Non-missing fraction:        "
            f"{row['non_missing_fraction']:.3f}"
        )

        print(
            "    Unique values:               "
            f"{row['unique_non_missing_values']}"
        )

        print(
            "    Has variation:               "
            f"{row['has_variation']}"
        )

        print(
            "    Endpoint ready:              "
            f"{row['endpoint_ready']}"
        )

        print(
            "    Status:                      "
            f"{row['endpoint_status']}"
        )

        print()


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M7.2B smoke test."""

    print()
    print("=" * 72)
    print("M7.2B CLINICAL ENDPOINT READINESS RESOLUTION")
    print("=" * 72)

    runner = M72BEndpointResolutionRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    summary = report[
        "summary"
    ]

    final_resolution = report[
        "final_endpoint_resolution"
    ]

    print()
    print("=" * 72)
    print("M7.2B RESULT SUMMARY")
    print("=" * 72)

    print(
        "Time-to-event endpoints assessed:           "
        f"{summary['time_to_event_endpoints_assessed']}"
    )

    print(
        "Time-to-event endpoints ready:              "
        f"{summary['time_to_event_endpoints_ready']}"
    )

    print(
        "Primary time-to-event endpoint available:   "
        f"{summary['primary_time_to_event_endpoint_available']}"
    )

    print()

    print(
        "Clinicopathologic endpoints assessed:       "
        f"{summary['clinicopathologic_endpoints_assessed']}"
    )

    print(
        "Clinicopathologic endpoints ready:          "
        f"{summary['clinicopathologic_endpoints_ready']}"
    )

    print(
        "Ready clinicopathologic endpoints:"
    )

    print(
        f"  {summary['ready_clinicopathologic_endpoints']}"
    )

    print()

    print(
        "Clinical validation route continues:        "
        f"{summary['clinical_validation_route_continues']}"
    )

    print(
        "Overall resolution status:"
    )

    print(
        f"  {summary['overall_resolution_status']}"
    )

    print_time_to_event(
        report
    )

    print_clinicopathologic(
        report
    )

    print()
    print("Final endpoint resolution:")
    print()

    print(
        "  Primary time-to-event endpoint:"
    )

    print(
        "   ",
        final_resolution[
            "primary_time_to_event_endpoint"
        ],
    )

    print(
        "  Time-to-event status:"
    )

    print(
        "   ",
        final_resolution[
            "time_to_event_resolution_status"
        ],
    )

    print(
        "  Clinicopathologic association allowed:"
    )

    print(
        "   ",
        final_resolution[
            "clinicopathologic_association_allowed"
        ],
    )

    print()

    print("Next stage:")
    print(
        f"  {summary['next_stage']}"
    )

    print()

    print("Downstream if tRF quantification becomes available:")
    print(
        f"  {summary['downstream_if_trf_quantifiable']}"
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
    print("M7.2B EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()