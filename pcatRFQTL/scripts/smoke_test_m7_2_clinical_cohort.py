"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_2_clinical_cohort.py

Description:
    Smoke test for M7.2 TCGA-PRAD Clinical Cohort Construction.

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

from pcatrfqtl.analysis.m7.runners.build_clinical_cohort import (
    M72ClinicalCohortRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m7_clinical_cohort.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "clinical_cohort"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


def print_endpoint_readiness(
    report: dict[str, Any],
) -> None:
    """Print endpoint-level readiness."""

    print()
    print("Clinical endpoint readiness:")
    print()

    for row in report[
        "endpoint_readiness"
    ]:

        print(
            f"  {row['endpoint_id']}"
        )

        print(
            "    Total cases:                 "
            f"{row['total_cases']}"
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
            f"{row['events']}"
        )

        print(
            "    Endpoint ready:              "
            f"{row['endpoint_ready']}"
        )

        print(
            "    Status:                      "
            f"{row['endpoint_status']}"
        )

        print(
            "    Primary endpoint selected:   "
            f"{row['primary_endpoint_selected']}"
        )

        print()


def print_selected_completeness(
    report: dict[str, Any],
) -> None:
    """Print important variable completeness."""

    important_fields = {
        "age_at_diagnosis_years",
        "vital_status",
        "tumor_grade",
        "ajcc_pathologic_stage",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "progression_or_recurrence",
        "days_to_recurrence",
        "max_follow_up_days",
    }

    print()
    print("Selected clinical variable completeness:")
    print()

    for row in report[
        "clinical_variable_completeness"
    ]:

        if row[
            "field"
        ] not in important_fields:
            continue

        print(
            f"  {row['field']:<32} "
            f"{row['non_missing_cases']:>4}/"
            f"{row['total_cases']:<4} "
            f"({row['non_missing_fraction']:.3f})"
        )


def main() -> None:
    """Execute M7.2 smoke test."""

    print()
    print("=" * 72)
    print("M7.2 TCGA-PRAD CLINICAL COHORT CONSTRUCTION")
    print("=" * 72)

    runner = M72ClinicalCohortRunner(
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
    print("M7.2 RESULT SUMMARY")
    print("=" * 72)

    print(
        "Project:                              "
        f"{summary['project_id']}"
    )

    print(
        "Cases acquired:                       "
        f"{summary['cases_acquired']}"
    )

    print(
        "Unique cases:                         "
        f"{summary['unique_cases']}"
    )

    print(
        "Clinically usable cases:              "
        f"{summary['clinically_usable_cases']}"
    )

    print()

    print(
        "OS usable cases:                      "
        f"{summary['os_usable_cases']}"
    )

    print(
        "OS events:                            "
        f"{summary['os_events']}"
    )

    print(
        "OS endpoint ready:                    "
        f"{summary['os_endpoint_ready']}"
    )

    print()

    print(
        "Recurrence usable cases:              "
        f"{summary['recurrence_usable_cases']}"
    )

    print(
        "Recurrence events:                    "
        f"{summary['recurrence_events']}"
    )

    print(
        "Recurrence endpoint ready:            "
        f"{summary['recurrence_endpoint_ready']}"
    )

    print()

    print(
        "Primary endpoint selected:            "
        f"{summary['primary_endpoint_selected']}"
    )

    print_endpoint_readiness(
        report
    )

    print_selected_completeness(
        report
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
    print("M7.2 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()