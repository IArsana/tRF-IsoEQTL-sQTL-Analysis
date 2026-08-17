"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4d_external_clinical_association_readiness.py

Description:
    Smoke test for M7.4D External Clinical Association Readiness.

    This test verifies that GSE80400 molecular evidence is retained while
    inferential clinical association is not performed when the external
    cohort lacks sufficiently supported clinical groups and verified
    run-to-sample mapping.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pcatrfqtl.analysis.m7.runners.assess_external_clinical_association_readiness import (
    M74DExternalClinicalAssociationReadinessRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


CONFIG = (
    ROOT
    /
    "configs"
    /
    "m7_external_clinical_association_readiness.yaml"
)


def validate_expected_state(
    report: dict[str, Any],
) -> None:
    """Validate expected conservative M7.4D state."""

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.4D"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "external_clinical_association_readiness"
    )

    summary = report[
        "summary"
    ]

    assert (
        summary[
            "dataset_key"
        ]
        ==
        "m7_gse80400"
    )

    assert (
        summary[
            "geo_accession"
        ]
        ==
        "GSE80400"
    )

    assert (
        summary[
            "source_samples"
        ]
        ==
        11
    )

    assert (
        summary[
            "ffpe_samples"
        ]
        ==
        1
    )

    assert (
        summary[
            "derived_or_nonindependent_samples"
        ]
        ==
        1
    )

    assert (
        summary[
            "run_sample_mapping_verified"
        ]
        is False
    )

    assert (
        summary[
            "endpoints_assessed"
        ]
        ==
        6
    )

    # ----------------------------------------------------------------------
    # No endpoint should currently be association-ready.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "endpoints_association_ready"
        ]
        ==
        0
    )

    assert (
        summary[
            "endpoints_not_ready"
        ]
        ==
        6
    )

    # ----------------------------------------------------------------------
    # No inferential analysis
    # ----------------------------------------------------------------------

    assert (
        summary[
            "association_model_fitted"
        ]
        is False
    )

    assert (
        summary[
            "hypothesis_test_performed"
        ]
        is False
    )

    assert (
        summary[
            "p_value_calculated"
        ]
        is False
    )

    assert (
        summary[
            "effect_size_estimated"
        ]
        is False
    )

    assert (
        summary[
            "survival_analysis_performed"
        ]
        is False
    )

    assert (
        summary[
            "recurrence_model_performed"
        ]
        is False
    )

    assert (
        summary[
            "clinical_validation_claimed"
        ]
        is False
    )

    assert (
        summary[
            "mature_24nt_trf_abundance_claimed"
        ]
        is False
    )

    # ----------------------------------------------------------------------
    # Molecular sequence evidence is retained.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "external_molecular_evidence_retained"
        ]
        is True
    )

    assert (
        summary[
            "candidate_abundance_interpretation"
        ]
        ==
        "EXACT_CANDIDATE_SEQUENCE_MOTIF_ABUNDANCE"
    )

    # ----------------------------------------------------------------------
    # Endpoint-level safeguards
    # ----------------------------------------------------------------------

    endpoints = report[
        "endpoint_readiness"
    ]

    assert (
        len(
            endpoints
        )
        ==
        6
    )

    for row in endpoints:

        assert (
            row[
                "association_ready"
            ]
            is False
        )

        assert (
            row[
                "association_model_fitted"
            ]
            is False
        )

        assert (
            row[
                "p_value_calculated"
            ]
            is False
        )

        assert (
            row[
                "effect_size_estimated"
            ]
            is False
        )

    # ----------------------------------------------------------------------
    # Overall
    # ----------------------------------------------------------------------

    assert (
        summary[
            "overall_status"
        ]
        ==
        "EXTERNAL_MOLECULAR_VALIDATION_RETAINED_CLINICAL_ASSOCIATION_NOT_JUSTIFIED"
    )

    assert (
        summary[
            "next_stage"
        ]
        ==
        "M7.5_VARIANT_LEVEL_VALIDATION_FEASIBILITY"
    )


def main() -> None:
    """Execute M7.4D smoke test."""

    print()
    print("=" * 72)
    print("M7.4D EXTERNAL CLINICAL ASSOCIATION READINESS")
    print("=" * 72)

    runner = M74DExternalClinicalAssociationReadinessRunner(
        project_root=ROOT,
        config_path=CONFIG,
    )

    report = runner.run()

    validate_expected_state(
        report
    )

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M7.4D RESULT SUMMARY")
    print("=" * 72)

    print(
        "Dataset:                          "
        f"{summary['geo_accession']}"
    )

    print(
        "Source samples:                   "
        f"{summary['source_samples']}"
    )

    print(
        "FFPE samples:                     "
        f"{summary['ffpe_samples']}"
    )

    print(
        "Derived/non-independent samples: "
        f"{summary['derived_or_nonindependent_samples']}"
    )

    print()

    print(
        "Run-sample mapping verified:      "
        f"{summary['run_sample_mapping_verified']}"
    )

    print(
        "Abundance interpretation:         "
        f"{summary['candidate_abundance_interpretation']}"
    )

    print()

    print(
        "Endpoints assessed:               "
        f"{summary['endpoints_assessed']}"
    )

    print(
        "Association-ready endpoints:      "
        f"{summary['endpoints_association_ready']}"
    )

    print(
        "Not-ready endpoints:              "
        f"{summary['endpoints_not_ready']}"
    )

    print()

    print("Endpoint readiness:")
    print()

    for row in report[
        "endpoint_readiness"
    ]:

        print(
            f"  {row['endpoint']}"
        )

        print(
            "    usable cases:               "
            f"{row['usable_cases']}"
        )

        print(
            "    groups:                     "
            f"{', '.join(row['group_counts'])}"
        )

        print(
            "    mapping verified:           "
            f"{row['run_sample_mapping_verified']}"
        )

        print(
            "    metadata complete:          "
            f"{row['metadata_complete']}"
        )

        print(
            "    minimum group support:      "
            f"{row['minimum_group_support_met']}"
        )

        print(
            "    association ready:          "
            f"{row['association_ready']}"
        )

        print(
            "    status:"
        )

        print(
            "      "
            f"{row['readiness_status']}"
        )

        print()

    print(
        "Association model fitted:         "
        f"{summary['association_model_fitted']}"
    )

    print(
        "P-value calculated:               "
        f"{summary['p_value_calculated']}"
    )

    print(
        "Effect size estimated:            "
        f"{summary['effect_size_estimated']}"
    )

    print(
        "Clinical validation claimed:      "
        f"{summary['clinical_validation_claimed']}"
    )

    print(
        "External molecular evidence kept: "
        f"{summary['external_molecular_evidence_retained']}"
    )

    print()

    print("Overall:")
    print(
        "  "
        f"{summary['overall_status']}"
    )

    print()

    print("Next:")
    print(
        "  "
        f"{summary['next_stage']}"
    )

    print()
    print("=" * 72)
    print("M7.4D SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()