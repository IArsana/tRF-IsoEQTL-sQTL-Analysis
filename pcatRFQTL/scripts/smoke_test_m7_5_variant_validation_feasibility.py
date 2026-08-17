"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_5_variant_validation_feasibility.py

Description:
    Smoke test for M7.5 Variant-Level Validation Feasibility.

    The test validates workflow invariants rather than hard-coding the number
    of curated associations because GWAS Catalog content may be updated.

    Required invariants:

        - exactly three prioritized candidates are assessed;
        - API success and zero results remain distinguishable;
        - API zero results are never negative evidence;
        - API failure is never negative evidence;
        - TCGA somatic calls are not used as germline;
        - independent replication remains unverified;
        - full summary-statistics audit remains the next stage.

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

from pcatrfqtl.analysis.m7.runners.assess_variant_validation_feasibility import (
    M75VariantValidationFeasibilityRunner,
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
    "m7_variant_validation_feasibility.yaml"
)


EXPECTED_RSIDS = {
    "rs10216902",
    "rs2328376",
    "rs1288100",
}


def validate_report(
    report: dict[str, Any],
) -> None:
    """Validate structural and scientific M7.5 invariants."""

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.5"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "variant_level_validation_feasibility"
    )

    # ----------------------------------------------------------------------
    # API identity
    # ----------------------------------------------------------------------

    assert (
        report[
            "api"
        ][
            "provider"
        ]
        ==
        "NHGRI_EBI_GWAS_CATALOG"
    )

    assert (
        report[
            "api"
        ][
            "version"
        ]
        ==
        "V2"
    )

    assert (
        report[
            "api"
        ][
            "endpoint"
        ]
        ==
        "/v2/associations"
    )

    assert (
        report[
            "api"
        ][
            "exact_variant_parameter"
        ]
        ==
        "rs_id"
    )

    # ----------------------------------------------------------------------
    # Candidate set
    # ----------------------------------------------------------------------

    summary = report[
        "summary"
    ]

    candidates = report[
        "candidate_feasibility"
    ]

    acquisitions = report[
        "acquisition_status"
    ]

    assert (
        summary[
            "candidates_assessed"
        ]
        ==
        3
    )

    assert (
        len(
            candidates
        )
        ==
        3
    )

    assert (
        len(
            acquisitions
        )
        ==
        3
    )

    assert {
        row[
            "rsid"
        ]
        for row in candidates
    } == EXPECTED_RSIDS

    assert {
        row[
            "rsid"
        ]
        for row in acquisitions
    } == EXPECTED_RSIDS

    # ----------------------------------------------------------------------
    # Acquisition safeguards
    # ----------------------------------------------------------------------

    for row in acquisitions:

        assert (
            row[
                "negative_evidence"
            ]
            is False
        )

        assert (
            row[
                "variant_absence_inferred"
            ]
            is False
        )

        assert (
            row[
                "records_retrieved"
            ]
            >=
            0
        )

        assert (
            row[
                "total_elements"
            ]
            >=
            0
        )

        if row[
            "request_success"
        ]:

            assert row[
                "acquisition_status"
            ] in {
                "REST_API_V2_QUERY_SUCCESSFUL_ASSOCIATIONS_RETURNED",
                "REST_API_V2_QUERY_SUCCESSFUL_NO_CURATED_ASSOCIATION",
            }

            assert (
                row[
                    "records_retrieved"
                ]
                ==
                row[
                    "total_elements"
                ]
            )

            if (
                row[
                    "total_elements"
                ]
                ==
                0
            ):

                assert (
                    row[
                        "acquisition_status"
                    ]
                    ==
                    "REST_API_V2_QUERY_SUCCESSFUL_NO_CURATED_ASSOCIATION"
                )

        else:

            assert row[
                "acquisition_status"
            ] in {
                "REST_API_V2_SERVICE_UNAVAILABLE",
                "REST_API_V2_RESPONSE_SCHEMA_INVALID",
            }

    # ----------------------------------------------------------------------
    # Candidate-level safeguards
    # ----------------------------------------------------------------------

    for row in candidates:

        assert (
            row[
                "study_independence_verified"
            ]
            is False
        )

        assert (
            row[
                "independent_replication_established"
            ]
            is False
        )

        assert (
            row[
                "full_summary_statistics_audit_required"
            ]
            is True
        )

        assert (
            row[
                "tcga_germline_route_available_now"
            ]
            is False
        )

        assert (
            row[
                "tcga_germline_deferred_controlled_access"
            ]
            is True
        )

        assert (
            row[
                "tcga_somatic_surrogate_allowed"
            ]
            is False
        )

        assert (
            row[
                "api_zero_result_interpreted_as_negative"
            ]
            is False
        )

        assert (
            row[
                "api_failure_interpreted_as_negative"
            ]
            is False
        )

        assert (
            row[
                "variant_absence_inferred"
            ]
            is False
        )

        assert row[
            "candidate_status"
        ] in {
            "CURATED_PROSTATE_ASSOCIATION_FOUND_INDEPENDENCE_UNRESOLVED",
            "CURATED_PROSTATE_ASSOCIATION_PREVIOUSLY_USED_NO_NEW_VALIDATION",
            "CURATED_ASSOCIATION_FOUND_NO_PROSTATE_TRAIT_MATCH",
            "NO_CURATED_ASSOCIATION_FULL_SUMMARY_STATISTICS_AUDIT_REQUIRED",
            "CURATED_ASSOCIATION_AUDIT_DEFERRED_API_UNAVAILABLE",
        }

    # ----------------------------------------------------------------------
    # Association-level safeguards
    # ----------------------------------------------------------------------

    for row in report[
        "association_inventory"
    ]:

        assert (
            row[
                "study_independence_verified"
            ]
            is False
        )

        assert (
            row[
                "independent_replication_claimed"
            ]
            is False
        )

        assert (
            row[
                "negative_evidence"
            ]
            is False
        )

    # ----------------------------------------------------------------------
    # Overall safeguards
    # ----------------------------------------------------------------------

    assert (
        summary[
            "independent_replications_verified"
        ]
        ==
        0
    )

    assert (
        summary[
            "independent_replication_claimed"
        ]
        is False
    )

    assert (
        summary[
            "study_independence_verified"
        ]
        is False
    )

    assert (
        summary[
            "phenotype_compatibility_verified"
        ]
        is False
    )

    assert (
        summary[
            "ancestry_compatibility_verified"
        ]
        is False
    )

    assert (
        summary[
            "effect_allele_harmonization_performed"
        ]
        is False
    )

    assert (
        summary[
            "effect_direction_concordance_claimed"
        ]
        is False
    )

    assert (
        summary[
            "full_summary_statistics_audit_performed"
        ]
        is False
    )

    assert (
        summary[
            "tcga_germline_route_available_now"
        ]
        is False
    )

    assert (
        summary[
            "tcga_germline_deferred_controlled_access"
        ]
        is True
    )

    assert (
        summary[
            "tcga_somatic_used_as_germline"
        ]
        is False
    )

    assert (
        summary[
            "api_zero_result_means_negative_evidence"
        ]
        is False
    )

    assert (
        summary[
            "api_failure_means_negative_evidence"
        ]
        is False
    )

    assert (
        summary[
            "api_zero_result_means_variant_absent"
        ]
        is False
    )

    assert (
        summary[
            "full_summary_statistics_audit_required_candidates"
        ]
        ==
        3
    )

    assert summary[
        "overall_status"
    ] in {
        "CURATED_VARIANT_EVIDENCE_AVAILABLE_INDEPENDENCE_UNRESOLVED",
        "CURATED_VARIANT_EVIDENCE_LIMITED_TO_PREVIOUSLY_USED_STUDIES",
        "NO_CURATED_VARIANT_EVIDENCE_FULL_SUMMARY_STATISTICS_AUDIT_REQUIRED",
        "CURATED_VARIANT_AUDIT_INCOMPLETE_API_UNAVAILABLE",
    }

    assert (
        summary[
            "next_stage"
        ]
        ==
        "M7.5B_FULL_SUMMARY_STATISTICS_AUDIT"
    )


def main() -> None:
    """Execute M7.5 smoke test."""

    print()
    print("=" * 72)
    print("M7.5 VARIANT-LEVEL VALIDATION FEASIBILITY")
    print("=" * 72)

    runner = M75VariantValidationFeasibilityRunner(
        project_root=ROOT,
        config_path=CONFIG,
    )

    report = runner.run()

    validate_report(
        report
    )

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M7.5 RESULT SUMMARY")
    print("=" * 72)

    print(
        "API provider:                         "
        f"{report['api']['provider']}"
    )

    print(
        "API version:                          "
        f"{report['api']['version']}"
    )

    print(
        "Endpoint:                             "
        f"{report['api']['endpoint']}"
    )

    print()

    print(
        "Candidates assessed:                  "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Successful API queries:               "
        f"{summary['api_queries_successful']}"
    )

    print(
        "Zero curated-result queries:          "
        f"{summary['api_queries_zero_curated_results']}"
    )

    print()

    print("Candidate feasibility:")
    print()

    for row in report[
        "candidate_feasibility"
    ]:

        print(
            f"  {row['rsid']} | "
            f"{row['priority_class']}"
        )

        print(
            "    API successful:                 "
            f"{row['api_query_successful']}"
        )

        print(
            "    curated records:                "
            f"{row['curated_records_returned']}"
        )

        print(
            "    exact candidate associations:   "
            f"{row['exact_candidate_associations']}"
        )

        print(
            "    prostate-trait candidates:      "
            f"{row['prostate_trait_candidate_associations']}"
        )

        print(
            "    reused prostate studies:        "
            f"{row['reused_prostate_studies']}"
        )

        print(
            "    non-reused prostate studies:    "
            f"{row['nonreused_prostate_studies']}"
        )

        print(
            "    independence verified:          "
            f"{row['study_independence_verified']}"
        )

        print(
            "    full sumstats audit required:   "
            f"{row['full_summary_statistics_audit_required']}"
        )

        print(
            "    status:"
        )

        print(
            "      "
            f"{row['candidate_status']}"
        )

        print()

    print("Acquisition:")
    print()

    for row in report[
        "acquisition_status"
    ]:

        print(
            f"  {row['rsid']}"
        )

        print(
            "    request successful:             "
            f"{row['request_success']}"
        )

        print(
            "    HTTP status:                    "
            f"{row['http_status']}"
        )

        print(
            "    total curated associations:     "
            f"{row['total_elements']}"
        )

        print(
            "    pages retrieved:                "
            f"{row['pages_retrieved']}"
        )

        print(
            "    status:"
        )

        print(
            "      "
            f"{row['acquisition_status']}"
        )

        print()

    print(
        "Independent replications verified:    "
        f"{summary['independent_replications_verified']}"
    )

    print(
        "Independent replication claimed:      "
        f"{summary['independent_replication_claimed']}"
    )

    print(
        "Full summary-statistics performed:    "
        f"{summary['full_summary_statistics_audit_performed']}"
    )

    print(
        "TCGA germline deferred:               "
        f"{summary['tcga_germline_deferred_controlled_access']}"
    )

    print(
        "TCGA somatic used as germline:        "
        f"{summary['tcga_somatic_used_as_germline']}"
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
    print("M7.5 SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()