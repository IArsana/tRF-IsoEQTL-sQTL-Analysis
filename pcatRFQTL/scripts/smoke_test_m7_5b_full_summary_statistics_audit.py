"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_5b_full_summary_statistics_audit.py

Description:
    Smoke test for M7.5B Full Summary Statistics Audit.

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

from pcatrfqtl.analysis.m7.runners.audit_full_summary_statistics import (
    M75BFullSummaryStatisticsAuditRunner,
)


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG = (
    ROOT
    /
    "configs"
    /
    "m7_full_summary_statistics_audit.yaml"
)


EXPECTED_RSIDS = {
    "rs10216902",
    "rs2328376",
    "rs1288100",
}


def validate_report(
    report: dict[str, Any],
) -> None:
    """Validate M7.5B invariants."""

    assert report[
        "milestone"
    ] == "M7.5B"

    assert report[
        "stage"
    ] == "full_summary_statistics_audit"

    summary = report[
        "summary"
    ]

    studies = report[
        "study_inventory"
    ]

    files = report[
        "file_manifest"
    ]

    candidates = report[
        "candidate_resolution"
    ]

    # ----------------------------------------------------------------------
    # Known study-discovery state
    # ----------------------------------------------------------------------

    assert (
        summary[
            "prostate_related_studies_discovered"
        ]
        ==
        35
    )

    assert (
        summary[
            "full_summary_statistics_available_studies"
        ]
        ==
        20
    )

    assert (
        summary[
            "eligible_full_sumstats_studies"
        ]
        ==
        13
    )

    assert (
        summary[
            "phenotype_excluded_or_unresolved_studies"
        ]
        ==
        7
    )

    # ----------------------------------------------------------------------
    # Candidate set
    # ----------------------------------------------------------------------

    assert (
        summary[
            "candidate_variants_assessed"
        ]
        ==
        3
    )

    assert {
        row["rsid"]
        for row in candidates
    } == EXPECTED_RSIDS

    # ----------------------------------------------------------------------
    # Eligibility
    # ----------------------------------------------------------------------

    eligible_accessions = {
        row["accession_id"]
        for row in studies
        if row[
            "eligible_for_exact_lookup"
        ]
    }

    assert (
        len(
            eligible_accessions
        )
        ==
        13
    )

    for row in studies:

        if row[
            "eligible_for_exact_lookup"
        ]:

            assert row[
                "full_summary_stats_available"
            ] is True

            assert row[
                "direct_prostate_phenotype"
            ] is True

            assert row[
                "phenotype_screening_status"
            ] == (
                "DIRECT_PROSTATE_CANCER_PHENOTYPE"
            )

        assert row[
            "phenotype_compatibility_verified"
        ] is False

        assert row[
            "ancestry_compatibility_verified"
        ] is False

        assert row[
            "study_independence_verified"
        ] is False

    # ----------------------------------------------------------------------
    # File resolution
    # ----------------------------------------------------------------------

    manifest_accessions = {
        row[
            "accession_id"
        ]
        for row in files
    }

    assert (
        manifest_accessions
        <=
        eligible_accessions
    )

    if summary[
        "full_summary_statistics_audit_complete"
    ]:

        assert (
            summary[
                "summary_statistics_files_selected"
            ]
            ==
            13
        )

        assert (
            summary[
                "summary_statistics_files_scanned_successfully"
            ]
            ==
            13
        )

    for row in files:

        if row[
            "scan_success"
        ]:

            assert row[
                "selected_for_scan"
            ] is True

            assert row[
                "eof_reached"
            ] is True

            assert isinstance(
                row[
                    "source_file_sha256"
                ],
                str,
            )

            assert (
                len(
                    row[
                        "source_file_sha256"
                    ]
                )
                ==
                64
            )

    # ----------------------------------------------------------------------
    # Specific file-selection invariants discovered in probe
    # ----------------------------------------------------------------------

    manifest_by_accession = {
        row[
            "accession_id"
        ]:
            row
        for row in files
    }

    if (
        "GCST90043894"
        in
        manifest_by_accession
    ):

        row = manifest_by_accession[
            "GCST90043894"
        ]

        assert row[
            "harmonised"
        ] is True

        assert row[
            "filename"
        ].endswith(
            ".h.tsv.gz"
        )

    if (
        "GCST90077635"
        in
        manifest_by_accession
    ):

        row = manifest_by_accession[
            "GCST90077635"
        ]

        assert row[
            "harmonised"
        ] is True

        assert row[
            "filename"
        ].endswith(
            ".h.tsv.gz"
        )

    for accession in (
        "GCST90296485",
        "GCST90296486",
    ):

        if accession in manifest_by_accession:

            row = manifest_by_accession[
                accession
            ]

            assert row[
                "compression"
            ] == "NONE"

            assert row[
                "filename"
            ].endswith(
                ".tsv"
            )

    # ----------------------------------------------------------------------
    # Exact hits
    # ----------------------------------------------------------------------

    for row in report[
        "exact_variant_hits"
    ]:

        assert row[
            "rs_id"
        ] in EXPECTED_RSIDS

        assert row[
            "exact_rsid_match"
        ] is True

    # ----------------------------------------------------------------------
    # Candidate scientific safeguards
    # ----------------------------------------------------------------------

    allowed_candidate_statuses = {
        "EXACT_VARIANT_FOUND_IN_NONREUSED_FULL_SUMMARY_STATISTICS_INDEPENDENCE_UNRESOLVED",
        "EXACT_VARIANT_FOUND_ONLY_IN_PREVIOUSLY_USED_FULL_SUMMARY_STATISTICS",
        "EXACT_VARIANT_NOT_FOUND_IN_AUDITED_FULL_SUMMARY_STATISTICS",
        "FULL_SUMMARY_STATISTICS_AUDIT_INCOMPLETE",
    }

    for row in candidates:

        assert (
            row[
                "candidate_status"
            ]
            in
            allowed_candidate_statuses
        )

        assert row[
            "independent_replication_established"
        ] is False

        assert row[
            "study_independence_verified"
        ] is False

        assert row[
            "phenotype_compatibility_verified"
        ] is False

        assert row[
            "ancestry_compatibility_verified"
        ] is False

        assert row[
            "effect_allele_harmonization_performed"
        ] is False

    # ----------------------------------------------------------------------
    # Global safeguards
    # ----------------------------------------------------------------------

    assert summary[
        "independent_replications_verified"
    ] == 0

    assert summary[
        "independent_replication_claimed"
    ] is False

    assert summary[
        "study_independence_verified"
    ] is False

    assert summary[
        "phenotype_compatibility_verified"
    ] is False

    assert summary[
        "ancestry_compatibility_verified"
    ] is False

    assert summary[
        "effect_allele_harmonization_performed"
    ] is False

    assert summary[
        "effect_direction_concordance_claimed"
    ] is False

    assert summary[
        "proxy_variant_lookup_performed"
    ] is False

    assert summary[
        "coordinate_substitution_performed"
    ] is False

    assert summary[
        "candidate_substitution_performed"
    ] is False

    assert summary[
        "excluded_phenotype_means_negative_evidence"
    ] is False

    assert summary[
        "candidate_not_found_means_no_association"
    ] is False

    assert summary[
        "overall_status"
    ] in {
        "FULL_SUMMARY_STATISTICS_EXACT_VARIANT_EVIDENCE_AVAILABLE_INDEPENDENCE_UNRESOLVED",
        "FULL_SUMMARY_STATISTICS_EVIDENCE_LIMITED_TO_PREVIOUSLY_USED_STUDIES",
        "NO_EXACT_VARIANT_FOUND_IN_AUDITED_FULL_SUMMARY_STATISTICS",
        "FULL_SUMMARY_STATISTICS_AUDIT_INCOMPLETE",
    }

    assert summary[
        "next_stage"
    ] == (
        "M7.5C_STUDY_INDEPENDENCE_AND_HARMONIZATION_AUDIT"
    )


def main() -> None:
    """Execute M7.5B."""

    print()
    print("=" * 76)
    print("M7.5B FULL SUMMARY STATISTICS AUDIT")
    print("=" * 76)

    runner = M75BFullSummaryStatisticsAuditRunner(
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
    print("=" * 76)
    print("M7.5B RESULT SUMMARY")
    print("=" * 76)

    print(
        "Prostate-related studies discovered:      "
        f"{summary['prostate_related_studies_discovered']}"
    )

    print(
        "Full summary statistics available:        "
        f"{summary['full_summary_statistics_available_studies']}"
    )

    print(
        "Direct phenotype studies:                 "
        f"{summary['direct_prostate_phenotype_studies']}"
    )

    print(
        "Phenotype excluded/unresolved:            "
        f"{summary['phenotype_excluded_or_unresolved_studies']}"
    )

    print(
        "Eligible studies:                         "
        f"{summary['eligible_full_sumstats_studies']}"
    )

    print(
        "Files selected:                           "
        f"{summary['summary_statistics_files_selected']}"
    )

    print(
        "Files scanned successfully:               "
        f"{summary['summary_statistics_files_scanned_successfully']}"
    )

    print(
        "Exact candidate rows found:               "
        f"{summary['exact_variant_rows_found']}"
    )

    print()

    print("Selected files:")
    print()

    for row in report[
        "file_manifest"
    ]:

        print(
            f"  {row['accession_id']}"
        )

        print(
            "    representation: "
            f"{row['file_representation']}"
        )

        print(
            "    compression:    "
            f"{row['compression']}"
        )

        print(
            "    harmonised:     "
            f"{row['harmonised']}"
        )

        print(
            "    filename:       "
            f"{row['filename']}"
        )

        print(
            "    scan success:   "
            f"{row['scan_success']}"
        )

        print(
            "    EOF reached:    "
            f"{row['eof_reached']}"
        )

        print(
            "    exact hits:     "
            f"{row['candidate_hits']}"
        )

        print()

    print("Candidate resolution:")
    print()

    for row in report[
        "candidate_resolution"
    ]:

        print(
            f"  {row['rsid']} | "
            f"{row['priority_class']}"
        )

        print(
            "    eligible studies:              "
            f"{row['eligible_studies_assessed']}"
        )

        print(
            "    auditable studies:              "
            f"{row['auditable_studies']}"
        )

        print(
            "    studies with exact variant:     "
            f"{row['studies_with_exact_variant']}"
        )

        print(
            "    reused studies:                 "
            f"{row['reused_studies_with_exact_variant']}"
        )

        print(
            "    non-reused studies:             "
            f"{row['nonreused_studies_with_exact_variant']}"
        )

        print(
            "    status:"
        )

        print(
            "      "
            f"{row['candidate_status']}"
        )

        print()

    print(
        "Independent replications verified:        "
        f"{summary['independent_replications_verified']}"
    )

    print(
        "Independent replication claimed:          "
        f"{summary['independent_replication_claimed']}"
    )

    print(
        "Full audit complete:                      "
        f"{summary['full_summary_statistics_audit_complete']}"
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
    print("=" * 76)
    print("M7.5B SMOKE TEST PASSED")
    print("=" * 76)


if __name__ == "__main__":
    main()