"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_3c_controlled_access_limitation.py

Description:
    Smoke test for M7.3C Controlled-Access Limitation Resolution.

    This smoke test verifies that:
        - TCGA aligned miRNA BAM resources remain technically available;
        - direct quantification is deferred because authorization is absent;
        - the access limitation is not interpreted as negative evidence;
        - unresolved tRF candidates remain sequence-limited;
        - all candidates remain retained;
        - no quantification or clinical association is falsely claimed.

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

from pcatrfqtl.analysis.m7.runners.resolve_controlled_access_limitation import (
    M73CControlledAccessLimitationRunner,
)


# ============================================================================
# Paths
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m7_controlled_access_limitation.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "controlled_access_limitation"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


# ============================================================================
# Candidate display
# ============================================================================


def print_candidate_resolution(
    report: dict[str, Any],
) -> None:
    """Display candidate-level M7.3C results."""

    print()
    print("Candidate controlled-access resolution:")
    print()

    candidates = report[
        "candidate_controlled_access_limitation"
    ]

    for row in candidates:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    tRF ID:                         "
            f"{row['trf_id']}"
        )

        print(
            "    Sequence:                       "
            f"{row['sequence_dna']}"
        )

        print(
            "    Sequence ready:                 "
            f"{row['sequence_ready_for_quantification']}"
        )

        print(
            "    Aligned miRNA BAM files:        "
            f"{row['aligned_mirna_bam_file_count']}"
        )

        print(
            "    Aligned miRNA BAM cases:        "
            f"{row['aligned_mirna_bam_case_count']}"
        )

        print(
            "    Controlled BAM files:           "
            f"{row['controlled_aligned_bam_count']}"
        )

        print(
            "    Resource technically available: "
            f"{row['resource_technically_available']}"
        )

        print(
            "    Sequence limited:               "
            f"{row['sequence_limited']}"
        )

        print(
            "    Access limited:                 "
            f"{row['access_limited']}"
        )

        print(
            "    Authorized access available:    "
            f"{row['authorized_access_available']}"
        )

        print(
            "    TCGA quantification deferred:   "
            f"{row['tcga_direct_quantification_deferred']}"
        )

        print(
            "    TCGA quantification performed:  "
            f"{row['tcga_quantification_performed']}"
        )

        print(
            "    TCGA expression negative:       "
            f"{row['tcga_expression_negative']}"
        )

        print(
            "    Candidate retained:             "
            f"{row['candidate_retained']}"
        )

        print(
            "    Final status:"
        )

        print(
            "      "
            f"{row['final_status']}"
        )

        print()


# ============================================================================
# Strict validation
# ============================================================================


def validate_expected_state(
    report: dict[str, Any],
) -> None:
    """
    Validate expected final M7.3C state.

    These assertions are intentionally strict because M7.3C is designed
    to become a locked milestone.
    """

    # ----------------------------------------------------------------------
    # Top-level identity
    # ----------------------------------------------------------------------

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.3C"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "controlled_access_limitation_resolution"
    )

    summary = report[
        "summary"
    ]

    candidates = report[
        "candidate_controlled_access_limitation"
    ]

    # ----------------------------------------------------------------------
    # Overall state
    # ----------------------------------------------------------------------

    assert (
        summary[
            "overall_resolution_status"
        ]
        ==
        "CONTROLLED_ACCESS_UNAVAILABLE_AT_ANALYSIS_TIME"
    )

    # ----------------------------------------------------------------------
    # Candidate counts
    # ----------------------------------------------------------------------

    assert (
        summary[
            "candidate_trfs_assessed"
        ]
        ==
        3
    )

    assert (
        summary[
            "candidates_access_limited"
        ]
        ==
        1
    )

    assert (
        summary[
            "candidates_sequence_limited"
        ]
        ==
        2
    )

    assert (
        summary[
            "candidates_tcga_quantification_deferred"
        ]
        ==
        1
    )

    # ----------------------------------------------------------------------
    # Authorization / quantification safeguards
    # ----------------------------------------------------------------------

    assert (
        summary[
            "controlled_access_authorization_available"
        ]
        is False
    )

    assert (
        summary[
            "tcga_direct_trf_quantification_performed"
        ]
        is False
    )

    assert (
        summary[
            "processed_mirna_substitution_performed"
        ]
        is False
    )

    assert (
        summary[
            "clinical_association_performed"
        ]
        is False
    )

    assert (
        summary[
            "clinical_association_deferred"
        ]
        is True
    )

    assert (
        summary[
            "access_limitation_interpreted_as_negative_evidence"
        ]
        is False
    )

    # ----------------------------------------------------------------------
    # Candidate identity
    # ----------------------------------------------------------------------

    by_rsid = {
        row[
            "lead_rsid"
        ]:
        row
        for row in candidates
    }

    assert set(
        by_rsid.keys()
    ) == {
        "rs10216902",
        "rs2328376",
        "rs1288100",
    }

    # ----------------------------------------------------------------------
    # rs10216902
    # ----------------------------------------------------------------------

    rs102 = by_rsid[
        "rs10216902"
    ]

    assert (
        rs102[
            "priority_rank"
        ]
        ==
        1
    )

    assert (
        rs102[
            "trf_id"
        ]
        ==
        "tRF-24-S3M8309N0Y"
    )

    assert (
        rs102[
            "sequence_dna"
        ]
        ==
        "GTAGTCGTGGCCGAGTGGTTAAGG"
    )

    assert (
        rs102[
            "sequence_ready_for_quantification"
        ]
        is True
    )

    assert (
        rs102[
            "resource_technically_available"
        ]
        is True
    )

    assert (
        rs102[
            "sequence_limited"
        ]
        is False
    )

    assert (
        rs102[
            "access_limited"
        ]
        is True
    )

    assert (
        rs102[
            "authorized_access_available"
        ]
        is False
    )

    assert (
        rs102[
            "tcga_direct_quantification_deferred"
        ]
        is True
    )

    assert (
        rs102[
            "tcga_quantification_performed"
        ]
        is False
    )

    assert (
        rs102[
            "tcga_expression_negative"
        ]
        is False
    )

    assert (
        rs102[
            "candidate_retained"
        ]
        is True
    )

    assert (
        rs102[
            "final_status"
        ]
        ==
        "DIRECT_TCGA_TRF_QUANTIFICATION_DEFERRED_CONTROLLED_ACCESS"
    )

    # ----------------------------------------------------------------------
    # rs2328376
    # ----------------------------------------------------------------------

    rs232 = by_rsid[
        "rs2328376"
    ]

    assert (
        rs232[
            "priority_rank"
        ]
        ==
        2
    )

    assert (
        rs232[
            "trf_id"
        ]
        ==
        "tRF-17-WS72092"
    )

    assert (
        rs232[
            "sequence_ready_for_quantification"
        ]
        is False
    )

    assert (
        rs232[
            "sequence_limited"
        ]
        is True
    )

    assert (
        rs232[
            "access_limited"
        ]
        is False
    )

    assert (
        rs232[
            "authorized_access_available"
        ]
        is False
    )

    assert (
        rs232[
            "tcga_direct_quantification_deferred"
        ]
        is False
    )

    assert (
        rs232[
            "tcga_quantification_performed"
        ]
        is False
    )

    assert (
        rs232[
            "tcga_expression_negative"
        ]
        is False
    )

    assert (
        rs232[
            "candidate_retained"
        ]
        is True
    )

    assert (
        rs232[
            "final_status"
        ]
        ==
        "DIRECT_TRF_QUANTIFICATION_BLOCKED_SEQUENCE_UNRESOLVED"
    )

    # ----------------------------------------------------------------------
    # rs1288100
    # ----------------------------------------------------------------------

    rs128 = by_rsid[
        "rs1288100"
    ]

    assert (
        rs128[
            "priority_rank"
        ]
        ==
        3
    )

    assert (
        rs128[
            "trf_id"
        ]
        ==
        "tRF-30-7EMQ18Y3E7QN"
    )

    assert (
        rs128[
            "sequence_ready_for_quantification"
        ]
        is False
    )

    assert (
        rs128[
            "sequence_limited"
        ]
        is True
    )

    assert (
        rs128[
            "access_limited"
        ]
        is False
    )

    assert (
        rs128[
            "authorized_access_available"
        ]
        is False
    )

    assert (
        rs128[
            "tcga_direct_quantification_deferred"
        ]
        is False
    )

    assert (
        rs128[
            "tcga_quantification_performed"
        ]
        is False
    )

    assert (
        rs128[
            "tcga_expression_negative"
        ]
        is False
    )

    assert (
        rs128[
            "candidate_retained"
        ]
        is True
    )

    assert (
        rs128[
            "final_status"
        ]
        ==
        "DIRECT_TRF_QUANTIFICATION_BLOCKED_SEQUENCE_UNRESOLVED"
    )

    # ----------------------------------------------------------------------
    # Continuation routes
    # ----------------------------------------------------------------------

    assert (
        summary[
            "next_external_validation_stage"
        ]
        ==
        "M7.4A_OPEN_SMALL_RNA_DATASET_DISCOVERY"
    )

    assert (
        summary[
            "next_variant_validation_stage"
        ]
        ==
        "M7.5_VARIANT_LEVEL_VALIDATION_FEASIBILITY"
    )

    assert (
        summary[
            "clinical_route"
        ]
        ==
        "RETAIN_CLINICAL_ENDPOINTS_ASSOCIATION_DEFERRED"
    )


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M7.3C smoke test."""

    print()
    print("=" * 72)
    print("M7.3C CONTROLLED-ACCESS LIMITATION RESOLUTION")
    print("=" * 72)

    runner = M73CControlledAccessLimitationRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    # ----------------------------------------------------------------------
    # Strict assertions
    # ----------------------------------------------------------------------

    validate_expected_state(
        report
    )

    summary = report[
        "summary"
    ]

    # ----------------------------------------------------------------------
    # Summary
    # ----------------------------------------------------------------------

    print()
    print("=" * 72)
    print("M7.3C RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidate tRFs assessed:                    "
        f"{summary['candidate_trfs_assessed']}"
    )

    print(
        "Candidates access-limited:                  "
        f"{summary['candidates_access_limited']}"
    )

    print(
        "Candidates sequence-limited:                "
        f"{summary['candidates_sequence_limited']}"
    )

    print(
        "TCGA quantification deferred:               "
        f"{summary['candidates_tcga_quantification_deferred']}"
    )

    print()

    print(
        "Controlled-access authorization available: "
        f"{summary['controlled_access_authorization_available']}"
    )

    print(
        "TCGA direct tRF quantification performed:  "
        f"{summary['tcga_direct_trf_quantification_performed']}"
    )

    print(
        "Processed miRNA substitution performed:    "
        f"{summary['processed_mirna_substitution_performed']}"
    )

    print(
        "Clinical association performed:            "
        f"{summary['clinical_association_performed']}"
    )

    print(
        "Clinical association deferred:             "
        f"{summary['clinical_association_deferred']}"
    )

    print(
        "Access limitation treated as negative:     "
        f"{summary['access_limitation_interpreted_as_negative_evidence']}"
    )

    print()

    print("Overall resolution:")
    print(
        f"  {summary['overall_resolution_status']}"
    )

    # ----------------------------------------------------------------------
    # Candidate details
    # ----------------------------------------------------------------------

    print_candidate_resolution(
        report
    )

    # ----------------------------------------------------------------------
    # Continuation
    # ----------------------------------------------------------------------

    print()
    print("Continuation routes:")
    print()

    print(
        "  External validation:"
    )
    print(
        "    "
        f"{summary['next_external_validation_stage']}"
    )

    print()

    print(
        "  Variant validation:"
    )
    print(
        "    "
        f"{summary['next_variant_validation_stage']}"
    )

    print()

    print(
        "  Clinical route:"
    )
    print(
        "    "
        f"{summary['clinical_route']}"
    )

    # ----------------------------------------------------------------------
    # QC
    # ----------------------------------------------------------------------

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
    print("M7.3C SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()