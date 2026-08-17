"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4b_external_small_rna_technical_eligibility.py

Description:
    Smoke test for M7.4B External Small-RNA Technical Eligibility.

    The smoke test verifies that:
        - all four M7.4A datasets are assessed;
        - at least one external raw-data route is technically eligible;
        - GSE80400 remains suitable for exact-sequence auditing;
        - GSE229904 and GSE89193 remain technically usable with caveats;
        - GSE290219 remains blocked by incomplete technical metadata;
        - no raw reads are downloaded;
        - no candidate sequence search or tRF quantification is performed.

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

from pcatrfqtl.analysis.m7.runners.assess_external_small_rna_technical_eligibility import (
    M74BExternalSmallRnaTechnicalEligibilityRunner,
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
    / "m7_external_small_rna_technical_eligibility.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "external_small_rna_technical_eligibility"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


# ============================================================================
# Display helper
# ============================================================================


def print_dataset_results(
    report: dict[str, Any],
) -> None:
    """Display dataset-level technical eligibility."""

    print()
    print("Dataset technical eligibility:")
    print()

    for row in report[
        "dataset_technical_eligibility"
    ]:

        print(
            f"  {row['dataset_key']}"
        )

        print(
            "    Raw public:                    "
            f"{row['raw_sequence_public']}"
        )

        print(
            "    Library strategy:              "
            f"{row['library_strategy']}"
        )

        print(
            "    Small-RNA library supported:   "
            f"{row['small_rna_library_supported']}"
        )

        print(
            "    Library preparation supported: "
            f"{row['library_preparation_supported']}"
        )

        print(
            "    Candidate length compatible:   "
            f"{row['candidate_length_compatible']}"
        )

        print(
            "    Read orientation resolvable:   "
            f"{row['read_orientation_resolvable']}"
        )

        print(
            "    Adapter handling resolvable:   "
            f"{row['adapter_handling_resolvable']}"
        )

        print(
            "    FFPE:                          "
            f"{row['ffpe']}"
        )

        print(
            "    Caveat count:                  "
            f"{row['technical_caveat_count']}"
        )

        print(
            "    Technical eligibility:         "
            f"{row['technical_eligibility']}"
        )

        print(
            "    Resolution class:              "
            f"{row['resolution_class']}"
        )

        print(
            "    Technical status:"
        )

        print(
            "      "
            f"{row['technical_status']}"
        )

        if row[
            "technical_caveats"
        ]:

            print(
                "    Caveats:"
            )

            for caveat in row[
                "technical_caveats"
            ]:

                print(
                    "      - "
                    f"{caveat}"
                )

        print()


# ============================================================================
# Strict expected-state validation
# ============================================================================


def validate_expected_state(
    report: dict[str, Any],
) -> None:
    """Validate expected locked M7.4B state."""

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.4B"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "external_small_rna_technical_eligibility"
    )

    summary = report[
        "summary"
    ]

    rows = report[
        "dataset_technical_eligibility"
    ]

    assert (
        summary[
            "datasets_assessed"
        ]
        ==
        4
    )

    assert (
        summary[
            "datasets_technically_eligible"
        ]
        ==
        3
    )

    # All three currently eligible routes carry documented caveats.
    assert (
        summary[
            "datasets_fully_eligible"
        ]
        ==
        0
    )

    assert (
        summary[
            "datasets_eligible_with_caveats"
        ]
        ==
        3
    )

    assert (
        summary[
            "datasets_metadata_incomplete"
        ]
        ==
        1
    )

    assert (
        summary[
            "datasets_not_eligible"
        ]
        ==
        0
    )

    assert (
        summary[
            "candidate_trf"
        ]
        ==
        "tRF-24-S3M8309N0Y"
    )

    assert (
        summary[
            "candidate_sequence"
        ]
        ==
        "GTAGTCGTGGCCGAGTGGTTAAGG"
    )

    assert (
        summary[
            "candidate_length_nt"
        ]
        ==
        24
    )

    # ----------------------------------------------------------------------
    # M7.4B safeguards
    # ----------------------------------------------------------------------

    assert (
        summary[
            "fastq_download_performed"
        ]
        is False
    )

    assert (
        summary[
            "candidate_sequence_search_performed"
        ]
        is False
    )

    assert (
        summary[
            "candidate_detection_claimed"
        ]
        is False
    )

    assert (
        summary[
            "candidate_absence_claimed"
        ]
        is False
    )

    assert (
        summary[
            "trf_quantification_performed"
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
            "overall_technical_status"
        ]
        ==
        "EXTERNAL_SMALL_RNA_QUANTIFICATION_ROUTE_AVAILABLE"
    )

    assert (
        summary[
            "next_stage"
        ]
        ==
        "M7.4C_CANDIDATE_TRF_EXTERNAL_QUANTIFICATION"
    )

    # ----------------------------------------------------------------------
    # Dataset-level
    # ----------------------------------------------------------------------

    by_key = {
        row[
            "dataset_key"
        ]:
        row
        for row in rows
    }

    assert set(
        by_key.keys()
    ) == {
        "m7_gse80400",
        "m7_gse229904",
        "m7_gse89193",
        "m7_gse290219",
    }

    # ----------------------------------------------------------------------
    # GSE80400
    # ----------------------------------------------------------------------

    row = by_key[
        "m7_gse80400"
    ]

    assert (
        row[
            "raw_sequence_public"
        ]
        is True
    )

    assert (
        row[
            "small_rna_library_supported"
        ]
        is True
    )

    assert (
        row[
            "candidate_length_compatible"
        ]
        is True
    )

    assert (
        row[
            "read_orientation_resolvable"
        ]
        is True
    )

    assert (
        row[
            "adapter_handling_resolvable"
        ]
        is True
    )

    assert (
        row[
            "technical_eligibility"
        ]
        is True
    )

    assert (
        row[
            "resolution_class"
        ]
        ==
        "TECHNICALLY_ELIGIBLE_WITH_CAVEATS"
    )

    assert (
        row[
            "technical_status"
        ]
        ==
        "TECHNICALLY_ELIGIBLE_WITH_CAVEATS"
    )

    # ----------------------------------------------------------------------
    # GSE229904
    # ----------------------------------------------------------------------

    row = by_key[
        "m7_gse229904"
    ]

    assert (
        row[
            "technical_eligibility"
        ]
        is True
    )

    assert (
        row[
            "candidate_length_compatible"
        ]
        is True
    )

    assert (
        row[
            "resolution_class"
        ]
        ==
        "TECHNICALLY_ELIGIBLE_WITH_CAVEATS"
    )

    # ----------------------------------------------------------------------
    # GSE89193
    # ----------------------------------------------------------------------

    row = by_key[
        "m7_gse89193"
    ]

    assert (
        row[
            "technical_eligibility"
        ]
        is True
    )

    assert (
        row[
            "ffpe"
        ]
        is True
    )

    assert (
        row[
            "resolution_class"
        ]
        ==
        "TECHNICALLY_ELIGIBLE_WITH_CAVEATS"
    )

    # ----------------------------------------------------------------------
    # GSE290219
    # ----------------------------------------------------------------------

    row = by_key[
        "m7_gse290219"
    ]

    assert (
        row[
            "technical_eligibility"
        ]
        is False
    )

    assert (
        row[
            "candidate_length_compatible"
        ]
        is None
    )

    assert (
        row[
            "read_orientation_resolvable"
        ]
        is False
    )

    assert (
        row[
            "resolution_class"
        ]
        ==
        "METADATA_INCOMPLETE"
    )

    assert (
        row[
            "technical_status"
        ]
        ==
        "TECHNICAL_ELIGIBILITY_INCOMPLETE_METADATA"
    )

    # ----------------------------------------------------------------------
    # No downstream analysis performed
    # ----------------------------------------------------------------------

    for row in rows:

        assert (
            row[
                "fastq_download_performed"
            ]
            is False
        )

        assert (
            row[
                "candidate_sequence_search_performed"
            ]
            is False
        )

        assert (
            row[
                "candidate_detected"
            ]
            is False
        )

        assert (
            row[
                "candidate_absent"
            ]
            is False
        )

        assert (
            row[
                "exact_match_counting_performed"
            ]
            is False
        )

        assert (
            row[
                "trf_quantification_performed"
            ]
            is False
        )

        assert (
            row[
                "clinical_association_performed"
            ]
            is False
        )


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute the M7.4B smoke test."""

    print()
    print("=" * 72)
    print("M7.4B EXTERNAL SMALL-RNA TECHNICAL ELIGIBILITY")
    print("=" * 72)

    runner = M74BExternalSmallRnaTechnicalEligibilityRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
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
    print("M7.4B RESULT SUMMARY")
    print("=" * 72)

    print(
        "Datasets assessed:                   "
        f"{summary['datasets_assessed']}"
    )

    print(
        "Technically eligible:               "
        f"{summary['datasets_technically_eligible']}"
    )

    print(
        "Fully eligible:                     "
        f"{summary['datasets_fully_eligible']}"
    )

    print(
        "Eligible with caveats:              "
        f"{summary['datasets_eligible_with_caveats']}"
    )

    print(
        "Metadata incomplete:                "
        f"{summary['datasets_metadata_incomplete']}"
    )

    print(
        "Not eligible:                       "
        f"{summary['datasets_not_eligible']}"
    )

    print()

    print(
        "Candidate tRF:                      "
        f"{summary['candidate_trf']}"
    )

    print(
        "Candidate sequence:                 "
        f"{summary['candidate_sequence']}"
    )

    print(
        "Candidate length:                   "
        f"{summary['candidate_length_nt']} nt"
    )

    print()

    print(
        "FASTQ download performed:           "
        f"{summary['fastq_download_performed']}"
    )

    print(
        "Candidate sequence search:          "
        f"{summary['candidate_sequence_search_performed']}"
    )

    print(
        "Candidate detection claimed:        "
        f"{summary['candidate_detection_claimed']}"
    )

    print(
        "Candidate absence claimed:          "
        f"{summary['candidate_absence_claimed']}"
    )

    print(
        "tRF quantification performed:       "
        f"{summary['trf_quantification_performed']}"
    )

    print(
        "Clinical association performed:     "
        f"{summary['clinical_association_performed']}"
    )

    print_dataset_results(
        report
    )

    print()
    print("Overall technical status:")
    print(
        "  "
        f"{summary['overall_technical_status']}"
    )

    print()

    print("Next stage:")
    print(
        "  "
        f"{summary['next_stage']}"
    )

    print()
    print("=" * 72)
    print("M7.4B SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()