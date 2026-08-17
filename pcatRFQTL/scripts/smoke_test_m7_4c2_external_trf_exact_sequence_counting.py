"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4c2_external_trf_exact_sequence_counting.py

Description:
    Smoke test for M7.4C.2 Exact Candidate tRF Sequence Counting.

    This smoke test validates the full 11-run GSE80400 analysis.

    It verifies:
        - all configured runs are processed;
        - FASTQ preprocessing accounting is internally consistent;
        - candidate-positive reads are decomposed into exact full-length
          and exact substring reads;
        - normalization metrics are non-negative;
        - exact sequence observation is supported by positive raw counts;
        - no approximate matching, parent-tRNA inference, genomic-origin
          inference, clinical association, or causal inference is performed.

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

from pcatrfqtl.analysis.m7.runners.count_external_trf_exact_sequence import (
    M74C2ExternalTrfExactSequenceRunner,
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
    /
    "configs"
    /
    "m7_external_trf_exact_sequence_counting.yaml"
)


EXPECTED_RUNS = {
    "SRR3400531",
    "SRR3400532",
    "SRR3400533",
    "SRR3400534",
    "SRR3400535",
    "SRR3400536",
    "SRR3400537",
    "SRR3400538",
    "SRR3400539",
    "SRR3400540",
    "SRR3400541",
}


# ============================================================================
# Validation
# ============================================================================


def validate_expected_state(
    report: dict[str, Any],
) -> None:
    """Validate expected full M7.4C.2 state."""

    # ----------------------------------------------------------------------
    # Identity
    # ----------------------------------------------------------------------

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.4C.2"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "external_trf_exact_candidate_sequence_counting"
    )

    summary = report[
        "summary"
    ]

    rows = report[
        "run_counts"
    ]

    # ----------------------------------------------------------------------
    # Full execution
    # ----------------------------------------------------------------------

    assert (
        summary[
            "execution_mode"
        ]
        ==
        "full"
    )

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
            "sra_study"
        ]
        ==
        "SRP073456"
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
    # All 11 runs
    # ----------------------------------------------------------------------

    assert (
        summary[
            "runs_processed"
        ]
        ==
        11
    )

    assert (
        len(
            rows
        )
        ==
        11
    )

    observed_runs = {
        row[
            "run_accession"
        ]
        for row in rows
    }

    assert (
        observed_runs
        ==
        EXPECTED_RUNS
    )

    assert (
        summary[
            "selected_run_set_complete"
        ]
        is True
    )

    assert (
        summary[
            "full_run_set_complete"
        ]
        is True
    )

    # ----------------------------------------------------------------------
    # Read-level totals
    # ----------------------------------------------------------------------

    assert (
        summary[
            "total_reads"
        ]
        >
        0
    )

    assert (
        summary[
            "usable_reads_after_preprocessing"
        ]
        >
        0
    )

    assert (
        0.0
        <=
        summary[
            "aggregate_adapter_trimmed_fraction"
        ]
        <=
        1.0
    )

    assert (
        0.0
        <=
        summary[
            "aggregate_usable_read_fraction"
        ]
        <=
        1.0
    )

    # ----------------------------------------------------------------------
    # Candidate counts
    # ----------------------------------------------------------------------

    assert (
        summary[
            "exact_candidate_positive_reads"
        ]
        >=
        0
    )

    assert (
        summary[
            "exact_full_length_candidate_reads"
        ]
        >=
        0
    )

    assert (
        summary[
            "exact_candidate_substring_reads"
        ]
        >=
        0
    )

    assert (
        summary[
            "exact_candidate_positive_reads"
        ]
        ==
        (
            summary[
                "exact_full_length_candidate_reads"
            ]
            +
            summary[
                "exact_candidate_substring_reads"
            ]
        )
    )

    # ----------------------------------------------------------------------
    # Pilot already demonstrated exact candidate sequence observation.
    #
    # Because SRR3400531 is part of the full run set and its FASTQ is reused,
    # the full analysis must retain at least that source-supported observation.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "exact_candidate_positive_reads"
        ]
        >
        0
    )

    assert (
        summary[
            "candidate_exact_match_observed"
        ]
        is True
    )

    assert (
        summary[
            "runs_with_exact_match"
        ]
        >=
        1
    )

    # ----------------------------------------------------------------------
    # Normalization
    # ----------------------------------------------------------------------

    assert (
        summary[
            "combined_candidate_cpm"
        ]
        >=
        0
    )

    assert (
        summary[
            "combined_full_length_candidate_cpm"
        ]
        >=
        0
    )

    assert (
        summary[
            "combined_substring_candidate_cpm"
        ]
        >=
        0
    )

    assert (
        summary[
            "combined_candidate_fraction"
        ]
        >=
        0
    )

    assert (
        summary[
            "run_candidate_cpm_min"
        ]
        >=
        0
    )

    assert (
        summary[
            "run_candidate_cpm_max"
        ]
        >=
        summary[
            "run_candidate_cpm_min"
        ]
    )

    assert (
        summary[
            "run_candidate_cpm_median"
        ]
        >=
        0
    )

    assert (
        summary[
            "run_candidate_cpm_mean"
        ]
        >=
        0
    )

    # ----------------------------------------------------------------------
    # Interpretation
    # ----------------------------------------------------------------------

    assert (
        summary[
            "exact_sequence_evidence_only"
        ]
        is True
    )

    assert (
        summary[
            "zero_count_interpreted_as_biological_absence"
        ]
        is False
    )

    # ----------------------------------------------------------------------
    # Safeguards
    # ----------------------------------------------------------------------

    assert (
        summary[
            "approximate_matching_performed"
        ]
        is False
    )

    assert (
        summary[
            "mismatches_allowed"
        ]
        is False
    )

    assert (
        summary[
            "indels_allowed"
        ]
        is False
    )

    assert (
        summary[
            "reverse_complement_search_performed"
        ]
        is False
    )

    assert (
        summary[
            "parent_trna_inference_performed"
        ]
        is False
    )

    assert (
        summary[
            "genomic_origin_inference_performed"
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
            "causal_inference_performed"
        ]
        is False
    )

    # ----------------------------------------------------------------------
    # Overall status
    # ----------------------------------------------------------------------

    assert (
        summary[
            "overall_status"
        ]
        ==
        "EXTERNAL_CANDIDATE_TRF_SEQUENCE_OBSERVED"
    )

    assert (
        summary[
            "next_stage"
        ]
        ==
        "M7.4D_EXTERNAL_CLINICAL_ASSOCIATION_READINESS"
    )

    # ----------------------------------------------------------------------
    # Run-level consistency
    # ----------------------------------------------------------------------

    for row in rows:

        assert (
            row[
                "total_reads"
            ]
            >
            0
        )

        assert (
            row[
                "usable_reads"
            ]
            >=
            0
        )

        accounted_reads = (
            row[
                "usable_reads"
            ]
            +
            row[
                "empty_after_trim"
            ]
            +
            row[
                "reads_shorter_than_min"
            ]
            +
            row[
                "reads_longer_than_max"
            ]
        )

        assert (
            accounted_reads
            ==
            row[
                "total_reads"
            ]
        )

        assert (
            row[
                "exact_candidate_positive_reads"
            ]
            ==
            (
                row[
                    "exact_full_length_candidate_reads"
                ]
                +
                row[
                    "exact_candidate_substring_reads"
                ]
            )
        )

        assert (
            0.0
            <=
            row[
                "adapter_trimmed_fraction"
            ]
            <=
            1.0
        )

        assert (
            0.0
            <=
            row[
                "usable_read_fraction"
            ]
            <=
            1.0
        )

        assert (
            row[
                "candidate_cpm"
            ]
            >=
            0
        )

        assert (
            row[
                "approximate_matching_performed"
            ]
            is False
        )

        assert (
            row[
                "mismatches_allowed"
            ]
            is False
        )

        assert (
            row[
                "indels_allowed"
            ]
            is False
        )

        assert (
            row[
                "reverse_complement_search_performed"
            ]
            is False
        )

        assert (
            row[
                "parent_trna_inference_performed"
            ]
            is False
        )

        assert (
            row[
                "genomic_origin_inference_performed"
            ]
            is False
        )

        assert (
            row[
                "clinical_association_performed"
            ]
            is False
        )

        assert (
            row[
                "causal_inference_performed"
            ]
            is False
        )


# ============================================================================
# Display
# ============================================================================


def print_run_results(
    report: dict[str, Any],
) -> None:
    """Print manuscript-useful per-run QC."""

    print()
    print("Per-run exact-sequence results:")
    print()

    for row in report[
        "run_counts"
    ]:

        print(
            f"  {row['run_accession']}"
        )

        print(
            "    Total reads:                "
            f"{row['total_reads']}"
        )

        print(
            "    Adapter-trimmed reads:      "
            f"{row['reads_adapter_trimmed']}"
        )

        print(
            "    Adapter-trimmed fraction:   "
            f"{row['adapter_trimmed_fraction']:.6f}"
        )

        print(
            "    Usable reads:               "
            f"{row['usable_reads']}"
        )

        print(
            "    Usable-read fraction:       "
            f"{row['usable_read_fraction']:.6f}"
        )

        print(
            "    Exact-positive reads:       "
            f"{row['exact_candidate_positive_reads']}"
        )

        print(
            "    Full-length exact reads:    "
            f"{row['exact_full_length_candidate_reads']}"
        )

        print(
            "    Exact substring reads:      "
            f"{row['exact_candidate_substring_reads']}"
        )

        print(
            "    Candidate CPM:              "
            f"{row['candidate_cpm']:.6f}"
        )

        print(
            "    Full-length CPM:            "
            f"{row['full_length_candidate_cpm']:.6f}"
        )

        print(
            "    Substring CPM:              "
            f"{row['substring_candidate_cpm']:.6f}"
        )

        print(
            "    Exact sequence observed:    "
            f"{row['candidate_exact_match_observed']}"
        )

        print(
            "    Status:"
        )

        print(
            "      "
            f"{row['run_status']}"
        )

        print()


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute full M7.4C.2 smoke test."""

    print()
    print("=" * 72)
    print("M7.4C.2 EXACT CANDIDATE tRF SEQUENCE COUNTING — FULL GSE80400")
    print("=" * 72)

    runner = M74C2ExternalTrfExactSequenceRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
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
    print("M7.4C.2 FULL RESULT")
    print("=" * 72)

    print(
        "Dataset:                            "
        f"{summary['geo_accession']}"
    )

    print(
        "SRA study:                          "
        f"{summary['sra_study']}"
    )

    print(
        "Candidate:                          "
        f"{summary['candidate_trf']}"
    )

    print(
        "Sequence:                           "
        f"{summary['candidate_sequence']}"
    )

    print()

    print(
        "Runs processed:                     "
        f"{summary['runs_processed']}"
    )

    print(
        "Runs with exact match:              "
        f"{summary['runs_with_exact_match']}"
    )

    print(
        "Runs without exact match:           "
        f"{summary['runs_without_exact_match']}"
    )

    print(
        "Full run set complete:              "
        f"{summary['full_run_set_complete']}"
    )

    print()

    print(
        "Total reads:                        "
        f"{summary['total_reads']}"
    )

    print(
        "Adapter-trimmed reads:              "
        f"{summary['reads_adapter_trimmed']}"
    )

    print(
        "Adapter-trimmed fraction:           "
        f"{summary['aggregate_adapter_trimmed_fraction']:.6f}"
    )

    print(
        "Reads shorter than minimum:         "
        f"{summary['reads_shorter_than_min']}"
    )

    print(
        "Reads longer than maximum:          "
        f"{summary['reads_longer_than_max']}"
    )

    print(
        "Usable reads:                       "
        f"{summary['usable_reads_after_preprocessing']}"
    )

    print(
        "Usable-read fraction:               "
        f"{summary['aggregate_usable_read_fraction']:.6f}"
    )

    print()

    print(
        "Exact candidate-positive reads:     "
        f"{summary['exact_candidate_positive_reads']}"
    )

    print(
        "Exact full-length candidate reads:  "
        f"{summary['exact_full_length_candidate_reads']}"
    )

    print(
        "Exact candidate substring reads:    "
        f"{summary['exact_candidate_substring_reads']}"
    )

    print()

    print(
        "Combined candidate CPM:             "
        f"{summary['combined_candidate_cpm']:.6f}"
    )

    print(
        "Combined full-length CPM:           "
        f"{summary['combined_full_length_candidate_cpm']:.6f}"
    )

    print(
        "Combined substring CPM:             "
        f"{summary['combined_substring_candidate_cpm']:.6f}"
    )

    print()

    print(
        "Run CPM minimum:                    "
        f"{summary['run_candidate_cpm_min']:.6f}"
    )

    print(
        "Run CPM median:                     "
        f"{summary['run_candidate_cpm_median']:.6f}"
    )

    print(
        "Run CPM mean:                       "
        f"{summary['run_candidate_cpm_mean']:.6f}"
    )

    print(
        "Run CPM maximum:                    "
        f"{summary['run_candidate_cpm_max']:.6f}"
    )

    print()

    print(
        "Candidate exact match observed:     "
        f"{summary['candidate_exact_match_observed']}"
    )

    print()

    print(
        "Approximate matching:               "
        f"{summary['approximate_matching_performed']}"
    )

    print(
        "Mismatches allowed:                 "
        f"{summary['mismatches_allowed']}"
    )

    print(
        "Indels allowed:                     "
        f"{summary['indels_allowed']}"
    )

    print(
        "Parent tRNA inference:              "
        f"{summary['parent_trna_inference_performed']}"
    )

    print(
        "Genomic-origin inference:           "
        f"{summary['genomic_origin_inference_performed']}"
    )

    print(
        "Clinical association:               "
        f"{summary['clinical_association_performed']}"
    )

    print(
        "Causal inference:                   "
        f"{summary['causal_inference_performed']}"
    )

    print_run_results(
        report
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
    print("M7.4C.2 FULL SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()