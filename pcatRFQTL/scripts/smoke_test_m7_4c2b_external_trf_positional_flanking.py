"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4c2b_external_trf_positional_flanking.py

Description:
    Smoke test for M7.4C.2B External tRF Positional and Flanking Audit.

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

from pcatrfqtl.analysis.m7.runners.audit_external_trf_positional_flanking import (
    M74C2BExternalTrfPositionalFlankingRunner,
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
    "m7_external_trf_positional_flanking_audit.yaml"
)


def main() -> None:

    print()
    print("=" * 72)
    print("M7.4C.2B POSITIONAL & FLANKING CONTEXT AUDIT")
    print("=" * 72)

    runner = M74C2BExternalTrfPositionalFlankingRunner(
        project_root=ROOT,
        config_path=CONFIG,
    )

    report = runner.run()

    summary = report[
        "summary"
    ]

    assert (
        summary[
            "runs_audited"
        ]
        ==
        11
    )

    assert (
        summary[
            "total_reads"
        ]
        ==
        160162133
    )

    # ----------------------------------------------------------------------
    # Every read must fall into exactly one length class.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "reads_lt15"
        ]
        +
        summary[
            "reads_15_to_35"
        ]
        +
        summary[
            "reads_gt35"
        ]
        ==
        summary[
            "total_reads"
        ]
    )

    # ----------------------------------------------------------------------
    # Candidate was previously observed.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "motif_positive_reads_all_lengths"
        ]
        >
        0
    )

    # ----------------------------------------------------------------------
    # Searching all lengths must preserve at least the canonical observations.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "motif_positive_reads_all_lengths"
        ]
        >=
        summary[
            "motif_positive_reads_15_to_35"
        ]
    )

    # ----------------------------------------------------------------------
    # Important methodological safeguard:
    # do NOT require short or long reads to exist.
    #
    # Their number is empirical, not something the pipeline should force.
    # ----------------------------------------------------------------------

    assert (
        summary[
            "reads_lt15"
        ]
        >=
        0
    )

    assert (
        summary[
            "reads_gt35"
        ]
        >=
        0
    )

    assert (
        summary[
            "motif_positive_reads_lt15"
        ]
        >=
        0
    )

    assert (
        summary[
            "motif_positive_reads_gt35"
        ]
        >=
        0
    )

    # ----------------------------------------------------------------------
    # Scientific safeguards
    # ----------------------------------------------------------------------

    assert (
        summary[
            "all_length_exact_search_performed"
        ]
        is True
    )

    assert (
        summary[
            "canonical_window_retained"
        ]
        is True
    )

    assert (
        summary[
            "artificial_read_extension_performed"
        ]
        is False
    )

    assert (
        summary[
            "artificial_read_truncation_performed"
        ]
        is False
    )

    assert (
        summary[
            "approximate_matching_performed"
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
            "overall_status"
        ]
        ==
        "EXTERNAL_TRF_POSITIONAL_CONTEXT_AUDIT_COMPLETE"
    )

    print()
    print("=" * 72)
    print("M7.4C.2B RESULT SUMMARY")
    print("=" * 72)

    print(
        f"Runs audited:                     "
        f"{summary['runs_audited']}"
    )

    print(
        f"Total reads:                      "
        f"{summary['total_reads']}"
    )

    print()

    print(
        f"Reads <15 nt:                     "
        f"{summary['reads_lt15']}"
    )

    print(
        f"Reads 15-35 nt:                   "
        f"{summary['reads_15_to_35']}"
    )

    print(
        f"Reads >35 nt:                     "
        f"{summary['reads_gt35']}"
    )

    print()

    print(
        f"Motif-positive all lengths:       "
        f"{summary['motif_positive_reads_all_lengths']}"
    )

    print(
        f"Motif-positive <15 nt:            "
        f"{summary['motif_positive_reads_lt15']}"
    )

    print(
        f"Motif-positive 15-35 nt:          "
        f"{summary['motif_positive_reads_15_to_35']}"
    )

    print(
        f"Motif-positive >35 nt:            "
        f"{summary['motif_positive_reads_gt35']}"
    )

    print()

    print(
        f"Full-length 24-nt reads:          "
        f"{summary['full_length_24nt_reads']}"
    )

    print(
        f"Candidate at 5' boundary:         "
        f"{summary['five_prime_boundary_reads']}"
    )

    print(
        f"Candidate at 3' boundary:         "
        f"{summary['three_prime_boundary_reads']}"
    )

    print(
        f"Candidate internal:               "
        f"{summary['internal_motif_reads']}"
    )

    print(
        f"Multiple occurrences/read:        "
        f"{summary['multiple_occurrence_reads']}"
    )

    print()

    print(
        f"Unique positional/flank contexts: "
        f"{summary['unique_positional_flanking_contexts']}"
    )

    print()

    print("Interpretation:")
    print(
        "  "
        f"{summary['interpretation_status']}"
    )

    print()

    print("Next:")
    print(
        "  "
        f"{summary['next_stage']}"
    )

    print()
    print("Top positional/flanking contexts:")
    print()

    for row in report[
        "top_motif_contexts"
    ][
        :15
    ]:

        print(
            f"  len={row['read_length_nt']} | "
            f"start={row['candidate_start_0based']} | "
            f"class={row['boundary_class']} | "
            f"left={row['left_flank_sequence']!r} | "
            f"right={row['right_flank_sequence']!r} | "
            f"n={row['occurrence_count']} | "
            f"runs={row['runs_observed']}"
        )

    print()
    print("=" * 72)
    print("M7.4C.2B SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()