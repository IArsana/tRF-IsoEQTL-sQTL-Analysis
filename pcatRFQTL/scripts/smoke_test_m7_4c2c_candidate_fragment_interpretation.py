"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4c2c_candidate_fragment_interpretation.py

Description:
    Smoke test for M7.4C.2C Candidate Fragment Interpretation.

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

from pcatrfqtl.analysis.m7.runners.interpret_candidate_fragment import (
    M74C2CCandidateFragmentInterpretationRunner,
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
    "m7_candidate_fragment_interpretation.yaml"
)


def validate_expected_state(
    report: dict[str, Any],
) -> None:
    """Validate expected M7.4C.2C interpretation."""

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.4C.2C"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "candidate_fragment_interpretation"
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
            "candidate_trf"
        ]
        ==
        "tRF-24-S3M8309N0Y"
    )

    assert (
        summary[
            "positive_reads"
        ]
        ==
        951
    )

    assert (
        summary[
            "positive_runs"
        ]
        ==
        11
    )

    assert (
        summary[
            "full_length_24nt_reads"
        ]
        ==
        0
    )

    assert (
        summary[
            "five_prime_boundary_reads"
        ]
        ==
        947
    )

    assert (
        summary[
            "three_prime_boundary_reads"
        ]
        ==
        0
    )

    assert (
        summary[
            "internal_motif_reads"
        ]
        ==
        4
    )

    expected_fraction = (
        947
        /
        951
    )

    assert abs(
        summary[
            "five_prime_boundary_fraction"
        ]
        -
        expected_fraction
    ) < 1e-12

    assert (
        summary[
            "direct_24nt_fragment_supported"
        ]
        is False
    )

    assert (
        summary[
            "longer_fragment_context_present"
        ]
        is True
    )

    assert (
        summary[
            "evidence_class"
        ]
        ==
        "CANDIDATE_5PRIME_SEQUENCE_WITH_LONGER_FRAGMENT_CONTEXT"
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
            "adapter_origin_assumed"
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

    assert (
        summary[
            "overall_status"
        ]
        ==
        "CANDIDATE_FRAGMENT_INTERPRETATION_COMPLETE"
    )

    assert (
        summary[
            "next_stage"
        ]
        ==
        "M7.4D_EXTERNAL_CLINICAL_ASSOCIATION_READINESS"
    )


def main() -> None:
    """Execute M7.4C.2C smoke test."""

    print()
    print("=" * 72)
    print("M7.4C.2C CANDIDATE FRAGMENT INTERPRETATION")
    print("=" * 72)

    runner = M74C2CCandidateFragmentInterpretationRunner(
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
    print("M7.4C.2C RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidate:                         "
        f"{summary['candidate_trf']}"
    )

    print(
        "Positive reads:                    "
        f"{summary['positive_reads']}"
    )

    print(
        "Positive runs:                     "
        f"{summary['positive_runs']}"
    )

    print()

    print(
        "Full-length 24-nt reads:           "
        f"{summary['full_length_24nt_reads']}"
    )

    print(
        "5' boundary reads:                 "
        f"{summary['five_prime_boundary_reads']}"
    )

    print(
        "3' boundary reads:                 "
        f"{summary['three_prime_boundary_reads']}"
    )

    print(
        "Internal motif reads:              "
        f"{summary['internal_motif_reads']}"
    )

    print()

    print(
        "5' boundary fraction:              "
        f"{summary['five_prime_boundary_fraction']:.6f}"
    )

    print(
        "Dominant context fraction:         "
        f"{summary['dominant_context_fraction']:.6f}"
    )

    print(
        "Dominant context runs:             "
        f"{summary['dominant_context_runs']}"
    )

    print()

    print(
        "Direct 24-nt fragment supported:   "
        f"{summary['direct_24nt_fragment_supported']}"
    )

    print(
        "Longer fragment context present:   "
        f"{summary['longer_fragment_context_present']}"
    )

    print()

    print("Evidence class:")
    print(
        "  "
        f"{summary['evidence_class']}"
    )

    print()

    print("Top flank architectures:")
    print()

    for row in report[
        "top_flank_architecture"
    ][
        :10
    ]:

        print(
            f"  rank={row['rank']} | "
            f"len={row['read_length_nt']} | "
            f"start={row['candidate_start_0based']} | "
            f"class={row['boundary_class']} | "
            f"left={row['left_flank_sequence']!r} | "
            f"right={row['right_flank_sequence']!r} | "
            f"n={row['occurrence_count']} | "
            f"fraction={row['occurrence_fraction']:.6f} | "
            f"runs={row['runs_observed']}"
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
    print("M7.4C.2C SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()