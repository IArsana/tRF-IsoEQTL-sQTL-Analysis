"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_3a_trf_sequence_resolution.py

Description:
    Smoke test for M7.3A Candidate tRF Sequence Resolution.

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

from pcatrfqtl.analysis.m7.runners.resolve_trf_sequences import (
    M73ATrfSequenceResolutionRunner,
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
    / "m7_trf_sequence_resolution.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "trf_sequence_resolution"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


# ============================================================================
# Display helpers
# ============================================================================


def _display_nullable_validation(
    value: Any,
) -> str:
    """Display nullable validation state clearly."""

    if value is None:

        return "NOT_ASSESSABLE"

    return str(
        value
    )


def print_candidate_resolution(
    report: dict[str, Any],
) -> None:
    """Print candidate-level tRF sequence resolution."""

    print()
    print("Candidate tRF sequence resolution:")
    print()

    for row in report[
        "candidate_trf_sequence_resolution"
    ]:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    tRF ID:                        "
            f"{row['trf_id']}"
        )

        print(
            "    Identifier source-confirmed:   "
            f"{row['exact_identifier_observed']}"
        )

        print(
            "    Exact sequence reported:       "
            f"{row['exact_sequence_reported']}"
        )

        print(
            "    Sequence:                      "
            f"{row['sequence_dna']}"
        )

        print(
            "    Identifier expected length:    "
            f"{row['identifier_expected_length']}"
        )

        print(
            "    Observed sequence length:       "
            f"{row['observed_sequence_length']}"
        )

        print(
            "    Sequence validation assessable:"
            f" {row['sequence_validation_assessable']}"
        )

        print(
            "    Alphabet valid:                "
            f"{_display_nullable_validation(row['sequence_alphabet_valid'])}"
        )

        print(
            "    Reported length consistent:     "
            f"{_display_nullable_validation(row['reported_length_consistent'])}"
        )

        print(
            "    Identifier length consistent:   "
            f"{_display_nullable_validation(row['identifier_length_consistent'])}"
        )

        print(
            "    Sequence validation passed:     "
            f"{row['sequence_validation_passed']}"
        )

        print(
            "    Sequence resolved:              "
            f"{row['sequence_resolved']}"
        )

        print(
            "    Quantification audit ready:     "
            f"{row['candidate_quantification_ready']}"
        )

        print(
            "    Resolution status:              "
            f"{row['sequence_resolution_status']}"
        )

        print()


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M7.3A smoke test."""

    print()
    print("=" * 72)
    print("M7.3A CANDIDATE tRF SEQUENCE RESOLUTION")
    print("=" * 72)

    runner = M73ATrfSequenceResolutionRunner(
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
    print("M7.3A RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidate tRFs assessed:                    "
        f"{summary['candidate_trfs_assessed']}"
    )

    print(
        "Identifiers source-confirmed:               "
        f"{summary['candidate_identifiers_source_confirmed']}"
    )

    print(
        "Sequences with assessable validation:       "
        f"{summary['candidate_sequences_with_assessable_validation']}"
    )

    print(
        "Candidate sequences resolved:               "
        f"{summary['candidate_sequences_resolved']}"
    )

    print(
        "Candidate sequences unresolved:             "
        f"{summary['candidate_sequences_unresolved']}"
    )

    print(
        "Candidates ready for quantification audit:  "
        f"{summary['candidates_ready_for_quantification_audit']}"
    )

    print()

    print(
        "Identifier decoding performed:              "
        f"{summary['identifier_decoding_performed']}"
    )

    print(
        "Sequence inference performed:               "
        f"{summary['sequence_inference_performed']}"
    )

    print(
        "Sequence substitution performed:            "
        f"{summary['sequence_substitution_performed']}"
    )

    print(
        "Unresolved interpreted as invalid:          "
        f"{summary['unresolved_sequence_interpreted_as_invalid']}"
    )

    print(
        "TCGA quantification performed:              "
        f"{summary['tcga_quantification_performed']}"
    )

    print()

    print(
        "Overall resolution:"
    )

    print(
        f"  {summary['overall_resolution_status']}"
    )

    print_candidate_resolution(
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
    print("M7.3A EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()