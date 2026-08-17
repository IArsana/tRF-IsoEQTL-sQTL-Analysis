"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_3b_trf_quantification_feasibility.py

Description:
    Smoke test for M7.3B tRF Quantification Feasibility Audit.

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

from pcatrfqtl.analysis.m7.runners.assess_trf_quantification_feasibility import (
    M73BTrfQuantificationFeasibilityRunner,
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
    / "m7_trf_quantification_feasibility.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "trf_quantification_feasibility"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


# ============================================================================
# Reporting
# ============================================================================


def print_candidate_feasibility(
    report: dict[str, Any],
) -> None:
    """Print candidate-specific quantification feasibility."""

    print()
    print("Candidate tRF quantification feasibility:")
    print()

    for row in report[
        "candidate_trf_quantification_feasibility"
    ]:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    tRF ID:                         "
            f"{row['trf_id']}"
        )

        print(
            "    Sequence resolved:              "
            f"{row['sequence_resolved']}"
        )

        print(
            "    Sequence:                       "
            f"{row['sequence_dna']}"
        )

        print(
            "    Aligned miRNA BAM files:        "
            f"{row['aligned_mirna_bam_file_count']}"
        )

        print(
            "    Cases with aligned BAM:         "
            f"{row['aligned_mirna_bam_case_count']}"
        )

        print(
            "    Controlled BAMs:                "
            f"{row['controlled_aligned_bam_count']}"
        )

        print(
            "    Open BAMs:                      "
            f"{row['open_aligned_bam_count']}"
        )

        print(
            "    Aligned resource available:     "
            f"{row['aligned_read_resource_available']}"
        )

        print(
            "    Controlled access required:     "
            f"{row['controlled_access_required']}"
        )

        print(
            "    Direct quantification ready:    "
            f"{row['direct_quantification_ready_now']}"
        )

        print(
            "    Quantification performed:       "
            f"{row['quantification_performed']}"
        )

        print(
            "    Status:"
        )

        print(
            "      "
            f"{row['candidate_quantification_feasibility_status']}"
        )

        print()


def main() -> None:
    """Execute M7.3B smoke test."""

    print()
    print("=" * 72)
    print("M7.3B tRF QUANTIFICATION FEASIBILITY AUDIT")
    print("=" * 72)

    runner = M73BTrfQuantificationFeasibilityRunner(
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
    print("M7.3B RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidate tRFs assessed:                    "
        f"{summary['candidate_trfs_assessed']}"
    )

    print(
        "Candidate sequences resolved:               "
        f"{summary['candidate_sequences_resolved']}"
    )

    print(
        "Candidates blocked by sequence resolution:  "
        f"{summary['candidates_blocked_by_sequence_resolution']}"
    )

    print()

    print(
        "GDC miRNA-Seq files observed:               "
        f"{summary['gdc_mirna_files_observed']}"
    )

    print(
        "Aligned miRNA BAM files:                    "
        f"{summary['aligned_mirna_bam_files']}"
    )

    print(
        "Processed miRNA expression files:           "
        f"{summary['processed_mirna_expression_files']}"
    )

    print(
        "Processed miRNA isoform files:              "
        f"{summary['processed_mirna_isoform_files']}"
    )

    print()

    print(
        "Candidates directly quantifiable now:       "
        f"{summary['candidates_directly_quantifiable_now']}"
    )

    print(
        "Candidates requiring controlled access:     "
        f"{summary['candidates_requiring_controlled_access']}"
    )

    print()

    print(
        "Processed miRNA used as tRF measurement:    "
        f"{summary['processed_mirna_used_as_trf_measurement']}"
    )

    print(
        "BAM download performed:                     "
        f"{summary['bam_download_performed']}"
    )

    print(
        "Read counting performed:                    "
        f"{summary['read_counting_performed']}"
    )

    print()

    print(
        "Overall feasibility:"
    )

    print(
        f"  {summary['overall_feasibility_status']}"
    )

    print_candidate_feasibility(
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
    print("M7.3B EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()