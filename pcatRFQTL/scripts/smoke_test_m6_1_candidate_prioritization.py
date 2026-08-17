"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_1_candidate_prioritization.py

Description:
    Smoke-test entry point for M6.1 final candidate prioritization.

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

from pcatrfqtl.analysis.m6.runners.prioritize_candidates import (
    M61CandidatePrioritizationRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m6_candidate_prioritization.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "candidate_prioritization"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "qc"
)


def validate_inputs() -> None:

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            f"M6.1 config not found: {CONFIG_PATH}"
        )


def print_header() -> None:

    print()
    print("=" * 72)
    print("M6.1 FINAL CANDIDATE PRIORITIZATION")
    print("=" * 72)

    print(
        f"Project root:       {ROOT}"
    )

    print(
        f"Config:             {CONFIG_PATH}"
    )

    print(
        f"Output directory:   {OUTPUT_DIRECTORY}"
    )

    print(
        f"QC directory:       {QC_DIRECTORY}"
    )

    print("=" * 72)


def print_candidate_priority(
    report: dict[str, Any],
) -> None:

    candidates = report.get(
        "candidate_priority",
        [],
    )

    print()
    print("Candidate priority:")
    print()

    for row in candidates:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    Priority class:          "
            f"{row['priority_class']}"
        )

        print(
            "    Disease evidence:        "
            f"{row['disease_evidence_class']}"
        )

        print(
            "    Regional genome-wide:    "
            f"{row['regional_genome_wide']}"
        )

        print(
            "    Recurrent suggestive:    "
            f"{row['recurrent_suggestive']}"
        )

        print(
            "    Moradi locus evidence:   "
            f"{row['moradi_locus_evidence']}"
        )

        print(
            "    Moradi proxy bridge:     "
            f"{row['moradi_direct_or_proxy_bridge']}"
        )

        print(
            "    Formal coloc ready:      "
            f"{row['formal_coloc_ready']}"
        )

        print(
            "    Continue to M6:          "
            f"{row['continue_to_m6_biological_interpretation']}"
        )

        print()


def print_summary(
    report: dict[str, Any],
) -> None:

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M6.1 RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidates assessed:               "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Candidates retained:               "
        f"{summary['candidates_retained']}"
    )

    print(
        "Highest-priority candidate:        "
        f"{summary['highest_priority_candidate']}"
    )

    print(
        "Highest-priority class:            "
        f"{summary['highest_priority_class']}"
    )

    print()

    print(
        "Formal coloc-ready candidates:     "
        f"{summary['formal_coloc_ready_candidates']}"
    )

    print(
        "Causal-claim-allowed candidates:   "
        f"{summary['causal_claim_allowed_candidates']}"
    )

    print()

    print(
        "Ranking method:"
    )

    print(
        f"  {summary['ranking_method']}"
    )

    print(
        "Additive score used:               "
        f"{summary['additive_score_used']}"
    )

    print(
        "Candidate exclusion performed:     "
        f"{summary['candidate_exclusion_performed']}"
    )

    print()

    print(
        "Next stage:"
    )

    print(
        f"  {summary['next_stage']}"
    )

    print_candidate_priority(
        report
    )

    print("=" * 72)


def main() -> None:

    validate_inputs()

    print_header()

    runner = M61CandidatePrioritizationRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print_summary(
        report
    )

    print()
    print("Full QC report:")
    print()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print()
    print("=" * 72)
    print("M6.1 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()