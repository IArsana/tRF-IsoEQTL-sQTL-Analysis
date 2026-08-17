"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_3_candidate_gene_integration.py

Description:
    Smoke test for M6.3 Candidate → Gene Integration.

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

from pcatrfqtl.analysis.m6.runners.integrate_candidate_genes import (
    M63CandidateGeneIntegrationRunner,
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
    / "m6_candidate_gene_integration.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "candidate_gene_integration"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "qc"
)


# ============================================================================
# Reporting
# ============================================================================


def print_candidates(
    report: dict[str, Any],
) -> None:
    """Print candidate-level gene integration."""

    rows = report[
        "final_candidate_gene_integration"
    ]

    print()
    print("Candidate → gene integration:")
    print()

    for row in rows:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    Priority class:                  "
            f"{row['priority_class']}"
        )

        print(
            "    Regulatory feature:              "
            f"{row['regulatory_feature_id']}"
        )

        print(
            "    Regulatory feature supported:    "
            f"{row['regulatory_feature_supported']}"
        )

        print(
            "    Integrated gene:                 "
            f"{row['integrated_gene']}"
        )

        print(
            "    Gene resolved:                   "
            f"{row['integrated_gene_resolved']}"
        )

        print(
            "    Assignment method:               "
            f"{row['gene_assignment_method']}"
        )

        print(
            "    Gene resolution status:          "
            f"{row['gene_resolution_status']}"
        )

        print(
            "    Candidate retained:              "
            f"{row['candidate_retained']}"
        )

        print(
            "    Biological interpretation:       "
            f"{row['biological_interpretation_status']}"
        )

        print(
            "    Causal claim allowed:            "
            f"{row['causal_claim_allowed']}"
        )

        print()


def print_summary(
    report: dict[str, Any],
) -> None:
    """Print M6.3 summary."""

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M6.3 RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidates assessed:                         "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Candidates with resolved gene:               "
        f"{summary['candidates_with_resolved_gene']}"
    )

    print(
        "Supported feature / unresolved gene:         "
        f"{summary['candidates_with_supported_feature_unresolved_gene']}"
    )

    print(
        "Candidates without defensible gene assignment:"
        f" {summary['candidates_without_defensible_gene_assignment']}"
    )

    print(
        "Candidates retained:                         "
        f"{summary['candidates_retained']}"
    )

    print()

    print(
        "Gene assignment created in M6.3:             "
        f"{summary['candidate_gene_assignment_created_in_m6_3']}"
    )

    print(
        "Nearest-gene mapping performed:              "
        f"{summary['nearest_gene_assignment_performed']}"
    )

    print(
        "SNP-proximity mapping performed:             "
        f"{summary['snp_proximity_gene_assignment_performed']}"
    )

    print(
        "Formal colocalization performed:             "
        f"{summary['formal_colocalization_performed']}"
    )

    print(
        "Causal inference performed:                  "
        f"{summary['causal_inference_performed']}"
    )

    print_candidates(
        report
    )

    print(
        "Next stage:"
    )

    print(
        f"  {summary['next_stage']}"
    )

    print()
    print("=" * 72)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M6.3 smoke test."""

    print()
    print("=" * 72)
    print("M6.3 CANDIDATE → GENE INTEGRATION")
    print("=" * 72)

    runner = M63CandidateGeneIntegrationRunner(
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
            allow_nan=False,
        )
    )

    print()
    print("=" * 72)
    print("M6.3 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()