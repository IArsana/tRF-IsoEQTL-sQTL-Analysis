"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_4b_pathway_readiness.py

Description:
    Smoke test for M6.4B Pathway Readiness Assessment.

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

from pcatrfqtl.analysis.m6.runners.assess_pathway_readiness import (
    M64BPathwayReadinessRunner,
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
    / "m6_pathway_readiness.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "pathway_readiness"
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


def print_candidate_readiness(
    report: dict[str, Any],
) -> None:
    """Print candidate-level pathway readiness."""

    rows = report.get(
        "candidate_pathway_readiness",
        [],
    )

    print()
    print("Candidate pathway readiness:")
    print()

    if not rows:

        print(
            "  No candidates available."
        )

        return

    for row in rows:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    Priority class:               "
            f"{row['priority_class']}"
        )

        print(
            "    Regulatory feature:           "
            f"{row['regulatory_feature_id']}"
        )

        print(
            "    Integrated gene:              "
            f"{row['integrated_gene']}"
        )

        print(
            "    Gene resolved:                "
            f"{row['integrated_gene_resolved']}"
        )

        print(
            "    Assignment method:            "
            f"{row['gene_assignment_method']}"
        )

        print(
            "    Gene provenance supported:    "
            f"{row['gene_provenance_supported']}"
        )

        print(
            "    Eligible for pathway set:      "
            f"{row['eligible_for_pathway_gene_set']}"
        )

        print(
            "    Exclusion reason:              "
            f"{row['pathway_exclusion_reason']}"
        )

        print()


def print_gene_set(
    report: dict[str, Any],
) -> None:
    """Print defensible pathway gene set."""

    rows = report.get(
        "defensible_gene_set",
        [],
    )

    print()
    print("Defensible pathway gene set:")
    print()

    if not rows:

        print(
            "  EMPTY"
        )

        return

    for row in rows:

        print(
            f"  {row['gene']}"
        )

        print(
            "    Candidates:          "
            f"{row['candidate_rsids']}"
        )

        print(
            "    Assignment methods:  "
            f"{row['assignment_methods']}"
        )

        print()


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M6.4B smoke test."""

    print()
    print("=" * 72)
    print("M6.4B PATHWAY READINESS ASSESSMENT")
    print("=" * 72)

    runner = M64BPathwayReadinessRunner(
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
    print("M6.4B RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidates assessed:                          "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Resolved candidate genes:                     "
        f"{summary['resolved_candidate_genes']}"
    )

    print(
        "Eligible candidate genes:                     "
        f"{summary['eligible_candidate_genes']}"
    )

    print(
        "Unique defensible genes:                      "
        f"{summary['unique_defensible_genes']}"
    )

    print(
        "Minimum unique genes required:                "
        f"{summary['minimum_unique_genes_required']}"
    )

    print()

    print(
        "Candidates with unresolved regulatory gene:   "
        f"{summary['candidates_with_unresolved_regulatory_feature_gene']}"
    )

    print(
        "Candidates without defensible gene assignment:"
        f" {summary['candidates_without_defensible_gene_assignment']}"
    )

    print()

    print(
        "Pathway enrichment ready:                     "
        f"{summary['pathway_enrichment_ready']}"
    )

    print(
        "Pathway readiness status:"
    )

    print(
        f"  {summary['pathway_readiness_status']}"
    )

    print()

    print(
        "Enrichment performed:                         "
        f"{summary['enrichment_performed']}"
    )

    print(
        "Gene set imputed:                             "
        f"{summary['gene_set_imputed']}"
    )

    print(
        "Nearest gene added:                           "
        f"{summary['nearest_gene_added']}"
    )

    print(
        "SNP-proximity gene added:                     "
        f"{summary['snp_proximity_gene_added']}"
    )

    print(
        "Cis-proximity gene added:                     "
        f"{summary['cis_proximity_gene_added']}"
    )

    print(
        "Unresolved feature converted to gene:         "
        f"{summary['unresolved_feature_converted_to_gene']}"
    )

    print(
        "Pathway claim generated:                      "
        f"{summary['pathway_claim_generated']}"
    )

    print_candidate_readiness(
        report
    )

    print_gene_set(
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
    print("M6.4B EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()