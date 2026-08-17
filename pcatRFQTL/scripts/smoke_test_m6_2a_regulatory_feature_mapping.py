"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_2a_regulatory_feature_mapping.py

Description:
    Smoke-test entry point for M6.2A Candidate → Regulatory Feature
    Resolution.

    The smoke test verifies that candidate-feature assignments inherited
    from M6.1 are reproducibly supported by the original Moradi S1, S2,
    and S7 supplementary datasets.

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

from pcatrfqtl.analysis.m6.runners.map_regulatory_features import (
    M62ARegulatoryFeatureMappingRunner,
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
    / "m6_regulatory_feature_mapping.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "regulatory_feature_mapping"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "qc"
)


# ============================================================================
# Validation
# ============================================================================


def validate_inputs() -> None:
    """Validate M6.2A configuration."""

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            "M6.2A configuration not found: "
            f"{CONFIG_PATH}"
        )


# ============================================================================
# Header
# ============================================================================


def print_header() -> None:
    """Print M6.2A header."""

    print()
    print("=" * 72)
    print("M6.2A CANDIDATE → REGULATORY FEATURE RESOLUTION")
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


# ============================================================================
# Feature evidence
# ============================================================================


def print_feature_evidence(
    report: dict[str, Any],
) -> None:
    """Print regulatory-feature evidence."""

    rows = report.get(
        "feature_evidence",
        [],
    )

    if not rows:

        print()
        print(
            "No assigned regulatory features were available."
        )

        return

    print()
    print("Regulatory feature evidence:")
    print()

    for row in rows:

        print(
            f"  {row['regulatory_feature_id']}"
        )

        print(
            "    Feature class:             "
            f"{row['feature_class']}"
        )

        print(
            "    Regulatory scope:          "
            f"{row['regulatory_scope']}"
        )

        print()

        print(
            "    S1 cis-QTL support:        "
            f"{row['s1_qtl_support']}"
        )

        print(
            "    S1 QTL rows:               "
            f"{row['s1_qtl_rows']}"
        )

        print(
            "    S1 unique QTL rsIDs:       "
            f"{row['s1_unique_qtl_rsids']}"
        )

        print(
            "    S1 minimum p-value:        "
            f"{row['s1_minimum_p_value']}"
        )

        print(
            "    S1 minimum FDR:            "
            f"{row['s1_minimum_fdr']}"
        )

        print()

        print(
            "    S2 GWAS-LD context:        "
            f"{row['s2_gwas_ld_context']}"
        )

        print(
            "    S2 GWAS-LD rows:           "
            f"{row['s2_gwas_ld_rows']}"
        )

        print(
            "    S2 unique QTL rsIDs:       "
            f"{row['s2_unique_qtl_rsids']}"
        )

        print(
            "    S2 unique tag rsIDs:       "
            f"{row['s2_unique_tag_rsids']}"
        )

        print(
            "    S2 maximum LD:             "
            f"{row['s2_maximum_ld']}"
        )

        print()

        print(
            "    S7 differential support:  "
            f"{row['s7_differential_support']}"
        )

        print(
            "    Differential log2FC:       "
            f"{row['s7_log2_fold_change']}"
        )

        print(
            "    Differential p-value:      "
            f"{row['s7_p_value']}"
        )

        print(
            "    Differential padj:         "
            f"{row['s7_adjusted_p_value']}"
        )

        print()

        print(
            "    Parent gene:               "
            f"{row['parent_gene']}"
        )

        print(
            "    Parent gene status:        "
            f"{row['parent_gene_status']}"
        )

        print(
            "    Confirmation:              "
            f"{row['overall_feature_confirmation']}"
        )

        print()


# ============================================================================
# Candidate mappings
# ============================================================================


def print_candidate_mapping(
    report: dict[str, Any],
) -> None:
    """Print candidate → regulatory-feature resolution."""

    rows = report.get(
        "candidate_mapping",
        [],
    )

    print()
    print("Candidate regulatory-feature mapping:")
    print()

    for row in rows:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    Priority class:           "
            f"{row['priority_class']}"
        )

        print(
            "    Regulatory feature:       "
            f"{row['regulatory_feature_id']}"
        )

        print(
            "    Feature class:            "
            f"{row['feature_class']}"
        )

        print(
            "    S1 sQTL support:          "
            f"{row['s1_qtl_support']}"
        )

        print(
            "    S2 GWAS-LD context:       "
            f"{row['s2_gwas_ld_context']}"
        )

        print(
            "    S7 differential support: "
            f"{row['s7_differential_support']}"
        )

        print(
            "    Parent gene status:       "
            f"{row['parent_gene_status']}"
        )

        print(
            "    Resolution status:        "
            f"{row['resolution_status']}"
        )

        print(
            "    Continue to M6.2B:        "
            f"{row['continue_to_gene_annotation']}"
        )

        print()


# ============================================================================
# Summary
# ============================================================================


def print_summary(
    report: dict[str, Any],
) -> None:
    """Print M6.2A summary."""

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M6.2A RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidates assessed:                  "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Candidates with regulatory feature:   "
        f"{summary['candidates_with_regulatory_feature']}"
    )

    print(
        "Candidates without regulatory feature:"
        f" {summary['candidates_without_regulatory_feature']}"
    )

    print()

    print(
        "Unique regulatory features assessed:  "
        f"{summary['unique_regulatory_features_assessed']}"
    )

    print(
        "Fully confirmed regulatory features:  "
        f"{summary['fully_confirmed_regulatory_features']}"
    )

    print()

    print(
        "Parent genes resolved:                "
        f"{summary['features_with_parent_gene_resolved']}"
    )

    print(
        "Parent genes pending:                 "
        f"{summary['features_with_parent_gene_pending']}"
    )

    print()

    print(
        "Parent-gene inference performed:      "
        f"{summary['parent_gene_inference_performed']}"
    )

    print(
        "Nearest-gene mapping performed:       "
        f"{summary['nearest_gene_mapping_performed']}"
    )

    print(
        "Formal colocalization performed:      "
        f"{summary['formal_colocalization_performed']}"
    )

    print(
        "Causal inference performed:           "
        f"{summary['causal_inference_performed']}"
    )

    print()

    print(
        "Next stage:"
    )

    print(
        f"  {summary['next_stage']}"
    )

    print_feature_evidence(
        report
    )

    print_candidate_mapping(
        report
    )

    print("=" * 72)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M6.2A smoke test."""

    validate_inputs()

    print_header()

    runner = M62ARegulatoryFeatureMappingRunner(
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
    print("M6.2A EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()