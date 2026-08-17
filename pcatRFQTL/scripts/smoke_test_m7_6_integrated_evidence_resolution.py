"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_6_integrated_evidence_resolution.py

Description:
    Smoke test for M7.6 Integrated Evidence Resolution.

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

from pcatrfqtl.analysis.m7.runners.resolve_integrated_evidence import (
    M76IntegratedEvidenceResolutionRunner,
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
    "m7_integrated_evidence_resolution.yaml"
)


EXPECTED_ORDER = [
    "rs10216902",
    "rs2328376",
    "rs1288100",
]


EXPECTED_CLASSES = {
    "rs10216902":
        "INTEGRATED_MULTI_DOMAIN_CONVERGENT_EVIDENCE_WITHOUT_CAUSAL_RESOLUTION",

    "rs2328376":
        "GENETIC_EVIDENCE_SUPPORTED_REGULATORY_AND_MOLECULAR_LINKS_UNRESOLVED",

    "rs1288100":
        "GENETIC_ASSOCIATION_CONTEXT_RETAINED_MECHANISTIC_LINK_UNRESOLVED",
}


def main() -> None:
    """Execute M7.6 smoke test."""

    print()
    print("=" * 96)
    print("M7.6 INTEGRATED EVIDENCE RESOLUTION")
    print("=" * 96)

    runner = M76IntegratedEvidenceResolutionRunner(
        project_root=ROOT,
        config_path=CONFIG,
    )

    report = runner.run()

    summary = report[
        "summary"
    ]

    candidates = report[
        "candidate_resolution"
    ]

    evidence = report[
        "candidate_evidence_matrix"
    ]

    domains = report[
        "evidence_domain_matrix"
    ]

    # ----------------------------------------------------------------------
    # Structure
    # ----------------------------------------------------------------------

    assert report[
        "milestone"
    ] == "M7.6"

    assert report[
        "stage"
    ] == (
        "integrated_evidence_resolution"
    )

    assert summary[
        "candidate_variants_assessed"
    ] == 3

    assert summary[
        "candidates_retained"
    ] == 3

    observed_order = [
        row[
            "rsid"
        ]
        for row in sorted(
            candidates,
            key=lambda row:
                row[
                    "priority_rank"
                ],
        )
    ]

    assert observed_order == EXPECTED_ORDER

    # ----------------------------------------------------------------------
    # Global safeguards
    # ----------------------------------------------------------------------

    assert summary[
        "independent_replications_verified"
    ] == 0

    assert summary[
        "causal_candidates"
    ] == 0

    assert summary[
        "candidate_priority_recalculated"
    ] is False

    assert summary[
        "arbitrary_additive_score_used"
    ] is False

    assert summary[
        "meta_analysis_performed"
    ] is False

    assert summary[
        "causal_inference_performed"
    ] is False

    assert summary[
        "candidates_with_tcga_direct_trf_quantification"
    ] == 0

    assert summary[
        "candidates_with_external_clinical_association"
    ] == 0

    assert summary[
        "candidates_with_formal_colocalization"
    ] == 0

    assert summary[
        "candidates_with_resolved_gene_assignment"
    ] == 0

    # ----------------------------------------------------------------------
    # Known domain state
    # ----------------------------------------------------------------------

    assert summary[
        "candidates_with_regulatory_feature_evidence"
    ] == 1

    assert summary[
        "candidates_with_resolved_trf_sequence"
    ] == 1

    assert summary[
        "candidates_with_external_molecular_evidence"
    ] == 1

    assert summary[
        "candidates_with_directional_heterogeneity"
    ] == 3

    # ----------------------------------------------------------------------
    # Candidate final classification
    # ----------------------------------------------------------------------

    candidate_map = {
        row[
            "rsid"
        ]:
            row
        for row in candidates
    }

    for rsid, expected in EXPECTED_CLASSES.items():

        assert candidate_map[
            rsid
        ][
            "integrated_evidence_class"
        ] == expected

        assert candidate_map[
            rsid
        ][
            "retain_candidate"
        ] is True

        assert candidate_map[
            rsid
        ][
            "independent_replication_verified"
        ] is False

        assert candidate_map[
            rsid
        ][
            "formal_colocalization_performed"
        ] is False

        assert candidate_map[
            rsid
        ][
            "clinical_validation_established"
        ] is False

        assert candidate_map[
            rsid
        ][
            "causal_relationship_established"
        ] is False

    # ----------------------------------------------------------------------
    # rs10216902 evidence architecture
    # ----------------------------------------------------------------------

    evidence_map = {
        row[
            "rsid"
        ]:
            row
        for row in evidence
    }

    rs102 = evidence_map[
        "rs10216902"
    ]

    assert rs102[
        "priority_class"
    ] == "HIGH_PRIORITY"

    assert rs102[
        "regulatory_feature_confirmed"
    ] is True

    assert rs102[
        "regulatory_feature_id"
    ] == "INT98200"

    assert rs102[
        "trf_exact_sequence_resolved"
    ] is True

    assert rs102[
        "trf_exact_sequence"
    ] == (
        "GTAGTCGTGGCCGAGTGGTTAAGG"
    )

    assert rs102[
        "external_molecular_evidence_available"
    ] is True

    assert rs102[
        "external_runs_assessed"
    ] == 11

    assert rs102[
        "external_runs_sequence_positive"
    ] == 11

    assert rs102[
        "external_sequence_positive_reads"
    ] == 951

    assert rs102[
        "external_five_prime_boundary_reads"
    ] == 947

    assert rs102[
        "standalone_mature_trf_confirmed"
    ] is False

    assert rs102[
        "cross_publication_directional_heterogeneity"
    ] is True

    # ----------------------------------------------------------------------
    # rs2328376
    # ----------------------------------------------------------------------

    rs232 = evidence_map[
        "rs2328376"
    ]

    assert rs232[
        "priority_class"
    ] == "INTERMEDIATE_PRIORITY"

    assert rs232[
        "trf_exact_sequence_resolved"
    ] is False

    assert rs232[
        "regulatory_feature_confirmed"
    ] is False

    assert rs232[
        "external_molecular_evidence_available"
    ] is False

    assert rs232[
        "within_publication_directional_heterogeneity"
    ] is True

    # ----------------------------------------------------------------------
    # rs1288100
    # ----------------------------------------------------------------------

    rs128 = evidence_map[
        "rs1288100"
    ]

    assert rs128[
        "priority_class"
    ] == "LOWER_PRIORITY"

    assert rs128[
        "trf_exact_sequence_resolved"
    ] is False

    assert rs128[
        "regulatory_feature_confirmed"
    ] is False

    assert rs128[
        "external_molecular_evidence_available"
    ] is False

    assert rs128[
        "within_publication_directional_heterogeneity"
    ] is True

    # ----------------------------------------------------------------------
    # Domain counts must remain descriptive
    # ----------------------------------------------------------------------

    for row in domains:

        assert row[
            "domain_count_is_descriptive_only"
        ] is True

        assert row[
            "domain_count_used_for_priority"
        ] is False

    # ----------------------------------------------------------------------
    # Overall
    # ----------------------------------------------------------------------

    assert summary[
        "overall_status"
    ] == (
        "INTEGRATED_EVIDENCE_RESOLUTION_COMPLETE_"
        "WITHOUT_CAUSAL_RESOLUTION"
    )

    assert summary[
        "next_stage"
    ] == (
        "MANUSCRIPT_INTEGRATED_RESULTS_SYNTHESIS"
    )

    # ----------------------------------------------------------------------
    # Print summary
    # ----------------------------------------------------------------------

    print()
    print("=" * 96)
    print("M7.6 RESULT SUMMARY")
    print("=" * 96)

    print(
        "Candidates assessed:                         "
        f"{summary['candidate_variants_assessed']}"
    )

    print(
        "Candidates retained:                         "
        f"{summary['candidates_retained']}"
    )

    print(
        "Regulatory-feature evidence:                 "
        f"{summary['candidates_with_regulatory_feature_evidence']}"
    )

    print(
        "Resolved tRF sequence:                       "
        f"{summary['candidates_with_resolved_trf_sequence']}"
    )

    print(
        "External molecular evidence:                 "
        f"{summary['candidates_with_external_molecular_evidence']}"
    )

    print(
        "Directional heterogeneity:                   "
        f"{summary['candidates_with_directional_heterogeneity']}"
    )

    print(
        "Direct TCGA tRF quantification:              "
        f"{summary['candidates_with_tcga_direct_trf_quantification']}"
    )

    print(
        "External clinical associations:              "
        f"{summary['candidates_with_external_clinical_association']}"
    )

    print(
        "Formal colocalization:                       "
        f"{summary['candidates_with_formal_colocalization']}"
    )

    print(
        "Resolved gene assignments:                   "
        f"{summary['candidates_with_resolved_gene_assignment']}"
    )

    print(
        "Independent replications verified:           "
        f"{summary['independent_replications_verified']}"
    )

    print(
        "Causal candidates:                           "
        f"{summary['causal_candidates']}"
    )

    print()
    print("Candidate integrated resolution:")
    print()

    for row in candidates:

        print(
            f"  {row['rsid']} | "
            f"{row['priority_class']}"
        )

        print(
            "    evidence domains: "
            f"{row['evidence_domains_supported']}"
        )

        print(
            "    class:"
        )

        print(
            "      "
            f"{row['integrated_evidence_class']}"
        )

        print(
            "    primary limitation:"
        )

        print(
            "      "
            f"{row['primary_limitation']}"
        )

        print(
            "    causal relationship established: "
            f"{row['causal_relationship_established']}"
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
    print("=" * 96)
    print("M7.6 SMOKE TEST PASSED")
    print("=" * 96)


if __name__ == "__main__":

    main()