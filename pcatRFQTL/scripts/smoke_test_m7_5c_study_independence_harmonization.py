"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_5c_study_independence_harmonization.py

Description:
    Smoke test for M7.5C Study Independence & Harmonization Audit.

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

from pcatrfqtl.analysis.m7.runners.audit_study_independence_harmonization import (
    M75CStudyIndependenceHarmonizationRunner,
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
    "m7_study_independence_harmonization.yaml"
)


EXPECTED_RSIDS = {
    "rs10216902",
    "rs2328376",
    "rs1288100",
}


def main() -> None:
    """Run M7.5C smoke test."""

    print()
    print("=" * 88)
    print("M7.5C STUDY INDEPENDENCE & HARMONIZATION AUDIT")
    print("=" * 88)

    runner = M75CStudyIndependenceHarmonizationRunner(
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

    family_direction = report[
        "publication_family_direction"
    ]

    # ----------------------------------------------------------------------
    # Basic structure
    # ----------------------------------------------------------------------

    assert report[
        "milestone"
    ] == "M7.5C"

    assert report[
        "stage"
    ] == (
        "study_independence_and_harmonization_audit"
    )

    assert summary[
        "candidate_variants_assessed"
    ] == 3

    assert {
        row[
            "rsid"
        ]
        for row in candidates
    } == EXPECTED_RSIDS

    assert summary[
        "effect_rows_assessed"
    ] == 18

    # ----------------------------------------------------------------------
    # Known current computational state
    # ----------------------------------------------------------------------

    assert summary[
        "publication_families_identified"
    ] == 7

    assert summary[
        "study_pairs_assessed"
    ] == 78

    assert summary[
        "study_pairs_with_explicit_relatedness"
    ] == 10

    assert summary[
        "directionally_usable_effect_rows"
    ] == 16

    assert summary[
        "directionally_unusable_effect_rows"
    ] == 2

    # ----------------------------------------------------------------------
    # Scientific safeguards
    # ----------------------------------------------------------------------

    assert summary[
        "independent_replications_verified"
    ] == 0

    assert summary[
        "independent_replication_claimed"
    ] is False

    assert summary[
        "study_independence_verified"
    ] is False

    assert summary[
        "different_accession_used_as_independence"
    ] is False

    assert summary[
        "different_pubmed_used_as_proof_of_independence"
    ] is False

    assert summary[
        "cohort_nonoverlap_used_as_proof_of_independence"
    ] is False

    assert summary[
        "nonharmonised_sources_used_for_directional_confirmation"
    ] is False

    assert summary[
        "p_value_used_to_select_direction"
    ] is False

    assert summary[
        "meta_analysis_performed"
    ] is False

    assert summary[
        "pooled_effect_estimated"
    ] is False

    # ----------------------------------------------------------------------
    # Effect harmonization
    # ----------------------------------------------------------------------

    effect_rows = report[
        "candidate_effect_harmonization"
    ]

    usable_rows = [
        row
        for row in effect_rows
        if row[
            "effect_direction_usable"
        ]
    ]

    assert len(
        usable_rows
    ) == 16

    for row in usable_rows:

        assert row[
            "harmonised_source"
        ] is True

        assert row[
            "allele_orientation"
        ] in {
            "SAME",
            "REVERSED",
        }

        assert row[
            "aligned_direction"
        ] in {
            "POSITIVE",
            "NEGATIVE",
            "ZERO",
        }

    # ----------------------------------------------------------------------
    # Publication-family direction
    # ----------------------------------------------------------------------

    allowed_family_directions = {
        "POSITIVE",
        "NEGATIVE",
        "MIXED",
        None,
    }

    for row in family_direction:

        assert row[
            "publication_family_direction"
        ] in allowed_family_directions

        assert row[
            "independent_replication_unit"
        ] is False

        assert row[
            "sample_independence_verified"
        ] is False

    # ----------------------------------------------------------------------
    # Candidate scientific safeguards
    # ----------------------------------------------------------------------

    allowed_candidate_statuses = {
        "CROSS_PUBLICATION_DIRECTIONAL_HETEROGENEITY_INDEPENDENCE_UNRESOLVED",
        "WITHIN_PUBLICATION_DIRECTIONAL_HETEROGENEITY_INDEPENDENCE_UNRESOLVED",
        "WITHIN_AND_CROSS_PUBLICATION_DIRECTIONAL_HETEROGENEITY_INDEPENDENCE_UNRESOLVED",
        "MULTI_PUBLICATION_DIRECTIONALLY_CONCORDANT_INDEPENDENCE_UNRESOLVED",
        "EXACT_VARIANT_EVIDENCE_LIMITED_TO_SINGLE_PUBLICATION_FAMILY",
        "MULTI_PUBLICATION_DIRECTION_UNRESOLVED",
    }

    for row in candidates:

        assert row[
            "candidate_status"
        ] in allowed_candidate_statuses

        assert row[
            "independent_replication_verified"
        ] is False

        assert row[
            "study_independence_verified"
        ] is False

    # ----------------------------------------------------------------------
    # Current expected directional pattern
    # ----------------------------------------------------------------------

    candidate_map = {
        row[
            "rsid"
        ]:
            row
        for row in candidates
    }

    # rs10216902:
    # PMID 34737426 positive
    # PMID 39024449 negative
    # PMID 39358599 positive
    assert candidate_map[
        "rs10216902"
    ][
        "within_publication_directional_heterogeneity"
    ] is False

    assert candidate_map[
        "rs10216902"
    ][
        "cross_publication_directional_heterogeneity"
    ] is True

    assert candidate_map[
        "rs10216902"
    ][
        "candidate_status"
    ] == (
        "CROSS_PUBLICATION_DIRECTIONAL_HETEROGENEITY_"
        "INDEPENDENCE_UNRESOLVED"
    )

    # rs2328376:
    # MVP publication family contains both negative and positive effects.
    assert candidate_map[
        "rs2328376"
    ][
        "within_publication_directional_heterogeneity"
    ] is True

    assert candidate_map[
        "rs2328376"
    ][
        "candidate_status"
    ] == (
        "WITHIN_PUBLICATION_DIRECTIONAL_HETEROGENEITY_"
        "INDEPENDENCE_UNRESOLVED"
    )

    # rs1288100:
    # MVP publication family also contains both positive and negative effects.
    assert candidate_map[
        "rs1288100"
    ][
        "within_publication_directional_heterogeneity"
    ] is True

    assert candidate_map[
        "rs1288100"
    ][
        "candidate_status"
    ] == (
        "WITHIN_PUBLICATION_DIRECTIONAL_HETEROGENEITY_"
        "INDEPENDENCE_UNRESOLVED"
    )

    # ----------------------------------------------------------------------
    # Overall
    # ----------------------------------------------------------------------

    assert summary[
        "overall_status"
    ] == (
        "MULTI_STUDY_VARIANT_EVIDENCE_DIRECTIONALLY_"
        "HETEROGENEOUS_INDEPENDENCE_UNRESOLVED"
    )

    assert summary[
        "next_stage"
    ] == (
        "M7.6_INTEGRATED_EVIDENCE_RESOLUTION"
    )

    # ----------------------------------------------------------------------
    # Print
    # ----------------------------------------------------------------------

    print()
    print("=" * 88)
    print("M7.5C RESULT SUMMARY")
    print("=" * 88)

    print(
        "Publication families identified:           "
        f"{summary['publication_families_identified']}"
    )

    print(
        "Study pairs assessed:                      "
        f"{summary['study_pairs_assessed']}"
    )

    print(
        "Pairs with explicit relatedness:           "
        f"{summary['study_pairs_with_explicit_relatedness']}"
    )

    print(
        "Effect rows assessed:                      "
        f"{summary['effect_rows_assessed']}"
    )

    print(
        "Directionally usable rows:                 "
        f"{summary['directionally_usable_effect_rows']}"
    )

    print(
        "Directionally unusable rows:               "
        f"{summary['directionally_unusable_effect_rows']}"
    )

    print(
        "Candidates with within-family heterogeneity:"
        f" {summary['candidates_with_within_publication_heterogeneity']}"
    )

    print(
        "Candidates with cross-family heterogeneity: "
        f"{summary['candidates_with_cross_publication_heterogeneity']}"
    )

    print()
    print("Publication-family direction:")
    print()

    for row in family_direction:

        print(
            f"  {row['rsid']} | "
            f"{row['publication_family_id']}"
        )

        print(
            "    studies:     "
            f"{row['study_accessions']}"
        )

        print(
            "    directions:  "
            f"{row['observed_directions']}"
        )

        print(
            "    family dir:  "
            f"{row['publication_family_direction']}"
        )

        print()

    print("Candidate resolution:")
    print()

    for row in candidates:

        print(
            f"  {row['rsid']} | "
            f"{row['priority_class']}"
        )

        print(
            "    studies:                    "
            f"{row['studies_with_exact_variant']}"
        )

        print(
            "    publication families:       "
            f"{row['publication_families_with_exact_variant']}"
        )

        print(
            "    usable effect rows:         "
            f"{row['directionally_usable_rows']}"
        )

        print(
            "    positive families:          "
            f"{row['publication_families_direction_positive']}"
        )

        print(
            "    negative families:          "
            f"{row['publication_families_direction_negative']}"
        )

        print(
            "    mixed families:             "
            f"{row['publication_families_direction_mixed']}"
        )

        print(
            "    within-publication hetero:  "
            f"{row['within_publication_directional_heterogeneity']}"
        )

        print(
            "    cross-publication hetero:   "
            f"{row['cross_publication_directional_heterogeneity']}"
        )

        print(
            "    independence verified:      "
            f"{row['independent_replication_verified']}"
        )

        print(
            "    status:"
        )

        print(
            "      "
            f"{row['candidate_status']}"
        )

        print()

    print(
        "Independent replications verified:         "
        f"{summary['independent_replications_verified']}"
    )

    print(
        "Independent replication claimed:           "
        f"{summary['independent_replication_claimed']}"
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
    print("=" * 88)
    print("M7.5C SMOKE TEST PASSED")
    print("=" * 88)


if __name__ == "__main__":

    main()