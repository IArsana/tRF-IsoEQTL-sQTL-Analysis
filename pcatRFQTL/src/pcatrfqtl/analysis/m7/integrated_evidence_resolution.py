"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/integrated_evidence_resolution.py

Description:
    Core logic for M7.6 Integrated Evidence Resolution.

    M7.6 performs structured synthesis of locked M6-M7 evidence.
    It performs no new statistical inference.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class M76Result:
    """Container for M7.6 outputs."""

    candidate_evidence_matrix: pd.DataFrame
    evidence_domain_matrix: pd.DataFrame
    candidate_limitation_matrix: pd.DataFrame
    candidate_resolution: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _candidate_map(
    config: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Return candidate config keyed by rsID."""

    return {
        str(
            candidate["rsid"]
        ):
            candidate
        for candidate in config[
            "candidates"
        ]
    }


def _ordered_candidates(
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    """Return candidates in locked M6 priority order."""

    return sorted(
        config[
            "candidates"
        ],
        key=lambda row:
            int(
                row[
                    "priority_rank"
                ]
            ),
    )


# ============================================================================
# Integrated evidence matrix
# ============================================================================


def build_candidate_evidence_matrix(
    *,
    m7_5b_resolution: pd.DataFrame,
    m7_5c_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one integrated evidence row per candidate."""

    config_map = _candidate_map(
        config
    )

    rows: list[dict[str, Any]] = []

    for candidate in _ordered_candidates(
        config
    ):

        rsid = str(
            candidate[
                "rsid"
            ]
        )

        if rsid not in config_map:

            raise RuntimeError(
                f"Candidate config missing: {rsid}"
            )

        m75b = m7_5b_resolution.loc[
            m7_5b_resolution[
                "rsid"
            ].astype(str)
            ==
            rsid
        ]

        m75c = m7_5c_resolution.loc[
            m7_5c_resolution[
                "rsid"
            ].astype(str)
            ==
            rsid
        ]

        if len(
            m75b
        ) != 1:

            raise RuntimeError(
                f"M7.5B candidate resolution not unique: {rsid}"
            )

        if len(
            m75c
        ) != 1:

            raise RuntimeError(
                f"M7.5C candidate resolution not unique: {rsid}"
            )

        m75b_row = m75b.iloc[
            0
        ]

        m75c_row = m75c.iloc[
            0
        ]

        disease = candidate[
            "disease_evidence"
        ]

        regulatory = candidate[
            "regulatory_evidence"
        ]

        trf = candidate[
            "trf_evidence"
        ]

        external = candidate[
            "external_molecular_evidence"
        ]

        tcga = candidate[
            "tcga_validation"
        ]

        clinical = candidate[
            "external_clinical_validation"
        ]

        coloc = candidate[
            "formal_colocalization"
        ]

        gene = candidate[
            "gene_assignment"
        ]

        rows.append(
            {
                "rsid":
                    rsid,

                "priority_rank":
                    int(
                        candidate[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        candidate[
                            "priority_class"
                        ]
                    ),

                "trf_id":
                    str(
                        candidate[
                            "trf_id"
                        ]
                    ),

                # ----------------------------------------------------------
                # Disease / genetic
                # ----------------------------------------------------------

                "m5_disease_class":
                    str(
                        disease[
                            "m5_class"
                        ]
                    ),

                "locus_level_disease_evidence":
                    bool(
                        disease[
                            "locus_level_disease_evidence"
                        ]
                    ),

                "exact_variant_full_sumstats_evidence":
                    bool(
                        disease[
                            "exact_variant_full_sumstats_evidence"
                        ]
                    ),

                "full_sumstats_studies_with_exact_variant":
                    int(
                        m75b_row[
                            "studies_with_exact_variant"
                        ]
                    ),

                "full_sumstats_nonreused_studies_with_exact_variant":
                    int(
                        m75b_row[
                            "nonreused_studies_with_exact_variant"
                        ]
                    ),

                "publication_families_with_exact_variant":
                    int(
                        m75c_row[
                            "publication_families_with_exact_variant"
                        ]
                    ),

                # ----------------------------------------------------------
                # Direction / independence
                # ----------------------------------------------------------

                "within_publication_directional_heterogeneity":
                    bool(
                        m75c_row[
                            "within_publication_directional_heterogeneity"
                        ]
                    ),

                "cross_publication_directional_heterogeneity":
                    bool(
                        m75c_row[
                            "cross_publication_directional_heterogeneity"
                        ]
                    ),

                "gwas_direction_resolution_status":
                    str(
                        m75c_row[
                            "candidate_status"
                        ]
                    ),

                "independent_replication_verified":
                    False,

                # ----------------------------------------------------------
                # Regulatory
                # ----------------------------------------------------------

                "retained_regulatory_feature":
                    bool(
                        regulatory[
                            "retained_regulatory_feature"
                        ]
                    ),

                "regulatory_feature_id":
                    regulatory[
                        "regulatory_feature_id"
                    ],

                "regulatory_feature_type":
                    regulatory[
                        "regulatory_feature_type"
                    ],

                "regulatory_feature_confirmed":
                    bool(
                        regulatory[
                            "regulatory_feature_confirmed"
                        ]
                    ),

                "regulatory_parent_gene_resolved":
                    bool(
                        regulatory[
                            "parent_gene_resolved"
                        ]
                    ),

                "regulatory_parent_gene_status":
                    str(
                        regulatory[
                            "parent_gene_status"
                        ]
                    ),

                # ----------------------------------------------------------
                # tRF
                # ----------------------------------------------------------

                "trf_identifier_resolved":
                    bool(
                        trf[
                            "identifier_resolved"
                        ]
                    ),

                "trf_exact_sequence_resolved":
                    bool(
                        trf[
                            "exact_sequence_resolved"
                        ]
                    ),

                "trf_exact_sequence":
                    trf[
                        "exact_sequence"
                    ],

                "trf_sequence_length_nt":
                    trf[
                        "sequence_length_nt"
                    ],

                "parent_trna_resolved":
                    bool(
                        trf[
                            "parent_trna_resolved"
                        ]
                    ),

                "trf_class_resolved":
                    bool(
                        trf[
                            "trf_class_resolved"
                        ]
                    ),

                # ----------------------------------------------------------
                # External molecular
                # ----------------------------------------------------------

                "external_molecular_evidence_available":
                    bool(
                        external[
                            "available"
                        ]
                    ),

                "external_molecular_dataset":
                    external[
                        "dataset"
                    ],

                "external_runs_assessed":
                    int(
                        external[
                            "runs_assessed"
                        ]
                    ),

                "external_runs_sequence_positive":
                    int(
                        external[
                            "runs_sequence_positive"
                        ]
                    ),

                "external_sequence_positive_reads":
                    int(
                        external[
                            "total_sequence_positive_reads"
                        ]
                    ),

                "external_five_prime_boundary_reads":
                    int(
                        external[
                            "five_prime_boundary_reads"
                        ]
                    ),

                "standalone_mature_trf_confirmed":
                    bool(
                        external[
                            "direct_mature_24nt_trf_confirmed"
                        ]
                    ),

                "external_molecular_evidence_class":
                    str(
                        external[
                            "evidence_class"
                        ]
                    ),

                # ----------------------------------------------------------
                # TCGA
                # ----------------------------------------------------------

                "tcga_exact_sequence_ready":
                    bool(
                        tcga[
                            "exact_sequence_ready"
                        ]
                    ),

                "tcga_direct_trf_quantification_performed":
                    bool(
                        tcga[
                            "direct_tcga_trf_quantification_performed"
                        ]
                    ),

                "tcga_validation_status":
                    str(
                        tcga[
                            "status"
                        ]
                    ),

                # ----------------------------------------------------------
                # Clinical
                # ----------------------------------------------------------

                "external_clinical_association_performed":
                    bool(
                        clinical[
                            "performed"
                        ]
                    ),

                "external_clinical_association_justified":
                    bool(
                        clinical[
                            "justified"
                        ]
                    ),

                "external_clinical_validation_status":
                    str(
                        clinical[
                            "status"
                        ]
                    ),

                # ----------------------------------------------------------
                # Coloc / gene / causality
                # ----------------------------------------------------------

                "formal_colocalization_justified":
                    bool(
                        coloc[
                            "justified"
                        ]
                    ),

                "formal_colocalization_performed":
                    bool(
                        coloc[
                            "performed"
                        ]
                    ),

                "formal_colocalization_status":
                    str(
                        coloc[
                            "status"
                        ]
                    ),

                "gene_assignment_resolved":
                    bool(
                        gene[
                            "resolved"
                        ]
                    ),

                "nearest_gene_substitution_used":
                    bool(
                        gene[
                            "nearest_gene_substitution_used"
                        ]
                    ),

                "causal_relationship_established":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Evidence-domain matrix
# ============================================================================


def build_evidence_domain_matrix(
    *,
    candidate_evidence_matrix: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build descriptive evidence-domain availability matrix."""

    rows: list[dict[str, Any]] = []

    for _, row in candidate_evidence_matrix.iterrows():

        domains = {
            "disease_genetic":
                bool(
                    row[
                        "locus_level_disease_evidence"
                    ]
                    and
                    row[
                        "exact_variant_full_sumstats_evidence"
                    ]
                ),

            "regulatory":
                bool(
                    row[
                        "regulatory_feature_confirmed"
                    ]
                ),

            "trf_identity":
                bool(
                    row[
                        "trf_identifier_resolved"
                    ]
                ),

            "sequence":
                bool(
                    row[
                        "trf_exact_sequence_resolved"
                    ]
                ),

            "external_molecular":
                bool(
                    row[
                        "external_molecular_evidence_available"
                    ]
                ),

            "tcga_direct":
                bool(
                    row[
                        "tcga_direct_trf_quantification_performed"
                    ]
                ),

            "clinical":
                bool(
                    row[
                        "external_clinical_association_performed"
                    ]
                ),

            "formal_coloc":
                bool(
                    row[
                        "formal_colocalization_performed"
                    ]
                ),

            "gene_assignment":
                bool(
                    row[
                        "gene_assignment_resolved"
                    ]
                ),
        }

        rows.append(
            {
                "rsid":
                    str(
                        row[
                            "rsid"
                        ]
                    ),

                "priority_rank":
                    int(
                        row[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        row[
                            "priority_class"
                        ]
                    ),

                "disease_genetic_evidence":
                    domains[
                        "disease_genetic"
                    ],

                "regulatory_feature_evidence":
                    domains[
                        "regulatory"
                    ],

                "trf_identity_evidence":
                    domains[
                        "trf_identity"
                    ],

                "trf_sequence_evidence":
                    domains[
                        "sequence"
                    ],

                "external_molecular_evidence":
                    domains[
                        "external_molecular"
                    ],

                "tcga_direct_trf_evidence":
                    domains[
                        "tcga_direct"
                    ],

                "clinical_association_evidence":
                    domains[
                        "clinical"
                    ],

                "formal_colocalization_evidence":
                    domains[
                        "formal_coloc"
                    ],

                "gene_assignment_evidence":
                    domains[
                        "gene_assignment"
                    ],

                "evidence_domains_supported":
                    int(
                        sum(
                            domains.values()
                        )
                    ),

                "domain_count_is_descriptive_only":
                    True,

                "domain_count_used_for_priority":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Limitation matrix
# ============================================================================


def build_candidate_limitation_matrix(
    *,
    candidate_evidence_matrix: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Generate explicit limitation list for each candidate."""

    vocab = config[
        "limitations"
    ]

    rows: list[dict[str, Any]] = []

    for _, row in candidate_evidence_matrix.iterrows():

        limitations: list[str] = []

        if not bool(
            row[
                "trf_exact_sequence_resolved"
            ]
        ):

            limitations.append(
                vocab[
                    "sequence_unresolved"
                ]
            )

        if not bool(
            row[
                "regulatory_feature_confirmed"
            ]
        ):

            limitations.append(
                vocab[
                    "regulatory_feature_unresolved"
                ]
            )

        if not bool(
            row[
                "gene_assignment_resolved"
            ]
        ):

            limitations.append(
                vocab[
                    "gene_annotation_unresolved"
                ]
            )

        if (
            str(
                row[
                    "tcga_validation_status"
                ]
            )
            ==
            "DIRECT_TCGA_TRF_QUANTIFICATION_DEFERRED_CONTROLLED_ACCESS"
        ):

            limitations.append(
                vocab[
                    "tcga_controlled_access"
                ]
            )

        if (
            bool(
                row[
                    "external_molecular_evidence_available"
                ]
            )
            and
            not bool(
                row[
                    "standalone_mature_trf_confirmed"
                ]
            )
        ):

            limitations.append(
                vocab[
                    "mature_fragment_not_confirmed"
                ]
            )

        if not bool(
            row[
                "external_molecular_evidence_available"
            ]
        ):

            limitations.append(
                vocab[
                    "external_molecular_unavailable"
                ]
            )

        if not bool(
            row[
                "external_clinical_association_justified"
            ]
        ):

            limitations.append(
                vocab[
                    "clinical_association_not_justified"
                ]
            )

        if not bool(
            row[
                "formal_colocalization_performed"
            ]
        ):

            limitations.append(
                vocab[
                    "formal_coloc_unavailable"
                ]
            )

        if (
            bool(
                row[
                    "within_publication_directional_heterogeneity"
                ]
            )
            or
            bool(
                row[
                    "cross_publication_directional_heterogeneity"
                ]
            )
        ):

            limitations.append(
                vocab[
                    "directional_heterogeneity"
                ]
            )

        limitations.append(
            vocab[
                "independence_unresolved"
            ]
        )

        limitations.append(
            vocab[
                "causal_unresolved"
            ]
        )

        rows.append(
            {
                "rsid":
                    str(
                        row[
                            "rsid"
                        ]
                    ),

                "priority_rank":
                    int(
                        row[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        row[
                            "priority_class"
                        ]
                    ),

                "limitation_count":
                    len(
                        limitations
                    ),

                "limitations":
                    limitations,

                "limitations_are_negative_evidence":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Integrated candidate classification
# ============================================================================


def _classify_candidate(
    *,
    row: pd.Series,
    config: dict[str, Any],
) -> str:
    """Assign evidence class without numerical scoring."""

    classes = config[
        "classification"
    ][
        "classes"
    ]

    rsid = str(
        row[
            "rsid"
        ]
    )

    # Highest priority candidate:
    # disease + regulatory feature + resolved sequence + external molecular
    if (
        rsid
        ==
        "rs10216902"
        and
        bool(
            row[
                "locus_level_disease_evidence"
            ]
        )
        and
        bool(
            row[
                "regulatory_feature_confirmed"
            ]
        )
        and
        bool(
            row[
                "trf_exact_sequence_resolved"
            ]
        )
        and
        bool(
            row[
                "external_molecular_evidence_available"
            ]
        )
    ):

        return classes[
            "high_convergent"
        ]

    # rs2328376 retains stronger recurrent disease evidence but lacks
    # resolved molecular/regulatory link.
    if (
        rsid
        ==
        "rs2328376"
        and
        bool(
            row[
                "locus_level_disease_evidence"
            ]
        )
        and
        bool(
            row[
                "exact_variant_full_sumstats_evidence"
            ]
        )
    ):

        return classes[
            "genetic_supported_mechanistic_unresolved"
        ]

    if (
        rsid
        ==
        "rs1288100"
        and
        bool(
            row[
                "locus_level_disease_evidence"
            ]
        )
        and
        bool(
            row[
                "exact_variant_full_sumstats_evidence"
            ]
        )
    ):

        return classes[
            "genetic_context_mechanistic_unresolved"
        ]

    return classes[
        "insufficient"
    ]


def build_candidate_resolution(
    *,
    candidate_evidence_matrix: pd.DataFrame,
    evidence_domain_matrix: pd.DataFrame,
    limitation_matrix: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build final integrated candidate resolution."""

    interpretations = config[
        "interpretation"
    ]

    rows: list[dict[str, Any]] = []

    for _, evidence_row in candidate_evidence_matrix.iterrows():

        rsid = str(
            evidence_row[
                "rsid"
            ]
        )

        domains = evidence_domain_matrix.loc[
            evidence_domain_matrix[
                "rsid"
            ]
            ==
            rsid
        ].iloc[
            0
        ]

        limitations = limitation_matrix.loc[
            limitation_matrix[
                "rsid"
            ]
            ==
            rsid
        ].iloc[
            0
        ]

        integrated_class = _classify_candidate(
            row=evidence_row,
            config=config,
        )

        limitation_values = list(
            limitations[
                "limitations"
            ]
        )

        primary_limitation = (
            limitation_values[
                0
            ]
            if limitation_values
            else None
        )

        rows.append(
            {
                "rsid":
                    rsid,

                "trf_id":
                    str(
                        evidence_row[
                            "trf_id"
                        ]
                    ),

                "priority_rank":
                    int(
                        evidence_row[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        evidence_row[
                            "priority_class"
                        ]
                    ),

                "evidence_domains_supported":
                    int(
                        domains[
                            "evidence_domains_supported"
                        ]
                    ),

                "integrated_evidence_class":
                    integrated_class,

                "primary_limitation":
                    primary_limitation,

                "limitation_count":
                    int(
                        limitations[
                            "limitation_count"
                        ]
                    ),

                "all_limitations":
                    limitation_values,

                "independent_replication_verified":
                    False,

                "formal_colocalization_performed":
                    False,

                "clinical_validation_established":
                    False,

                "causal_relationship_established":
                    False,

                "recommended_interpretation":
                    str(
                        interpretations[
                            rsid
                        ]
                    ),

                "retain_candidate":
                    True,
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            "priority_rank",
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    candidate_evidence_matrix: pd.DataFrame,
    candidate_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.6 summary."""

    candidates_with_regulatory = int(
        candidate_evidence_matrix[
            "regulatory_feature_confirmed"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    candidates_with_sequence = int(
        candidate_evidence_matrix[
            "trf_exact_sequence_resolved"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    candidates_with_external_molecular = int(
        candidate_evidence_matrix[
            "external_molecular_evidence_available"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    candidates_with_directional_heterogeneity = int(
        (
            candidate_evidence_matrix[
                "within_publication_directional_heterogeneity"
            ]
            .fillna(False)
            .astype(bool)
            |
            candidate_evidence_matrix[
                "cross_publication_directional_heterogeneity"
            ]
            .fillna(False)
            .astype(bool)
        ).sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.6",

                "candidate_variants_assessed":
                    int(
                        len(
                            candidate_resolution
                        )
                    ),

                "candidates_retained":
                    int(
                        candidate_resolution[
                            "retain_candidate"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "candidates_with_regulatory_feature_evidence":
                    candidates_with_regulatory,

                "candidates_with_resolved_trf_sequence":
                    candidates_with_sequence,

                "candidates_with_external_molecular_evidence":
                    candidates_with_external_molecular,

                "candidates_with_directional_heterogeneity":
                    candidates_with_directional_heterogeneity,

                "candidates_with_tcga_direct_trf_quantification":
                    int(
                        candidate_evidence_matrix[
                            "tcga_direct_trf_quantification_performed"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "candidates_with_external_clinical_association":
                    int(
                        candidate_evidence_matrix[
                            "external_clinical_association_performed"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "candidates_with_formal_colocalization":
                    int(
                        candidate_evidence_matrix[
                            "formal_colocalization_performed"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "candidates_with_resolved_gene_assignment":
                    int(
                        candidate_evidence_matrix[
                            "gene_assignment_resolved"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "independent_replications_verified":
                    0,

                "causal_candidates":
                    0,

                "candidate_priority_recalculated":
                    False,

                "arbitrary_additive_score_used":
                    False,

                "meta_analysis_performed":
                    False,

                "causal_inference_performed":
                    False,

                "overall_status":
                    config[
                        "overall_status"
                    ],

                "next_stage":
                    config[
                        "next_stage"
                    ],
            }
        ]
    )


# ============================================================================
# Validation
# ============================================================================


def _validate_result(
    *,
    candidate_evidence_matrix: pd.DataFrame,
    evidence_domain_matrix: pd.DataFrame,
    candidate_resolution: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict[str, Any],
) -> None:
    """Enforce M7.6 structural and scientific invariants."""

    expected_rsids = {
        str(
            value
        )
        for value in config[
            "validation"
        ][
            "expected_rsids"
        ]
    }

    observed_rsids = set(
        candidate_resolution[
            "rsid"
        ].astype(str)
    )

    if observed_rsids != expected_rsids:

        raise RuntimeError(
            "M7.6 candidate set mismatch."
        )

    if len(
        candidate_resolution
    ) != int(
        config[
            "validation"
        ][
            "expected_candidate_count"
        ]
    ):

        raise RuntimeError(
            "M7.6 candidate count mismatch."
        )

    observed_order = (
        candidate_resolution
        .sort_values(
            "priority_rank"
        )[
            "rsid"
        ]
        .astype(str)
        .tolist()
    )

    expected_order = [
        str(
            value
        )
        for value in config[
            "validation"
        ][
            "expected_priority_order"
        ]
    ]

    if observed_order != expected_order:

        raise RuntimeError(
            "M7.6 changed locked candidate priority."
        )

    if bool(
        candidate_resolution[
            "causal_relationship_established"
        ]
        .fillna(False)
        .astype(bool)
        .any()
    ):

        raise RuntimeError(
            "M7.6 causal safeguard violated."
        )

    if bool(
        candidate_resolution[
            "formal_colocalization_performed"
        ]
        .fillna(False)
        .astype(bool)
        .any()
    ):

        raise RuntimeError(
            "M7.6 colocalization safeguard violated."
        )

    if bool(
        candidate_resolution[
            "independent_replication_verified"
        ]
        .fillna(False)
        .astype(bool)
        .any()
    ):

        raise RuntimeError(
            "M7.6 independent replication safeguard violated."
        )

    if bool(
        evidence_domain_matrix[
            "domain_count_used_for_priority"
        ]
        .fillna(False)
        .astype(bool)
        .any()
    ):

        raise RuntimeError(
            "M7.6 domain-count scoring safeguard violated."
        )

    if len(
        summary
    ) != 1:

        raise RuntimeError(
            "M7.6 summary must contain exactly one row."
        )


# ============================================================================
# Public API
# ============================================================================


def assess_integrated_evidence_resolution(
    *,
    m7_5b_resolution: pd.DataFrame,
    m7_5c_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> M76Result:
    """Execute M7.6 integrated evidence synthesis."""

    candidate_evidence_matrix = (
        build_candidate_evidence_matrix(
            m7_5b_resolution=m7_5b_resolution,
            m7_5c_resolution=m7_5c_resolution,
            config=config,
        )
    )

    evidence_domain_matrix = (
        build_evidence_domain_matrix(
            candidate_evidence_matrix=(
                candidate_evidence_matrix
            ),
            config=config,
        )
    )

    candidate_limitation_matrix = (
        build_candidate_limitation_matrix(
            candidate_evidence_matrix=(
                candidate_evidence_matrix
            ),
            config=config,
        )
    )

    candidate_resolution = (
        build_candidate_resolution(
            candidate_evidence_matrix=(
                candidate_evidence_matrix
            ),
            evidence_domain_matrix=(
                evidence_domain_matrix
            ),
            limitation_matrix=(
                candidate_limitation_matrix
            ),
            config=config,
        )
    )

    summary = build_summary(
        candidate_evidence_matrix=(
            candidate_evidence_matrix
        ),
        candidate_resolution=(
            candidate_resolution
        ),
        config=config,
    )

    _validate_result(
        candidate_evidence_matrix=(
            candidate_evidence_matrix
        ),
        evidence_domain_matrix=(
            evidence_domain_matrix
        ),
        candidate_resolution=(
            candidate_resolution
        ),
        summary=summary,
        config=config,
    )

    return M76Result(
        candidate_evidence_matrix=(
            candidate_evidence_matrix
        ),
        evidence_domain_matrix=(
            evidence_domain_matrix
        ),
        candidate_limitation_matrix=(
            candidate_limitation_matrix
        ),
        candidate_resolution=(
            candidate_resolution
        ),
        summary=summary,
    )