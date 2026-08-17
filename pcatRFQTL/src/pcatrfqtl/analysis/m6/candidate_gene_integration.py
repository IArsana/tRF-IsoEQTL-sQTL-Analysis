"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/candidate_gene_integration.py

Description:
    Core logic for M6.3 Candidate → Gene Integration.

    M6.3 integrates candidate-level gene evidence that has already been
    established by upstream analyses.

    The stage does NOT create new gene assignments.

    In particular:

        rs10216902 → INT98200

    is preserved as a supported regulatory-feature relationship even when
    the parent gene of INT98200 remains unresolved.

    An unresolved parent gene is treated as an annotation limitation and
    not as biological evidence that the feature lacks a parent gene.

    M6.3 does NOT:
        - assign nearest genes;
        - map genes from SNP proximity;
        - map genes from cis proximity;
        - decode regulatory feature numeric identifiers;
        - infer genes from candidate chromosome position;
        - perform colocalization;
        - perform fine-mapping;
        - make causal claims.

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


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class CandidateGeneIntegrationResult:
    """Container for M6.3 outputs."""

    evidence_matrix: pd.DataFrame
    final_integration: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional scalar text."""

    if value is None:
        return None

    try:

        if pd.isna(
            value
        ):
            return None

    except (
        TypeError,
        ValueError,
    ):

        pass

    normalized = str(
        value
    ).strip()

    return normalized or None


def _normalize_rsid(
    value: Any,
) -> str | None:
    """Normalize canonical rsID."""

    value = _normalize_text(
        value
    )

    if value is None:
        return None

    value = value.lower()

    if not value.startswith(
        "rs"
    ):
        return None

    if not value[
        2:
    ].isdigit():
        return None

    return value


# ============================================================================
# Candidate base
# ============================================================================


def extract_candidate_base(
    *,
    m6_1_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract locked candidate prioritization."""

    rows: list[
        dict[str, Any]
    ] = []

    candidates = m6_1_qc.get(
        "candidate_priority",
        [],
    )

    if not isinstance(
        candidates,
        list,
    ):

        raise RuntimeError(
            "M6.1 candidate_priority must be a list."
        )

    for candidate in candidates:

        lead_rsid = _normalize_rsid(
            candidate.get(
                "lead_rsid"
            )
        )

        if lead_rsid is None:

            raise RuntimeError(
                "Invalid lead_rsid in M6.1 candidate priority."
            )

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

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

                "continue_biological_interpretation":
                    bool(
                        candidate.get(
                            "continue_m6",
                            True,
                        )
                    ),
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No M6.1 candidates available for M6.3."
        )

    return result.sort_values(
        [
            "priority_rank",
            "lead_rsid",
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Regulatory feature evidence
# ============================================================================


def extract_regulatory_feature_evidence(
    *,
    m6_2a_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract candidate → regulatory-feature relationships from M6.2A."""

    rows: list[
        dict[str, Any]
    ] = []

    for record in m6_2a_qc.get(
        "candidate_mapping",
        [],
    ):

        lead_rsid = _normalize_rsid(
            record.get(
                "lead_rsid"
            )
        )

        if lead_rsid is None:
            continue

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "regulatory_feature_id":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_id"
                        )
                    ),

                "regulatory_feature_class":
                    _normalize_text(
                        record.get(
                            "feature_class"
                        )
                    ),

                "regulatory_scope":
                    _normalize_text(
                        record.get(
                            "regulatory_scope"
                        )
                    ),

                "regulatory_s1_qtl_support":
                    bool(
                        record.get(
                            "s1_qtl_support",
                            False,
                        )
                    ),

                "regulatory_s2_gwas_ld_context":
                    bool(
                        record.get(
                            "s2_gwas_ld_context",
                            False,
                        )
                    ),

                "regulatory_s7_differential_support":
                    bool(
                        record.get(
                            "s7_differential_support",
                            False,
                        )
                    ),

                "regulatory_feature_resolution_status":
                    _normalize_text(
                        record.get(
                            "resolution_status"
                        )
                    ),

                "regulatory_feature_present":
                    (
                        _normalize_text(
                            record.get(
                                "regulatory_feature_id"
                            )
                        )
                        is not None
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Source gene resolution
# ============================================================================


def extract_source_gene_resolution(
    *,
    m6_2b_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract source-native gene annotation results."""

    rows: list[
        dict[str, Any]
    ] = []

    for record in m6_2b_qc.get(
        "feature_resolution",
        [],
    ):

        lead_rsid = _normalize_rsid(
            record.get(
                "lead_rsid"
            )
        )

        if lead_rsid is None:
            continue

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "regulatory_feature_id":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_id"
                        )
                    ),

                "source_parent_gene":
                    _normalize_text(
                        record.get(
                            "parent_gene"
                        )
                    ),

                "source_parent_gene_resolved":
                    bool(
                        record.get(
                            "parent_gene_resolved",
                            False,
                        )
                    ),

                "source_gene_annotation_status":
                    _normalize_text(
                        record.get(
                            "gene_annotation_status"
                        )
                    ),

                "source_gene_assignment_method":
                    _normalize_text(
                        record.get(
                            "gene_assignment_method"
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Reconstruction resolution
# ============================================================================


def extract_reconstruction_resolution(
    *,
    m6_2c_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract event reconstruction result from M6.2C.1."""

    rows: list[
        dict[str, Any]
    ] = []

    for record in m6_2c_qc.get(
        "feature_readiness",
        [],
    ):

        lead_rsid = _normalize_rsid(
            record.get(
                "lead_rsid"
            )
        )

        if lead_rsid is None:
            continue

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "regulatory_feature_id":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_id"
                        )
                    ),

                "reconstruction_classification":
                    _normalize_text(
                        record.get(
                            "reconstruction_classification"
                        )
                    ),

                "exact_reconstruction_allowed":
                    bool(
                        record.get(
                            "exact_reconstruction_allowed",
                            False,
                        )
                    ),

                "source_compatible_reconstruction_allowed":
                    bool(
                        record.get(
                            "source_compatible_reconstruction_allowed",
                            False,
                        )
                    ),

                "reconstructed_parent_gene":
                    _normalize_text(
                        record.get(
                            "parent_gene"
                        )
                    ),

                "reconstructed_parent_gene_assigned":
                    bool(
                        record.get(
                            "parent_gene_assigned",
                            False,
                        )
                    ),

                "reconstruction_gene_assignment_allowed":
                    bool(
                        record.get(
                            "parent_gene_assignment_from_reconstruction_allowed",
                            False,
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Evidence matrix
# ============================================================================


def build_candidate_gene_evidence_matrix(
    *,
    candidates: pd.DataFrame,
    regulatory_features: pd.DataFrame,
    source_gene_resolution: pd.DataFrame,
    reconstruction_resolution: pd.DataFrame,
) -> pd.DataFrame:
    """Build candidate-level gene evidence matrix."""

    result = candidates.copy()

    # ----------------------------------------------------------------------
    # Regulatory-feature branch
    # ----------------------------------------------------------------------

    if not regulatory_features.empty:

        result = result.merge(
            regulatory_features,
            on="lead_rsid",
            how="left",
            validate="one_to_one",
        )

    else:

        result[
            "regulatory_feature_id"
        ] = None

    # ----------------------------------------------------------------------
    # Source gene branch
    # ----------------------------------------------------------------------

    if not source_gene_resolution.empty:

        source = source_gene_resolution.copy()

        result = result.merge(
            source,
            on=[
                "lead_rsid",
                "regulatory_feature_id",
            ],
            how="left",
            validate="one_to_one",
        )

    else:

        result[
            "source_parent_gene"
        ] = None

        result[
            "source_parent_gene_resolved"
        ] = False

    # ----------------------------------------------------------------------
    # Event reconstruction branch
    # ----------------------------------------------------------------------

    if not reconstruction_resolution.empty:

        reconstruction = (
            reconstruction_resolution.copy()
        )

        result = result.merge(
            reconstruction,
            on=[
                "lead_rsid",
                "regulatory_feature_id",
            ],
            how="left",
            validate="one_to_one",
        )

    else:

        result[
            "reconstructed_parent_gene"
        ] = None

        result[
            "reconstructed_parent_gene_assigned"
        ] = False

    # ----------------------------------------------------------------------
    # Normalize boolean columns after left joins.
    # ----------------------------------------------------------------------

    boolean_columns = [
        "regulatory_feature_present",
        "regulatory_s1_qtl_support",
        "regulatory_s2_gwas_ld_context",
        "regulatory_s7_differential_support",
        "source_parent_gene_resolved",
        "exact_reconstruction_allowed",
        "source_compatible_reconstruction_allowed",
        "reconstructed_parent_gene_assigned",
        "reconstruction_gene_assignment_allowed",
    ]

    for column in boolean_columns:

        if column not in result.columns:

            result[
                column
            ] = False

        result[
            column
        ] = (
            result[
                column
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
        )

    # ----------------------------------------------------------------------
    # Resolve any defensible gene already established upstream.
    # ----------------------------------------------------------------------

    final_genes: list[
        str | None
    ] = []

    gene_methods: list[
        str | None
    ] = []

    for row in result.to_dict(
        orient="records"
    ):

        source_gene = _normalize_text(
            row.get(
                "source_parent_gene"
            )
        )

        reconstructed_gene = _normalize_text(
            row.get(
                "reconstructed_parent_gene"
            )
        )

        if (
            source_gene is not None
            and
            bool(
                row.get(
                    "source_parent_gene_resolved",
                    False,
                )
            )
        ):

            final_genes.append(
                source_gene
            )

            gene_methods.append(
                "SOURCE_REPORTED_GENE"
            )

        elif (
            reconstructed_gene is not None
            and
            bool(
                row.get(
                    "reconstructed_parent_gene_assigned",
                    False,
                )
            )
            and
            bool(
                row.get(
                    "reconstruction_gene_assignment_allowed",
                    False,
                )
            )
        ):

            final_genes.append(
                reconstructed_gene
            )

            gene_methods.append(
                "SOURCE_REPRODUCIBLE_EVENT_GENE"
            )

        else:

            final_genes.append(
                None
            )

            gene_methods.append(
                None
            )

    result[
        "integrated_gene"
    ] = final_genes

    result[
        "integrated_gene_assignment_method"
    ] = gene_methods

    result[
        "integrated_gene_resolved"
    ] = (
        result[
            "integrated_gene"
        ]
        .notna()
    )

    result[
        "gene_assignment_created_in_m6_3"
    ] = False

    result[
        "nearest_gene_assignment_performed"
    ] = False

    result[
        "snp_proximity_gene_assignment_performed"
    ] = False

    result[
        "causal_claim_allowed"
    ] = False

    return result


# ============================================================================
# Final candidate-gene integration
# ============================================================================


def build_final_candidate_gene_integration(
    *,
    evidence_matrix: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Create final candidate-level M6.3 interpretation."""

    statuses = config[
        "resolution_status"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for record in evidence_matrix.to_dict(
        orient="records"
    ):

        lead_rsid = str(
            record[
                "lead_rsid"
            ]
        )

        feature_id = _normalize_text(
            record.get(
                "regulatory_feature_id"
            )
        )

        gene = _normalize_text(
            record.get(
                "integrated_gene"
            )
        )

        feature_present = bool(
            record.get(
                "regulatory_feature_present",
                False,
            )
        )

        # ------------------------------------------------------------------
        # Defensible gene assignment already exists.
        # ------------------------------------------------------------------

        if gene is not None:

            resolution = str(
                statuses[
                    "resolved_single_gene"
                ]
            )

            biological_status = (
                "CANDIDATE_WITH_SOURCE_SUPPORTED_GENE_ASSIGNMENT"
            )

        # ------------------------------------------------------------------
        # Supported regulatory feature but unresolved parent gene.
        # ------------------------------------------------------------------

        elif feature_present:

            resolution = str(
                statuses[
                    "regulatory_feature_gene_unresolved"
                ]
            )

            biological_status = (
                "CANDIDATE_WITH_SUPPORTED_REGULATORY_FEATURE_"
                "AND_UNRESOLVED_PARENT_GENE"
            )

        # ------------------------------------------------------------------
        # Candidate retained, but no defensible gene route yet.
        # ------------------------------------------------------------------

        else:

            resolution = str(
                statuses[
                    "no_gene_assignment"
                ]
            )

            biological_status = str(
                statuses[
                    "candidate_retained_without_gene"
                ]
            )

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "priority_rank":
                    int(
                        record[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        record[
                            "priority_class"
                        ]
                    ),

                "regulatory_feature_id":
                    feature_id,

                "regulatory_feature_class":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_class"
                        )
                    ),

                "regulatory_feature_supported":
                    feature_present,

                "integrated_gene":
                    gene,

                "integrated_gene_resolved":
                    bool(
                        gene is not None
                    ),

                "gene_assignment_method":
                    _normalize_text(
                        record.get(
                            "integrated_gene_assignment_method"
                        )
                    ),

                "gene_resolution_status":
                    resolution,

                "biological_interpretation_status":
                    biological_status,

                "candidate_retained":
                    True,

                "continue_biological_interpretation":
                    bool(
                        record.get(
                            "continue_biological_interpretation",
                            True,
                        )
                    ),

                "nearest_gene_assignment_performed":
                    False,

                "snp_proximity_gene_assignment_performed":
                    False,

                "causal_claim_allowed":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    return result.sort_values(
        [
            "priority_rank",
            "lead_rsid",
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    final_integration: pd.DataFrame,
) -> pd.DataFrame:
    """Build one-row M6.3 summary."""

    candidates = int(
        final_integration[
            "lead_rsid"
        ].nunique()
    )

    resolved = int(
        final_integration[
            "integrated_gene_resolved"
        ]
        .astype(
            bool
        )
        .sum()
    )

    with_feature_unresolved = int(
        (
            final_integration[
                "regulatory_feature_supported"
            ].astype(
                bool
            )
            &
            ~final_integration[
                "integrated_gene_resolved"
            ].astype(
                bool
            )
        ).sum()
    )

    without_gene_route = int(
        (
            ~final_integration[
                "regulatory_feature_supported"
            ].astype(
                bool
            )
            &
            ~final_integration[
                "integrated_gene_resolved"
            ].astype(
                bool
            )
        ).sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.3",

                "candidates_assessed":
                    candidates,

                "candidates_with_resolved_gene":
                    resolved,

                "candidates_with_supported_feature_unresolved_gene":
                    with_feature_unresolved,

                "candidates_without_defensible_gene_assignment":
                    without_gene_route,

                "candidates_retained":
                    candidates,

                "candidate_gene_assignment_created_in_m6_3":
                    False,

                "nearest_gene_assignment_performed":
                    False,

                "snp_proximity_gene_assignment_performed":
                    False,

                "formal_colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,

                "next_stage":
                    "M6.4_BIOLOGICAL_GENE_AND_PATHWAY_ANNOTATION",
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def integrate_candidate_genes(
    *,
    m6_1_qc: dict[str, Any],
    m6_2a_qc: dict[str, Any],
    m6_2b_qc: dict[str, Any],
    m6_2c_qc: dict[str, Any],
    config: dict[str, Any],
) -> CandidateGeneIntegrationResult:
    """Execute M6.3 Candidate → Gene Integration."""

    candidates = extract_candidate_base(
        m6_1_qc=m6_1_qc,
    )

    regulatory_features = (
        extract_regulatory_feature_evidence(
            m6_2a_qc=m6_2a_qc,
        )
    )

    source_gene_resolution = (
        extract_source_gene_resolution(
            m6_2b_qc=m6_2b_qc,
        )
    )

    reconstruction_resolution = (
        extract_reconstruction_resolution(
            m6_2c_qc=m6_2c_qc,
        )
    )

    evidence_matrix = (
        build_candidate_gene_evidence_matrix(
            candidates=candidates,
            regulatory_features=regulatory_features,
            source_gene_resolution=source_gene_resolution,
            reconstruction_resolution=reconstruction_resolution,
        )
    )

    final_integration = (
        build_final_candidate_gene_integration(
            evidence_matrix=evidence_matrix,
            config=config,
        )
    )

    summary = build_summary(
        final_integration=final_integration,
    )

    return CandidateGeneIntegrationResult(
        evidence_matrix=evidence_matrix,
        final_integration=final_integration,
        summary=summary,
    )