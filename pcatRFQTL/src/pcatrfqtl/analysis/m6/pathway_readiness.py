"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/pathway_readiness.py

Description:
    Core logic for M6.4B Pathway Readiness Assessment.

    This stage determines whether the current candidate-level gene evidence
    provides a defensible gene set for pathway enrichment.

    M6.4B does NOT perform pathway enrichment.

    It does NOT:
        - create new gene assignments;
        - add nearest genes;
        - infer genes from SNP proximity;
        - infer genes from cis proximity;
        - convert unresolved regulatory features into genes;
        - treat tRF identifiers as genes;
        - perform colocalization;
        - perform fine-mapping;
        - make causal pathway claims.

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
class PathwayReadinessResult:
    """Container for M6.4B outputs."""

    defensible_gene_set: pd.DataFrame
    candidate_readiness: pd.DataFrame
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


# ============================================================================
# Candidate extraction
# ============================================================================


def extract_candidate_gene_state(
    *,
    m6_3_qc: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract final candidate-gene state from locked M6.3 QC.
    """

    records = m6_3_qc.get(
        "final_candidate_gene_integration",
        [],
    )

    if not isinstance(
        records,
        list,
    ):

        raise RuntimeError(
            "M6.3 final_candidate_gene_integration must be a list."
        )

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        rows.append(
            {
                "lead_rsid":
                    _normalize_text(
                        record.get(
                            "lead_rsid"
                        )
                    ),

                "priority_rank":
                    record.get(
                        "priority_rank"
                    ),

                "priority_class":
                    _normalize_text(
                        record.get(
                            "priority_class"
                        )
                    ),

                "regulatory_feature_id":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_id"
                        )
                    ),

                "regulatory_feature_supported":
                    bool(
                        record.get(
                            "regulatory_feature_supported",
                            False,
                        )
                    ),

                "integrated_gene":
                    _normalize_text(
                        record.get(
                            "integrated_gene"
                        )
                    ),

                "integrated_gene_resolved":
                    bool(
                        record.get(
                            "integrated_gene_resolved",
                            False,
                        )
                    ),

                "gene_assignment_method":
                    _normalize_text(
                        record.get(
                            "gene_assignment_method"
                        )
                    ),

                "gene_resolution_status":
                    _normalize_text(
                        record.get(
                            "gene_resolution_status"
                        )
                    ),

                "candidate_retained":
                    bool(
                        record.get(
                            "candidate_retained",
                            True,
                        )
                    ),

                "continue_biological_interpretation":
                    bool(
                        record.get(
                            "continue_biological_interpretation",
                            True,
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# M6.4A consistency audit
# ============================================================================


def extract_trf_annotation_state(
    *,
    m6_4a_qc: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract candidate biological annotation state from M6.4A.

    M6.4A is used only for consistency/context validation.
    It does not supply new pathway genes.
    """

    records = m6_4a_qc.get(
        "candidate_trf_biological_annotation",
        [],
    )

    if not isinstance(
        records,
        list,
    ):

        raise RuntimeError(
            "M6.4A candidate_trf_biological_annotation must be a list."
        )

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        rows.append(
            {
                "lead_rsid":
                    _normalize_text(
                        record.get(
                            "lead_rsid"
                        )
                    ),

                "trf_count":
                    int(
                        record.get(
                            "trf_count",
                            0,
                        )
                    ),

                "trf_ids":
                    record.get(
                        "trf_ids",
                        [],
                    ),

                "rna_processing_context":
                    _normalize_text(
                        record.get(
                            "rna_processing_context"
                        )
                    ),

                "trf_annotation_status":
                    _normalize_text(
                        record.get(
                            "annotation_status"
                        )
                    ),

                "m6_4a_integrated_gene":
                    _normalize_text(
                        record.get(
                            "integrated_gene"
                        )
                    ),

                "m6_4a_gene_resolved":
                    bool(
                        record.get(
                            "integrated_gene_resolved",
                            False,
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Candidate readiness
# ============================================================================


def build_candidate_readiness(
    *,
    candidates: pd.DataFrame,
    trf_annotations: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Determine whether each candidate contributes a defensible pathway gene.
    """

    result = candidates.copy()

    if not trf_annotations.empty:

        result = result.merge(
            trf_annotations,
            on="lead_rsid",
            how="left",
            validate="one_to_one",
        )

    allowed_methods = {
        str(value)
        for value
        in config[
            "allowed_gene_assignment_methods"
        ]
    }

    exclusion_reasons = config[
        "candidate_exclusion_reason"
    ]

    pathway_eligible: list[
        bool
    ] = []

    exclusion_reason: list[
        str | None
    ] = []

    provenance_supported: list[
        bool
    ] = []

    for row in result.to_dict(
        orient="records"
    ):

        gene = _normalize_text(
            row.get(
                "integrated_gene"
            )
        )

        method = _normalize_text(
            row.get(
                "gene_assignment_method"
            )
        )

        resolved = bool(
            row.get(
                "integrated_gene_resolved",
                False,
            )
        )

        regulatory_feature = _normalize_text(
            row.get(
                "regulatory_feature_id"
            )
        )

        # ------------------------------------------------------------------
        # No resolved gene
        # ------------------------------------------------------------------

        if not resolved:

            pathway_eligible.append(
                False
            )

            provenance_supported.append(
                False
            )

            if regulatory_feature is not None:

                exclusion_reason.append(
                    str(
                        exclusion_reasons[
                            "regulatory_feature_gene_unresolved"
                        ]
                    )
                )

            else:

                exclusion_reason.append(
                    str(
                        exclusion_reasons[
                            "no_defensible_gene_assignment"
                        ]
                    )
                )

            continue

        # ------------------------------------------------------------------
        # Resolved flag but gene missing
        # ------------------------------------------------------------------

        if gene is None:

            pathway_eligible.append(
                False
            )

            provenance_supported.append(
                False
            )

            exclusion_reason.append(
                str(
                    exclusion_reasons[
                        "gene_identity_missing"
                    ]
                )
            )

            continue

        # ------------------------------------------------------------------
        # Assignment provenance not allowed
        # ------------------------------------------------------------------

        if method not in allowed_methods:

            pathway_eligible.append(
                False
            )

            provenance_supported.append(
                False
            )

            exclusion_reason.append(
                str(
                    exclusion_reasons[
                        "unsupported_gene_assignment_method"
                    ]
                )
            )

            continue

        # ------------------------------------------------------------------
        # Defensible candidate gene
        # ------------------------------------------------------------------

        pathway_eligible.append(
            True
        )

        provenance_supported.append(
            True
        )

        exclusion_reason.append(
            None
        )

    result[
        "gene_provenance_supported"
    ] = provenance_supported

    result[
        "eligible_for_pathway_gene_set"
    ] = pathway_eligible

    result[
        "pathway_exclusion_reason"
    ] = exclusion_reason

    result[
        "gene_added_in_m6_4b"
    ] = False

    result[
        "nearest_gene_added"
    ] = False

    result[
        "unresolved_feature_converted_to_gene"
    ] = False

    return result


# ============================================================================
# Gene-set construction
# ============================================================================


def build_defensible_gene_set(
    *,
    candidate_readiness: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build unique defensible pathway gene set.

    Only candidates already passing upstream gene resolution and provenance
    checks are included.
    """

    eligible = candidate_readiness.loc[
        candidate_readiness[
            "eligible_for_pathway_gene_set"
        ].astype(
            bool
        )
    ].copy()

    if eligible.empty:

        return pd.DataFrame(
            columns=[
                "gene",
                "candidate_count",
                "candidate_rsids",
                "priority_classes",
                "assignment_methods",
            ]
        )

    eligible[
        "gene"
    ] = eligible[
        "integrated_gene"
    ].astype(
        str
    )

    rows: list[
        dict[str, Any]
    ] = []

    for gene, group in eligible.groupby(
        "gene",
        sort=True,
    ):

        rows.append(
            {
                "gene":
                    str(
                        gene
                    ),

                "candidate_count":
                    int(
                        group[
                            "lead_rsid"
                        ].nunique()
                    ),

                "candidate_rsids":
                    sorted(
                        {
                            str(value)
                            for value
                            in group[
                                "lead_rsid"
                            ]
                            .dropna()
                            .tolist()
                        }
                    ),

                "priority_classes":
                    sorted(
                        {
                            str(value)
                            for value
                            in group[
                                "priority_class"
                            ]
                            .dropna()
                            .tolist()
                        }
                    ),

                "assignment_methods":
                    sorted(
                        {
                            str(value)
                            for value
                            in group[
                                "gene_assignment_method"
                            ]
                            .dropna()
                            .tolist()
                        }
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Final readiness decision
# ============================================================================


def determine_readiness(
    *,
    candidate_readiness: pd.DataFrame,
    gene_set: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[
    str,
    bool,
]:
    """
    Determine whether pathway enrichment is justified.
    """

    statuses = config[
        "status"
    ]

    minimum_genes = int(
        config[
            "readiness_policy"
        ][
            "minimum_unique_genes"
        ]
    )

    if candidate_readiness.empty:

        return (
            str(
                statuses[
                    "empty_candidate_set"
                ]
            ),
            False,
        )

    unique_genes = int(
        len(
            gene_set
        )
    )

    resolved_count = int(
        candidate_readiness[
            "integrated_gene_resolved"
        ]
        .astype(
            bool
        )
        .sum()
    )

    eligible_count = int(
        candidate_readiness[
            "eligible_for_pathway_gene_set"
        ]
        .astype(
            bool
        )
        .sum()
    )

    # ----------------------------------------------------------------------
    # No resolved genes at all
    # ----------------------------------------------------------------------

    if resolved_count == 0:

        return (
            str(
                statuses[
                    "no_resolved_genes"
                ]
            ),
            False,
        )

    # ----------------------------------------------------------------------
    # Resolved genes exist but none passes provenance requirements
    # ----------------------------------------------------------------------

    if eligible_count == 0:

        return (
            str(
                statuses[
                    "unsupported_gene_provenance"
                ]
            ),
            False,
        )

    # ----------------------------------------------------------------------
    # Valid genes exist but too few for enrichment
    # ----------------------------------------------------------------------

    if unique_genes < minimum_genes:

        return (
            str(
                statuses[
                    "insufficient_gene_count"
                ]
            ),
            False,
        )

    return (
        str(
            statuses[
                "ready"
            ]
        ),
        True,
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    candidate_readiness: pd.DataFrame,
    gene_set: pd.DataFrame,
    readiness_status: str,
    enrichment_ready: bool,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M6.4B summary."""

    resolved_candidate_genes = int(
        candidate_readiness[
            "integrated_gene_resolved"
        ]
        .astype(
            bool
        )
        .sum()
        if not candidate_readiness.empty
        else 0
    )

    eligible_candidate_genes = int(
        candidate_readiness[
            "eligible_for_pathway_gene_set"
        ]
        .astype(
            bool
        )
        .sum()
        if not candidate_readiness.empty
        else 0
    )

    unresolved_regulatory_features = int(
        (
            candidate_readiness[
                "regulatory_feature_id"
            ].notna()
            &
            ~candidate_readiness[
                "integrated_gene_resolved"
            ].astype(
                bool
            )
        ).sum()
        if not candidate_readiness.empty
        else 0
    )

    no_gene_assignment = int(
        (
            candidate_readiness[
                "integrated_gene"
            ].isna()
            &
            candidate_readiness[
                "regulatory_feature_id"
            ].isna()
        ).sum()
        if not candidate_readiness.empty
        else 0
    )

    next_stage = (
        config[
            "next_stage"
        ][
            "if_ready"
        ]
        if enrichment_ready
        else
        config[
            "next_stage"
        ][
            "if_blocked"
        ]
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.4B",

                "candidates_assessed":
                    int(
                        len(
                            candidate_readiness
                        )
                    ),

                "resolved_candidate_genes":
                    resolved_candidate_genes,

                "eligible_candidate_genes":
                    eligible_candidate_genes,

                "unique_defensible_genes":
                    int(
                        len(
                            gene_set
                        )
                    ),

                "minimum_unique_genes_required":
                    int(
                        config[
                            "readiness_policy"
                        ][
                            "minimum_unique_genes"
                        ]
                    ),

                "candidates_with_unresolved_regulatory_feature_gene":
                    unresolved_regulatory_features,

                "candidates_without_defensible_gene_assignment":
                    no_gene_assignment,

                "pathway_enrichment_ready":
                    enrichment_ready,

                "pathway_readiness_status":
                    readiness_status,

                "enrichment_performed":
                    False,

                "gene_set_imputed":
                    False,

                "nearest_gene_added":
                    False,

                "snp_proximity_gene_added":
                    False,

                "cis_proximity_gene_added":
                    False,

                "unresolved_feature_converted_to_gene":
                    False,

                "pathway_claim_generated":
                    False,

                "causal_inference_performed":
                    False,

                "next_stage":
                    str(
                        next_stage
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def assess_pathway_readiness(
    *,
    m6_3_qc: dict[str, Any],
    m6_4a_qc: dict[str, Any],
    config: dict[str, Any],
) -> PathwayReadinessResult:
    """Execute M6.4B Pathway Readiness Assessment."""

    candidates = extract_candidate_gene_state(
        m6_3_qc=m6_3_qc,
    )

    trf_annotations = extract_trf_annotation_state(
        m6_4a_qc=m6_4a_qc,
    )

    candidate_readiness = build_candidate_readiness(
        candidates=candidates,
        trf_annotations=trf_annotations,
        config=config,
    )

    gene_set = build_defensible_gene_set(
        candidate_readiness=candidate_readiness,
    )

    (
        readiness_status,
        enrichment_ready,
    ) = determine_readiness(
        candidate_readiness=candidate_readiness,
        gene_set=gene_set,
        config=config,
    )

    summary = build_summary(
        candidate_readiness=candidate_readiness,
        gene_set=gene_set,
        readiness_status=readiness_status,
        enrichment_ready=enrichment_ready,
        config=config,
    )

    return PathwayReadinessResult(
        defensible_gene_set=gene_set,
        candidate_readiness=candidate_readiness,
        summary=summary,
    )