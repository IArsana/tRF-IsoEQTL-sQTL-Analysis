"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/validation_readiness.py

Description:
    Core logic for M7.1 TCGA-PRAD Clinical Validation Readiness Audit.

    This stage evaluates whether each validation route is currently
    methodologically available.

    M7.1 does NOT:
        - download TCGA data;
        - perform clinical association analysis;
        - quantify candidate tRFs;
        - perform survival analysis;
        - analyze germline genotypes;
        - reconstruct INT98200;
        - substitute miRNA expression for tRF expression;
        - substitute somatic variants for germline variants;
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
class ValidationReadinessResult:
    """Container for M7.1 outputs."""

    resource_readiness: pd.DataFrame
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


def extract_locked_candidates(
    *,
    m6_4b_qc: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract retained candidate state from locked M6.4B.

    Candidate and tRF identities are inherited rather than manually
    recreated in M7.1.
    """

    records = m6_4b_qc.get(
        "candidate_pathway_readiness",
        [],
    )

    if not isinstance(
        records,
        list,
    ):

        raise RuntimeError(
            "M6.4B candidate_pathway_readiness must be a list."
        )

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        if not bool(
            record.get(
                "candidate_retained",
                True,
            )
        ):
            continue

        trf_ids_raw = record.get(
            "trf_ids",
            [],
        )

        if not isinstance(
            trf_ids_raw,
            list,
        ):
            trf_ids_raw = []

        trf_ids = sorted(
            {
                str(value).strip()
                for value in trf_ids_raw
                if _normalize_text(value) is not None
            }
        )

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

                "trf_count":
                    int(
                        record.get(
                            "trf_count",
                            len(trf_ids),
                        )
                    ),

                "trf_ids":
                    trf_ids,

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

                "rna_processing_context":
                    _normalize_text(
                        record.get(
                            "rna_processing_context"
                        )
                    ),

                "trf_annotation_status":
                    _normalize_text(
                        record.get(
                            "trf_annotation_status"
                        )
                    ),
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No retained candidates available for M7.1."
        )

    if result[
        "lead_rsid"
    ].isna().any():

        raise RuntimeError(
            "M7.1 encountered candidate with missing rsID."
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
# Resource readiness
# ============================================================================


def build_resource_readiness(
    *,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build validation-domain capability audit."""

    capabilities = config[
        "gdc_capabilities"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for resource_id, resource in capabilities.items():

        rows.append(
            {
                "resource_id":
                    str(
                        resource_id
                    ),

                "domain":
                    str(
                        resource[
                            "domain"
                        ]
                    ),

                "public_resource_available":
                    bool(
                        resource[
                            "public_resource_available"
                        ]
                    ),

                "controlled_access_required":
                    bool(
                        resource[
                            "controlled_access_required"
                        ]
                    ),

                "direct_candidate_measurement":
                    bool(
                        resource[
                            "direct_candidate_measurement"
                        ]
                    ),

                "capability_status":
                    str(
                        resource[
                            "capability_status"
                        ]
                    ),

                "interpretation":
                    str(
                        resource[
                            "interpretation"
                        ]
                    ),

                "next_action":
                    str(
                        resource[
                            "next_action"
                        ]
                    ),

                "data_downloaded_in_m7_1":
                    False,

                "analysis_performed_in_m7_1":
                    False,
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
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Build candidate-specific validation-route readiness.

    Current M6.4A annotations provide tRF identifiers, but no exact sequence
    resolution. Therefore direct tRF quantification is not yet allowed.
    """

    statuses = config[
        "statuses"
    ]

    target_regulatory_feature = str(
        config[
            "regulatory_event_policy"
        ][
            "target_feature"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for record in candidates.to_dict(
        orient="records"
    ):

        lead_rsid = str(
            record[
                "lead_rsid"
            ]
        )

        trf_ids = record.get(
            "trf_ids",
            [],
        )

        has_trf_identity = bool(
            trf_ids
        )

        # ------------------------------------------------------------------
        # tRF validation route
        # ------------------------------------------------------------------

        if has_trf_identity:

            trf_status = str(
                statuses[
                    "trf"
                ][
                    "identity_available_sequence_unresolved"
                ]
            )

            trf_sequence_resolved = False
            trf_quantification_ready = False

        else:

            trf_status = str(
                statuses[
                    "trf"
                ][
                    "identity_missing"
                ]
            )

            trf_sequence_resolved = False
            trf_quantification_ready = False

        # ------------------------------------------------------------------
        # Variant validation route
        # ------------------------------------------------------------------

        variant_status = str(
            statuses[
                "variant"
            ][
                "controlled_access_required"
            ]
        )

        variant_direct_validation_ready = False

        # ------------------------------------------------------------------
        # Regulatory-event route
        # ------------------------------------------------------------------

        regulatory_feature = _normalize_text(
            record.get(
                "regulatory_feature_id"
            )
        )

        if (
            regulatory_feature
            ==
            target_regulatory_feature
        ):

            regulatory_event_status = str(
                statuses[
                    "regulatory_event"
                ][
                    "not_reproducible"
                ]
            )

            regulatory_event_direct_validation_ready = False

        elif regulatory_feature is not None:

            regulatory_event_status = str(
                statuses[
                    "regulatory_event"
                ][
                    "not_reproducible"
                ]
            )

            regulatory_event_direct_validation_ready = False

        else:

            regulatory_event_status = None
            regulatory_event_direct_validation_ready = False

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

                "trf_ids":
                    trf_ids,

                "trf_identity_available":
                    has_trf_identity,

                "trf_sequence_resolved":
                    trf_sequence_resolved,

                "trf_validation_status":
                    trf_status,

                "direct_trf_quantification_ready":
                    trf_quantification_ready,

                "processed_mirna_accepted_as_trf_measurement":
                    False,

                "variant_validation_status":
                    variant_status,

                "individual_germline_genotype_ready":
                    False,

                "direct_variant_validation_ready":
                    variant_direct_validation_ready,

                "regulatory_feature_id":
                    regulatory_feature,

                "regulatory_feature_supported":
                    bool(
                        record.get(
                            "regulatory_feature_supported",
                            False,
                        )
                    ),

                "regulatory_event_validation_status":
                    regulatory_event_status,

                "direct_regulatory_event_validation_ready":
                    regulatory_event_direct_validation_ready,

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

                "clinical_metadata_route_ready":
                    True,

                "candidate_retained_for_m7":
                    True,

                "candidate_substitution_performed":
                    False,

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    ).sort_values(
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
    resource_readiness: pd.DataFrame,
    candidate_readiness: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.1 summary."""

    candidates_assessed = int(
        len(
            candidate_readiness
        )
    )

    candidates_with_trf_identity = int(
        candidate_readiness[
            "trf_identity_available"
        ]
        .astype(bool)
        .sum()
    )

    candidates_trf_quantification_ready = int(
        candidate_readiness[
            "direct_trf_quantification_ready"
        ]
        .astype(bool)
        .sum()
    )

    candidates_variant_validation_ready = int(
        candidate_readiness[
            "direct_variant_validation_ready"
        ]
        .astype(bool)
        .sum()
    )

    regulatory_features_present = int(
        candidate_readiness[
            "regulatory_feature_id"
        ]
        .notna()
        .sum()
    )

    regulatory_events_ready = int(
        candidate_readiness[
            "direct_regulatory_event_validation_ready"
        ]
        .astype(bool)
        .sum()
    )

    controlled_domains = int(
        resource_readiness[
            "controlled_access_required"
        ]
        .astype(bool)
        .sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.1",

                "target_project":
                    str(
                        config[
                            "target_cohort"
                        ][
                            "project_id"
                        ]
                    ),

                "candidates_assessed":
                    candidates_assessed,

                "candidates_with_trf_identity":
                    candidates_with_trf_identity,

                "candidates_with_resolved_trf_sequence":
                    0,

                "candidates_ready_for_direct_trf_quantification":
                    candidates_trf_quantification_ready,

                "candidates_ready_for_direct_variant_validation":
                    candidates_variant_validation_ready,

                "candidates_with_regulatory_feature":
                    regulatory_features_present,

                "regulatory_events_ready_for_direct_validation":
                    regulatory_events_ready,

                "clinical_metadata_route_ready":
                    True,

                "controlled_access_domains":
                    controlled_domains,

                "clinical_association_performed":
                    False,

                "trf_quantification_performed":
                    False,

                "variant_validation_performed":
                    False,

                "regulatory_event_validation_performed":
                    False,

                "overall_readiness_status":
                    "PARTIAL_VALIDATION_ROUTES_AVAILABLE",

                "clinical_next_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "clinical"
                        ]
                    ),

                "trf_next_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "trf"
                        ]
                    ),

                "variant_next_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "variant"
                        ]
                    ),

                "regulatory_event_next_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "regulatory_event"
                        ]
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def assess_validation_readiness(
    *,
    m6_4b_qc: dict[str, Any],
    config: dict[str, Any],
) -> ValidationReadinessResult:
    """Execute M7.1 TCGA-PRAD validation readiness audit."""

    candidates = extract_locked_candidates(
        m6_4b_qc=m6_4b_qc,
    )

    resource_readiness = build_resource_readiness(
        config=config,
    )

    candidate_readiness = build_candidate_readiness(
        candidates=candidates,
        config=config,
    )

    summary = build_summary(
        resource_readiness=resource_readiness,
        candidate_readiness=candidate_readiness,
        config=config,
    )

    return ValidationReadinessResult(
        resource_readiness=resource_readiness,
        candidate_readiness=candidate_readiness,
        summary=summary,
    )