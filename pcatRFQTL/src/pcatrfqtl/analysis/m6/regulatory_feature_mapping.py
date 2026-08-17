"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/regulatory_feature_mapping.py

Description:
    Core logic for M6.2A Candidate → Regulatory Feature Resolution.

    M6.2A inherits candidate-to-regulatory-feature assignments from the
    locked M6.1 candidate prioritization output and validates those features
    against the original Moradi supplementary datasets.

    For retained intron-regulatory features, the current implementation
    evaluates three evidence layers:

        S1:
            cis intron-retention sQTL evidence.

        S2:
            Moradi GWAS-LD context for prostate cancer.

        S7:
            differential intron-retention evidence.

    M6.2A intentionally does NOT assign a parent gene.

    A feature such as INT98200 may be biologically supported while its
    parent gene remains unresolved in the available Moradi supplementary
    tables. Such features are carried forward with:

        PENDING_SOURCE_ANNOTATION

    Important safeguards:
        - Candidate-feature relationships are inherited from M6.1.
        - Raw Moradi data are not modified.
        - Feature IDs are matched exactly.
        - Candidate absence is not interpreted as biological absence.
        - Nearest-gene assignment is prohibited.
        - Positional gene inference is prohibited.
        - Regulatory support is not interpreted as causality.
        - No formal colocalization is performed.
        - No fine-mapping is performed.
        - No causal inference is performed.

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

import numpy as np
import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class RegulatoryFeatureMappingResult:
    """Container for M6.2A analytical outputs."""

    feature_evidence: pd.DataFrame

    candidate_mapping: pd.DataFrame

    summary: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize simple textual identifiers."""

    if value is None:

        return None

    if pd.isna(
        value
    ):

        return None

    normalized = str(
        value
    ).strip()

    if not normalized:

        return None

    return normalized


def _normalize_rsid(
    value: Any,
) -> str | None:
    """Normalize canonical rsID when possible."""

    normalized = _normalize_text(
        value
    )

    if normalized is None:

        return None

    normalized = normalized.lower()

    if not normalized.startswith(
        "rs"
    ):

        return None

    if not normalized[
        2:
    ].isdigit():

        return None

    return normalized


def _numeric(
    series: pd.Series,
) -> pd.Series:
    """Convert series to numeric conservatively."""

    return pd.to_numeric(
        series,
        errors="coerce",
    )


def _feature_class(
    feature_id: str,
    *,
    config: dict[str, Any],
) -> str:
    """Resolve feature class from configured identifier prefixes."""

    prefixes = config[
        "feature_policy"
    ][
        "supported_feature_prefixes"
    ]

    normalized = str(
        feature_id
    ).strip().upper()

    for prefix, properties in prefixes.items():

        if normalized.startswith(
            str(
                prefix
            ).upper()
        ):

            return str(
                properties[
                    "feature_class"
                ]
            )

    return (
        "UNKNOWN_REGULATORY_FEATURE"
    )


def _disease_matches(
    series: pd.Series,
    *,
    config: dict[str, Any],
) -> pd.Series:
    """Match configured prostate-cancer disease aliases."""

    aliases = {
        str(
            value
        )
        .strip()
        .lower()
        for value
        in config[
            "disease_context"
        ][
            "disease_aliases"
        ]
    }

    normalized = (
        series
        .astype("string")
        .str.strip()
        .str.lower()
    )

    return normalized.isin(
        aliases
    )


# ============================================================================
# M6.1 candidate-feature extraction
# ============================================================================


def extract_candidate_feature_assignments(
    *,
    m6_1_qc: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract candidate → regulatory-feature assignments from locked M6.1.

    M6.2A does not discover new candidate-feature relationships.
    """

    candidate_rows = m6_1_qc.get(
        "candidate_priority",
        [],
    )

    output_rows: list[
        dict[str, Any]
    ] = []

    for candidate in candidate_rows:

        lead_rsid = _normalize_rsid(
            candidate.get(
                "lead_rsid"
            )
        )

        if lead_rsid is None:

            raise ValueError(
                "M6.1 candidate contains invalid lead_rsid."
            )

        raw_features = candidate.get(
            "regulatory_feature_ids",
            [],
        )

        if raw_features is None:

            raw_features = []

        normalized_features = [
            feature
            for feature
            in (
                _normalize_text(
                    value
                )
                for value
                in raw_features
            )
            if feature is not None
        ]

        # ------------------------------------------------------------------
        # Candidate without retained Moradi feature.
        # ------------------------------------------------------------------

        if not normalized_features:

            output_rows.append(
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

                    "regulatory_feature_id":
                        None,

                    "candidate_feature_assignment_present":
                        False,

                    "assignment_source":
                        "M6.1_LOCKED_CANDIDATE_PRIORITIZATION",
                }
            )

            continue

        # ------------------------------------------------------------------
        # Candidate with retained regulatory feature(s).
        # ------------------------------------------------------------------

        for feature_id in normalized_features:

            output_rows.append(
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

                    "regulatory_feature_id":
                        feature_id,

                    "candidate_feature_assignment_present":
                        True,

                    "assignment_source":
                        "M6.1_LOCKED_CANDIDATE_PRIORITIZATION",
                }
            )

    return pd.DataFrame(
        output_rows
    )


# ============================================================================
# Raw Moradi evidence extraction
# ============================================================================


def build_feature_evidence(
    *,
    assignments: pd.DataFrame,
    s1: pd.DataFrame,
    s2: pd.DataFrame,
    s7: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve evidence for every assigned regulatory feature."""

    assigned_features = (
        assignments.loc[
            assignments[
                "candidate_feature_assignment_present"
            ].astype(
                bool
            ),
            "regulatory_feature_id",
        ]
        .dropna()
        .astype(
            str
        )
        .drop_duplicates()
        .tolist()
    )

    if not assigned_features:

        return pd.DataFrame(
            columns=[
                "regulatory_feature_id",
                "feature_class",
                "regulatory_scope",
                "s1_qtl_support",
                "s1_qtl_rows",
                "s1_unique_qtl_rsids",
                "s1_minimum_p_value",
                "s1_minimum_fdr",
                "s2_gwas_ld_context",
                "s2_gwas_ld_rows",
                "s2_unique_qtl_rsids",
                "s2_unique_tag_rsids",
                "s2_maximum_ld",
                "s7_differential_support",
                "s7_differential_rows",
                "s7_log2_fold_change",
                "s7_p_value",
                "s7_adjusted_p_value",
                "parent_gene",
                "parent_gene_status",
                "overall_feature_confirmation",
            ]
        )

    s1_config = config[
        "moradi_sources"
    ][
        "s1_cis_qtl"
    ]

    s2_config = config[
        "moradi_sources"
    ][
        "s2_gwas_ld"
    ]

    s7_config = config[
        "moradi_sources"
    ][
        "s7_differential_intron"
    ]

    s1_feature_column = str(
        s1_config[
            "feature_column"
        ]
    )

    s2_feature_column = str(
        s2_config[
            "feature_column"
        ]
    )

    s7_feature_column = str(
        s7_config[
            "feature_column"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for feature_id in assigned_features:

        # ==================================================================
        # S1 — cis intron QTL
        # ==================================================================

        s1_mask = (
            s1[
                s1_feature_column
            ]
            .astype("string")
            .str.strip()
            .eq(
                feature_id
            )
        )

        s1_hits = s1.loc[
            s1_mask
        ].copy()

        s1_variant_column = str(
            s1_config[
                "variant_column"
            ]
        )

        s1_p_column = str(
            s1_config[
                "p_value_column"
            ]
        )

        s1_fdr_column = str(
            s1_config[
                "fdr_column"
            ]
        )

        if not s1_hits.empty:

            s1_rsids = (
                s1_hits[
                    s1_variant_column
                ]
                .map(
                    _normalize_rsid
                )
                .dropna()
            )

            s1_p_values = (
                _numeric(
                    s1_hits[
                        s1_p_column
                    ]
                )
                .dropna()
            )

            s1_fdr_values = (
                _numeric(
                    s1_hits[
                        s1_fdr_column
                    ]
                )
                .dropna()
            )

            s1_unique_rsids = int(
                s1_rsids.nunique()
            )

            s1_minimum_p = (
                None
                if s1_p_values.empty
                else float(
                    s1_p_values.min()
                )
            )

            s1_minimum_fdr = (
                None
                if s1_fdr_values.empty
                else float(
                    s1_fdr_values.min()
                )
            )

        else:

            s1_unique_rsids = 0
            s1_minimum_p = None
            s1_minimum_fdr = None

        s1_support = bool(
            len(
                s1_hits
            )
            >=
            int(
                config[
                    "evidence_policy"
                ][
                    "s1_qtl_support"
                ][
                    "minimum_rows"
                ]
            )
        )

        # ==================================================================
        # S2 — GWAS-LD disease context
        #
        # NOTE:
        # Actual Moradi S2 headers are counterintuitive:
        #
        #   sQTL_SNP      → coordinate
        #   sQTL_SNP-pos  → rsID
        #   tag_SNP       → coordinate
        #   tag_SNP_pos   → rsID
        #
        # The configuration intentionally follows the observed data values.
        # ==================================================================

        s2_feature_mask = (
            s2[
                s2_feature_column
            ]
            .astype("string")
            .str.strip()
            .eq(
                feature_id
            )
        )

        disease_mask = _disease_matches(
            s2[
                str(
                    s2_config[
                        "disease_column"
                    ]
                )
            ],
            config=config,
        )

        s2_hits = s2.loc[
            s2_feature_mask
            &
            disease_mask
        ].copy()

        qtl_variant_column = str(
            s2_config[
                "qtl_variant_column"
            ]
        )

        tag_variant_column = str(
            s2_config[
                "tag_variant_column"
            ]
        )

        ld_column = str(
            s2_config[
                "ld_column"
            ]
        )

        if not s2_hits.empty:

            s2_qtl_rsids = (
                s2_hits[
                    qtl_variant_column
                ]
                .map(
                    _normalize_rsid
                )
                .dropna()
            )

            s2_tag_rsids = (
                s2_hits[
                    tag_variant_column
                ]
                .map(
                    _normalize_rsid
                )
                .dropna()
            )

            s2_ld_values = (
                _numeric(
                    s2_hits[
                        ld_column
                    ]
                )
                .dropna()
            )

            s2_unique_qtl_rsids = int(
                s2_qtl_rsids.nunique()
            )

            s2_unique_tag_rsids = int(
                s2_tag_rsids.nunique()
            )

            s2_maximum_ld = (
                None
                if s2_ld_values.empty
                else float(
                    s2_ld_values.max()
                )
            )

        else:

            s2_unique_qtl_rsids = 0
            s2_unique_tag_rsids = 0
            s2_maximum_ld = None

        s2_support = bool(
            len(
                s2_hits
            )
            >=
            int(
                config[
                    "evidence_policy"
                ][
                    "s2_gwas_ld_context"
                ][
                    "minimum_rows"
                ]
            )
        )

        # ==================================================================
        # S7 — differential intron event
        # ==================================================================

        s7_mask = (
            s7[
                s7_feature_column
            ]
            .astype("string")
            .str.strip()
            .eq(
                feature_id
            )
        )

        s7_hits = s7.loc[
            s7_mask
        ].copy()

        if not s7_hits.empty:

            log2fc_values = (
                _numeric(
                    s7_hits[
                        str(
                            s7_config[
                                "log2_fold_change_column"
                            ]
                        )
                    ]
                )
                .dropna()
            )

            p_values = (
                _numeric(
                    s7_hits[
                        str(
                            s7_config[
                                "p_value_column"
                            ]
                        )
                    ]
                )
                .dropna()
            )

            padj_values = (
                _numeric(
                    s7_hits[
                        str(
                            s7_config[
                                "adjusted_p_value_column"
                            ]
                        )
                    ]
                )
                .dropna()
            )

            differential_log2fc = (
                None
                if log2fc_values.empty
                else float(
                    log2fc_values.iloc[
                        0
                    ]
                )
            )

            differential_p = (
                None
                if p_values.empty
                else float(
                    p_values.iloc[
                        0
                    ]
                )
            )

            differential_padj = (
                None
                if padj_values.empty
                else float(
                    padj_values.iloc[
                        0
                    ]
                )
            )

        else:

            differential_log2fc = None
            differential_p = None
            differential_padj = None

        padj_threshold = float(
            config[
                "evidence_policy"
            ][
                "s7_differential_support"
            ][
                "adjusted_p_value_threshold"
            ]
        )

        s7_support = bool(
            differential_padj is not None
            and
            differential_padj
            <=
            padj_threshold
        )

        # ==================================================================
        # Overall feature confirmation
        # ==================================================================

        evidence_layers = int(
            s1_support
        ) + int(
            s2_support
        ) + int(
            s7_support
        )

        if evidence_layers == 3:

            confirmation = (
                "CONFIRMED_SQTL_GWAS_LD_AND_DIFFERENTIAL_FEATURE"
            )

        elif evidence_layers >= 1:

            confirmation = (
                "PARTIALLY_CONFIRMED_REGULATORY_FEATURE"
            )

        else:

            confirmation = (
                "REGULATORY_FEATURE_NOT_FOUND_IN_EXPECTED_MORADI_SOURCES"
            )

        rows.append(
            {
                "regulatory_feature_id":
                    feature_id,

                "feature_class":
                    _feature_class(
                        feature_id,
                        config=config,
                    ),

                "regulatory_scope":
                    (
                        "CIS"
                    ),

                # ----------------------------------------------------------
                # S1
                # ----------------------------------------------------------

                "s1_qtl_support":
                    s1_support,

                "s1_qtl_rows":
                    int(
                        len(
                            s1_hits
                        )
                    ),

                "s1_unique_qtl_rsids":
                    s1_unique_rsids,

                "s1_minimum_p_value":
                    s1_minimum_p,

                "s1_minimum_fdr":
                    s1_minimum_fdr,

                # ----------------------------------------------------------
                # S2
                # ----------------------------------------------------------

                "s2_gwas_ld_context":
                    s2_support,

                "s2_gwas_ld_rows":
                    int(
                        len(
                            s2_hits
                        )
                    ),

                "s2_unique_qtl_rsids":
                    s2_unique_qtl_rsids,

                "s2_unique_tag_rsids":
                    s2_unique_tag_rsids,

                "s2_maximum_ld":
                    s2_maximum_ld,

                # ----------------------------------------------------------
                # S7
                # ----------------------------------------------------------

                "s7_differential_support":
                    s7_support,

                "s7_differential_rows":
                    int(
                        len(
                            s7_hits
                        )
                    ),

                "s7_log2_fold_change":
                    differential_log2fc,

                "s7_p_value":
                    differential_p,

                "s7_adjusted_p_value":
                    differential_padj,

                # ----------------------------------------------------------
                # Annotation status
                # ----------------------------------------------------------

                "parent_gene":
                    None,

                "parent_gene_status":
                    str(
                        config[
                            "feature_policy"
                        ][
                            "unresolved_parent_gene_status"
                        ]
                    ),

                "parent_gene_inferred":
                    False,

                "nearest_gene_mapping_performed":
                    False,

                # ----------------------------------------------------------
                # Overall
                # ----------------------------------------------------------

                "evidence_layer_count":
                    evidence_layers,

                "overall_feature_confirmation":
                    confirmation,

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Candidate-level mapping
# ============================================================================


def build_candidate_mapping(
    *,
    assignments: pd.DataFrame,
    feature_evidence: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build final M6.2A candidate → regulatory-feature mapping."""

    feature_lookup: dict[
        str,
        dict[str, Any],
    ] = {}

    if not feature_evidence.empty:

        feature_lookup = {
            str(
                row[
                    "regulatory_feature_id"
                ]
            ):
                row
            for row
            in feature_evidence.to_dict(
                orient="records"
            )
        }

    statuses = config[
        "resolution_status"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for assignment in assignments.to_dict(
        orient="records"
    ):

        lead_rsid = str(
            assignment[
                "lead_rsid"
            ]
        )

        feature_present = bool(
            assignment[
                "candidate_feature_assignment_present"
            ]
        )

        feature_id = assignment.get(
            "regulatory_feature_id"
        )

        # ------------------------------------------------------------------
        # No retained Moradi feature
        # ------------------------------------------------------------------

        if not feature_present:

            rows.append(
                {
                    "lead_rsid":
                        lead_rsid,

                    "priority_rank":
                        int(
                            assignment[
                                "priority_rank"
                            ]
                        ),

                    "priority_class":
                        str(
                            assignment[
                                "priority_class"
                            ]
                        ),

                    "regulatory_feature_id":
                        None,

                    "feature_class":
                        None,

                    "regulatory_scope":
                        None,

                    "s1_qtl_support":
                        False,

                    "s2_gwas_ld_context":
                        False,

                    "s7_differential_support":
                        False,

                    "differential_log2_fold_change":
                        None,

                    "differential_adjusted_p_value":
                        None,

                    "parent_gene":
                        None,

                    "parent_gene_status":
                        "NOT_APPLICABLE_NO_RETAINED_FEATURE",

                    "resolution_status":
                        str(
                            statuses[
                                "no_feature_assigned"
                            ]
                        ),

                    "continue_to_gene_annotation":
                        False,

                    "biological_interpretation_status":
                        (
                            "CANDIDATE_RETAINED_WITHOUT_MORADI_FEATURE_MAPPING"
                        ),

                    "causal_claim_allowed":
                        False,
                }
            )

            continue

        # ------------------------------------------------------------------
        # Retained feature
        # ------------------------------------------------------------------

        feature = feature_lookup.get(
            str(
                feature_id
            )
        )

        if feature is None:

            rows.append(
                {
                    "lead_rsid":
                        lead_rsid,

                    "priority_rank":
                        int(
                            assignment[
                                "priority_rank"
                            ]
                        ),

                    "priority_class":
                        str(
                            assignment[
                                "priority_class"
                            ]
                        ),

                    "regulatory_feature_id":
                        feature_id,

                    "feature_class":
                        None,

                    "regulatory_scope":
                        None,

                    "s1_qtl_support":
                        False,

                    "s2_gwas_ld_context":
                        False,

                    "s7_differential_support":
                        False,

                    "differential_log2_fold_change":
                        None,

                    "differential_adjusted_p_value":
                        None,

                    "parent_gene":
                        None,

                    "parent_gene_status":
                        str(
                            statuses[
                                "feature_annotation_pending"
                            ]
                        ),

                    "resolution_status":
                        str(
                            statuses[
                                "feature_not_found"
                            ]
                        ),

                    "continue_to_gene_annotation":
                        False,

                    "biological_interpretation_status":
                        (
                            "ASSIGNED_FEATURE_NOT_CONFIRMED_IN_RAW_SOURCE"
                        ),

                    "causal_claim_allowed":
                        False,
                }
            )

            continue

        confirmation = str(
            feature[
                "overall_feature_confirmation"
            ]
        )

        fully_confirmed = (
            confirmation
            ==
            "CONFIRMED_SQTL_GWAS_LD_AND_DIFFERENTIAL_FEATURE"
        )

        resolution_status = (
            str(
                statuses[
                    "feature_confirmed"
                ]
            )
            if fully_confirmed
            else str(
                statuses[
                    "feature_partially_confirmed"
                ]
            )
        )

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "priority_rank":
                    int(
                        assignment[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        assignment[
                            "priority_class"
                        ]
                    ),

                "regulatory_feature_id":
                    str(
                        feature_id
                    ),

                "feature_class":
                    str(
                        feature[
                            "feature_class"
                        ]
                    ),

                "regulatory_scope":
                    str(
                        feature[
                            "regulatory_scope"
                        ]
                    ),

                "s1_qtl_support":
                    bool(
                        feature[
                            "s1_qtl_support"
                        ]
                    ),

                "s2_gwas_ld_context":
                    bool(
                        feature[
                            "s2_gwas_ld_context"
                        ]
                    ),

                "s7_differential_support":
                    bool(
                        feature[
                            "s7_differential_support"
                        ]
                    ),

                "differential_log2_fold_change":
                    feature[
                        "s7_log2_fold_change"
                    ],

                "differential_adjusted_p_value":
                    feature[
                        "s7_adjusted_p_value"
                    ],

                "parent_gene":
                    None,

                "parent_gene_status":
                    str(
                        feature[
                            "parent_gene_status"
                        ]
                    ),

                "resolution_status":
                    resolution_status,

                "continue_to_gene_annotation":
                    True,

                "biological_interpretation_status":
                    (
                        "REGULATORY_FEATURE_SUPPORTED_GENE_IDENTITY_PENDING"
                    ),

                "causal_claim_allowed":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if not result.empty:

        result = result.sort_values(
            [
                "priority_rank",
                "lead_rsid",
                "regulatory_feature_id",
            ],
            na_position="last",
            kind="stable",
        ).reset_index(
            drop=True
        )

    return result


# ============================================================================
# Summary
# ============================================================================


def build_mapping_summary(
    *,
    candidate_mapping: pd.DataFrame,
    feature_evidence: pd.DataFrame,
) -> pd.DataFrame:
    """Build one-row M6.2A summary."""

    candidates = int(
        candidate_mapping[
            "lead_rsid"
        ].nunique()
    )

    assigned = int(
        candidate_mapping[
            "regulatory_feature_id"
        ]
        .notna()
        .groupby(
            candidate_mapping[
                "lead_rsid"
            ]
        )
        .any()
        .sum()
    )

    no_feature = int(
        candidates
        -
        assigned
    )

    confirmed_features = int(
        (
            feature_evidence[
                "overall_feature_confirmation"
            ]
            ==
            "CONFIRMED_SQTL_GWAS_LD_AND_DIFFERENTIAL_FEATURE"
        ).sum()
        if not feature_evidence.empty
        else 0
    )

    unresolved_gene_features = int(
        feature_evidence[
            "parent_gene"
        ]
        .isna()
        .sum()
        if not feature_evidence.empty
        else 0
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.2A",

                "candidates_assessed":
                    candidates,

                "candidates_with_regulatory_feature":
                    assigned,

                "candidates_without_regulatory_feature":
                    no_feature,

                "unique_regulatory_features_assessed":
                    int(
                        feature_evidence[
                            "regulatory_feature_id"
                        ].nunique()
                        if not feature_evidence.empty
                        else 0
                    ),

                "fully_confirmed_regulatory_features":
                    confirmed_features,

                "features_with_parent_gene_resolved":
                    int(
                        0
                    ),

                "features_with_parent_gene_pending":
                    unresolved_gene_features,

                "parent_gene_inference_performed":
                    False,

                "nearest_gene_mapping_performed":
                    False,

                "formal_colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,

                "next_stage":
                    (
                        "M6.2B_REGULATORY_FEATURE_TO_GENE_ANNOTATION"
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def resolve_regulatory_features(
    *,
    m6_1_qc: dict[str, Any],
    s1: pd.DataFrame,
    s2: pd.DataFrame,
    s7: pd.DataFrame,
    config: dict[str, Any],
) -> RegulatoryFeatureMappingResult:
    """Execute complete M6.2A candidate-regulatory feature resolution."""

    assignments = extract_candidate_feature_assignments(
        m6_1_qc=m6_1_qc,
    )

    feature_evidence = build_feature_evidence(
        assignments=assignments,
        s1=s1,
        s2=s2,
        s7=s7,
        config=config,
    )

    candidate_mapping = build_candidate_mapping(
        assignments=assignments,
        feature_evidence=feature_evidence,
        config=config,
    )

    summary = build_mapping_summary(
        candidate_mapping=candidate_mapping,
        feature_evidence=feature_evidence,
    )

    return RegulatoryFeatureMappingResult(
        feature_evidence=feature_evidence,
        candidate_mapping=candidate_mapping,
        summary=summary,
    )