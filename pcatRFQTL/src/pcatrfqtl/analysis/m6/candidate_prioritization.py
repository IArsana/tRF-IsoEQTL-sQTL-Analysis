"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/candidate_prioritization.py

Description:
    Core logic for M6.1 final candidate prioritization.

    M6.1 converts the locked M3-M5 evidence into a deterministic biological
    interpretation priority order.

    The ranking is lexicographic rather than based on arbitrary additive
    point scoring.

    Primary hierarchy:
        1. Regional genome-wide disease evidence.
        2. Multi-study suggestive disease evidence.
        3. Single-study suggestive disease evidence.
        4. No regional disease evidence.

    Regulatory evidence is used only as a secondary evidence descriptor /
    tie-breaker and must not be interpreted as proof of causality.

    All retained tRF-QTL candidates remain eligible for subsequent M6
    biological interpretation regardless of priority rank.

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
class CandidatePrioritizationResult:
    """Container for M6.1 outputs."""

    evidence_matrix: pd.DataFrame

    priorities: pd.DataFrame

    summary: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _as_bool(
    value: Any,
) -> bool:
    """Conservatively normalize boolean-like values."""

    if isinstance(
        value,
        bool,
    ):

        return value

    if value is None:

        return False

    if isinstance(
        value,
        str,
    ):

        normalized = (
            value
            .strip()
            .lower()
        )

        if normalized in {
            "true",
            "1",
            "yes",
        }:

            return True

        if normalized in {
            "false",
            "0",
            "no",
            "",
        }:

            return False

    return bool(
        value
    )


def _validate_rsid(
    rsid: str,
) -> str:
    """Normalize and validate candidate rsID."""

    canonical = (
        str(
            rsid
        )
        .strip()
        .lower()
    )

    if not canonical.startswith(
        "rs"
    ):

        raise ValueError(
            f"Invalid candidate rsID: {rsid}"
        )

    if not canonical[
        2:
    ].isdigit():

        raise ValueError(
            f"Invalid candidate rsID: {rsid}"
        )

    return canonical


# ============================================================================
# M5.6D validation
# ============================================================================


def validate_final_coloc_candidates(
    *,
    config: dict[str, Any],
    final_coloc_qc: dict[str, Any],
) -> None:
    """
    Verify that M6.1 candidate universe matches the locked M5.6D candidates.
    """

    configured = {
        _validate_rsid(
            rsid
        )
        for rsid
        in config[
            "candidate_leads"
        ]
    }

    upstream_rows = final_coloc_qc.get(
        "candidate_status",
        [],
    )

    upstream = {
        _validate_rsid(
            row[
                "lead_rsid"
            ]
        )
        for row
        in upstream_rows
        if row.get(
            "lead_rsid"
        )
        is not None
    }

    if configured != upstream:

        raise RuntimeError(
            "M6.1 candidate universe does not match locked M5.6D. "
            f"Configured={sorted(configured)} | "
            f"M5.6D={sorted(upstream)}"
        )

    for row in upstream_rows:

        final_status = str(
            row.get(
                "final_coloc_status",
                ""
            )
        )

        if final_status != (
            "FORMAL_COLOC_NOT_JUSTIFIED_DATA_LIMITATION"
        ):

            raise RuntimeError(
                "Unexpected M5.6D colocalization status for "
                f"{row.get('lead_rsid')}: {final_status}"
            )


# ============================================================================
# Evidence matrix
# ============================================================================


def build_candidate_evidence_matrix(
    *,
    config: dict[str, Any],
    final_coloc_qc: dict[str, Any],
) -> pd.DataFrame:
    """Build standardized evidence matrix for all final candidates."""

    upstream_candidates = {
        _validate_rsid(
            row[
                "lead_rsid"
            ]
        ):
            row
        for row
        in final_coloc_qc.get(
            "candidate_status",
            [],
        )
        if row.get(
            "lead_rsid"
        )
        is not None
    }

    rows: list[
        dict[str, Any]
    ] = []

    for raw_rsid, evidence in config[
        "candidate_leads"
    ].items():

        rsid = _validate_rsid(
            raw_rsid
        )

        disease = evidence[
            "disease_evidence"
        ]

        regulatory = evidence[
            "regulatory_evidence"
        ]

        ld = evidence[
            "ld_evidence"
        ]

        trfqtl = evidence[
            "trfqtl_candidate"
        ]

        coloc = upstream_candidates[
            rsid
        ]

        rows.append(
            {
                "lead_rsid":
                    rsid,

                # ----------------------------------------------------------
                # tRF-QTL
                # ----------------------------------------------------------

                "trfqtl_candidate_retained":
                    _as_bool(
                        trfqtl.get(
                            "retained"
                        )
                    ),

                # ----------------------------------------------------------
                # Disease evidence
                # ----------------------------------------------------------

                "disease_evidence_class":
                    str(
                        disease[
                            "class"
                        ]
                    ),

                "disease_architecture":
                    str(
                        disease[
                            "architecture"
                        ]
                    ),

                "regional_genome_wide":
                    _as_bool(
                        disease.get(
                            "regional_genome_wide"
                        )
                    ),

                "recurrent_suggestive":
                    _as_bool(
                        disease.get(
                            "recurrent_suggestive"
                        )
                    ),

                "single_study_suggestive":
                    _as_bool(
                        disease.get(
                            "single_study_suggestive"
                        )
                    ),

                "direct_lead_suggestive":
                    _as_bool(
                        disease.get(
                            "direct_lead_disease_association",
                            {},
                        ).get(
                            "suggestive"
                        )
                    ),

                # ----------------------------------------------------------
                # Regulatory evidence
                # ----------------------------------------------------------

                "moradi_locus_evidence":
                    _as_bool(
                        regulatory.get(
                            "moradi_locus_evidence"
                        )
                    ),

                "moradi_direct_or_proxy_bridge":
                    _as_bool(
                        regulatory.get(
                            "moradi_direct_or_proxy_bridge"
                        )
                    ),

                "regulatory_feature_ids":
                    list(
                        regulatory.get(
                            "regulatory_feature_ids",
                            [],
                        )
                    ),

                # ----------------------------------------------------------
                # LD
                # ----------------------------------------------------------

                "candidate_reference_ld_available":
                    _as_bool(
                        ld.get(
                            "candidate_reference_ld_available"
                        )
                    ),

                "disease_to_regulatory_proxy_ld_bridge":
                    _as_bool(
                        ld.get(
                            "disease_to_regulatory_proxy_ld_bridge"
                        )
                    ),

                # ----------------------------------------------------------
                # Final M5.6D state
                # ----------------------------------------------------------

                "final_coloc_status":
                    str(
                        coloc.get(
                            "final_coloc_status"
                        )
                    ),

                "formal_coloc_ready":
                    _as_bool(
                        coloc.get(
                            "formal_coloc_ready"
                        )
                    ),

                "formal_colocalization_performed":
                    _as_bool(
                        coloc.get(
                            "formal_colocalization_performed"
                        )
                    ),

                "regulatory_evidence_retained":
                    _as_bool(
                        coloc.get(
                            "regulatory_evidence_should_be_retained"
                        )
                    ),

                "regional_disease_evidence_retained":
                    _as_bool(
                        coloc.get(
                            "regional_disease_evidence_should_be_retained"
                        )
                    ),

                "ld_evidence_retained":
                    _as_bool(
                        coloc.get(
                            "ld_evidence_should_be_retained"
                        )
                    ),

                "causal_claim_allowed":
                    _as_bool(
                        coloc.get(
                            "causal_claim_allowed"
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Deterministic prioritization
# ============================================================================


def prioritize_candidates(
    evidence_matrix: pd.DataFrame,
    *,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Rank candidates using deterministic lexicographic evidence hierarchy.
    """

    if evidence_matrix.empty:

        raise ValueError(
            "Candidate evidence matrix is empty."
        )

    hierarchy = config[
        "priority_hierarchy"
    ][
        "disease_evidence_order"
    ]

    working = evidence_matrix.copy()

    working[
        "_disease_order"
    ] = (
        working[
            "disease_evidence_class"
        ]
        .map(
            hierarchy
        )
    )

    if working[
        "_disease_order"
    ].isna().any():

        unknown = (
            working.loc[
                working[
                    "_disease_order"
                ].isna(),
                "disease_evidence_class",
            ]
            .astype(
                str
            )
            .unique()
            .tolist()
        )

        raise ValueError(
            "Unknown disease evidence classes: "
            f"{unknown}"
        )

    retained_required = _as_bool(
        config[
            "priority_hierarchy"
        ][
            "trfqtl_requirement"
        ][
            "retained_candidate_required"
        ]
    )

    if retained_required:

        nonretained = working.loc[
            ~working[
                "trfqtl_candidate_retained"
            ].astype(
                bool
            ),
            "lead_rsid",
        ].tolist()

        if nonretained:

            raise RuntimeError(
                "M6.1 contains candidate(s) not retained by "
                f"tRF-QTL selection: {nonretained}"
            )

    # ----------------------------------------------------------------------
    # Lexicographic ranking
    #
    # 1. Disease evidence hierarchy.
    # 2. Regulatory direct/proxy bridge.
    # 3. Regulatory locus evidence.
    # 4. rsID only as deterministic final tie-break.
    #
    # No additive evidence score is constructed.
    # ----------------------------------------------------------------------

    working = working.sort_values(
        by=[
            "_disease_order",
            "moradi_direct_or_proxy_bridge",
            "moradi_locus_evidence",
            "lead_rsid",
        ],
        ascending=[
            True,
            False,
            False,
            True,
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )

    working[
        "priority_rank"
    ] = range(
        1,
        len(
            working
        )
        + 1,
    )

    class_config = config[
        "priority_classes"
    ]

    rank_to_class = {
        1:
            str(
                class_config[
                    "rank_1"
                ][
                    "label"
                ]
            ),

        2:
            str(
                class_config[
                    "rank_2"
                ][
                    "label"
                ]
            ),

        3:
            str(
                class_config[
                    "rank_3"
                ][
                    "label"
                ]
            ),
    }

    working[
        "priority_class"
    ] = working[
        "priority_rank"
    ].map(
        rank_to_class
    )

    working[
        "ranking_method"
    ] = (
        "DETERMINISTIC_LEXICOGRAPHIC_EVIDENCE_HIERARCHY"
    )

    working[
        "additive_score_used"
    ] = False

    working[
        "continue_to_m6_biological_interpretation"
    ] = True

    working[
        "deep_annotation_priority"
    ] = (
        working[
            "priority_rank"
        ]
        ==
        1
    )

    working[
        "priority_interpretation"
    ] = working.apply(
        _build_priority_interpretation,
        axis=1,
    )

    return working.drop(
        columns=[
            "_disease_order",
        ]
    )


# ============================================================================
# Candidate interpretation
# ============================================================================


def _build_priority_interpretation(
    row: pd.Series,
) -> str:
    """Build conservative candidate-level interpretation."""

    evidence_class = str(
        row[
            "disease_evidence_class"
        ]
    )

    regulatory_bridge = bool(
        row[
            "moradi_direct_or_proxy_bridge"
        ]
    )

    if (
        evidence_class
        ==
        "REGIONAL_GENOME_WIDE"
        and
        regulatory_bridge
    ):

        return (
            "Highest-priority locus based on regional genome-wide disease "
            "evidence with additional regulatory convergence. The evidence "
            "does not establish formal colocalization or causality."
        )

    if (
        evidence_class
        ==
        "MULTI_STUDY_SUGGESTIVE"
    ):

        return (
            "Intermediate-priority locus supported by recurrent suggestive "
            "regional disease evidence. Formal colocalization and causal "
            "regulatory inference are not supported."
        )

    if (
        evidence_class
        ==
        "SINGLE_STUDY_SUGGESTIVE"
    ):

        return (
            "Lower-priority retained locus supported by single-study "
            "suggestive regional disease evidence. It remains eligible for "
            "biological annotation but requires stronger validation."
        )

    return (
        "Retained candidate for biological interpretation under the current "
        "integrated evidence framework."
    )


# ============================================================================
# Summary
# ============================================================================


def build_priority_summary(
    priorities: pd.DataFrame,
) -> pd.DataFrame:
    """Build one-row M6.1 summary."""

    if priorities.empty:

        raise ValueError(
            "Prioritized candidate table is empty."
        )

    top = priorities.iloc[
        0
    ]

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.1",

                "candidates_assessed":
                    int(
                        len(
                            priorities
                        )
                    ),

                "candidates_retained":
                    int(
                        priorities[
                            "continue_to_m6_biological_interpretation"
                        ]
                        .astype(
                            bool
                        )
                        .sum()
                    ),

                "highest_priority_candidate":
                    str(
                        top[
                            "lead_rsid"
                        ]
                    ),

                "highest_priority_class":
                    str(
                        top[
                            "priority_class"
                        ]
                    ),

                "formal_coloc_ready_candidates":
                    int(
                        priorities[
                            "formal_coloc_ready"
                        ]
                        .astype(
                            bool
                        )
                        .sum()
                    ),

                "causal_claim_allowed_candidates":
                    int(
                        priorities[
                            "causal_claim_allowed"
                        ]
                        .astype(
                            bool
                        )
                        .sum()
                    ),

                "ranking_method":
                    (
                        "DETERMINISTIC_LEXICOGRAPHIC_EVIDENCE_HIERARCHY"
                    ),

                "additive_score_used":
                    False,

                "candidate_exclusion_performed":
                    False,

                "next_stage":
                    (
                        "M6.2_CANDIDATE_TO_REGULATORY_FEATURE_MAPPING"
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def run_candidate_prioritization(
    *,
    config: dict[str, Any],
    final_coloc_qc: dict[str, Any],
) -> CandidatePrioritizationResult:
    """Execute M6.1 final candidate prioritization."""

    validate_final_coloc_candidates(
        config=config,
        final_coloc_qc=final_coloc_qc,
    )

    evidence_matrix = (
        build_candidate_evidence_matrix(
            config=config,
            final_coloc_qc=final_coloc_qc,
        )
    )

    priorities = (
        prioritize_candidates(
            evidence_matrix,
            config=config,
        )
    )

    summary = (
        build_priority_summary(
            priorities
        )
    )

    return CandidatePrioritizationResult(
        evidence_matrix=evidence_matrix,
        priorities=priorities,
        summary=summary,
    )