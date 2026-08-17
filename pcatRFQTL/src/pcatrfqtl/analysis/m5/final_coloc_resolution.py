"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/final_coloc_resolution.py

Description:
    Core logic for M5.6D final colocalization resolution.

    M5.6D integrates locked upstream evidence from:

        M5.5
            Colocalization readiness assessment.

        M5.6A
            Moradi public-data gap resolution.

        M5.6B
            Moradi reconstruction feasibility audit.

        M5.6C.1
            Direct CancerSplicingQTL PRAD dataset inspection.

    The purpose of this stage is to determine whether any valid route to
    formal colocalization remains after all data-gap investigations.

    Expected final state under the current evidence:

        FORMAL_COLOC_NOT_JUSTIFIED_DATA_LIMITATION

    This is a methodological resolution, not a negative biological finding.

    Failure to perform formal colocalization must NOT be interpreted as:
        - evidence against regulatory involvement;
        - evidence against shared locus biology;
        - evidence that the candidate has no QTL effect;
        - evidence that disease and regulatory signals are independent.

    Scientific safeguards:
        - Significant-only QTL resources cannot resolve the formal coloc gap.
        - Approximate reconstruction cannot substitute for the original dense
          QTL association matrix.
        - Missing candidate rows are not interpreted as null effects.
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

import pandas as pd


# ============================================================================
# Final status constants
# ============================================================================


FORMAL_COLOC_NOT_JUSTIFIED_DATA_LIMITATION = (
    "FORMAL_COLOC_NOT_JUSTIFIED_DATA_LIMITATION"
)

FORMAL_COLOC_INPUT_RESOLVED_REASSESSMENT_REQUIRED = (
    "FORMAL_COLOC_INPUT_RESOLVED_REASSESSMENT_REQUIRED"
)


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class FinalColocResolutionResult:
    """Container for M5.6D outputs."""

    evidence_chain: pd.DataFrame

    candidates: pd.DataFrame

    summary: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _nested_get(
    payload: dict[str, Any],
    *keys: str,
    default: Any = None,
) -> Any:
    """Safely retrieve nested JSON values."""

    current: Any = payload

    for key in keys:

        if not isinstance(
            current,
            dict,
        ):

            return default

        if key not in current:

            return default

        current = current[
            key
        ]

    return current


def _as_bool(
    value: Any,
) -> bool:
    """Normalize boolean-like values conservatively."""

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


# ============================================================================
# Evidence-chain construction
# ============================================================================


def build_final_evidence_chain(
    *,
    m5_5: dict[str, Any],
    m5_6a: dict[str, Any],
    m5_6b: dict[str, Any],
    m5_6c_1: dict[str, Any],
) -> pd.DataFrame:
    """Build one row for each locked colocalization decision stage."""

    rows: list[
        dict[str, Any]
    ] = []

    # ----------------------------------------------------------------------
    # M5.5
    # ----------------------------------------------------------------------

    rows.append(
        {
            "stage":
                "M5.5",

            "evidence_source":
                "Moradi_QTL_plus_GWAS_readiness",

            "question":
                "Are current GWAS and QTL inputs suitable for formal coloc?",

            "observed_state":
                (
                    "QTL_FULL_LOCUS_UNAVAILABLE"
                ),

            "formal_coloc_route_available":
                bool(
                    int(
                        _nested_get(
                            m5_5,
                            "summary",
                            "formal_coloc_ready_pairs",
                            default=0,
                        )
                    )
                    >
                    0
                ),

            "blocking_reason":
                (
                    "Moradi QTL resource lacks dense full-locus "
                    "summary statistics."
                ),

            "biological_negative_evidence":
                False,
        }
    )

    # ----------------------------------------------------------------------
    # M5.6A
    # ----------------------------------------------------------------------

    moradi_state = str(
        _nested_get(
            m5_6a,
            "summary",
            "moradi_resolution_state",
            default="UNKNOWN",
        )
    )

    rows.append(
        {
            "stage":
                "M5.6A",

            "evidence_source":
                "Moradi_public_data_resolution",

            "question":
                (
                    "Can the missing full Moradi QTL statistics be obtained "
                    "from a public release?"
                ),

            "observed_state":
                moradi_state,

            "formal_coloc_route_available":
                bool(
                    int(
                        _nested_get(
                            m5_6a,
                            "summary",
                            "formal_coloc_ready_candidates",
                            default=0,
                        )
                    )
                    >
                    0
                ),

            "blocking_reason":
                (
                    "Public full unfiltered Moradi QTL summary statistics "
                    "were not identified."
                ),

            "biological_negative_evidence":
                False,
        }
    )

    # ----------------------------------------------------------------------
    # M5.6B
    # ----------------------------------------------------------------------

    reconstruction_state = str(
        _nested_get(
            m5_6b,
            "summary",
            "reconstruction_readiness",
            default="UNKNOWN",
        )
    )

    rows.append(
        {
            "stage":
                "M5.6B",

            "evidence_source":
                "Moradi_reconstruction_audit",

            "question":
                (
                    "Can the original Moradi QTL analysis be reproduced "
                    "with sufficient fidelity?"
                ),

            "observed_state":
                reconstruction_state,

            "formal_coloc_route_available":
                _as_bool(
                    _nested_get(
                        m5_6b,
                        "summary",
                        "reconstruction_can_proceed",
                        default=False,
                    )
                ),

            "blocking_reason":
                (
                    "Original QTL analysis is theoretically reconstructable "
                    "but not reproducible with sufficient fidelity."
                ),

            "biological_negative_evidence":
                False,
        }
    )

    # ----------------------------------------------------------------------
    # M5.6C.1
    # ----------------------------------------------------------------------

    alternative_state = str(
        _nested_get(
            m5_6c_1,
            "summary",
            "classification",
            default="UNKNOWN",
        )
    )

    rows.append(
        {
            "stage":
                "M5.6C.1",

            "evidence_source":
                "CancerSplicingQTL_PRAD",

            "question":
                (
                    "Can CancerSplicingQTL provide a valid alternative "
                    "dense PRAD QTL input?"
                ),

            "observed_state":
                alternative_state,

            "formal_coloc_route_available":
                _as_bool(
                    _nested_get(
                        m5_6c_1,
                        "summary",
                        "formal_coloc_ready",
                        default=False,
                    )
                ),

            "blocking_reason":
                (
                    "CancerSplicingQTL PRAD is significance-filtered and "
                    "does not provide the complete tested association universe."
                ),

            "biological_negative_evidence":
                False,
        }
    )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Final route resolution
# ============================================================================


def determine_final_resolution(
    *,
    evidence_chain: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[
    str,
    bool,
]:
    """Resolve whether any valid formal-coloc route remains."""

    if evidence_chain.empty:

        raise ValueError(
            "Final colocalization evidence chain is empty."
        )

    route_available = bool(
        evidence_chain[
            "formal_coloc_route_available"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .any()
    )

    policy = config[
        "resolution_policy"
    ]

    if route_available:

        return (
            str(
                policy[
                    "final_status_if_valid_qtl_input_available"
                ]
            ),
            True,
        )

    return (
        str(
            policy[
                "final_status_if_no_valid_qtl_input"
            ]
        ),
        False,
    )


# ============================================================================
# Candidate-level resolution
# ============================================================================


def build_candidate_final_resolution(
    *,
    config: dict[str, Any],
    final_status: str,
    valid_qtl_route_available: bool,
    m5_6c_1: dict[str, Any],
) -> pd.DataFrame:
    """Build final colocalization status for each candidate lead."""

    alternative_candidate_status = {
        str(
            row.get(
                "lead_rsid"
            )
        ).strip().lower():
            row
        for row
        in m5_6c_1.get(
            "candidate_status",
            [],
        )
        if row.get(
            "lead_rsid"
        )
        is not None
    }

    retained_evidence = list(
        config[
            "evidence_preservation"
        ][
            "retain"
        ]
    )

    prohibited_claims = list(
        config[
            "evidence_preservation"
        ][
            "prohibit_claims"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for lead_rsid in config[
        "candidate_leads"
    ]:

        canonical = (
            str(
                lead_rsid
            )
            .strip()
            .lower()
        )

        alternative_row = (
            alternative_candidate_status.get(
                canonical,
                {},
            )
        )

        alternative_direct_present = _as_bool(
            alternative_row.get(
                "direct_rsid_present",
                False,
            )
        )

        rows.append(
            {
                "lead_rsid":
                    canonical,

                "final_coloc_status":
                    final_status,

                "valid_dense_qtl_route_available":
                    valid_qtl_route_available,

                "formal_coloc_ready":
                    False,

                "coloc_abf_ready":
                    False,

                "coloc_susie_ready":
                    False,

                "formal_colocalization_performed":
                    False,

                "posterior_h4_available":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,

                "cancersplicingqtl_direct_record_present":
                    alternative_direct_present,

                "cancersplicingqtl_absence_interpreted_as_null":
                    False,

                "regulatory_evidence_should_be_retained":
                    True,

                "regional_disease_evidence_should_be_retained":
                    True,

                "ld_evidence_should_be_retained":
                    True,

                "retained_evidence_classes":
                    retained_evidence,

                "prohibited_claim_classes":
                    prohibited_claims,

                "recommended_analysis_framework":
                    (
                        "LOCUS_LEVEL_CONVERGENT_EVIDENCE_WITHOUT_FORMAL_COLOC"
                    ),

                "recommended_manuscript_interpretation":
                    (
                        "Formal colocalization was not justified because "
                        "suitable dense feature-specific QTL summary statistics "
                        "were unavailable. Regional disease, LD, and regulatory "
                        "evidence should therefore be interpreted as convergent "
                        "locus-level evidence rather than evidence of a shared "
                        "causal variant."
                    ),

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Final summary
# ============================================================================


def build_final_summary(
    *,
    evidence_chain: pd.DataFrame,
    candidates: pd.DataFrame,
    final_status: str,
    valid_qtl_route_available: bool,
) -> pd.DataFrame:
    """Build one-row M5.6D summary."""

    routes_assessed = int(
        len(
            evidence_chain
        )
    )

    routes_available = int(
        evidence_chain[
            "formal_coloc_route_available"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M5.6D",

                "candidate_leads_assessed":
                    int(
                        len(
                            candidates
                        )
                    ),

                "formal_coloc_routes_assessed":
                    routes_assessed,

                "formal_coloc_routes_available":
                    routes_available,

                "valid_dense_qtl_route_available":
                    valid_qtl_route_available,

                "final_coloc_status":
                    final_status,

                "formal_coloc_ready_candidates":
                    int(
                        candidates[
                            "formal_coloc_ready"
                        ]
                        .fillna(
                            False
                        )
                        .sum()
                    ),

                "formal_coloc_blocked_candidates":
                    int(
                        len(
                            candidates
                        )
                    ),

                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,

                "final_interpretation":
                    (
                        "Formal colocalization is not justified under the "
                        "available QTL evidence. Existing disease, LD, and "
                        "regulatory evidence remains valid as locus-level "
                        "convergent evidence."
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def resolve_final_coloc_status(
    *,
    config: dict[str, Any],
    m5_5: dict[str, Any],
    m5_6a: dict[str, Any],
    m5_6b: dict[str, Any],
    m5_6c_1: dict[str, Any],
) -> FinalColocResolutionResult:
    """Execute complete M5.6D final colocalization resolution."""

    evidence_chain = (
        build_final_evidence_chain(
            m5_5=m5_5,
            m5_6a=m5_6a,
            m5_6b=m5_6b,
            m5_6c_1=m5_6c_1,
        )
    )

    (
        final_status,
        valid_qtl_route_available,
    ) = determine_final_resolution(
        evidence_chain=evidence_chain,
        config=config,
    )

    candidates = (
        build_candidate_final_resolution(
            config=config,
            final_status=final_status,
            valid_qtl_route_available=valid_qtl_route_available,
            m5_6c_1=m5_6c_1,
        )
    )

    summary = (
        build_final_summary(
            evidence_chain=evidence_chain,
            candidates=candidates,
            final_status=final_status,
            valid_qtl_route_available=valid_qtl_route_available,
        )
    )

    return FinalColocResolutionResult(
        evidence_chain=evidence_chain,
        candidates=candidates,
        summary=summary,
    )