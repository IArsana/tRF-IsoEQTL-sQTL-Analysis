"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/coloc_data_gap.py

Description:
    Core logic for M5.6 colocalization data-gap resolution.

    M5.6 evaluates whether the regulatory-QTL limitations identified during
    M5.5 can be resolved using:

        1. the original Moradi public resource;
        2. reconstruction from underlying source data;
        3. an alternative prostate regulatory-QTL resource.

    This stage does not perform colocalization.

    Resolution states:
        MORADI_FULL_QTL_PUBLICLY_AVAILABLE
        MORADI_FULL_QTL_RECONSTRUCTABLE
        MORADI_FULL_QTL_NOT_IDENTIFIED_PUBLICLY
        ALTERNATIVE_QTL_RESOURCE_AVAILABLE
        NO_VALID_COLOC_INPUT_RESOLUTION

    Scientific safeguards:
        - Significant-only QTL tables are not promoted to full-locus data.
        - Absence of a QTL row is not interpreted as a null association.
        - No synthetic standard errors are created.
        - No missing QTL associations are imputed.
        - No formal colocalization is performed.
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
# Resolution states
# ============================================================================


MORADI_FULL_QTL_PUBLICLY_AVAILABLE = (
    "MORADI_FULL_QTL_PUBLICLY_AVAILABLE"
)

MORADI_FULL_QTL_RECONSTRUCTABLE = (
    "MORADI_FULL_QTL_RECONSTRUCTABLE"
)

MORADI_FULL_QTL_NOT_IDENTIFIED_PUBLICLY = (
    "MORADI_FULL_QTL_NOT_IDENTIFIED_PUBLICLY"
)

ALTERNATIVE_QTL_RESOURCE_AVAILABLE = (
    "ALTERNATIVE_QTL_RESOURCE_AVAILABLE"
)

NO_VALID_COLOC_INPUT_RESOLUTION = (
    "NO_VALID_COLOC_INPUT_RESOLUTION"
)


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class DataGapResolutionResult:
    """Container for M5.6 resolution output."""

    sources: pd.DataFrame

    resolution: pd.DataFrame

    candidates: pd.DataFrame


# ============================================================================
# Source audit
# ============================================================================


def build_source_audit(
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build deterministic source-level M5.6 audit."""

    moradi = config[
        "moradi"
    ]

    reconstruction = config[
        "reconstruction"
    ]

    alternative = config[
        "alternative_resources"
    ][
        "cancer_splicing_qtl"
    ]

    rows = [
        {
            "source_id":
                "MORADI_SUPPLEMENT",

            "source_type":
                "ORIGINAL_QTL_SUPPLEMENT",

            "cancer":
                "PRAD",

            "genome_build":
                moradi[
                    "genome_build"
                ],

            "qtl_scope":
                moradi[
                    "available_public_artifact"
                ][
                    "scope"
                ],

            "full_locus_statistics":
                bool(
                    moradi[
                        "available_public_artifact"
                    ][
                        "full_unfiltered_association_statistics"
                    ]
                ),

            "feature_specific":
                True,

            "formal_coloc_ready":
                False,

            "role":
                "CURRENT_REGULATORY_EVIDENCE",

            "blocking_reason":
                (
                    "Public artifact contains significant QTL associations "
                    "rather than complete unfiltered feature-specific "
                    "locus-wide association statistics."
                ),
        },

        {
            "source_id":
                "MORADI_RECONSTRUCTION",

            "source_type":
                "RECONSTRUCTION_PATH",

            "cancer":
                "PRAD",

            "genome_build":
                moradi[
                    "genome_build"
                ],

            "qtl_scope":
                "RECOMPUTED_FULL_ASSOCIATION_SCAN",

            "full_locus_statistics":
                bool(
                    reconstruction[
                        "theoretically_possible"
                    ]
                ),

            "feature_specific":
                True,

            "formal_coloc_ready":
                False,

            "role":
                "POTENTIAL_FUTURE_RESOLUTION",

            "blocking_reason":
                (
                    "Exact reconstruction requires underlying genotype, "
                    "expression/splicing phenotypes, covariates, QC filters, "
                    "and reproduction of the original MatrixEQTL model."
                ),
        },

        {
            "source_id":
                "CANCER_SPLICING_QTL_PRAD",

            "source_type":
                "ALTERNATIVE_QTL_RESOURCE",

            "cancer":
                alternative[
                    "cancer"
                ],

            "genome_build":
                None,

            "qtl_scope":
                alternative[
                    "data_type"
                ],

            "full_locus_statistics":
                False,

            "feature_specific":
                True,

            "formal_coloc_ready":
                bool(
                    alternative[
                        "formal_coloc_status"
                    ][
                        "ready"
                    ]
                ),

            "role":
                alternative[
                    "role"
                ],

            "blocking_reason":
                alternative[
                    "formal_coloc_status"
                ][
                    "reason"
                ],
        },
    ]

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Resolution
# ============================================================================


def resolve_moradi_gap(
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve the current Moradi colocalization-data gap."""

    moradi = config[
        "moradi"
    ]

    public_status = (
        str(
            moradi[
                "public_full_summary_statistics"
            ][
                "status"
            ]
        )
        .strip()
        .lower()
    )

    reconstruction_enabled = bool(
        config[
            "reconstruction"
        ][
            "enable_in_current_pipeline"
        ]
    )

    reconstruction_possible = bool(
        config[
            "reconstruction"
        ][
            "theoretically_possible"
        ]
    )

    if public_status == "available":

        state = (
            MORADI_FULL_QTL_PUBLICLY_AVAILABLE
        )

        formal_coloc_next = True

    elif (
        reconstruction_enabled
        and
        reconstruction_possible
    ):

        state = (
            MORADI_FULL_QTL_RECONSTRUCTABLE
        )

        formal_coloc_next = False

    else:

        state = (
            MORADI_FULL_QTL_NOT_IDENTIFIED_PUBLICLY
        )

        formal_coloc_next = False

    return pd.DataFrame(
        [
            {
                "resource":
                    "Moradi_2022",

                "resolution_state":
                    state,

                "public_full_qtl_statistics_identified":
                    bool(
                        public_status
                        ==
                        "available"
                    ),

                "reconstruction_theoretically_possible":
                    reconstruction_possible,

                "reconstruction_enabled":
                    reconstruction_enabled,

                "formal_coloc_can_proceed":
                    formal_coloc_next,

                "recommended_action":
                    (
                        "HARMONIZE_AND_REASSESS_COLOC"
                        if formal_coloc_next
                        else
                        "DOCUMENT_QTL_DATA_LIMITATION"
                    ),
            }
        ]
    )


# ============================================================================
# Candidate-level resolution
# ============================================================================


def build_candidate_resolution(
    config: dict[str, Any],
    *,
    moradi_resolution: pd.DataFrame,
) -> pd.DataFrame:
    """Propagate source-level resolution to candidate leads."""

    if moradi_resolution.empty:

        raise ValueError(
            "Moradi resolution table is empty."
        )

    source_state = str(
        moradi_resolution.iloc[
            0
        ][
            "resolution_state"
        ]
    )

    formal_coloc = bool(
        moradi_resolution.iloc[
            0
        ][
            "formal_coloc_can_proceed"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for lead in config[
        "candidate_leads"
    ]:

        rows.append(
            {
                "lead_rsid":
                    str(
                        lead
                    ),

                "moradi_resolution_state":
                    source_state,

                "formal_coloc_ready":
                    formal_coloc,

                "coloc_abf_ready":
                    formal_coloc,

                "coloc_susie_ready":
                    False,

                "current_resolution":
                    (
                        "QTL_DATA_GAP_RESOLVED"
                        if formal_coloc
                        else
                        "QTL_DATA_GAP_UNRESOLVED"
                    ),

                "recommended_analysis":
                    (
                        "FORMAL_COLOCATION_READINESS_REASSESSMENT"
                        if formal_coloc
                        else
                        "LOCUS_LEVEL_CONVERGENCE_ONLY"
                    ),

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Public API
# ============================================================================


def resolve_coloc_data_gap(
    config: dict[str, Any],
) -> DataGapResolutionResult:
    """Run deterministic M5.6 resolution logic."""

    sources = build_source_audit(
        config
    )

    resolution = resolve_moradi_gap(
        config
    )

    candidates = build_candidate_resolution(
        config,
        moradi_resolution=resolution,
    )

    return DataGapResolutionResult(
        sources=sources,
        resolution=resolution,
        candidates=candidates,
    )