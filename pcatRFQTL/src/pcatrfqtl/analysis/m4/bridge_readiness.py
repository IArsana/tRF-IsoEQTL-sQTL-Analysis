"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/bridge_readiness.py

Description:
    M4.4 Regulatory Bridge Readiness assessment.

    Determines whether each SNP-tRF candidate already has sufficient
    exact multi-omic evidence or requires locus-aware evaluation in M5.

    Bridge-readiness classification is workflow routing only.
    It is not a biological candidate ranking system.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from enum import StrEnum

import pandas as pd


class RegulatoryBridgeState(StrEnum):
    """Workflow readiness states."""

    EXACT_BRIDGE_AVAILABLE = (
        "EXACT_BRIDGE_AVAILABLE"
    )

    PARTIAL_EXACT_BRIDGE = (
        "PARTIAL_EXACT_BRIDGE"
    )

    LOCUS_BRIDGE_REQUIRED = (
        "LOCUS_BRIDGE_REQUIRED"
    )


class M44RegulatoryBridgeReadinessBuilder:
    """Assess readiness for M5 locus-aware integration."""

    REQUIRED_COLUMNS = frozenset(
        {
            "m4_candidate_id",
            "trfqtl_rsid",
            "trf_id",
            "direct_gwas_match",
            "moradi_direct_qtl_match",
            "moradi_source_tag_match",
            "cis_exact_support",
            "trans_exact_support",
            "exact_multiomic_support",
        }
    )

    @classmethod
    def build(
        cls,
        evidence: pd.DataFrame,
    ) -> pd.DataFrame:
        """Build candidate-level regulatory bridge assessment."""

        missing = (
            cls.REQUIRED_COLUMNS
            - set(
                evidence.columns
            )
        )

        if missing:
            raise ValueError(
                "M4.4 input is missing required columns: "
                f"{sorted(missing)}"
            )

        result = evidence.copy()

        result[
            "regulatory_exact_support"
        ] = (
            result[
                "moradi_direct_qtl_match"
            ]
            | result[
                "moradi_source_tag_match"
            ]
        ).astype(
            "boolean"
        )

        result[
            "bridge_has_gwas_anchor"
        ] = result[
            "direct_gwas_match"
        ].astype(
            "boolean"
        )

        result[
            "bridge_has_regulatory_anchor"
        ] = result[
            "regulatory_exact_support"
        ].astype(
            "boolean"
        )

        result[
            "bridge_has_cis_anchor"
        ] = result[
            "cis_exact_support"
        ].astype(
            "boolean"
        )

        result[
            "bridge_has_trans_anchor"
        ] = result[
            "trans_exact_support"
        ].astype(
            "boolean"
        )

        result[
            "bridge_exact_layer_count"
        ] = (
            result[
                "bridge_has_gwas_anchor"
            ].astype(
                "Int64"
            )
            + result[
                "bridge_has_regulatory_anchor"
            ].astype(
                "Int64"
            )
        )

        result[
            "regulatory_bridge_state"
        ] = (
            RegulatoryBridgeState
            .LOCUS_BRIDGE_REQUIRED
            .value
        )

        partial = (
            result[
                "bridge_exact_layer_count"
            ]
            == 1
        )

        result.loc[
            partial,
            "regulatory_bridge_state",
        ] = (
            RegulatoryBridgeState
            .PARTIAL_EXACT_BRIDGE
            .value
        )

        exact = (
            result[
                "exact_multiomic_support"
            ]
            .fillna(False)
            .eq(True)
        )

        result.loc[
            exact,
            "regulatory_bridge_state",
        ] = (
            RegulatoryBridgeState
            .EXACT_BRIDGE_AVAILABLE
            .value
        )

        # --------------------------------------------------------------
        # M5 workflow routing
        # --------------------------------------------------------------

        result[
            "m5_requires_ld_lookup"
        ] = (
            result[
                "regulatory_bridge_state"
            ]
            != (
                RegulatoryBridgeState
                .EXACT_BRIDGE_AVAILABLE
                .value
            )
        ).astype(
            "boolean"
        )

        result[
            "m5_requires_locus_mapping"
        ] = result[
            "m5_requires_ld_lookup"
        ].astype(
            "boolean"
        )

        result[
            "m5_requires_colocalization_assessment"
        ] = result[
            "m5_requires_locus_mapping"
        ].astype(
            "boolean"
        )

        result[
            "m5_search_qtl_neighbors"
        ] = (
            ~result[
                "bridge_has_regulatory_anchor"
            ]
        ).astype(
            "boolean"
        )

        result[
            "m5_search_gwas_neighbors"
        ] = (
            ~result[
                "bridge_has_gwas_anchor"
            ]
        ).astype(
            "boolean"
        )

        result[
            "m5_search_cis_regulatory_bridge"
        ] = (
            ~result[
                "bridge_has_cis_anchor"
            ]
        ).astype(
            "boolean"
        )

        result[
            "m5_search_trans_regulatory_bridge"
        ] = (
            ~result[
                "bridge_has_trans_anchor"
            ]
        ).astype(
            "boolean"
        )

        # --------------------------------------------------------------
        # Explicit interpretation safeguards
        # --------------------------------------------------------------

        result[
            "bridge_ld_already_computed"
        ] = False

        result[
            "bridge_colocalization_already_computed"
        ] = False

        result[
            "bridge_causal_relationship_established"
        ] = False

        result[
            "bridge_ranked"
        ] = False

        result[
            "bridge_readiness_reason"
        ] = (
            "Exact GWAS-to-tRF-QTL-to-regulatory bridge is absent; "
            "locus-aware LD and colocalization assessment is required."
        )

        result.loc[
            partial,
            "bridge_readiness_reason",
        ] = (
            "One exact external evidence layer is available, but the "
            "complete multi-omic bridge remains incomplete."
        )

        result.loc[
            exact,
            "bridge_readiness_reason",
        ] = (
            "Exact GWAS and regulatory evidence are both available "
            "for the tRF-QTL variant."
        )

        return result