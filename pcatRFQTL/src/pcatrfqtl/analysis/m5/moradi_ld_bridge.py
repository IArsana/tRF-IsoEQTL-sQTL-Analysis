"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/moradi_ld_bridge.py

Description:
    Builds M5.4 LD-aware regulatory bridges between tRF-QTL lead
    variants and the Unified Moradi QTL Index.

    Two regulatory bridge targets are kept distinct:

        LD_TO_MORADI_QTL
            The external LD proxy matches a Moradi direct QTL variant.

        LD_TO_MORADI_SOURCE_TAG
            The external LD proxy matches a tag variant reported by
            Moradi.

    The second category must not be interpreted as newly calculated
    lead-to-QTL LD through the Moradi tag relationship. It represents
    two separately recorded evidence relationships.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import pandas as pd


class M54MoradiLDBridgeBuilder:
    """Build LD-aware Moradi regulatory evidence."""

    REQUIRED_MORADI_COLUMNS = frozenset(
        {
            "m4_source_table",
            "m4_source_sheet",
            "m4_regulatory_scope",
            "m4_feature_class",
            "qtl_rsid",
            "qtl_rsid_usable",
            "tag_rsid",
            "tag_rsid_usable",
            "feature_id",
            "feature_identity_key",
        }
    )

    @staticmethod
    def _require_columns(
        dataframe: pd.DataFrame,
        columns: frozenset[str],
    ) -> None:
        """Require common Moradi index schema."""

        missing = (
            columns
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                "M5.4 Moradi index missing columns: "
                f"{sorted(missing)}"
            )

    @staticmethod
    def _candidate_ld(
        loci: pd.DataFrame,
        ld_evidence: pd.DataFrame,
    ) -> pd.DataFrame:
        """Attach imported LD proxies to candidate identities."""

        candidate_columns = [
            column
            for column
            in (
                "m5_locus_id",
                "m4_candidate_id",
                "trfqtl_rsid",
                "trf_id",
                "trf_identity_key",
                "lead_chromosome",
                "lead_position",
            )
            if column
            in loci.columns
        ]

        candidates = loci[
            candidate_columns
        ].copy()

        candidates[
            "trfqtl_rsid"
        ] = (
            candidates[
                "trfqtl_rsid"
            ]
            .astype(
                "string"
            )
            .str.lower()
        )

        usable_ld = (
            ld_evidence.loc[
                ld_evidence[
                    "m5_bridge_eligible"
                ]
                .fillna(False)
                .eq(True)
            ]
            .copy()
        )

        if usable_ld.empty:
            return pd.DataFrame()

        return candidates.merge(
            usable_ld,
            how="inner",
            left_on="trfqtl_rsid",
            right_on="lead_rsid",
            validate="one_to_many",
        )

    @staticmethod
    def _moradi_projection(
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """Project regulatory evidence fields."""

        preferred = (
            "m4_source_dataset",
            "m4_source_table",
            "m4_source_sheet",
            "m4_regulatory_scope",
            "m4_feature_class",
            "m4_source_row",
            "qtl_rsid",
            "qtl_rsid_usable",
            "qtl_variant_identity_key",
            "qtl_variant_status",
            "tag_rsid",
            "tag_rsid_usable",
            "tag_variant_identity_key",
            "tag_variant_status",
            "feature_id",
            "feature_identity_key",
            "feature_type",
            "feature_usable",
            "feature_status",
            "disease_id",
            "is_primary_disease",
            "p_value",
            "fdr",
            "beta",
            "ld_r2",
        )

        return moradi[
            [
                column
                for column
                in preferred
                if column
                in moradi.columns
            ]
        ].copy()

    # ------------------------------------------------------------------
    # Direct QTL bridge
    # ------------------------------------------------------------------

    @classmethod
    def _direct_qtl_matches(
        cls,
        candidate_ld: pd.DataFrame,
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """Match LD proxies against direct Moradi QTL rsIDs."""

        evidence = moradi.loc[
            moradi[
                "qtl_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            & moradi[
                "qtl_rsid"
            ]
            .notna()
        ].copy()

        if evidence.empty:
            return pd.DataFrame()

        evidence = cls._moradi_projection(
            evidence
        )

        matched = candidate_ld.merge(
            evidence,
            how="inner",
            left_on="neighbor_rsid",
            right_on="qtl_rsid",
            validate="many_to_many",
        )

        if matched.empty:
            return matched

        matched[
            "m5_bridge_type"
        ] = "LD_TO_MORADI_QTL"

        matched[
            "m5_bridge_variant_rsid"
        ] = matched[
            "qtl_rsid"
        ]

        matched[
            "m5_moradi_target_role"
        ] = "DIRECT_QTL_VARIANT"

        matched[
            "m5_moradi_source_tag_relationship_used"
        ] = False

        return matched

    # ------------------------------------------------------------------
    # Source tag bridge
    # ------------------------------------------------------------------

    @classmethod
    def _source_tag_matches(
        cls,
        candidate_ld: pd.DataFrame,
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """Match LD proxies against Moradi source tag variants."""

        evidence = moradi.loc[
            moradi[
                "tag_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            & moradi[
                "tag_rsid"
            ]
            .notna()
        ].copy()

        if evidence.empty:
            return pd.DataFrame()

        evidence = cls._moradi_projection(
            evidence
        )

        matched = candidate_ld.merge(
            evidence,
            how="inner",
            left_on="neighbor_rsid",
            right_on="tag_rsid",
            validate="many_to_many",
        )

        if matched.empty:
            return matched

        matched[
            "m5_bridge_type"
        ] = (
            "LD_TO_MORADI_SOURCE_TAG"
        )

        matched[
            "m5_bridge_variant_rsid"
        ] = matched[
            "tag_rsid"
        ]

        matched[
            "m5_moradi_target_role"
        ] = "SOURCE_REPORTED_TAG_VARIANT"

        matched[
            "m5_moradi_source_tag_relationship_used"
        ] = True

        return matched

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def build_table(
        cls,
        loci: pd.DataFrame,
        ld_evidence: pd.DataFrame,
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """Build regulatory bridges for one Moradi index part."""

        cls._require_columns(
            moradi,
            cls.REQUIRED_MORADI_COLUMNS,
        )

        candidate_ld = cls._candidate_ld(
            loci,
            ld_evidence,
        )

        if candidate_ld.empty:
            return pd.DataFrame()

        direct = cls._direct_qtl_matches(
            candidate_ld,
            moradi,
        )

        tag = cls._source_tag_matches(
            candidate_ld,
            moradi,
        )

        frames = [
            frame
            for frame
            in (
                direct,
                tag,
            )
            if not frame.empty
        ]

        if not frames:
            return pd.DataFrame()

        result = pd.concat(
            frames,
            ignore_index=True,
            sort=False,
        )

        result[
            "m5_bridge_supported_by_ld"
        ] = True

        result[
            "m5_primary_regulatory_ld_bridge"
        ] = result[
            "m5_primary_ld_support"
        ].astype(
            "boolean"
        )

        result[
            "m5_source_reported_ld_recomputed"
        ] = False

        result[
            "m5_colocalization_performed"
        ] = False

        result[
            "m5_causal_inference_performed"
        ] = False

        return result