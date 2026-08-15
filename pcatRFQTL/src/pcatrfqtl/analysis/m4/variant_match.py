"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/variant_match.py

Description:
    Exact variant-level integration between prostate cancer tRF-QTL
    associations and the M4.1 Unified Moradi QTL Index.

    Two evidence modes are preserved independently:

        QTL_RSID_EXACT
            The tRF-QTL variant rsID exactly matches a Moradi QTL
            variant rsID.

        SOURCE_TAG_RSID_EXACT
            The tRF-QTL variant rsID exactly matches a tag variant
            reported by the Moradi source.

    SOURCE_TAG_RSID_EXACT must not be interpreted as newly calculated
    linkage disequilibrium. No LD inference is performed in this stage.

    This stage does not perform:

        - coordinate matching
        - genome liftover
        - LD calculation or expansion
        - colocalization
        - statistical filtering
        - candidate ranking
        - causal inference

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


class M4VariantMatchMode(StrEnum):
    """Exact variant match modes used by M4.2."""

    QTL_RSID_EXACT = "QTL_RSID_EXACT"

    SOURCE_TAG_RSID_EXACT = (
        "SOURCE_TAG_RSID_EXACT"
    )


class TRFQTLMoradiVariantMatcher:
    """
    Match PRAD tRF-QTL variants against one Moradi index part.

    Matching is performed only using canonical rsID equality.
    """

    REQUIRED_TRFQTL_COLUMNS = frozenset(
        {
            "harm_variant_rsid",
            "harm_variant_rsid_usable",
            "harm_feature_id",
            "harm_feature_identity_key",
            "harm_feature_usable",
        }
    )

    REQUIRED_MORADI_COLUMNS = frozenset(
        {
            "m4_source_table",
            "m4_source_sheet",
            "m4_regulatory_scope",
            "m4_feature_class",
            "m4_source_row",
            "qtl_rsid",
            "qtl_rsid_usable",
            "tag_rsid",
            "tag_rsid_usable",
            "feature_id",
            "feature_identity_key",
            "feature_type",
        }
    )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    @classmethod
    def _require_columns(
        cls,
        dataframe: pd.DataFrame,
        required: frozenset[str],
        *,
        label: str,
    ) -> None:
        """Require a dataframe schema."""

        missing = (
            required
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                f"{label} is missing required columns: "
                f"{sorted(missing)}"
            )

    # ------------------------------------------------------------------
    # tRF-QTL preparation
    # ------------------------------------------------------------------

    @classmethod
    def prepare_trfqtl(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Prepare exact-match-eligible PRAD tRF-QTL rows.

        No deduplication is performed.
        """

        cls._require_columns(
            dataframe,
            cls.REQUIRED_TRFQTL_COLUMNS,
            label="M4.2 tRF-QTL input",
        )

        result = dataframe.copy()

        eligible = (
            result[
                "harm_variant_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            &
            result[
                "harm_feature_usable"
            ]
            .fillna(False)
            .eq(True)
            &
            result[
                "harm_variant_rsid"
            ]
            .notna()
        )

        result = (
            result.loc[
                eligible
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        result[
            "m4_trfqtl_input_row"
        ] = pd.Series(
            range(
                1,
                len(
                    result
                )
                + 1,
            ),
            dtype="int64",
        )

        return result

    # ------------------------------------------------------------------
    # Output projection
    # ------------------------------------------------------------------

    @staticmethod
    def _candidate_projection(
        trfqtl: pd.DataFrame,
    ) -> pd.DataFrame:
        """Create compact tRF-QTL candidate representation."""

        columns = {
            "m4_trfqtl_input_row":
                "m4_trfqtl_input_row",

            "harm_variant_rsid":
                "trfqtl_rsid",

            "harm_variant_identity_key":
                "trfqtl_variant_identity_key",

            "harm_feature_id":
                "trf_id",

            "harm_feature_identity_key":
                "trf_identity_key",

            "harm_feature_type":
                "trf_feature_type",

            "harm_disease_id":
                "trfqtl_disease_id",

            "harm_is_primary_disease":
                "trfqtl_is_primary_disease",
        }

        available = {
            source:
                target
            for (
                source,
                target,
            ) in columns.items()
            if source
            in trfqtl.columns
        }

        return (
            trfqtl[
                list(
                    available
                )
            ]
            .rename(
                columns=available
            )
            .copy()
        )

    @staticmethod
    def _moradi_projection(
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """Create compact regulatory-evidence representation."""

        preferred_columns = (
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
            "disease_status",
            "p_value",
            "fdr",
            "beta",
            "ld_r2",
        )

        columns = [
            column
            for column
            in preferred_columns
            if column
            in moradi.columns
        ]

        return moradi[
            columns
        ].copy()

    # ------------------------------------------------------------------
    # Direct QTL match
    # ------------------------------------------------------------------

    @classmethod
    def _match_qtl_variant(
        cls,
        trfqtl: pd.DataFrame,
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """Match tRF-QTL rsIDs to Moradi direct QTL rsIDs."""

        eligible = (
            moradi[
                "qtl_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            &
            moradi[
                "qtl_rsid"
            ]
            .notna()
        )

        evidence = (
            moradi.loc[
                eligible
            ]
            .copy()
        )

        if evidence.empty:
            return pd.DataFrame()

        candidates = (
            cls._candidate_projection(
                trfqtl
            )
        )

        evidence = (
            cls._moradi_projection(
                evidence
            )
        )

        matched = candidates.merge(
            evidence,
            how="inner",
            left_on="trfqtl_rsid",
            right_on="qtl_rsid",
            validate="many_to_many",
        )

        if matched.empty:
            return matched

        matched[
            "m4_match_mode"
        ] = (
            M4VariantMatchMode
            .QTL_RSID_EXACT
            .value
        )

        matched[
            "m4_matched_moradi_rsid"
        ] = matched[
            "qtl_rsid"
        ]

        matched[
            "m4_source_reported_tag_evidence"
        ] = False

        matched[
            "m4_ld_inference_performed"
        ] = False

        matched[
            "m4_colocalization_performed"
        ] = False

        return matched

    # ------------------------------------------------------------------
    # Source tag match
    # ------------------------------------------------------------------

    @classmethod
    def _match_tag_variant(
        cls,
        trfqtl: pd.DataFrame,
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Match tRF-QTL rsIDs to source-reported Moradi tag rsIDs.

        This represents source-reported tag evidence only. It is not
        newly inferred LD evidence.
        """

        eligible = (
            moradi[
                "tag_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            &
            moradi[
                "tag_rsid"
            ]
            .notna()
        )

        evidence = (
            moradi.loc[
                eligible
            ]
            .copy()
        )

        if evidence.empty:
            return pd.DataFrame()

        candidates = (
            cls._candidate_projection(
                trfqtl
            )
        )

        evidence = (
            cls._moradi_projection(
                evidence
            )
        )

        matched = candidates.merge(
            evidence,
            how="inner",
            left_on="trfqtl_rsid",
            right_on="tag_rsid",
            validate="many_to_many",
        )

        if matched.empty:
            return matched

        matched[
            "m4_match_mode"
        ] = (
            M4VariantMatchMode
            .SOURCE_TAG_RSID_EXACT
            .value
        )

        matched[
            "m4_matched_moradi_rsid"
        ] = matched[
            "tag_rsid"
        ]

        matched[
            "m4_source_reported_tag_evidence"
        ] = True

        matched[
            "m4_ld_inference_performed"
        ] = False

        matched[
            "m4_colocalization_performed"
        ] = False

        return matched

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def match_table(
        cls,
        trfqtl: pd.DataFrame,
        moradi: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Match eligible tRF-QTL variants against one Moradi index part.

        The returned representation is long-form evidence. A single
        Moradi row may therefore generate separate direct-QTL and
        source-tag evidence records when applicable.
        """

        cls._require_columns(
            moradi,
            cls.REQUIRED_MORADI_COLUMNS,
            label="M4.1 Moradi QTL index part",
        )

        prepared_trfqtl = (
            cls.prepare_trfqtl(
                trfqtl
            )
        )

        direct_matches = (
            cls._match_qtl_variant(
                prepared_trfqtl,
                moradi,
            )
        )

        tag_matches = (
            cls._match_tag_variant(
                prepared_trfqtl,
                moradi,
            )
        )

        frames = [
            dataframe
            for dataframe
            in (
                direct_matches,
                tag_matches,
            )
            if not dataframe.empty
        ]

        if not frames:
            return pd.DataFrame()

        result = pd.concat(
            frames,
            ignore_index=True,
            sort=False,
        )

        result[
            "m4_exact_rsid_match"
        ] = True

        result[
            "m4_coordinate_match_performed"
        ] = False

        result[
            "m4_candidate_ranked"
        ] = False

        return result