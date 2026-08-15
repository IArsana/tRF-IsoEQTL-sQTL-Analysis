"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/evidence_summary.py

Description:
    Builds the M4.3 candidate-level multi-omic exact-evidence summary.

    The summary preserves one row per prostate cancer SNP-tRF pair and
    consolidates exact evidence obtained from previous milestones:

        - M3 exact prostate cancer GWAS support
        - M4.2 exact Moradi direct-QTL support
        - M4.2 exact Moradi source-reported tag support
        - cis regulatory support
        - trans regulatory support
        - intron, exon, and transcript-isoform evidence

    This stage summarizes evidence only.

    It does not perform:

        - genome-coordinate matching
        - genome liftover
        - LD calculation
        - LD expansion
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

from typing import Any

import pandas as pd


class M43ExactEvidenceSummaryBuilder:
    """
    Build one-row-per-SNP-tRF exact evidence summaries.

    The builder assumes the prostate tRF-QTL table has already been
    restricted to the primary disease during M3.
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

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    @staticmethod
    def _require_columns(
        dataframe: pd.DataFrame,
        required: frozenset[str],
        *,
        label: str,
    ) -> None:
        """Require expected dataframe columns."""

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
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_string(
        dataframe: pd.DataFrame,
        column: str,
    ) -> pd.Series:
        """Return nullable string column or missing values."""

        if column not in dataframe.columns:
            return pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype="string",
            )

        return dataframe[
            column
        ].astype(
            "string"
        )

    @staticmethod
    def _safe_boolean(
        dataframe: pd.DataFrame,
        column: str,
        *,
        default: bool = False,
    ) -> pd.Series:
        """Return nullable boolean column or a constant default."""

        if column not in dataframe.columns:
            return pd.Series(
                default,
                index=dataframe.index,
                dtype="boolean",
            )

        return dataframe[
            column
        ].fillna(
            default
        ).astype(
            "boolean"
        )

    @staticmethod
    def _empty_boolean(
        index: pd.Index,
    ) -> pd.Series:
        """Return a False boolean series."""

        return pd.Series(
            False,
            index=index,
            dtype="boolean",
        )

    @staticmethod
    def _empty_integer(
        index: pd.Index,
    ) -> pd.Series:
        """Return a zero integer series."""

        return pd.Series(
            0,
            index=index,
            dtype="Int64",
        )

    # ------------------------------------------------------------------
    # Candidate preparation
    # ------------------------------------------------------------------

    @classmethod
    def prepare_candidates(
        cls,
        trfqtl: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Build stable one-row-per-SNP-tRF candidate scaffold.

        Association rows are not biologically deduplicated. Only exact
        duplicate SNP-feature identities are collapsed for the
        candidate-level summary representation.
        """

        cls._require_columns(
            trfqtl,
            cls.REQUIRED_TRFQTL_COLUMNS,
            label="M4.3 prostate tRF-QTL input",
        )

        eligible = (
            trfqtl[
                "harm_variant_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            &
            trfqtl[
                "harm_feature_usable"
            ]
            .fillna(False)
            .eq(True)
            &
            trfqtl[
                "harm_variant_rsid"
            ]
            .notna()
            &
            trfqtl[
                "harm_feature_identity_key"
            ]
            .notna()
        )

        source = (
            trfqtl.loc[
                eligible
            ]
            .copy()
        )

        source[
            "m4_candidate_key"
        ] = (
            source[
                "harm_variant_rsid"
            ]
            .astype(str)
            + "|"
            + source[
                "harm_feature_identity_key"
            ]
            .astype(str)
        )

        source = (
            source.sort_values(
                by=[
                    "harm_variant_rsid",
                    "harm_feature_identity_key",
                ],
                kind="stable",
            )
            .drop_duplicates(
                subset=[
                    "m4_candidate_key",
                ],
                keep="first",
            )
            .reset_index(
                drop=True
            )
        )

        output = pd.DataFrame(
            index=source.index
        )

        output[
            "m4_candidate_key"
        ] = source[
            "m4_candidate_key"
        ].astype(
            "string"
        )

        output[
            "trfqtl_rsid"
        ] = source[
            "harm_variant_rsid"
        ].astype(
            "string"
        )

        output[
            "trfqtl_variant_identity_key"
        ] = cls._safe_string(
            source,
            "harm_variant_identity_key",
        )

        output[
            "trf_id"
        ] = source[
            "harm_feature_id"
        ].astype(
            "string"
        )

        output[
            "trf_identity_key"
        ] = source[
            "harm_feature_identity_key"
        ].astype(
            "string"
        )

        output[
            "trf_feature_type"
        ] = cls._safe_string(
            source,
            "harm_feature_type",
        )

        output[
            "disease_id"
        ] = cls._safe_string(
            source,
            "harm_disease_id",
        )

        output[
            "is_primary_disease"
        ] = cls._safe_boolean(
            source,
            "harm_is_primary_disease",
            default=True,
        )

        output[
            "m4_candidate_id"
        ] = [
            f"M4C{index:03d}"
            for index
            in range(
                1,
                len(
                    output
                )
                + 1,
            )
        ]

        return output

    # ------------------------------------------------------------------
    # M3 support
    # ------------------------------------------------------------------

    @classmethod
    def add_m3_exact_gwas_support(
        cls,
        candidates: pd.DataFrame,
        *,
        m3_candidate_scaffold: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """
        Add direct GWAS evidence.

        When an M3 candidate scaffold is available, compatible direct
        GWAS fields are consumed.

        Otherwise the locked M3 result is represented as no exact GWAS
        support for the current three PRAD tRF-QTL candidates.
        """

        result = candidates.copy()

        result[
            "direct_gwas_match"
        ] = False

        result[
            "direct_gwas_evidence_rows"
        ] = 0

        result[
            "direct_gwas_support_source"
        ] = "M3_LOCKED_EXACT_OVERLAP_RESULT"

        if (
            m3_candidate_scaffold is None
            or m3_candidate_scaffold.empty
        ):
            return result

        rsid_column = None

        for candidate in (
            "trfqtl_rsid",
            "harm_variant_rsid",
            "variant_rsid",
        ):
            if candidate in m3_candidate_scaffold.columns:
                rsid_column = candidate
                break

        support_column = None

        for candidate in (
            "direct_gwas_match",
            "has_direct_gwas_match",
            "exact_gwas_support",
        ):
            if candidate in m3_candidate_scaffold.columns:
                support_column = candidate
                break

        if (
            rsid_column is None
            or support_column is None
        ):
            return result

        scaffold = (
            m3_candidate_scaffold[
                [
                    rsid_column,
                    support_column,
                ]
            ]
            .copy()
        )

        scaffold[
            support_column
        ] = (
            scaffold[
                support_column
            ]
            .fillna(False)
            .astype(bool)
        )

        summary = (
            scaffold.groupby(
                rsid_column,
                dropna=True,
            )[
                support_column
            ]
            .agg(
                [
                    "any",
                    "sum",
                ]
            )
            .reset_index()
            .rename(
                columns={
                    "any":
                        "_direct_gwas_match",
                    "sum":
                        "_direct_gwas_rows",
                }
            )
        )

        result = result.merge(
            summary,
            how="left",
            left_on="trfqtl_rsid",
            right_on=rsid_column,
        )

        result[
            "direct_gwas_match"
        ] = (
            result[
                "_direct_gwas_match"
            ]
            .fillna(False)
            .astype(
                "boolean"
            )
        )

        result[
            "direct_gwas_evidence_rows"
        ] = (
            result[
                "_direct_gwas_rows"
            ]
            .fillna(0)
            .astype(
                "Int64"
            )
        )

        result[
            "direct_gwas_support_source"
        ] = "M3_CANDIDATE_SCAFFOLD"

        drop_columns = [
            column
            for column
            in (
                rsid_column,
                "_direct_gwas_match",
                "_direct_gwas_rows",
            )
            if column
            in result.columns
            and column
            != "trfqtl_rsid"
        ]

        return result.drop(
            columns=drop_columns
        )

    # ------------------------------------------------------------------
    # Moradi evidence aggregation
    # ------------------------------------------------------------------

    @staticmethod
    def _aggregate_moradi_matches(
        matches: pd.DataFrame,
    ) -> pd.DataFrame:
        """Aggregate long-form M4.2 evidence to variant level."""

        if matches.empty:
            return pd.DataFrame()

        required = {
            "trfqtl_rsid",
            "m4_match_mode",
        }

        missing = (
            required
            - set(
                matches.columns
            )
        )

        if missing:
            raise ValueError(
                "M4.2 match table is missing required columns: "
                f"{sorted(missing)}"
            )

        source = matches.copy()

        source[
            "_direct_qtl"
        ] = (
            source[
                "m4_match_mode"
            ]
            == "QTL_RSID_EXACT"
        )

        source[
            "_source_tag"
        ] = (
            source[
                "m4_match_mode"
            ]
            == "SOURCE_TAG_RSID_EXACT"
        )

        if (
            "m4_regulatory_scope"
            in source.columns
        ):
            source[
                "_cis"
            ] = (
                source[
                    "m4_regulatory_scope"
                ]
                == "cis"
            )

            source[
                "_trans"
            ] = (
                source[
                    "m4_regulatory_scope"
                ]
                == "trans"
            )
        else:
            source[
                "_cis"
            ] = False

            source[
                "_trans"
            ] = False

        if (
            "m4_feature_class"
            in source.columns
        ):
            source[
                "_intron"
            ] = (
                source[
                    "m4_feature_class"
                ]
                == "INTRON"
            )

            source[
                "_exon"
            ] = (
                source[
                    "m4_feature_class"
                ]
                == "EXON"
            )

            source[
                "_isoform"
            ] = (
                source[
                    "m4_feature_class"
                ]
                == "TRANSCRIPT_ISOFORM"
            )
        else:
            source[
                "_intron"
            ] = False

            source[
                "_exon"
            ] = False

            source[
                "_isoform"
            ] = False

        source[
            "_feature_key"
        ] = (
            source[
                "feature_identity_key"
            ]
            if "feature_identity_key"
            in source.columns
            else pd.NA
        )

        grouped = (
            source.groupby(
                "trfqtl_rsid",
                dropna=True,
            )
            .agg(
                moradi_exact_evidence_rows=(
                    "m4_match_mode",
                    "size",
                ),
                moradi_direct_qtl_match=(
                    "_direct_qtl",
                    "any",
                ),
                moradi_source_tag_match=(
                    "_source_tag",
                    "any",
                ),
                cis_exact_support=(
                    "_cis",
                    "any",
                ),
                trans_exact_support=(
                    "_trans",
                    "any",
                ),
                intron_exact_support=(
                    "_intron",
                    "any",
                ),
                exon_exact_support=(
                    "_exon",
                    "any",
                ),
                isoform_exact_support=(
                    "_isoform",
                    "any",
                ),
                unique_regulatory_features=(
                    "_feature_key",
                    "nunique",
                ),
            )
            .reset_index()
        )

        return grouped

    @classmethod
    def add_moradi_exact_support(
        cls,
        candidates: pd.DataFrame,
        matches: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add exact Moradi regulatory evidence."""

        result = candidates.copy()

        defaults: dict[str, Any] = {
            "moradi_exact_evidence_rows":
                cls._empty_integer(
                    result.index
                ),

            "moradi_direct_qtl_match":
                cls._empty_boolean(
                    result.index
                ),

            "moradi_source_tag_match":
                cls._empty_boolean(
                    result.index
                ),

            "cis_exact_support":
                cls._empty_boolean(
                    result.index
                ),

            "trans_exact_support":
                cls._empty_boolean(
                    result.index
                ),

            "intron_exact_support":
                cls._empty_boolean(
                    result.index
                ),

            "exon_exact_support":
                cls._empty_boolean(
                    result.index
                ),

            "isoform_exact_support":
                cls._empty_boolean(
                    result.index
                ),

            "unique_regulatory_features":
                cls._empty_integer(
                    result.index
                ),
        }

        if matches.empty:
            for (
                column,
                values,
            ) in defaults.items():
                result[
                    column
                ] = values

            return result

        aggregated = (
            cls._aggregate_moradi_matches(
                matches
            )
        )

        result = result.merge(
            aggregated,
            how="left",
            on="trfqtl_rsid",
        )

        boolean_columns = (
            "moradi_direct_qtl_match",
            "moradi_source_tag_match",
            "cis_exact_support",
            "trans_exact_support",
            "intron_exact_support",
            "exon_exact_support",
            "isoform_exact_support",
        )

        integer_columns = (
            "moradi_exact_evidence_rows",
            "unique_regulatory_features",
        )

        for column in boolean_columns:
            result[
                column
            ] = (
                result[
                    column
                ]
                .fillna(False)
                .astype(
                    "boolean"
                )
            )

        for column in integer_columns:
            result[
                column
            ] = (
                result[
                    column
                ]
                .fillna(0)
                .astype(
                    "Int64"
                )
            )

        return result

    # ------------------------------------------------------------------
    # Exact multi-omic state
    # ------------------------------------------------------------------

    @staticmethod
    def finalize(
        candidates: pd.DataFrame,
    ) -> pd.DataFrame:
        """Calculate exact evidence state and M5 transition flags."""

        result = candidates.copy()

        result[
            "moradi_any_exact_support"
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
            "regulatory_any_exact_support"
        ] = (
            result[
                "cis_exact_support"
            ]
            | result[
                "trans_exact_support"
            ]
        ).astype(
            "boolean"
        )

        result[
            "exact_multiomic_support"
        ] = (
            result[
                "direct_gwas_match"
            ]
            & result[
                "moradi_any_exact_support"
            ]
        ).astype(
            "boolean"
        )

        result[
            "exact_evidence_layer_count"
        ] = (
            result[
                "direct_gwas_match"
            ].astype(
                "Int64"
            )
            + result[
                "moradi_any_exact_support"
            ].astype(
                "Int64"
            )
        )

        result[
            "m5_ld_required"
        ] = (
            ~result[
                "exact_multiomic_support"
            ]
        ).astype(
            "boolean"
        )

        result[
            "m5_colocalization_candidate"
        ] = result[
            "m5_ld_required"
        ].astype(
            "boolean"
        )

        result[
            "m4_exact_evidence_status"
        ] = "NO_EXACT_MULTIOMIC_BRIDGE"

        result.loc[
            result[
                "moradi_any_exact_support"
            ],
            "m4_exact_evidence_status",
        ] = "REGULATORY_EXACT_SUPPORT_ONLY"

        result.loc[
            result[
                "direct_gwas_match"
            ],
            "m4_exact_evidence_status",
        ] = "GWAS_EXACT_SUPPORT_ONLY"

        result.loc[
            result[
                "exact_multiomic_support"
            ],
            "m4_exact_evidence_status",
        ] = "EXACT_MULTIOMIC_BRIDGE"

        result[
            "m4_ld_inference_performed"
        ] = False

        result[
            "m4_colocalization_performed"
        ] = False

        result[
            "m4_candidate_ranked"
        ] = False

        return result

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def build(
        cls,
        trfqtl: pd.DataFrame,
        moradi_matches: pd.DataFrame,
        *,
        m3_candidate_scaffold: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Build complete M4.3 exact-evidence summary."""

        candidates = cls.prepare_candidates(
            trfqtl
        )

        candidates = (
            cls.add_m3_exact_gwas_support(
                candidates,
                m3_candidate_scaffold=(
                    m3_candidate_scaffold
                ),
            )
        )

        candidates = (
            cls.add_moradi_exact_support(
                candidates,
                moradi_matches,
            )
        )

        return cls.finalize(
            candidates
        )