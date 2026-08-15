"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/moradi_index.py

Description:
    Builds a unified regulatory index from harmonized Moradi QTL tables.

    The index preserves the distinction between:

        - direct QTL variants
        - GWAS/tag variants
        - cis and trans regulatory scope
        - intron, exon, and transcript-isoform features

    The unified index is designed for downstream M4 variant-level
    integration with prostate cancer tRF-QTL candidates.

    This stage does not perform:

        - tRF-QTL matching
        - LD expansion
        - genome liftover
        - colocalization
        - significance filtering
        - candidate ranking
        - deduplication

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


@dataclass(frozen=True)
class MoradiIndexSpec:
    """Description of one harmonized Moradi QTL source table."""

    table_id: str
    source_sheet: str
    regulatory_scope: str
    feature_class: str
    has_tag_variant: bool


MORADI_INDEX_SPECS: dict[str, MoradiIndexSpec] = {
    "s1_cis_intron": MoradiIndexSpec(
        table_id="s1_cis_intron",
        source_sheet="S1",
        regulatory_scope="cis",
        feature_class="INTRON",
        has_tag_variant=False,
    ),
    "s1_cis_exon": MoradiIndexSpec(
        table_id="s1_cis_exon",
        source_sheet="S1",
        regulatory_scope="cis",
        feature_class="EXON",
        has_tag_variant=False,
    ),
    "s1_cis_isoform": MoradiIndexSpec(
        table_id="s1_cis_isoform",
        source_sheet="S1",
        regulatory_scope="cis",
        feature_class="TRANSCRIPT_ISOFORM",
        has_tag_variant=False,
    ),
    "s2_cis_intron": MoradiIndexSpec(
        table_id="s2_cis_intron",
        source_sheet="S2",
        regulatory_scope="cis",
        feature_class="INTRON",
        has_tag_variant=True,
    ),
    "s2_cis_exon": MoradiIndexSpec(
        table_id="s2_cis_exon",
        source_sheet="S2",
        regulatory_scope="cis",
        feature_class="EXON",
        has_tag_variant=True,
    ),
    "s2_cis_isoform": MoradiIndexSpec(
        table_id="s2_cis_isoform",
        source_sheet="S2",
        regulatory_scope="cis",
        feature_class="TRANSCRIPT_ISOFORM",
        has_tag_variant=True,
    ),
    "s3_trans_intron": MoradiIndexSpec(
        table_id="s3_trans_intron",
        source_sheet="S3",
        regulatory_scope="trans",
        feature_class="INTRON",
        has_tag_variant=False,
    ),
    "s3_trans_exon": MoradiIndexSpec(
        table_id="s3_trans_exon",
        source_sheet="S3",
        regulatory_scope="trans",
        feature_class="EXON",
        has_tag_variant=False,
    ),
    "s3_trans_isoform": MoradiIndexSpec(
        table_id="s3_trans_isoform",
        source_sheet="S3",
        regulatory_scope="trans",
        feature_class="TRANSCRIPT_ISOFORM",
        has_tag_variant=False,
    ),
    "s4_trans_intron": MoradiIndexSpec(
        table_id="s4_trans_intron",
        source_sheet="S4",
        regulatory_scope="trans",
        feature_class="INTRON",
        has_tag_variant=True,
    ),
    "s4_trans_exon": MoradiIndexSpec(
        table_id="s4_trans_exon",
        source_sheet="S4",
        regulatory_scope="trans",
        feature_class="EXON",
        has_tag_variant=True,
    ),
    "s4_trans_isoform": MoradiIndexSpec(
        table_id="s4_trans_isoform",
        source_sheet="S4",
        regulatory_scope="trans",
        feature_class="TRANSCRIPT_ISOFORM",
        has_tag_variant=True,
    ),
}


class MoradiQTLIndexBuilder:
    """
    Convert harmonized Moradi QTL tables into one common index schema.

    Schema aliases are resolved conservatively to support the existing
    locked harmonization outputs without requiring them to be rewritten.
    """

    # ------------------------------------------------------------------
    # Column resolution
    # ------------------------------------------------------------------

    @staticmethod
    def _first_existing_column(
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> str | None:
        """Return the first matching source column."""

        for column in candidates:
            if column in dataframe.columns:
                return column

        return None

    @classmethod
    def _series_or_na(
        cls,
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
        *,
        dtype: str = "string",
    ) -> pd.Series:
        """Return a schema-compatible column or an NA series."""

        column = cls._first_existing_column(
            dataframe,
            candidates,
        )

        if column is None:
            return pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype=dtype,
            )

        if dtype == "string":
            return dataframe[
                column
            ].astype(
                "string"
            )

        return dataframe[
            column
        ].astype(
            dtype
        )

    @classmethod
    def _boolean_or_false(
        cls,
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> pd.Series:
        """Return a nullable boolean column or False."""

        column = cls._first_existing_column(
            dataframe,
            candidates,
        )

        if column is None:
            return pd.Series(
                False,
                index=dataframe.index,
                dtype="boolean",
            )

        return dataframe[
            column
        ].astype(
            "boolean"
        )

    @classmethod
    def _numeric_or_na(
        cls,
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> pd.Series:
        """Return optional numeric evidence."""

        column = cls._first_existing_column(
            dataframe,
            candidates,
        )

        if column is None:
            return pd.Series(
                float("nan"),
                index=dataframe.index,
                dtype="float64",
            )

        return pd.to_numeric(
            dataframe[
                column
            ],
            errors="coerce",
        )

    # ------------------------------------------------------------------
    # QTL variant
    # ------------------------------------------------------------------

    @classmethod
    def _add_qtl_variant(
        cls,
        output: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Add direct QTL-variant identity."""

        output[
            "qtl_rsid"
        ] = cls._series_or_na(
            source,
            (
                "harm_qtl_variant_rsid",
                "harm_variant_rsid",
                "qtl_variant_rsid",
            ),
        )

        output[
            "qtl_rsid_usable"
        ] = cls._boolean_or_false(
            source,
            (
                "harm_qtl_variant_rsid_usable",
                "harm_variant_rsid_usable",
                "qtl_variant_rsid_usable",
            ),
        )

        output[
            "qtl_variant_identity_key"
        ] = cls._series_or_na(
            source,
            (
                "harm_qtl_variant_identity_key",
                "harm_variant_identity_key",
                "qtl_variant_identity_key",
            ),
        )

        output[
            "qtl_variant_status"
        ] = cls._series_or_na(
            source,
            (
                "harm_qtl_variant_status",
                "harm_variant_status",
                "qtl_variant_status",
            ),
        )

        output[
            "qtl_coordinate_join_allowed"
        ] = cls._boolean_or_false(
            source,
            (
                "harm_qtl_coordinate_join_allowed",
                "harm_coordinate_join_allowed",
                "harm_variant_coordinate_join_allowed",
            ),
        )

    # ------------------------------------------------------------------
    # Tag variant
    # ------------------------------------------------------------------

    @classmethod
    def _add_tag_variant(
        cls,
        output: pd.DataFrame,
        source: pd.DataFrame,
        *,
        has_tag_variant: bool,
    ) -> None:
        """Add GWAS/tag-variant identity when supplied by source."""

        if not has_tag_variant:

            output[
                "tag_rsid"
            ] = pd.Series(
                pd.NA,
                index=source.index,
                dtype="string",
            )

            output[
                "tag_rsid_usable"
            ] = pd.Series(
                False,
                index=source.index,
                dtype="boolean",
            )

            output[
                "tag_variant_identity_key"
            ] = pd.Series(
                pd.NA,
                index=source.index,
                dtype="string",
            )

            output[
                "tag_variant_status"
            ] = pd.Series(
                pd.NA,
                index=source.index,
                dtype="string",
            )

            output[
                "tag_coordinate_join_allowed"
            ] = pd.Series(
                False,
                index=source.index,
                dtype="boolean",
            )

            return

        output[
            "tag_rsid"
        ] = cls._series_or_na(
            source,
            (
                "harm_tag_variant_rsid",
                "tag_variant_rsid",
                "harm_tag_rsid",
            ),
        )

        output[
            "tag_rsid_usable"
        ] = cls._boolean_or_false(
            source,
            (
                "harm_tag_variant_rsid_usable",
                "tag_variant_rsid_usable",
                "harm_tag_rsid_usable",
            ),
        )

        output[
            "tag_variant_identity_key"
        ] = cls._series_or_na(
            source,
            (
                "harm_tag_variant_identity_key",
                "tag_variant_identity_key",
            ),
        )

        output[
            "tag_variant_status"
        ] = cls._series_or_na(
            source,
            (
                "harm_tag_variant_status",
                "tag_variant_status",
            ),
        )

        output[
            "tag_coordinate_join_allowed"
        ] = cls._boolean_or_false(
            source,
            (
                "harm_tag_coordinate_join_allowed",
                "tag_coordinate_join_allowed",
            ),
        )

    # ------------------------------------------------------------------
    # Feature
    # ------------------------------------------------------------------

    @classmethod
    def _add_feature(
        cls,
        output: pd.DataFrame,
        source: pd.DataFrame,
        *,
        feature_class: str,
    ) -> None:
        """Add regulatory feature identity."""

        output[
            "feature_id"
        ] = cls._series_or_na(
            source,
            (
                "harm_feature_id",
                "feature_id",
            ),
        )

        output[
            "feature_identity_key"
        ] = cls._series_or_na(
            source,
            (
                "harm_feature_identity_key",
                "feature_identity_key",
            ),
        )

        output[
            "feature_type"
        ] = cls._series_or_na(
            source,
            (
                "harm_feature_type",
                "feature_type",
            ),
        )

        missing_type = (
            output[
                "feature_type"
            ]
            .isna()
        )

        output.loc[
            missing_type,
            "feature_type",
        ] = feature_class

        output[
            "feature_usable"
        ] = cls._boolean_or_false(
            source,
            (
                "harm_feature_usable",
                "feature_usable",
            ),
        )

        output[
            "feature_status"
        ] = cls._series_or_na(
            source,
            (
                "harm_feature_status",
                "feature_status",
            ),
        )

    # ------------------------------------------------------------------
    # Disease
    # ------------------------------------------------------------------

    @classmethod
    def _add_disease(
        cls,
        output: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Preserve harmonized disease context."""

        output[
            "disease_id"
        ] = cls._series_or_na(
            source,
            (
                "harm_disease_id",
                "disease_id",
            ),
        )

        output[
            "is_primary_disease"
        ] = cls._boolean_or_false(
            source,
            (
                "harm_is_primary_disease",
                "is_primary_disease",
            ),
        )

        output[
            "disease_status"
        ] = cls._series_or_na(
            source,
            (
                "harm_disease_status",
                "disease_status",
            ),
        )

    # ------------------------------------------------------------------
    # Statistical evidence
    # ------------------------------------------------------------------

    @classmethod
    def _add_statistics(
        cls,
        output: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Preserve available association statistics."""

        output[
            "p_value"
        ] = cls._numeric_or_na(
            source,
            (
                "p_value",
                "qtl_p_value",
                "source_p_value",
            ),
        )

        output[
            "fdr"
        ] = cls._numeric_or_na(
            source,
            (
                "fdr",
                "FDR",
                "source_fdr",
            ),
        )

        output[
            "beta"
        ] = cls._numeric_or_na(
            source,
            (
                "beta",
                "qtl_beta",
                "source_beta",
            ),
        )

        output[
            "ld_r2"
        ] = cls._numeric_or_na(
            source,
            (
                "ld_r2",
                "LD",
                "source_ld",
            ),
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def build_table(
        cls,
        dataframe: pd.DataFrame,
        *,
        spec: MoradiIndexSpec,
    ) -> pd.DataFrame:
        """
        Convert one harmonized Moradi table into common M4 index schema.
        """

        source = dataframe.reset_index(
            drop=True
        )

        output = pd.DataFrame(
            index=source.index
        )

        output[
            "m4_source_dataset"
        ] = "moradi"

        output[
            "m4_source_table"
        ] = spec.table_id

        output[
            "m4_source_sheet"
        ] = spec.source_sheet

        output[
            "m4_regulatory_scope"
        ] = spec.regulatory_scope

        output[
            "m4_feature_class"
        ] = spec.feature_class

        output[
            "m4_has_tag_variant"
        ] = spec.has_tag_variant

        output[
            "m4_source_row"
        ] = pd.Series(
            range(
                1,
                len(
                    source
                )
                + 1,
            ),
            dtype="int64",
        )

        cls._add_qtl_variant(
            output,
            source,
        )

        cls._add_tag_variant(
            output,
            source,
            has_tag_variant=(
                spec.has_tag_variant
            ),
        )

        cls._add_feature(
            output,
            source,
            feature_class=(
                spec.feature_class
            ),
        )

        cls._add_disease(
            output,
            source,
        )

        cls._add_statistics(
            output,
            source,
        )

        output[
            "m4_direct_variant_match_ready"
        ] = (
            output[
                "qtl_rsid_usable"
            ]
            & output[
                "feature_usable"
            ]
        ).astype(
            "boolean"
        )

        output[
            "m4_tag_variant_match_ready"
        ] = (
            output[
                "tag_rsid_usable"
            ]
            & output[
                "feature_usable"
            ]
        ).astype(
            "boolean"
        )

        output[
            "m4_ld_inference_performed"
        ] = False

        output[
            "m4_colocalization_performed"
        ] = False

        if (
            len(
                output
            )
            != len(
                source
            )
        ):
            raise RuntimeError(
                f"Moradi index cardinality changed for "
                f"{spec.table_id}: "
                f"{len(source)} -> {len(output)}"
            )

        return output