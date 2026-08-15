"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/overlap.py

Description:
    Core scientific logic for M3.3 direct GWAS × tRF-QTL overlap.

    M3.3 identifies direct shared variants between:

        - prostate-cancer GWAS associations from M3.1
        - prostate-cancer tRF-QTL associations from M3.2

    Matching is performed exclusively by exact canonical rsID identity.

    The stage intentionally does not:

        - use LD proxy variants
        - perform coordinate-based matching
        - perform genome liftover
        - collapse GWAS studies
        - rank biological candidates
        - perform colocalization
        - infer causal relationships

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


class DirectOverlapAnalyzer:
    """Identify direct canonical-rsID overlap between GWAS and tRF-QTL."""

    GWAS_REQUIRED_COLUMNS = {
        "harm_variant_rsid",
        "m3_direct_overlap_eligible",
    }

    TRFQTL_REQUIRED_COLUMNS = {
        "harm_variant_rsid",
        "harm_feature_id",
        "harm_feature_identity_key",
        "m3_direct_overlap_eligible",
    }

    @classmethod
    def validate_gwas_schema(
        cls,
        dataframe: pd.DataFrame,
    ) -> None:
        """Require M3.1 fields."""

        missing = (
            cls.GWAS_REQUIRED_COLUMNS
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                "M3.1 GWAS input is missing required M3.3 columns: "
                f"{sorted(missing)}"
            )

    @classmethod
    def validate_trfqtl_schema(
        cls,
        dataframe: pd.DataFrame,
    ) -> None:
        """Require M3.2 fields."""

        missing = (
            cls.TRFQTL_REQUIRED_COLUMNS
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                "M3.2 tRF-QTL input is missing required M3.3 columns: "
                f"{sorted(missing)}"
            )

    @classmethod
    def eligible_gwas(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return GWAS rows eligible for exact rsID matching."""

        cls.validate_gwas_schema(
            dataframe
        )

        mask = (
            dataframe[
                "m3_direct_overlap_eligible"
            ]
            .fillna(False)
            .eq(True)
            & dataframe[
                "harm_variant_rsid"
            ]
            .notna()
        )

        return (
            dataframe.loc[
                mask
            ]
            .copy()
        )

    @classmethod
    def eligible_trfqtl(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return tRF-QTL rows eligible for exact rsID matching."""

        cls.validate_trfqtl_schema(
            dataframe
        )

        mask = (
            dataframe[
                "m3_direct_overlap_eligible"
            ]
            .fillna(False)
            .eq(True)
            & dataframe[
                "harm_variant_rsid"
            ]
            .notna()
            & dataframe[
                "harm_feature_identity_key"
            ]
            .notna()
        )

        return (
            dataframe.loc[
                mask
            ]
            .copy()
        )

    @classmethod
    def shared_rsids(
        cls,
        gwas: pd.DataFrame,
        trfqtl: pd.DataFrame,
    ) -> set[str]:
        """Return canonical rsIDs shared by eligible GWAS and tRF-QTL."""

        gwas_eligible = (
            cls.eligible_gwas(
                gwas
            )
        )

        trfqtl_eligible = (
            cls.eligible_trfqtl(
                trfqtl
            )
        )

        gwas_rsids = set(
            gwas_eligible[
                "harm_variant_rsid"
            ]
            .dropna()
            .astype(str)
        )

        trfqtl_rsids = set(
            trfqtl_eligible[
                "harm_variant_rsid"
            ]
            .dropna()
            .astype(str)
        )

        return (
            gwas_rsids
            & trfqtl_rsids
        )

    @classmethod
    def build_overlap_table(
        cls,
        gwas: pd.DataFrame,
        trfqtl: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Build row-level direct GWAS × tRF-QTL overlap table.

        GWAS association multiplicity is intentionally preserved.

        Therefore one shared rsID may produce multiple output rows when
        the GWAS Catalog contains several prostate-cancer associations
        for the same canonical variant.
        """

        gwas_eligible = (
            cls.eligible_gwas(
                gwas
            )
        )

        trfqtl_eligible = (
            cls.eligible_trfqtl(
                trfqtl
            )
        )

        shared = (
            cls.shared_rsids(
                gwas_eligible,
                trfqtl_eligible,
            )
        )

        if not shared:
            return pd.DataFrame()

        gwas_shared = (
            gwas_eligible[
                gwas_eligible[
                    "harm_variant_rsid"
                ]
                .isin(
                    shared
                )
            ]
            .copy()
        )

        trfqtl_shared = (
            trfqtl_eligible[
                trfqtl_eligible[
                    "harm_variant_rsid"
                ]
                .isin(
                    shared
                )
            ]
            .copy()
        )

        overlap = (
            gwas_shared.merge(
                trfqtl_shared,
                on="harm_variant_rsid",
                how="inner",
                suffixes=(
                    "_gwas",
                    "_trfqtl",
                ),
                validate="many_to_many",
            )
        )

        overlap[
            "m3_overlap_method"
        ] = "canonical_rsid_exact"

        overlap[
            "m3_direct_variant_overlap"
        ] = True

        return overlap