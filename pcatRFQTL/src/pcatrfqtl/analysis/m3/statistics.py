"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/statistics.py

Description:
    Statistical and descriptive summary logic for M3.5.

    M3.5 summarizes the first-stage prostate-cancer GWAS × tRF-QTL
    integration without introducing LD, colocalization, causal
    inference, or biological candidate ranking.

    The stage quantifies:

        - GWAS variant representation
        - PRAD tRF-QTL representation
        - direct canonical-rsID overlap
        - candidate-level direct evidence status
        - direct overlap proportions

    Because M3.3 may yield zero exact variant overlap, M3.5 treats
    absence of direct overlap as a valid scientific result rather than
    assigning unsupported candidate rankings.

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


class M3StatisticsAnalyzer:
    """Calculate descriptive statistics for M3 integration."""

    REQUIRED_GWAS_COLUMNS = {
        "harm_variant_rsid",
        "m3_direct_overlap_eligible",
    }

    REQUIRED_TRFQTL_COLUMNS = {
        "harm_variant_rsid",
        "harm_feature_id",
        "harm_feature_identity_key",
        "m3_direct_overlap_eligible",
    }

    REQUIRED_CANDIDATE_COLUMNS = {
        "m3_candidate_id",
        "harm_variant_rsid",
        "harm_feature_id",
        "harm_feature_identity_key",
        "m3_direct_gwas_match",
        "m3_gwas_association_count",
        "m3_candidate_direct_evidence",
    }

    @staticmethod
    def _require_columns(
        dataframe: pd.DataFrame,
        required: set[str],
        dataset_name: str,
    ) -> None:
        """Require columns needed for M3.5."""

        missing = (
            required
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                f"{dataset_name} is missing required M3.5 columns: "
                f"{sorted(missing)}"
            )

    @staticmethod
    def _eligible_rows(
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return direct-overlap eligible records."""

        return (
            dataframe.loc[
                dataframe[
                    "m3_direct_overlap_eligible"
                ]
                .fillna(False)
                .eq(True)
            ]
            .copy()
        )

    @staticmethod
    def _percentage(
        numerator: int,
        denominator: int,
    ) -> float:
        """Return percentage with zero-denominator protection."""

        if denominator == 0:
            return 0.0

        return float(
            numerator
            / denominator
            * 100.0
        )

    @classmethod
    def summarize(
        cls,
        *,
        gwas: pd.DataFrame,
        trfqtl: pd.DataFrame,
        candidates: pd.DataFrame,
    ) -> dict[str, Any]:
        """Calculate M3.5 descriptive statistics."""

        cls._require_columns(
            gwas,
            cls.REQUIRED_GWAS_COLUMNS,
            "M3.1 GWAS",
        )

        cls._require_columns(
            trfqtl,
            cls.REQUIRED_TRFQTL_COLUMNS,
            "M3.2 tRF-QTL",
        )

        cls._require_columns(
            candidates,
            cls.REQUIRED_CANDIDATE_COLUMNS,
            "M3.4 candidate table",
        )

        gwas_eligible = (
            cls._eligible_rows(
                gwas
            )
        )

        trfqtl_eligible = (
            cls._eligible_rows(
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

        shared_rsids = (
            gwas_rsids
            & trfqtl_rsids
        )

        direct_candidate_rows = int(
            candidates[
                "m3_direct_gwas_match"
            ]
            .fillna(False)
            .eq(True)
            .sum()
        )

        candidate_rows = len(
            candidates
        )

        unique_candidate_rsids = int(
            candidates[
                "harm_variant_rsid"
            ]
            .dropna()
            .nunique()
        )

        unique_candidate_trfs = int(
            candidates[
                "harm_feature_id"
            ]
            .dropna()
            .nunique()
        )

        unique_variant_trf_pairs = int(
            candidates[
                [
                    "harm_variant_rsid",
                    "harm_feature_identity_key",
                ]
            ]
            .dropna()
            .drop_duplicates()
            .shape[
                0
            ]
        )

        statistics: dict[
            str,
            Any,
        ] = {
            "gwas": {
                "association_rows":
                    int(
                        len(
                            gwas
                        )
                    ),

                "eligible_association_rows":
                    int(
                        len(
                            gwas_eligible
                        )
                    ),

                "unique_eligible_rsids":
                    int(
                        len(
                            gwas_rsids
                        )
                    ),
            },

            "trfqtl": {
                "association_rows":
                    int(
                        len(
                            trfqtl
                        )
                    ),

                "eligible_association_rows":
                    int(
                        len(
                            trfqtl_eligible
                        )
                    ),

                "unique_eligible_rsids":
                    int(
                        len(
                            trfqtl_rsids
                        )
                    ),

                "unique_trfs":
                    int(
                        trfqtl[
                            "harm_feature_id"
                        ]
                        .dropna()
                        .nunique()
                    ),
            },

            "candidate_table": {
                "rows":
                    int(
                        candidate_rows
                    ),

                "unique_rsids":
                    unique_candidate_rsids,

                "unique_trfs":
                    unique_candidate_trfs,

                "unique_variant_trf_pairs":
                    unique_variant_trf_pairs,

                "direct_match_rows":
                    direct_candidate_rows,

                "no_direct_match_rows":
                    int(
                        candidate_rows
                        - direct_candidate_rows
                    ),
            },

            "direct_overlap": {
                "shared_rsids":
                    sorted(
                        shared_rsids
                    ),

                "unique_shared_rsids":
                    int(
                        len(
                            shared_rsids
                        )
                    ),

                "trfqtl_variant_overlap_rate_percent":
                    cls._percentage(
                        len(
                            shared_rsids
                        ),
                        len(
                            trfqtl_rsids
                        ),
                    ),

                "gwas_variant_overlap_rate_percent":
                    cls._percentage(
                        len(
                            shared_rsids
                        ),
                        len(
                            gwas_rsids
                        ),
                    ),

                "candidate_direct_match_rate_percent":
                    cls._percentage(
                        direct_candidate_rows,
                        candidate_rows,
                    ),
            },
        }

        return statistics

    @classmethod
    def build_candidate_statistics_table(
        cls,
        candidates: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Build concise M3.5 candidate statistics table.

        This table does not rank candidates.
        """

        cls._require_columns(
            candidates,
            cls.REQUIRED_CANDIDATE_COLUMNS,
            "M3.4 candidate table",
        )

        preferred_columns = [
            "m3_candidate_id",
            "harm_variant_rsid",
            "harm_variant_identity_key",
            "harm_feature_id",
            "harm_feature_identity_key",
            "m3_direct_overlap_eligible",
            "m3_direct_gwas_match",
            "m3_gwas_association_count",
            "m3_candidate_direct_evidence",
        ]

        available_columns = [
            column
            for column in preferred_columns
            if column in candidates.columns
        ]

        result = (
            candidates[
                available_columns
            ]
            .copy()
        )

        result[
            "m3_evidence_level"
        ] = "TRFQTL_ONLY"

        result.loc[
            result[
                "m3_direct_gwas_match"
            ]
            .fillna(False)
            .eq(True),
            "m3_evidence_level",
        ] = "DIRECT_GWAS_TRFQTL"

        result[
            "m3_rank"
        ] = pd.Series(
            [
                pd.NA
                for _ in range(
                    len(
                        result
                    )
                )
            ],
            dtype="Int64",
        )

        result[
            "m3_ranking_performed"
        ] = False

        return result