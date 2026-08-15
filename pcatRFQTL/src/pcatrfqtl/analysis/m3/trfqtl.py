"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/trfqtl.py

Description:
    Core scientific logic for M3.2 prostate-cancer tRF-QTL selection.

    M3.2 selects harmonized Cancer-tRFQTL S2 records whose disease
    context corresponds to prostate cancer.

    All selected prostate-cancer tRF-QTL association rows are retained.

    Records are annotated for eligibility in the later M3.3 direct
    GWAS × tRF-QTL overlap analysis.

    Direct overlap eligibility requires:

        - prostate-cancer disease context
        - usable canonical rsID
        - usable tRF feature identity

    M3.2 intentionally does not:

        - join tRF-QTL records with GWAS
        - deduplicate variants
        - deduplicate tRFs
        - collapse multiple SNP-tRF associations
        - perform coordinate matching
        - perform genome liftover
        - rank candidate tRFs
        - apply additional significance thresholds

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


class ProstateTRFQTLSelector:
    """Select prostate-cancer Cancer-tRFQTL S2 association records."""

    PRIMARY_DISEASE_ID = "prostate_cancer"

    REQUIRED_COLUMNS = {
        "harm_is_primary_disease",
        "harm_disease_id",
        "harm_variant_rsid",
        "harm_variant_rsid_usable",
        "harm_variant_identity_key",
        "harm_feature_id",
        "harm_feature_identity_key",
        "harm_feature_usable",
        "harm_feature_type",
    }

    @classmethod
    def validate_schema(
        cls,
        dataframe: pd.DataFrame,
    ) -> None:
        """Require harmonized Cancer-tRFQTL fields needed by M3.2."""

        missing = (
            cls.REQUIRED_COLUMNS
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                "Harmonized Cancer-tRFQTL S2 input is missing "
                f"required M3.2 columns: {sorted(missing)}"
            )

    @classmethod
    def select_primary_disease(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Select prostate-cancer tRF-QTL rows."""

        cls.validate_schema(
            dataframe
        )

        mask = (
            dataframe[
                "harm_is_primary_disease"
            ]
            .fillna(False)
            .eq(True)
            & dataframe[
                "harm_disease_id"
            ]
            .eq(
                cls.PRIMARY_DISEASE_ID
            )
        )

        return (
            dataframe.loc[
                mask
            ]
            .copy()
        )

    @staticmethod
    def add_overlap_eligibility(
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Annotate eligibility for M3.3 direct rsID overlap.

        A row is eligible when both its canonical rsID and tRF feature
        identity are usable.
        """

        result = (
            dataframe.copy()
        )

        result[
            "m3_direct_overlap_eligible"
        ] = (
            result[
                "harm_variant_rsid_usable"
            ]
            .fillna(False)
            .eq(True)
            & result[
                "harm_variant_rsid"
            ]
            .notna()
            & result[
                "harm_feature_usable"
            ]
            .fillna(False)
            .eq(True)
            & result[
                "harm_feature_identity_key"
            ]
            .notna()
        )

        result[
            "m3_selection_reason"
        ] = (
            "harmonized_primary_disease"
        )

        return result

    @classmethod
    def prepare(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Produce the M3.2 prostate-cancer tRF-QTL subset."""

        selected = (
            cls.select_primary_disease(
                dataframe
            )
        )

        selected = (
            cls.add_overlap_eligibility(
                selected
            )
        )

        return selected