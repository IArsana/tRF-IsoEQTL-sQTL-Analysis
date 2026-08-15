"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/candidate_loci.py

Description:
    Builds M5.2 candidate locus definitions around prostate cancer
    tRF-QTL lead variants.

    Candidate loci are physical screening regions only.

    They must not be interpreted as:

        - LD blocks
        - credible sets
        - causal loci
        - colocalized regions

    Actual LD relationships are evaluated in later M5 stages.

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

from pcatrfqtl.analysis.m5.ld_strategy import (
    LDReferenceStrategy,
)


class M52CandidateLocusBuilder:
    """Build one physical locus per M4 SNP-tRF candidate."""

    REQUIRED_BRIDGE_COLUMNS = frozenset(
        {
            "m4_candidate_id",
            "trfqtl_rsid",
            "trf_id",
            "regulatory_bridge_state",
            "m5_requires_ld_lookup",
        }
    )

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _require_columns(
        dataframe: pd.DataFrame,
        required: frozenset[str],
        *,
        label: str,
    ) -> None:
        """Require dataframe columns."""

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

    @staticmethod
    def _first_existing_column(
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> str | None:
        """Return first available schema-compatible column."""

        for column in candidates:
            if column in dataframe.columns:
                return column

        return None

    @staticmethod
    def _normalize_chromosome(
        value: Any,
    ) -> str | None:
        """Normalize chromosome to chrN notation."""

        if value is None:
            return None

        try:
            if pd.isna(
                value
            ):
                return None
        except (
            TypeError,
            ValueError,
        ):
            pass

        chromosome = str(
            value
        ).strip()

        if not chromosome:
            return None

        if chromosome.lower().startswith(
            "chr"
        ):
            chromosome = chromosome[
                3:
            ]

        chromosome = chromosome.upper()

        if chromosome == "M":
            chromosome = "MT"

        valid = {
            str(index)
            for index
            in range(
                1,
                23,
            )
        } | {
            "X",
            "Y",
            "MT",
        }

        if chromosome not in valid:
            return None

        return f"chr{chromosome}"

    @staticmethod
    def _normalize_position(
        value: Any,
    ) -> int | None:
        """Normalize genomic position."""

        if value is None:
            return None

        try:
            if pd.isna(
                value
            ):
                return None
        except (
            TypeError,
            ValueError,
        ):
            pass

        try:
            position = int(
                value
            )
        except (
            TypeError,
            ValueError,
        ):
            return None

        if position <= 0:
            return None

        return position

    # ------------------------------------------------------------------
    # Coordinate extraction
    # ------------------------------------------------------------------

    @classmethod
    def prepare_trfqtl_coordinates(
        cls,
        trfqtl: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Extract canonical lead variant coordinates from harmonized
        Cancer-tRFQTL data.
        """

        rsid_column = (
            cls._first_existing_column(
                trfqtl,
                (
                    "harm_variant_rsid",
                    "variant_rsid",
                    "snp_rsid",
                ),
            )
        )

        chromosome_column = (
            cls._first_existing_column(
                trfqtl,
                (
                    "harm_variant_chromosome",
                    "harm_chromosome",
                    "chromosome",
                ),
            )
        )

        position_column = (
            cls._first_existing_column(
                trfqtl,
                (
                    "harm_variant_position",
                    "harm_position",
                    "position",
                ),
            )
        )

        if rsid_column is None:
            raise ValueError(
                "M5.2 could not resolve tRF-QTL rsID column."
            )

        if chromosome_column is None:
            raise ValueError(
                "M5.2 could not resolve tRF-QTL chromosome column."
            )

        if position_column is None:
            raise ValueError(
                "M5.2 could not resolve tRF-QTL position column."
            )

        source = trfqtl[
            [
                rsid_column,
                chromosome_column,
                position_column,
            ]
        ].copy()

        source = source.rename(
            columns={
                rsid_column:
                    "lead_rsid",

                chromosome_column:
                    "_source_chromosome",

                position_column:
                    "_source_position",
            }
        )

        source[
            "lead_chromosome"
        ] = [
            cls._normalize_chromosome(
                value
            )
            for value
            in source[
                "_source_chromosome"
            ]
        ]

        source[
            "lead_position"
        ] = [
            cls._normalize_position(
                value
            )
            for value
            in source[
                "_source_position"
            ]
        ]

        source[
            "lead_coordinate_usable"
        ] = (
            source[
                "lead_chromosome"
            ].notna()
            &
            source[
                "lead_position"
            ].notna()
        )

        source = (
            source[
                [
                    "lead_rsid",
                    "lead_chromosome",
                    "lead_position",
                    "lead_coordinate_usable",
                ]
            ]
            .drop_duplicates(
                subset=[
                    "lead_rsid",
                    "lead_chromosome",
                    "lead_position",
                ]
            )
            .reset_index(
                drop=True
            )
        )

        conflicting = (
            source.loc[
                source[
                    "lead_coordinate_usable"
                ],
                "lead_rsid",
            ]
            .value_counts()
        )

        conflicting = conflicting[
            conflicting
            > 1
        ]

        if not conflicting.empty:
            raise RuntimeError(
                "M5.2 detected conflicting coordinates for rsIDs: "
                f"{sorted(conflicting.index.tolist())}"
            )

        return source

    # ------------------------------------------------------------------
    # Locus construction
    # ------------------------------------------------------------------

    @classmethod
    def build(
        cls,
        bridge_readiness: pd.DataFrame,
        trfqtl: pd.DataFrame,
        *,
        strategy: LDReferenceStrategy,
    ) -> pd.DataFrame:
        """Build candidate screening loci."""

        strategy.validate()

        cls._require_columns(
            bridge_readiness,
            cls.REQUIRED_BRIDGE_COLUMNS,
            label="M5.2 bridge readiness",
        )

        candidates = bridge_readiness.copy()

        coordinates = (
            cls.prepare_trfqtl_coordinates(
                trfqtl
            )
        )

        candidates = candidates.merge(
            coordinates,
            how="left",
            left_on="trfqtl_rsid",
            right_on="lead_rsid",
            validate="many_to_one",
        )

        if (
            "lead_rsid"
            in candidates.columns
        ):
            candidates = candidates.drop(
                columns=[
                    "lead_rsid",
                ]
            )

        candidates[
            "m5_locus_id"
        ] = [
            f"M5L{index:03d}"
            for index
            in range(
                1,
                len(
                    candidates
                )
                + 1,
            )
        ]

        candidates[
            "m5_reference_assembly"
        ] = strategy.reference_assembly

        candidates[
            "m5_source_build"
        ] = strategy.source_build

        candidates[
            "m5_screening_window_bp"
        ] = strategy.screening_window_bp

        candidates[
            "locus_start"
        ] = pd.Series(
            pd.NA,
            index=candidates.index,
            dtype="Int64",
        )

        candidates[
            "locus_end"
        ] = pd.Series(
            pd.NA,
            index=candidates.index,
            dtype="Int64",
        )

        usable = (
            candidates[
                "lead_coordinate_usable"
            ]
            .fillna(False)
            .eq(True)
        )

        candidates.loc[
            usable,
            "locus_start",
        ] = (
            candidates.loc[
                usable,
                "lead_position",
            ]
            .astype(
                "int64"
            )
            .sub(
                strategy.screening_window_bp
            )
            .clip(
                lower=1
            )
            .astype(
                "Int64"
            )
        )

        candidates.loc[
            usable,
            "locus_end",
        ] = (
            candidates.loc[
                usable,
                "lead_position",
            ]
            .astype(
                "int64"
            )
            .add(
                strategy.screening_window_bp
            )
            .astype(
                "Int64"
            )
        )

        candidates[
            "m5_locus_definition"
        ] = "PHYSICAL_SCREENING_WINDOW"

        candidates[
            "m5_locus_is_ld_block"
        ] = False

        candidates[
            "m5_locus_is_credible_set"
        ] = False

        candidates[
            "m5_ld_calculated"
        ] = False

        candidates[
            "m5_ld_reference_panel_available"
        ] = (
            strategy.reference_panel_name
            is not None
        )

        candidates[
            "m5_primary_r2_threshold"
        ] = strategy.primary_r2_threshold

        candidates[
            "m5_secondary_r2_threshold"
        ] = strategy.secondary_r2_threshold

        candidates[
            "m5_ld_query_ready"
        ] = (
            usable
            &
            candidates[
                "m5_requires_ld_lookup"
            ]
            .fillna(False)
            .eq(True)
        ).astype(
            "boolean"
        )

        candidates[
            "m5_locus_status"
        ] = "READY_FOR_LD_QUERY"

        candidates.loc[
            ~usable,
            "m5_locus_status",
        ] = "LEAD_COORDINATE_UNRESOLVED"

        candidates.loc[
            (
                usable
                &
                ~candidates[
                    "m5_requires_ld_lookup"
                ]
                .fillna(False)
                .eq(True)
            ),
            "m5_locus_status",
        ] = "LD_QUERY_NOT_REQUIRED"

        candidates[
            "m5_colocalization_performed"
        ] = False

        candidates[
            "m5_candidate_ranked"
        ] = False

        candidates[
            "m5_causal_inference_performed"
        ] = False

        return candidates