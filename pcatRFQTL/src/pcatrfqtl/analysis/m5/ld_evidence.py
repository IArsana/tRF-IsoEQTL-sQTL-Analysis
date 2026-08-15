"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/ld_evidence.py

Description:
    Provider-neutral normalization of externally calculated linkage
    disequilibrium evidence for M5.

    This module does not calculate LD.

    Expected evidence represents relationships between a tRF-QTL lead
    variant and neighboring variants obtained from an external LD
    reference panel or locally calculated reference dataset.

    Required biological distinction:

        independently calculated/reference LD
            !=
        Moradi source-reported tag-SNP LD

    Only independently supplied M5 LD evidence is normalized here.

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
from typing import Any

import pandas as pd


class LDEvidenceClass(StrEnum):
    """Operational LD evidence categories."""

    HIGH_LD = "HIGH_LD"
    MODERATE_LD = "MODERATE_LD"
    BELOW_THRESHOLD = "BELOW_THRESHOLD"
    INVALID = "INVALID"


class M5LDEvidenceNormalizer:
    """
    Normalize provider-specific LD tables into one M5 schema.

    Input tables may use aliases such as:

        lead_rsid / query_rsid / lead_snp
        neighbor_rsid / proxy_rsid / variant_rsid
        r2 / R2 / r_squared
        d_prime / Dprime
        population / pop
        reference_panel / panel

    At minimum, lead rsID, neighbor rsID, and r² must be resolvable.
    """

    RSID_PATTERN = r"^rs\d+$"

    # ------------------------------------------------------------------
    # Column resolution
    # ------------------------------------------------------------------

    @staticmethod
    def _first_existing_column(
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> str | None:
        """Return first existing compatible column."""

        for column in candidates:
            if column in dataframe.columns:
                return column

        return None

    @classmethod
    def _required_column(
        cls,
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
        *,
        label: str,
    ) -> str:
        """Resolve a required LD evidence column."""

        column = cls._first_existing_column(
            dataframe,
            candidates,
        )

        if column is None:
            raise ValueError(
                f"LD evidence is missing required {label} column. "
                f"Tried aliases: {candidates}"
            )

        return column

    @classmethod
    def _optional_string(
        cls,
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> pd.Series:
        """Return optional string field."""

        column = cls._first_existing_column(
            dataframe,
            candidates,
        )

        if column is None:
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

    @classmethod
    def _optional_numeric(
        cls,
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> pd.Series:
        """Return optional numeric field."""

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
    # Normalization
    # ------------------------------------------------------------------

    @classmethod
    def normalize(
        cls,
        dataframe: pd.DataFrame,
        *,
        primary_r2_threshold: float,
        secondary_r2_threshold: float,
        source_file: str | None = None,
    ) -> pd.DataFrame:
        """Normalize one external/precomputed LD evidence table."""

        if not (
            0.0
            <= secondary_r2_threshold
            <= primary_r2_threshold
            <= 1.0
        ):
            raise ValueError(
                "LD thresholds must satisfy "
                "0 <= secondary <= primary <= 1."
            )

        source = dataframe.reset_index(
            drop=True
        )

        lead_column = cls._required_column(
            source,
            (
                "lead_rsid",
                "query_rsid",
                "index_rsid",
                "lead_snp",
                "query_snp",
            ),
            label="lead rsID",
        )

        neighbor_column = cls._required_column(
            source,
            (
                "neighbor_rsid",
                "proxy_rsid",
                "variant_rsid",
                "target_rsid",
                "linked_rsid",
                "proxy_snp",
            ),
            label="neighbor rsID",
        )

        r2_column = cls._required_column(
            source,
            (
                "r2",
                "R2",
                "r_squared",
                "ld_r2",
            ),
            label="r²",
        )

        result = pd.DataFrame(
            index=source.index
        )

        result[
            "lead_rsid"
        ] = (
            source[
                lead_column
            ]
            .astype(
                "string"
            )
            .str.strip()
            .str.lower()
        )

        result[
            "neighbor_rsid"
        ] = (
            source[
                neighbor_column
            ]
            .astype(
                "string"
            )
            .str.strip()
            .str.lower()
        )

        result[
            "r2"
        ] = pd.to_numeric(
            source[
                r2_column
            ],
            errors="coerce",
        )

        result[
            "d_prime"
        ] = cls._optional_numeric(
            source,
            (
                "d_prime",
                "Dprime",
                "D_prime",
                "dprime",
            ),
        )

        result[
            "ld_population"
        ] = cls._optional_string(
            source,
            (
                "population",
                "pop",
                "ld_population",
                "reference_population",
            ),
        )

        result[
            "ld_reference_panel"
        ] = cls._optional_string(
            source,
            (
                "reference_panel",
                "panel",
                "ld_reference_panel",
                "reference_dataset",
            ),
        )

        result[
            "ld_reference_build"
        ] = cls._optional_string(
            source,
            (
                "reference_build",
                "genome_build",
                "build",
                "assembly",
            ),
        )

        result[
            "ld_provider"
        ] = cls._optional_string(
            source,
            (
                "provider",
                "ld_provider",
                "source",
            ),
        )

        result[
            "ld_source_file"
        ] = (
            source_file
            if source_file is not None
            else pd.NA
        )

        result[
            "lead_rsid_valid"
        ] = (
            result[
                "lead_rsid"
            ]
            .str.match(
                cls.RSID_PATTERN,
                na=False,
            )
            .astype(
                "boolean"
            )
        )

        result[
            "neighbor_rsid_valid"
        ] = (
            result[
                "neighbor_rsid"
            ]
            .str.match(
                cls.RSID_PATTERN,
                na=False,
            )
            .astype(
                "boolean"
            )
        )

        result[
            "r2_valid"
        ] = (
            result[
                "r2"
            ]
            .between(
                0.0,
                1.0,
                inclusive="both",
            )
            .fillna(False)
            .astype(
                "boolean"
            )
        )

        result[
            "is_self_pair"
        ] = (
            result[
                "lead_rsid"
            ]
            == result[
                "neighbor_rsid"
            ]
        ).astype(
            "boolean"
        )

        result[
            "ld_pair_usable"
        ] = (
            result[
                "lead_rsid_valid"
            ]
            & result[
                "neighbor_rsid_valid"
            ]
            & result[
                "r2_valid"
            ]
            & ~result[
                "is_self_pair"
            ]
        ).astype(
            "boolean"
        )

        result[
            "ld_evidence_class"
        ] = (
            LDEvidenceClass.INVALID.value
        )

        valid = result[
            "ld_pair_usable"
        ].fillna(
            False
        )

        result.loc[
            valid,
            "ld_evidence_class",
        ] = (
            LDEvidenceClass
            .BELOW_THRESHOLD
            .value
        )

        moderate = (
            valid
            & (
                result[
                    "r2"
                ]
                >= secondary_r2_threshold
            )
        )

        result.loc[
            moderate,
            "ld_evidence_class",
        ] = (
            LDEvidenceClass
            .MODERATE_LD
            .value
        )

        high = (
            valid
            & (
                result[
                    "r2"
                ]
                >= primary_r2_threshold
            )
        )

        result.loc[
            high,
            "ld_evidence_class",
        ] = (
            LDEvidenceClass
            .HIGH_LD
            .value
        )

        result[
            "m5_bridge_eligible"
        ] = moderate.astype(
            "boolean"
        )

        result[
            "m5_primary_ld_support"
        ] = high.astype(
            "boolean"
        )

        result[
            "m5_secondary_ld_support"
        ] = moderate.astype(
            "boolean"
        )

        result[
            "m5_ld_calculated_by_pipeline"
        ] = False

        result[
            "m5_ld_imported_evidence"
        ] = True

        result[
            "m5_colocalization_performed"
        ] = False

        return result

    # ------------------------------------------------------------------
    # Multi-file normalization
    # ------------------------------------------------------------------

    @classmethod
    def combine(
        cls,
        tables: list[
            tuple[
                pd.DataFrame,
                str,
            ]
        ],
        *,
        primary_r2_threshold: float,
        secondary_r2_threshold: float,
    ) -> pd.DataFrame:
        """Normalize and combine multiple external LD tables."""

        frames: list[
            pd.DataFrame
        ] = []

        for (
            dataframe,
            source_file,
        ) in tables:

            normalized = cls.normalize(
                dataframe,
                primary_r2_threshold=(
                    primary_r2_threshold
                ),
                secondary_r2_threshold=(
                    secondary_r2_threshold
                ),
                source_file=source_file,
            )

            frames.append(
                normalized
            )

        if not frames:
            return pd.DataFrame()

        result = pd.concat(
            frames,
            ignore_index=True,
            sort=False,
        )

        return result