"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/trfqtl_supplementary.py

Description:
    Harmonization utilities for supplementary Cancer-tRFQTL tables.

    Supported backfill sheets:

        S1  - survival-associated tRFQTLs
        S3  - cross-ancestry GWAS metadata
        S4  - tRF immune-infiltration associations
        S5  - tRF drug-response associations
        S6  - Chinese cohort metadata
        S7  - experimental sequence metadata
        S8  - UK Biobank cohort metadata
        S9  - PLCO cohort metadata
        S11 - POU2F1-regulated genes

    S2 and S10 are intentionally excluded because their harmonized
    outputs already belong to the locked core pipeline.

    Harmonization is annotation-first and preserves row cardinality.
    No disease filtering, genome liftover, statistical filtering,
    deduplication, cross-source joining, or candidate ranking occurs.

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
import re
from typing import Any

import pandas as pd


DATASET_ID = "cancer_trfqtl"
SOURCE_BUILD = "hg19"
CANONICAL_BUILD = "GRCh37"


# ---------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class SupplementaryHarmonizationSpec:
    """Harmonization policy for one supplementary sheet."""

    sheet: str
    harmonize_disease: bool = False
    harmonize_variant: bool = False
    harmonize_feature: bool = False
    feature_type: str | None = None
    metadata_only: bool = False


SHEET_REGISTRY: dict[str, SupplementaryHarmonizationSpec] = {
    "S1": SupplementaryHarmonizationSpec(
        sheet="S1",
        harmonize_disease=True,
        harmonize_variant=True,
        harmonize_feature=True,
        feature_type="TRF",
    ),
    "S3": SupplementaryHarmonizationSpec(
        sheet="S3",
        harmonize_disease=True,
        metadata_only=True,
    ),
    "S4": SupplementaryHarmonizationSpec(
        sheet="S4",
        harmonize_disease=True,
        harmonize_feature=True,
        feature_type="TRF",
    ),
    "S5": SupplementaryHarmonizationSpec(
        sheet="S5",
        harmonize_disease=True,
        harmonize_feature=True,
        feature_type="TRF",
    ),
    "S6": SupplementaryHarmonizationSpec(
        sheet="S6",
        metadata_only=True,
    ),
    "S7": SupplementaryHarmonizationSpec(
        sheet="S7",
        metadata_only=True,
    ),
    "S8": SupplementaryHarmonizationSpec(
        sheet="S8",
        metadata_only=True,
    ),
    "S9": SupplementaryHarmonizationSpec(
        sheet="S9",
        metadata_only=True,
    ),
    "S11": SupplementaryHarmonizationSpec(
        sheet="S11",
        harmonize_feature=True,
        feature_type="GENE",
    ),
}


# ---------------------------------------------------------------------
# Harmonizer
# ---------------------------------------------------------------------


class CancerTRFQTLSupplementaryHarmonizer:
    """
    Harmonize non-core Cancer-tRFQTL supplementary tables.

    Harmonization is conservative:

        - source identifiers are preserved
        - PRAD/prostate cancer labels are resolved to prostate_cancer
        - other disease labels remain source-preserved
        - tRF identifiers remain source-namespace identifiers
        - gene symbols and Ensembl gene IDs are retained
        - valid rsIDs are canonicalized
        - hg19 coordinates are parsed without liftover
    """

    SUPPORTED_SHEETS = frozenset(
        SHEET_REGISTRY
    )

    LOCKED_CORE_SHEETS = frozenset(
        {
            "S2",
            "S10",
        }
    )

    RSID_PATTERN = re.compile(
        r"^rs\d+$",
        flags=re.IGNORECASE,
    )

    COORDINATE_PATTERN = re.compile(
        r"^(?:chr)?(?P<chromosome>[0-9]{1,2}|X|Y|MT|M):(?P<position>\d+)$",
        flags=re.IGNORECASE,
    )

    ENSEMBL_GENE_PATTERN = re.compile(
        r"^ENSG\d+(?:\.\d+)?$",
        flags=re.IGNORECASE,
    )

    PROSTATE_ALIASES = frozenset(
        {
            "prad",
            "prostate cancer",
            "prostate_cancer",
            "prostate carcinoma",
            "prostate_carcinoma",
            "prostate adenocarcinoma",
            "prostate_adenocarcinoma",
        }
    )

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------

    @classmethod
    def require_supported_sheet(
        cls,
        sheet: str,
    ) -> None:
        """Require a sheet handled by supplementary harmonization."""

        if sheet in cls.LOCKED_CORE_SHEETS:
            raise ValueError(
                f"{sheet} belongs to the locked Cancer-tRFQTL "
                "core harmonization pipeline."
            )

        if sheet not in cls.SUPPORTED_SHEETS:
            raise ValueError(
                f"Unsupported Cancer-tRFQTL supplementary sheet: {sheet}"
            )

    @staticmethod
    def _clean_text(
        value: Any,
    ) -> str | None:
        """Convert a non-missing scalar to stripped text."""

        if value is None:
            return None

        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            pass

        text = str(
            value
        ).strip()

        if not text:
            return None

        return text

    @staticmethod
    def _series(
        values: list[Any],
        *,
        dtype: str = "object",
    ) -> pd.Series:
        """Create a reset-index series."""

        return pd.Series(
            values,
            dtype=dtype,
        )

    # ------------------------------------------------------------------
    # Disease harmonization
    # ------------------------------------------------------------------

    @classmethod
    def _harmonize_disease_value(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Harmonize disease context.

        Only explicitly recognized prostate-cancer aliases are mapped to
        the internal prostate_cancer identity. Other disease labels are
        preserved without ontology inference.
        """

        raw = cls._clean_text(
            value
        )

        if raw is None:
            return {
                "harm_disease_id":
                    pd.NA,
                "harm_disease_label":
                    pd.NA,
                "harm_disease_status":
                    "UNRESOLVED",
                "harm_disease_method":
                    "NONE",
                "harm_is_primary_disease":
                    False,
            }

        normalized = (
            raw.lower()
            .strip()
            .replace("-", " ")
        )

        normalized = re.sub(
            r"\s+",
            " ",
            normalized,
        )

        if normalized in cls.PROSTATE_ALIASES:
            return {
                "harm_disease_id":
                    "prostate_cancer",
                "harm_disease_label":
                    "Prostate cancer",
                "harm_disease_status":
                    "RESOLVED",
                "harm_disease_method":
                    "EXACT_ALIAS",
                "harm_is_primary_disease":
                    True,
            }

        return {
            "harm_disease_id":
                raw,
            "harm_disease_label":
                raw,
            "harm_disease_status":
                "PRESERVED",
            "harm_disease_method":
                "SOURCE_LABEL",
            "harm_is_primary_disease":
                False,
        }

    @classmethod
    def _add_disease_harmonization(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add harmonized disease columns."""

        result = dataframe.copy()

        if "disease_raw" not in result.columns:
            raise ValueError(
                "Disease harmonization requires 'disease_raw'."
            )

        harmonized = [
            cls._harmonize_disease_value(
                value
            )
            for value in result[
                "disease_raw"
            ]
        ]

        result[
            "harm_disease_id"
        ] = cls._series(
            [
                item[
                    "harm_disease_id"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_disease_label"
        ] = cls._series(
            [
                item[
                    "harm_disease_label"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_disease_status"
        ] = cls._series(
            [
                item[
                    "harm_disease_status"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_disease_method"
        ] = cls._series(
            [
                item[
                    "harm_disease_method"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_is_primary_disease"
        ] = cls._series(
            [
                bool(
                    item[
                        "harm_is_primary_disease"
                    ]
                )
                for item in harmonized
            ],
            dtype="boolean",
        )

        return result

    # ------------------------------------------------------------------
    # Variant harmonization
    # ------------------------------------------------------------------

    @classmethod
    def _canonical_rsid(
        cls,
        value: Any,
    ) -> str | None:
        """Return a canonical lowercase rsID when valid."""

        text = cls._clean_text(
            value
        )

        if text is None:
            return None

        if not cls.RSID_PATTERN.fullmatch(
            text
        ):
            return None

        return text.lower()

    @classmethod
    def _parse_coordinate(
        cls,
        value: Any,
    ) -> tuple[str | None, int | None]:
        """Parse chromosome:position without performing liftover."""

        text = cls._clean_text(
            value
        )

        if text is None:
            return (
                None,
                None,
            )

        match = cls.COORDINATE_PATTERN.fullmatch(
            text
        )

        if match is None:
            return (
                None,
                None,
            )

        chromosome = (
            match.group(
                "chromosome"
            )
            .upper()
        )

        if chromosome == "M":
            chromosome = "MT"

        position = int(
            match.group(
                "position"
            )
        )

        if position <= 0:
            return (
                None,
                None,
            )

        return (
            chromosome,
            position,
        )

    @classmethod
    def _add_variant_harmonization(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Add hg19/GRCh37 variant identity annotations."""

        result = dataframe.copy()

        if "variant_raw" not in result.columns:
            raise ValueError(
                "Variant harmonization requires 'variant_raw'."
            )

        canonical_rsids: list[Any] = []
        rsid_usable: list[bool] = []

        chromosomes: list[Any] = []
        positions: list[Any] = []
        coordinate_usable: list[bool] = []

        identity_keys: list[Any] = []
        statuses: list[str] = []
        methods: list[str] = []

        coordinate_values = (
            result[
                "coordinate_raw"
            ]
            if "coordinate_raw"
            in result.columns
            else pd.Series(
                [
                    pd.NA
                    for _ in range(
                        len(
                            result
                        )
                    )
                ]
            )
        )

        for (
            variant_raw,
            coordinate_raw,
        ) in zip(
            result[
                "variant_raw"
            ],
            coordinate_values,
            strict=True,
        ):

            rsid = cls._canonical_rsid(
                variant_raw
            )

            (
                chromosome,
                position,
            ) = cls._parse_coordinate(
                coordinate_raw
            )

            has_rsid = (
                rsid is not None
            )

            has_coordinate = (
                chromosome is not None
                and position is not None
            )

            canonical_rsids.append(
                rsid
                if rsid is not None
                else pd.NA
            )

            rsid_usable.append(
                has_rsid
            )

            chromosomes.append(
                chromosome
                if chromosome is not None
                else pd.NA
            )

            positions.append(
                position
                if position is not None
                else pd.NA
            )

            coordinate_usable.append(
                has_coordinate
            )

            if has_rsid:
                identity_keys.append(
                    f"rsid:{rsid}"
                )

                statuses.append(
                    "RESOLVED"
                )

                methods.append(
                    "RSID"
                )

            elif has_coordinate:
                identity_keys.append(
                    (
                        f"coord:{SOURCE_BUILD}:"
                        f"{chromosome}:{position}"
                    )
                )

                statuses.append(
                    "PARTIAL"
                )

                methods.append(
                    "COORDINATE"
                )

            else:
                identity_keys.append(
                    pd.NA
                )

                statuses.append(
                    "UNRESOLVED"
                )

                methods.append(
                    "NONE"
                )

        result[
            "harm_variant_rsid"
        ] = cls._series(
            canonical_rsids,
            dtype="string",
        )

        result[
            "harm_variant_rsid_usable"
        ] = cls._series(
            rsid_usable,
            dtype="boolean",
        )

        result[
            "harm_variant_chromosome"
        ] = cls._series(
            chromosomes,
            dtype="string",
        )

        result[
            "harm_variant_position"
        ] = cls._series(
            positions,
            dtype="Int64",
        )

        result[
            "harm_variant_source_build"
        ] = SOURCE_BUILD

        result[
            "harm_variant_reference_assembly"
        ] = CANONICAL_BUILD

        result[
            "harm_variant_coordinate_usable"
        ] = cls._series(
            coordinate_usable,
            dtype="boolean",
        )

        result[
            "harm_variant_coordinate_join_allowed"
        ] = cls._series(
            coordinate_usable,
            dtype="boolean",
        )

        result[
            "harm_variant_identity_key"
        ] = cls._series(
            identity_keys,
            dtype="string",
        )

        result[
            "harm_variant_status"
        ] = cls._series(
            statuses,
            dtype="string",
        )

        result[
            "harm_variant_method"
        ] = cls._series(
            methods,
            dtype="string",
        )

        return result

    # ------------------------------------------------------------------
    # Feature harmonization
    # ------------------------------------------------------------------

    @classmethod
    def _harmonize_trf_value(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """Preserve a valid source tRF identifier."""

        raw = cls._clean_text(
            value
        )

        if raw is None:
            return {
                "id":
                    pd.NA,
                "key":
                    pd.NA,
                "status":
                    "UNRESOLVED",
                "usable":
                    False,
                "namespace":
                    pd.NA,
            }

        return {
            "id":
                raw,
            "key":
                f"trf:source:{raw}",
            "status":
                "PRESERVED",
            "usable":
                True,
            "namespace":
                "Cancer-tRFQTL",
        }

    @classmethod
    def _harmonize_gene_value(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """Preserve gene symbol or Ensembl gene identity."""

        raw = cls._clean_text(
            value
        )

        if raw is None:
            return {
                "id":
                    pd.NA,
                "key":
                    pd.NA,
                "status":
                    "UNRESOLVED",
                "usable":
                    False,
                "namespace":
                    pd.NA,
            }

        if cls.ENSEMBL_GENE_PATTERN.fullmatch(
            raw
        ):
            base = raw.split(
                ".",
                maxsplit=1,
            )[0]

            return {
                "id":
                    raw,
                "key":
                    f"gene:ensembl:{base}",
                "status":
                    "RESOLVED",
                "usable":
                    True,
                "namespace":
                    "ENSEMBL",
            }

        return {
            "id":
                raw,
            "key":
                f"gene:symbol:{raw.upper()}",
            "status":
                "PRESERVED",
            "usable":
                True,
            "namespace":
                "SYMBOL",
        }

    @classmethod
    def _add_feature_harmonization(
        cls,
        dataframe: pd.DataFrame,
        *,
        feature_type: str,
    ) -> pd.DataFrame:
        """Add harmonized tRF or gene identity."""

        result = dataframe.copy()

        if "feature_raw" not in result.columns:
            raise ValueError(
                "Feature harmonization requires 'feature_raw'."
            )

        if feature_type == "TRF":
            harmonizer = (
                cls._harmonize_trf_value
            )

        elif feature_type == "GENE":
            harmonizer = (
                cls._harmonize_gene_value
            )

        else:
            raise ValueError(
                "Unsupported supplementary feature type: "
                f"{feature_type}"
            )

        harmonized = [
            harmonizer(
                value
            )
            for value in result[
                "feature_raw"
            ]
        ]

        result[
            "harm_feature_id"
        ] = cls._series(
            [
                item[
                    "id"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_feature_identity_key"
        ] = cls._series(
            [
                item[
                    "key"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_feature_type"
        ] = feature_type

        result[
            "harm_feature_namespace"
        ] = cls._series(
            [
                item[
                    "namespace"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_feature_status"
        ] = cls._series(
            [
                item[
                    "status"
                ]
                for item in harmonized
            ],
            dtype="string",
        )

        result[
            "harm_feature_usable"
        ] = cls._series(
            [
                bool(
                    item[
                        "usable"
                    ]
                )
                for item in harmonized
            ],
            dtype="boolean",
        )

        result[
            "harm_feature_source_anomaly"
        ] = False

        return result

    # ------------------------------------------------------------------
    # Metadata annotation
    # ------------------------------------------------------------------

    @staticmethod
    def _add_metadata_annotation(
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Explicitly mark a source as metadata-only."""

        result = dataframe.copy()

        result[
            "harm_metadata_only"
        ] = True

        result[
            "harm_cross_source_join_ready"
        ] = False

        return result

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def harmonize(
        cls,
        dataframe: pd.DataFrame,
        *,
        sheet: str,
    ) -> pd.DataFrame:
        """
        Harmonize one supplementary Cancer-tRFQTL table.

        Input and output cardinality must remain identical.
        """

        cls.require_supported_sheet(
            sheet
        )

        spec = SHEET_REGISTRY[
            sheet
        ]

        result = dataframe.copy()

        input_rows = len(
            result
        )

        result[
            "harm_dataset_id"
        ] = DATASET_ID

        result[
            "harm_source_sheet"
        ] = sheet

        result[
            "harmonization_stage"
        ] = "M2-backfill"

        result[
            "harm_liftover_performed"
        ] = False

        if spec.harmonize_disease:
            result = (
                cls._add_disease_harmonization(
                    result
                )
            )

        if spec.harmonize_variant:
            result = (
                cls._add_variant_harmonization(
                    result
                )
            )

        if spec.harmonize_feature:
            if spec.feature_type is None:
                raise RuntimeError(
                    f"{sheet} requires feature harmonization "
                    "but no feature type was configured."
                )

            result = (
                cls._add_feature_harmonization(
                    result,
                    feature_type=spec.feature_type,
                )
            )

        if spec.metadata_only:
            result = (
                cls._add_metadata_annotation(
                    result
                )
            )
        else:
            result[
                "harm_metadata_only"
            ] = False

        if (
            len(
                result
            )
            != input_rows
        ):
            raise RuntimeError(
                f"Cancer-tRFQTL {sheet} harmonization changed "
                f"cardinality: {input_rows} -> {len(result)}"
            )

        return result