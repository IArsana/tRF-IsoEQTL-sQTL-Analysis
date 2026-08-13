"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/moradi_de.py

Description:
    Source-specific standardization utilities for Moradi differential
    expression supplementary datasets S5-S10.

    Supported datasets:

        S5:
            Exon-level differential expression, cis-related set.

        S6:
            Exon-level differential expression, trans-related set.

        S7:
            Intron-level differential expression, cis-related set.

        S8:
            Intron-level differential expression, trans-related set.

        S9:
            Transcript/isoform-level differential expression,
            cis-related set.

        S10:
            Transcript/isoform-level differential expression,
            trans-related set.

    Responsibilities:

        - canonical feature identifier standardization
        - DESeq2-like numerical field standardization
        - probability-field standardization
        - explicit feature type and analysis scope
        - feature usability flags
        - row-level standardization status
        - preservation of raw source identifiers

    Important:
        Known malformed source identifiers such as "INT1e+05" are
        retained as raw values and are not converted to canonical
        INT identifiers.

        These datasets represent differential-expression outputs and
        must not be treated as QTL datasets.

    This module performs no file I/O and never modifies raw data.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import re
from typing import Any, Literal

import pandas as pd

from pcatrfqtl.standardization.common import (
    StandardizationUtils,
)


FeatureType = Literal[
    "exon",
    "intron",
    "transcript",
]

AnalysisScope = Literal[
    "cis",
    "trans",
]


class MoradiDEStandardizer:
    """Standardize validated Moradi differential-expression records."""

    EXON_PATTERN = re.compile(
        r"^EX\d+$",
        re.IGNORECASE,
    )

    INTRON_PATTERN = re.compile(
        r"^INT\d+$",
        re.IGNORECASE,
    )

    TRANSCRIPT_PATTERN = re.compile(
        r"^ENST\d+(?:\.\d+)?$",
        re.IGNORECASE,
    )

    SCIENTIFIC_INTRON_PATTERN = re.compile(
        r"^INT\d+(?:\.\d+)?[eE][+-]?\d+$",
        re.IGNORECASE,
    )

    # ==================================================================
    # Feature identifiers
    # ==================================================================

    @classmethod
    def is_scientific_notation_like_intron(
        cls,
        value: Any,
    ) -> bool:
        """
        Detect malformed intron identifiers resembling scientific
        notation.

        Example:
            INT1e+05
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return False

        text = str(
            value
        ).strip()

        return (
            cls.SCIENTIFIC_INTRON_PATTERN
            .fullmatch(
                text
            )
            is not None
        )

    @classmethod
    def standardize_feature(
        cls,
        value: Any,
        feature_type: FeatureType,
    ) -> dict[str, Any]:
        """
        Standardize a Moradi DE feature identifier.

        Canonical examples:
            EX123
            INT123
            ENST00000318325
            ENST00000318325.6

        Malformed identifiers such as INT1e+05 are preserved only in
        feature_raw by the record-level standardizer and are not
        converted.
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "feature_id":
                    None,

                "feature_type":
                    feature_type,

                "feature_usable":
                    False,

                "feature_source":
                    "missing",

                "source_feature_anomaly":
                    False,
            }

        text = str(
            value
        ).strip()

        if feature_type == "exon":
            pattern = (
                cls.EXON_PATTERN
            )

        elif feature_type == "intron":
            pattern = (
                cls.INTRON_PATTERN
            )

        elif feature_type == "transcript":
            pattern = (
                cls.TRANSCRIPT_PATTERN
            )

        else:
            raise ValueError(
                "Unsupported Moradi DE feature type: "
                f"{feature_type!r}"
            )

        if (
            pattern.fullmatch(
                text
            )
            is None
        ):
            source_anomaly = (
                feature_type
                == "intron"
                and cls
                .is_scientific_notation_like_intron(
                    value
                )
            )

            return {
                "feature_id":
                    None,

                "feature_type":
                    feature_type,

                "feature_usable":
                    False,

                "feature_source":
                    (
                        "scientific_notation_like_source_identifier"
                        if source_anomaly
                        else "non_canonical_feature_id"
                    ),

                "source_feature_anomaly":
                    bool(
                        source_anomaly
                    ),
            }

        if (
            feature_type
            == "transcript"
        ):
            standardized = (
                text.upper()
            )

        else:
            prefix = (
                "EX"
                if feature_type
                == "exon"
                else "INT"
            )

            suffix = text[
                len(
                    prefix
                ):
            ]

            standardized = (
                prefix
                + suffix
            )

        return {
            "feature_id":
                standardized,

            "feature_type":
                feature_type,

            "feature_usable":
                True,

            "feature_source":
                "canonical_source_identifier",

            "source_feature_anomaly":
                False,
        }

    # ==================================================================
    # Numeric helpers
    # ==================================================================

    @classmethod
    def standardize_numeric(
        cls,
        value: Any,
        *,
        allow_missing: bool = True,
    ) -> float | None:
        """Return a finite numeric representation or None."""

        result = (
            StandardizationUtils
            .standardize_numeric(
                value,
                allow_missing=(
                    allow_missing
                ),
            )
        )

        if (
            result.standardized_value
            is None
        ):
            return None

        if not result.usable:
            return None

        return float(
            result.standardized_value
        )

    @classmethod
    def standardize_nonnegative_numeric(
        cls,
        value: Any,
    ) -> float | None:
        """
        Standardize a finite non-negative numeric field.

        Used for baseMean.
        """

        numeric = (
            cls.standardize_numeric(
                value
            )
        )

        if numeric is None:
            return None

        if numeric < 0:
            return None

        return numeric

    @classmethod
    def standardize_probability(
        cls,
        value: Any,
    ) -> float | None:
        """Return a probability in the inclusive range 0-1 or None."""

        result = (
            StandardizationUtils
            .standardize_probability(
                value,
                allow_zero=True,
                allow_missing=True,
            )
        )

        if (
            result.standardized_value
            is None
        ):
            return None

        if not result.usable:
            return None

        return float(
            result.standardized_value
        )

    # ==================================================================
    # Record
    # ==================================================================

    @classmethod
    def standardize_record(
        cls,
        record: dict[str, Any],
        *,
        source_table: str,
        feature_type: FeatureType,
        analysis_scope: AnalysisScope,
    ) -> dict[str, Any]:
        """
        Standardize one Moradi differential-expression record.
        """

        feature_raw = (
            record.get(
                "Unnamed: 0"
            )
        )

        feature = (
            cls.standardize_feature(
                feature_raw,
                feature_type,
            )
        )

        base_mean = (
            cls.standardize_nonnegative_numeric(
                record.get(
                    "baseMean"
                )
            )
        )

        log2_fold_change = (
            cls.standardize_numeric(
                record.get(
                    "log2FoldChange"
                )
            )
        )

        lfc_se = (
            cls.standardize_numeric(
                record.get(
                    "lfcSE"
                )
            )
        )

        statistic = (
            cls.standardize_numeric(
                record.get(
                    "stat"
                )
            )
        )

        p_value = (
            cls.standardize_probability(
                record.get(
                    "pvalue"
                )
            )
        )

        adjusted_p_value = (
            cls.standardize_probability(
                record.get(
                    "padj"
                )
            )
        )

        de_statistics_usable = (
            base_mean is not None
            and log2_fold_change is not None
            and lfc_se is not None
            and statistic is not None
            and p_value is not None
        )

        adjusted_p_value_available = (
            adjusted_p_value
            is not None
        )

        de_record_usable = (
            feature[
                "feature_usable"
            ]
            and de_statistics_usable
        )

        if de_record_usable:
            standardization_status = (
                "STANDARDIZED"
            )

        elif (
            feature[
                "feature_usable"
            ]
            or de_statistics_usable
        ):
            standardization_status = (
                "PARTIAL"
            )

        else:
            standardization_status = (
                "PRESERVED"
            )

        return {
            "source_table":
                source_table,

            "analysis_scope":
                analysis_scope,

            "data_type":
                "differential_expression",

            "feature_raw":
                feature_raw,

            **feature,

            "base_mean_raw":
                record.get(
                    "baseMean"
                ),

            "base_mean":
                base_mean,

            "log2_fold_change_raw":
                record.get(
                    "log2FoldChange"
                ),

            "log2_fold_change":
                log2_fold_change,

            "lfc_se_raw":
                record.get(
                    "lfcSE"
                ),

            "lfc_se":
                lfc_se,

            "statistic_raw":
                record.get(
                    "stat"
                ),

            "statistic":
                statistic,

            "p_value_raw":
                record.get(
                    "pvalue"
                ),

            "p_value":
                p_value,

            "adjusted_p_value_raw":
                record.get(
                    "padj"
                ),

            "adjusted_p_value":
                adjusted_p_value,

            "de_statistics_usable":
                bool(
                    de_statistics_usable
                ),

            "adjusted_p_value_available":
                bool(
                    adjusted_p_value_available
                ),

            "de_record_usable":
                bool(
                    de_record_usable
                ),

            "standardization_status":
                standardization_status,
        }

    # ==================================================================
    # DataFrame
    # ==================================================================

    @classmethod
    def standardize_dataframe(
        cls,
        dataframe: pd.DataFrame,
        *,
        source_table: str,
        feature_type: FeatureType,
        analysis_scope: AnalysisScope,
    ) -> pd.DataFrame:
        """
        Standardize one validated Moradi DE DataFrame.

        Row cardinality is preserved.
        """

        records = (
            dataframe
            .to_dict(
                orient="records"
            )
        )

        standardized = [
            cls.standardize_record(
                record,
                source_table=(
                    source_table
                ),
                feature_type=(
                    feature_type
                ),
                analysis_scope=(
                    analysis_scope
                ),
            )
            for record
            in records
        ]

        return pd.DataFrame(
            standardized
        )