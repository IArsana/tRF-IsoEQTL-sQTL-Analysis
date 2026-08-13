"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/validators/moradi_de.py

Description:
    Field-level validation utilities for Moradi differential-expression
    supplementary datasets S5-S10.

    Supported source feature types:
        - exon events: EX*
        - intron events: INT*
        - transcript isoforms: ENST*

    Supported statistical fields:
        - baseMean
        - log2FoldChange
        - lfcSE
        - stat
        - pvalue
        - padj

    The validator performs structural and numeric validation only.
    Raw source data are never modified.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import math
import re
from typing import Any

import pandas as pd


class MoradiDEValidator:
    """Field-level validator for Moradi DE datasets S5-S10."""

    ANOMALY_IDS = {
        "feature": "MORADI-DE-FEATURE-001",
        "base_mean": "MORADI-DE-BASEMEAN-001",
        "numeric": "MORADI-DE-NUMERIC-001",
        "probability": "MORADI-DE-PROBABILITY-001",
        "missing": "MORADI-DE-MISSING-001",
    }

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

    @staticmethod
    def _is_missing(
        value: Any,
    ) -> bool:
        """Return True when a value should be treated as missing."""

        if value is None:
            return True

        try:
            if pd.isna(value):
                return True
        except (
            TypeError,
            ValueError,
        ):
            pass

        if isinstance(
            value,
            str,
        ):
            return (
                value.strip().lower()
                in {
                    "",
                    "na",
                    "nan",
                    "none",
                    "null",
                }
            )

        return False

    @classmethod
    def _format_error(
        cls,
        rule_id: str,
        field_name: str,
        index: int,
        value: Any,
        message: str,
    ) -> str:
        """Create a standardized validation message."""

        return (
            f"[{rule_id}] "
            f"{field_name}[{index}] "
            f"{message}: {value!r}"
        )

    # ==================================================================
    # Feature identifiers
    # ==================================================================

    @classmethod
    def validate_features(
        cls,
        values: list[Any],
        feature_type: str,
        field_name: str = "feature_id",
        *,
        allow_missing: bool = False,
    ) -> list[str]:
        """
        Validate Moradi feature identifiers.

        Parameters
        ----------
        feature_type:
            One of:
                exon
                intron
                transcript
        """

        patterns = {
            "exon":
                cls.EXON_PATTERN,
            "intron":
                cls.INTRON_PATTERN,
            "transcript":
                cls.TRANSCRIPT_PATTERN,
        }

        if feature_type not in patterns:
            raise ValueError(
                "Unsupported Moradi DE feature type: "
                f"{feature_type!r}"
            )

        pattern = patterns[
            feature_type
        ]

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):

            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            (
                                "missing feature "
                                "identifier"
                            ),
                        )
                    )

                continue

            text = str(
                value
            ).strip()

            if (
                pattern.fullmatch(
                    text
                )
                is None
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "feature"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "invalid "
                            f"{feature_type} "
                            "identifier"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # baseMean
    # ==================================================================

    @classmethod
    def validate_base_mean(
        cls,
        values: list[Any],
        field_name: str = "baseMean",
        *,
        allow_missing: bool = False,
    ) -> list[str]:
        """Validate non-negative finite baseMean values."""

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):

            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            "missing baseMean",
                        )
                    )

                continue

            try:
                numeric = float(
                    value
                )
            except (
                TypeError,
                ValueError,
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "base_mean"
                        ],
                        field_name,
                        index,
                        value,
                        "non-numeric baseMean",
                    )
                )

                continue

            if (
                not math.isfinite(
                    numeric
                )
                or numeric < 0
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "base_mean"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "baseMean must be "
                            "finite and >= 0"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # Generic numeric fields
    # ==================================================================

    @classmethod
    def validate_numeric(
        cls,
        values: list[Any],
        field_name: str,
        *,
        allow_missing: bool = False,
    ) -> list[str]:
        """Validate generic finite numeric fields."""

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):

            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            (
                                "missing numeric "
                                "value"
                            ),
                        )
                    )

                continue

            try:
                numeric = float(
                    value
                )
            except (
                TypeError,
                ValueError,
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "numeric"
                        ],
                        field_name,
                        index,
                        value,
                        "non-numeric value",
                    )
                )

                continue

            if not math.isfinite(
                numeric
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "numeric"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "non-finite numeric "
                            "value"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # Probability fields
    # ==================================================================

    @classmethod
    def validate_probability(
        cls,
        values: list[Any],
        field_name: str,
        *,
        allow_missing: bool = False,
    ) -> list[str]:
        """
        Validate probability-like values.

        Valid range:
            0 <= value <= 1
        """

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):

            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            (
                                "missing probability "
                                "value"
                            ),
                        )
                    )

                continue

            try:
                numeric = float(
                    value
                )
            except (
                TypeError,
                ValueError,
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "probability"
                        ],
                        field_name,
                        index,
                        value,
                        "non-numeric probability",
                    )
                )

                continue

            if (
                not math.isfinite(
                    numeric
                )
                or numeric < 0
                or numeric > 1
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "probability"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "value outside valid "
                            "range 0 <= value <= 1"
                        ),
                    )
                )

        return errors