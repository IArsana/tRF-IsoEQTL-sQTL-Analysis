"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/common.py

Description:
    Source-independent standardization primitives for the pcatRFQTL
    research pipeline.

    Utilities include:

        - missing-value detection
        - whitespace-normalized strings
        - chromosome standardization
        - genomic-position standardization
        - finite numeric conversion
        - probability conversion
        - boolean conversion
        - identifier normalization

    These utilities derive standardized representations only.
    Raw source values are never modified.

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

from pcatrfqtl.standardization.models import (
    StandardizationResult,
    StandardizationStatus,
)


class StandardizationUtils:
    """Common source-independent standardization utilities."""

    MISSING_STRINGS = {
        "",
        "na",
        "nan",
        "none",
        "null",
    }

    CHROMOSOME_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)$",
        re.IGNORECASE,
    )

    RSID_PATTERN = re.compile(
        r"^rs\d+$",
        re.IGNORECASE,
    )

    ENST_PATTERN = re.compile(
        r"^ENST\d+(?:\.\d+)?$",
        re.IGNORECASE,
    )

    # ==================================================================
    # Missing values
    # ==================================================================

    @classmethod
    def is_missing(
        cls,
        value: Any,
    ) -> bool:
        """Return True when a source value should be treated as missing."""

        if value is None:
            return True

        try:
            if pd.isna(
                value
            ):
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
                value
                .strip()
                .lower()
                in cls.MISSING_STRINGS
            )

        return False

    # ==================================================================
    # Strings
    # ==================================================================

    @classmethod
    def standardize_string(
        cls,
        value: Any,
        *,
        collapse_whitespace: bool = True,
        empty_as_missing: bool = True,
    ) -> StandardizationResult:
        """Standardize a generic string while preserving source provenance."""

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=False,
            )

        text = str(
            value
        ).strip()

        if collapse_whitespace:
            text = re.sub(
                r"\s+",
                " ",
                text,
            )

        if (
            empty_as_missing
            and not text
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="empty_string",
                usable=False,
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=text,
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source="string_cleanup",
            usable=True,
        )

    # ==================================================================
    # Chromosomes
    # ==================================================================

    @classmethod
    def standardize_chromosome(
        cls,
        value: Any,
    ) -> StandardizationResult:
        """
        Standardize chromosome labels.

        Examples:
            chr1 -> 1
            chrX -> X
            M    -> MT
            chrM -> MT
            MT   -> MT
        """

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=False,
            )

        text = str(
            value
        ).strip()

        if (
            cls.CHROMOSOME_PATTERN
            .fullmatch(
                text
            )
            is None
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="unrecognized_chromosome",
                usable=False,
                note=(
                    "Value does not match a canonical "
                    "human chromosome representation."
                ),
            )

        if text.lower().startswith(
            "chr"
        ):
            text = text[
                3:
            ]

        text = text.upper()

        if text == "M":
            text = "MT"

        if text.isdigit():
            text = str(
                int(
                    text
                )
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=text,
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source="chromosome_normalization",
            usable=True,
        )

    # ==================================================================
    # Genomic positions
    # ==================================================================

    @classmethod
    def standardize_position(
        cls,
        value: Any,
    ) -> StandardizationResult:
        """Standardize a positive integer genomic position."""

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=False,
            )

        try:
            numeric = float(
                value
            )
        except (
            TypeError,
            ValueError,
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_numeric_position",
                usable=False,
            )

        if (
            not math.isfinite(
                numeric
            )
            or numeric <= 0
            or not numeric.is_integer()
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="invalid_position",
                usable=False,
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=int(
                numeric
            ),
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source="integer_position",
            usable=True,
        )

    # ==================================================================
    # Numeric values
    # ==================================================================

    @classmethod
    def standardize_numeric(
        cls,
        value: Any,
        *,
        allow_missing: bool = True,
    ) -> StandardizationResult:
        """Convert a source value to a finite float."""

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=allow_missing,
            )

        try:
            numeric = float(
                value
            )
        except (
            TypeError,
            ValueError,
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_numeric",
                usable=False,
            )

        if not math.isfinite(
            numeric
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_finite_numeric",
                usable=False,
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=numeric,
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source="numeric_conversion",
            usable=True,
        )

    # ==================================================================
    # Probabilities
    # ==================================================================

    @classmethod
    def standardize_probability(
        cls,
        value: Any,
        *,
        allow_zero: bool = True,
        allow_missing: bool = True,
    ) -> StandardizationResult:
        """
        Standardize probability-like values.

        Default valid range:
            0 <= value <= 1
        """

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=allow_missing,
            )

        try:
            numeric = float(
                value
            )
        except (
            TypeError,
            ValueError,
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_numeric_probability",
                usable=False,
            )

        if not math.isfinite(
            numeric
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_finite_probability",
                usable=False,
            )

        if (
            numeric < 0
            or numeric > 1
            or (
                not allow_zero
                and numeric == 0
            )
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="out_of_range_probability",
                usable=False,
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=numeric,
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source="probability_conversion",
            usable=True,
        )

    # ==================================================================
    # Boolean values
    # ==================================================================

    @classmethod
    def standardize_boolean(
        cls,
        value: Any,
    ) -> StandardizationResult:
        """Standardize common boolean representations."""

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=False,
            )

        if isinstance(
            value,
            bool,
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=value,
                status=(
                    StandardizationStatus
                    .PRESERVED
                ),
                source="boolean",
                usable=True,
            )

        if isinstance(
            value,
            int,
        ) and value in {
            0,
            1,
        }:
            return StandardizationResult(
                raw_value=value,
                standardized_value=bool(
                    value
                ),
                status=(
                    StandardizationStatus
                    .STANDARDIZED
                ),
                source="integer_boolean",
                usable=True,
            )

        text = str(
            value
        ).strip().lower()

        truthy = {
            "1",
            "true",
            "yes",
            "y",
        }

        falsy = {
            "0",
            "false",
            "no",
            "n",
        }

        if text in truthy:
            return StandardizationResult(
                raw_value=value,
                standardized_value=True,
                status=(
                    StandardizationStatus
                    .STANDARDIZED
                ),
                source="string_boolean",
                usable=True,
            )

        if text in falsy:
            return StandardizationResult(
                raw_value=value,
                standardized_value=False,
                status=(
                    StandardizationStatus
                    .STANDARDIZED
                ),
                source="string_boolean",
                usable=True,
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=None,
            status=(
                StandardizationStatus
                .UNSUPPORTED
            ),
            source="unrecognized_boolean",
            usable=False,
        )

    # ==================================================================
    # Canonical identifiers
    # ==================================================================

    @classmethod
    def standardize_rsid(
        cls,
        value: Any,
    ) -> StandardizationResult:
        """Standardize a canonical dbSNP rsID."""

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=False,
            )

        text = str(
            value
        ).strip()

        if (
            cls.RSID_PATTERN
            .fullmatch(
                text
            )
            is None
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_canonical_rsid",
                usable=False,
            )

        standardized = (
            "rs"
            + text[
                2:
            ]
        )

        return StandardizationResult(
            raw_value=value,
            standardized_value=standardized,
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source="rsid_normalization",
            usable=True,
        )

    @classmethod
    def standardize_transcript_id(
        cls,
        value: Any,
        *,
        strip_version: bool = False,
    ) -> StandardizationResult:
        """
        Standardize an Ensembl transcript identifier.

        Version suffixes are preserved by default.
        """

        if cls.is_missing(
            value
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .MISSING
                ),
                source="missing",
                usable=False,
            )

        text = str(
            value
        ).strip().upper()

        if (
            cls.ENST_PATTERN
            .fullmatch(
                text
            )
            is None
        ):
            return StandardizationResult(
                raw_value=value,
                standardized_value=None,
                status=(
                    StandardizationStatus
                    .UNSUPPORTED
                ),
                source="non_canonical_transcript_id",
                usable=False,
            )

        if (
            strip_version
            and "."
            in text
        ):
            text = text.split(
                ".",
                maxsplit=1,
            )[0]

            source = (
                "ensembl_transcript_version_removed"
            )

        else:
            source = (
                "ensembl_transcript_normalization"
            )

        return StandardizationResult(
            raw_value=value,
            standardized_value=text,
            status=(
                StandardizationStatus
                .STANDARDIZED
            ),
            source=source,
            usable=True,
        )