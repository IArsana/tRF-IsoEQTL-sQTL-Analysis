"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/runners/gwas_catalog.py

Description:
    Dataset-level validation runner for the GWAS Catalog association
    export used by the pcatRFQTL research pipeline.

    The runner validates the complete GWAS Catalog TSV using chunked
    processing to reduce memory usage.

    Validation includes:

        - required schema
        - heterogeneous GWAS variant descriptors
        - chromosome identifiers
        - interaction-aware chromosome representations
        - genomic positions
        - interaction-aware position representations
        - association P-values
        - missing genomic mappings
        - missing P-values
        - zero P-values
        - exact duplicate rows
        - variant descriptor classification statistics

    Important design decisions:

        - GWAS Catalog "SNPS" is not treated as a strict rsID field.
        - Missing genomic mappings are retained and reported as WARNING.
        - Missing P-values are retained and reported as WARNING.
        - P-VALUE = 0 is retained and reported as WARNING.
        - Source descriptors classified as "other" are retained.
        - Raw source data are never modified.

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
from pathlib import Path
from typing import Any

import pandas as pd

from pcatrfqtl.validation.anomaly import (
    AnomalyFactory,
)
from pcatrfqtl.validation.validators.gwas import (
    GWASValidator,
)


class GWASCatalogValidator:
    """
    Validate a GWAS Catalog association TSV export.

    The source file is processed in chunks so that large GWAS Catalog
    exports can be validated without loading the entire dataset into
    memory.
    """

    DATASET_ID = "gwas_catalog"

    REQUIRED_COLUMNS = [
        "SNPS",
        "CHR_ID",
        "CHR_POS",
        "P-VALUE",
    ]

    ERROR_PATTERN = re.compile(
        r"^(?:\[[^\]]+\]\s*)?"
        r"(?P<field>.+?)"
        r"\[(?P<index>\d+)\]"
        r".*?:\s*"
        r"(?P<value>.+)$"
    )

    def __init__(
        self,
        file_path: str | Path,
        *,
        chunk_size: int = 100_000,
        anomaly_sample_limit: int = 1_000,
    ) -> None:
        """Initialize GWAS Catalog validation."""

        self.file_path = Path(
            file_path
        )

        self.chunk_size = int(
            chunk_size
        )

        self.anomaly_sample_limit = int(
            anomaly_sample_limit
        )

        self.validator = (
            GWASValidator()
        )

        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies: list[
            dict[str, Any]
        ] = []

        self.total_anomaly_occurrences = 0

    # ==================================================================
    # Status
    # ==================================================================

    @staticmethod
    def _status_from_severity_counts(
        severity_counts: dict[
            str,
            int,
        ],
    ) -> str:
        """Determine final validation status."""

        if (
            severity_counts.get(
                "ERROR",
                0,
            )
            > 0
        ):
            return "FAIL"

        if (
            severity_counts.get(
                "WARNING",
                0,
            )
            > 0
        ):
            return "WARNING"

        return "PASS"

    # ==================================================================
    # Anomaly handling
    # ==================================================================

    def _add_anomaly(
        self,
        *,
        row: int | None,
        column: str | None,
        value: Any,
        anomaly_type: str,
        category: str,
        severity: str,
        message: str,
        action: str,
    ) -> None:
        """
        Register a validation anomaly.

        All anomaly occurrences contribute to the total count, while
        only a limited number are stored in the JSON report.
        """

        self.total_anomaly_occurrences += 1

        if (
            len(self.anomalies)
            >= self.anomaly_sample_limit
        ):
            return

        anomaly = (
            self.anomaly_factory
            .create(
                type=anomaly_type,
                category=category,
                severity=severity,
                dataset=self.DATASET_ID,
                file=str(
                    self.file_path
                ),
                sheet=None,
                row=row,
                column=column,
                value=value,
                message=message,
                action=action,
            )
        )

        self.anomalies.append(
            anomaly.to_dict()
        )

    @classmethod
    def _extract_error_metadata(
        cls,
        error: str,
    ) -> tuple[
        int | None,
        Any,
    ]:
        """
        Extract a chunk-local index and value from a validation error.
        """

        match = (
            cls.ERROR_PATTERN.search(
                error
            )
        )

        if match is None:
            return None, None

        index = int(
            match.group(
                "index"
            )
        )

        value = (
            match.group(
                "value"
            )
            .strip()
            .rstrip(".")
        )

        if (
            len(value) >= 2
            and value[0]
            == value[-1]
            and value[0]
            in {
                "'",
                '"',
            }
        ):
            value = value[
                1:-1
            ]

        return index, value

    def _add_validation_errors(
        self,
        *,
        errors: list[str],
        chunk: pd.DataFrame,
        chunk_start_row: int,
        column: str,
        anomaly_type: str,
        category: str,
        severity: str,
        action: str,
    ) -> None:
        """
        Convert field-level validator messages into anomaly records.
        """

        for error in errors:

            (
                local_index,
                parsed_value,
            ) = (
                self
                ._extract_error_metadata(
                    error
                )
            )

            row: int | None = None
            value: Any = (
                parsed_value
            )

            if local_index is not None:

                # TSV line 1 contains the header.
                row = (
                    chunk_start_row
                    + local_index
                    + 2
                )

                if (
                    0 <= local_index
                    < len(chunk)
                    and column
                    in chunk.columns
                ):
                    value = (
                        chunk.iloc[
                            local_index
                        ][column]
                    )

            self._add_anomaly(
                row=row,
                column=column,
                value=value,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity=severity,
                message=error,
                action=action,
            )

    # ==================================================================
    # Schema
    # ==================================================================

    def _read_columns(
        self,
    ) -> list[str]:
        """Read the TSV header only."""

        frame = pd.read_csv(
            self.file_path,
            sep="\t",
            nrows=0,
        )

        return [
            str(column).strip()
            for column
            in frame.columns
        ]

    def _validate_schema(
        self,
        columns: list[str],
    ) -> list[str]:
        """Return missing required columns."""

        available = set(
            columns
        )

        return [
            column
            for column
            in self.REQUIRED_COLUMNS
            if column not in available
        ]

    # ==================================================================
    # Missing values
    # ==================================================================

    @staticmethod
    def _missing_count(
        series: pd.Series,
    ) -> int:
        """Count values considered missing by GWASValidator."""

        return int(
            series
            .map(
                GWASValidator
                ._is_missing
            )
            .sum()
        )

    # ==================================================================
    # Validation
    # ==================================================================

    def validate(
        self,
    ) -> dict[str, Any]:
        """Validate the complete GWAS Catalog association TSV."""

        # --------------------------------------------------------------
        # Source file
        # --------------------------------------------------------------

        if not self.file_path.exists():
            raise FileNotFoundError(
                "GWAS Catalog file not found: "
                f"{self.file_path}"
            )

        if not self.file_path.is_file():
            raise ValueError(
                "GWAS Catalog path is not a file: "
                f"{self.file_path}"
            )

        # Make validator instance reusable.
        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies = []

        self.total_anomaly_occurrences = 0

        # --------------------------------------------------------------
        # Schema
        # --------------------------------------------------------------

        columns = (
            self._read_columns()
        )

        missing_columns = (
            self._validate_schema(
                columns
            )
        )

        if missing_columns:

            for column in missing_columns:

                self._add_anomaly(
                    row=None,
                    column=column,
                    value=None,
                    anomaly_type=(
                        "MISSING_REQUIRED_COLUMN"
                    ),
                    category="SCHEMA",
                    severity="ERROR",
                    message=(
                        "Required GWAS Catalog column "
                        f"{column!r} was not found."
                    ),
                    action=(
                        "STOP_DATASET_VALIDATION"
                    ),
                )

            return {
                "dataset_id":
                    self.DATASET_ID,

                "file":
                    str(
                        self.file_path
                    ),

                "status":
                    "FAIL",

                "rows":
                    0,

                "columns":
                    len(columns),

                "chunks_processed":
                    0,

                "chunk_size":
                    self.chunk_size,

                "schema": {
                    "available_columns":
                        columns,
                    "required_columns":
                        list(
                            self.REQUIRED_COLUMNS
                        ),
                    "missing_required_columns":
                        missing_columns,
                },

                "validation_error_counts":
                    {},

                "variant_descriptor_statistics":
                    {},

                "missing_values":
                    {},

                "mapping_statistics":
                    {},

                "pvalue_statistics":
                    {},

                "duplicates":
                    {},

                "anomalies":
                    list(
                        self.anomalies
                    ),

                "validation_summary": {
                    "severity_counts": {
                        "WARNING": 0,
                        "ERROR":
                            len(
                                missing_columns
                            ),
                    },

                    "total_anomaly_occurrences":
                        int(
                            self
                            .total_anomaly_occurrences
                        ),

                    "stored_anomaly_samples":
                        len(
                            self.anomalies
                        ),

                    "anomaly_sample_limit":
                        int(
                            self
                            .anomaly_sample_limit
                        ),
                },
            }

        # --------------------------------------------------------------
        # Counters
        # --------------------------------------------------------------

        total_rows = 0
        chunk_count = 0

        validation_error_counts = {
            "variant_descriptors": 0,
            "chromosomes": 0,
            "positions": 0,
            "p_values": 0,
        }

        variant_descriptor_statistics: dict[
            str,
            int,
        ] = {}

        missing_values = {
            column: 0
            for column in columns
        }

        mapping_statistics = {
            "missing_snps": 0,
            "missing_chromosomes": 0,
            "missing_positions": 0,
            "complete_genomic_mapping": 0,
        }

        pvalue_statistics = {
            "missing_p_values": 0,
            "zero_p_values": 0,
        }

        severity_counts = {
            "WARNING": 0,
            "ERROR": 0,
        }

        seen_hashes: set[int] = (
            set()
        )

        duplicate_excess_rows = 0

        # --------------------------------------------------------------
        # Chunk reader
        # --------------------------------------------------------------

        reader = pd.read_csv(
            self.file_path,
            sep="\t",
            chunksize=self.chunk_size,
            low_memory=False,
        )

        for chunk in reader:

            chunk_count += 1

            chunk = chunk.copy()

            chunk.columns = [
                str(column).strip()
                for column
                in chunk.columns
            ]

            chunk_start_row = (
                total_rows
            )

            total_rows += len(
                chunk
            )

            # ==========================================================
            # Missing statistics
            # ==========================================================

            for column in columns:

                if column not in chunk.columns:
                    continue

                missing_values[
                    column
                ] += (
                    self
                    ._missing_count(
                        chunk[column]
                    )
                )

            snp_missing = (
                chunk["SNPS"]
                .map(
                    self.validator
                    ._is_missing
                )
            )

            chromosome_missing = (
                chunk["CHR_ID"]
                .map(
                    self.validator
                    ._is_missing
                )
            )

            position_missing = (
                chunk["CHR_POS"]
                .map(
                    self.validator
                    ._is_missing
                )
            )

            pvalue_missing = (
                chunk["P-VALUE"]
                .map(
                    self.validator
                    ._is_missing
                )
            )

            mapping_statistics[
                "missing_snps"
            ] += int(
                snp_missing.sum()
            )

            mapping_statistics[
                "missing_chromosomes"
            ] += int(
                chromosome_missing.sum()
            )

            mapping_statistics[
                "missing_positions"
            ] += int(
                position_missing.sum()
            )

            pvalue_statistics[
                "missing_p_values"
            ] += int(
                pvalue_missing.sum()
            )

            complete_mapping = (
                ~snp_missing
                & ~chromosome_missing
                & ~position_missing
            )

            mapping_statistics[
                "complete_genomic_mapping"
            ] += int(
                complete_mapping.sum()
            )

            # ==========================================================
            # Variant descriptor classification
            # ==========================================================

            classifications = (
                chunk["SNPS"]
                .map(
                    self.validator
                    .classify_variant_descriptor
                )
            )

            for (
                classification,
                count,
            ) in (
                classifications
                .value_counts(
                    dropna=False
                )
                .items()
            ):

                key = str(
                    classification
                )

                variant_descriptor_statistics[
                    key
                ] = (
                    variant_descriptor_statistics
                    .get(
                        key,
                        0,
                    )
                    + int(count)
                )

            # ==========================================================
            # Field validation
            # ==========================================================

            variant_errors = (
                self.validator
                .validate_variant_descriptors(
                    chunk[
                        "SNPS"
                    ].tolist(),
                    "SNPS",
                    allow_missing=True,
                )
            )

            chromosome_errors = (
                self.validator
                .validate_chromosomes(
                    chunk[
                        "CHR_ID"
                    ].tolist(),
                    "CHR_ID",
                    allow_missing=True,
                )
            )

            position_errors = (
                self.validator
                .validate_positions(
                    chunk[
                        "CHR_POS"
                    ].tolist(),
                    "CHR_POS",
                    allow_missing=True,
                )
            )

            pvalue_errors = (
                self.validator
                .validate_p_values(
                    chunk[
                        "P-VALUE"
                    ].tolist(),
                    "P-VALUE",
                    allow_zero=True,
                    allow_missing=True,
                )
            )

            validation_error_counts[
                "variant_descriptors"
            ] += len(
                variant_errors
            )

            validation_error_counts[
                "chromosomes"
            ] += len(
                chromosome_errors
            )

            validation_error_counts[
                "positions"
            ] += len(
                position_errors
            )

            validation_error_counts[
                "p_values"
            ] += len(
                pvalue_errors
            )

            # ==========================================================
            # P-value statistics
            # ==========================================================

            numeric_p = pd.to_numeric(
                chunk["P-VALUE"],
                errors="coerce",
            )

            pvalue_statistics[
                "zero_p_values"
            ] += int(
                (
                    numeric_p
                    == 0
                ).sum()
            )

            # ==========================================================
            # Convert genuine validation errors to anomalies
            # ==========================================================

            field_error_groups = [
                (
                    variant_errors,
                    "SNPS",
                    "INVALID_VARIANT_DESCRIPTOR",
                    "FORMAT",
                ),
                (
                    chromosome_errors,
                    "CHR_ID",
                    "INVALID_CHROMOSOME",
                    "COORDINATE",
                ),
                (
                    position_errors,
                    "CHR_POS",
                    "INVALID_POSITION",
                    "COORDINATE",
                ),
                (
                    pvalue_errors,
                    "P-VALUE",
                    "INVALID_P_VALUE",
                    "RANGE",
                ),
            ]

            for (
                errors,
                column,
                anomaly_type,
                category,
            ) in field_error_groups:

                if not errors:
                    continue

                self._add_validation_errors(
                    errors=errors,
                    chunk=chunk,
                    chunk_start_row=(
                        chunk_start_row
                    ),
                    column=column,
                    anomaly_type=(
                        anomaly_type
                    ),
                    category=category,
                    severity="ERROR",
                    action=(
                        "FLAG_FOR_REVIEW"
                    ),
                )

                severity_counts[
                    "ERROR"
                ] += len(
                    errors
                )

            # ==========================================================
            # Exact duplicate detection
            # ==========================================================

            row_hashes = (
                pd.util
                .hash_pandas_object(
                    chunk,
                    index=False,
                )
                .astype(
                    "uint64"
                )
            )

            for hash_value in row_hashes:

                numeric_hash = int(
                    hash_value
                )

                if (
                    numeric_hash
                    in seen_hashes
                ):
                    duplicate_excess_rows += 1

                else:
                    seen_hashes.add(
                        numeric_hash
                    )

        # ==================================================================
        # Dataset-level warnings
        # ==================================================================

        # --------------------------------------------------------------
        # Missing genomic mapping
        # --------------------------------------------------------------

        missing_mapping_fields = [
            (
                "SNPS",
                mapping_statistics[
                    "missing_snps"
                ],
            ),
            (
                "CHR_ID",
                mapping_statistics[
                    "missing_chromosomes"
                ],
            ),
            (
                "CHR_POS",
                mapping_statistics[
                    "missing_positions"
                ],
            ),
        ]

        for (
            column,
            count,
        ) in missing_mapping_fields:

            if count <= 0:
                continue

            severity_counts[
                "WARNING"
            ] += 1

            self._add_anomaly(
                row=None,
                column=column,
                value=None,
                anomaly_type=(
                    "MISSING_GENOMIC_MAPPING"
                ),
                category="MISSING",
                severity="WARNING",
                message=(
                    f"{count:,} GWAS Catalog "
                    f"records have missing "
                    f"{column} values."
                ),
                action=(
                    "RETAIN_RECORD_AND_EXCLUDE_FROM_"
                    "COORDINATE_DEPENDENT_ANALYSIS"
                ),
            )

        # --------------------------------------------------------------
        # Missing P-values
        # --------------------------------------------------------------

        missing_p_values = int(
            pvalue_statistics[
                "missing_p_values"
            ]
        )

        if missing_p_values > 0:

            severity_counts[
                "WARNING"
            ] += 1

            self._add_anomaly(
                row=None,
                column="P-VALUE",
                value=None,
                anomaly_type=(
                    "MISSING_P_VALUE"
                ),
                category="MISSING",
                severity="WARNING",
                message=(
                    f"{missing_p_values:,} GWAS Catalog "
                    "records have missing P-VALUE values."
                ),
                action=(
                    "RETAIN_RAW_RECORD_AND_EXCLUDE_FROM_"
                    "P_VALUE_DEPENDENT_ANALYSIS"
                ),
            )

        # --------------------------------------------------------------
        # Zero P-values
        # --------------------------------------------------------------

        zero_p_values = int(
            pvalue_statistics[
                "zero_p_values"
            ]
        )

        if zero_p_values > 0:

            severity_counts[
                "WARNING"
            ] += 1

            self._add_anomaly(
                row=None,
                column="P-VALUE",
                value=None,
                anomaly_type=(
                    "ZERO_P_VALUE"
                ),
                category=(
                    "NUMERIC_PRECISION"
                ),
                severity="WARNING",
                message=(
                    f"{zero_p_values:,} GWAS Catalog "
                    "records contain P-VALUE = 0. "
                    "Source values are retained without "
                    "replacement or imputation."
                ),
                action=(
                    "RETAIN_RAW_VALUE_AND_HANDLE_DURING_"
                    "NORMALIZATION"
                ),
            )

        # --------------------------------------------------------------
        # Unclassified variant descriptors
        # --------------------------------------------------------------

        other_count = int(
            variant_descriptor_statistics
            .get(
                "other",
                0,
            )
        )

        if other_count > 0:

            severity_counts[
                "WARNING"
            ] += 1

            self._add_anomaly(
                row=None,
                column="SNPS",
                value=None,
                anomaly_type=(
                    "UNCLASSIFIED_VARIANT_DESCRIPTOR"
                ),
                category="SOURCE",
                severity="WARNING",
                message=(
                    f"{other_count:,} GWAS Catalog "
                    "records contain non-empty source "
                    "variant descriptors outside the "
                    "currently classified representation "
                    "set."
                ),
                action=(
                    "RETAIN_SOURCE_DESCRIPTOR_AND_REVIEW_"
                    "DURING_NORMALIZATION"
                ),
            )

        # --------------------------------------------------------------
        # Duplicate associations
        # --------------------------------------------------------------

        if duplicate_excess_rows > 0:

            severity_counts[
                "WARNING"
            ] += 1

            self._add_anomaly(
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "EXACT_DUPLICATE_ASSOCIATION"
                ),
                category="DUPLICATE",
                severity="WARNING",
                message=(
                    f"{duplicate_excess_rows:,} exact "
                    "duplicate excess association rows "
                    "were detected."
                ),
                action=(
                    "RETAIN_RAW_AND_REVIEW_DURING_"
                    "NORMALIZATION"
                ),
            )

        # ==================================================================
        # Final status
        # ==================================================================

        status = (
            self
            ._status_from_severity_counts(
                severity_counts
            )
        )

        # ==================================================================
        # Report
        # ==================================================================

        return {
            "dataset_id":
                self.DATASET_ID,

            "file":
                str(
                    self.file_path
                ),

            "status":
                status,

            "rows":
                int(
                    total_rows
                ),

            "columns":
                int(
                    len(columns)
                ),

            "chunks_processed":
                int(
                    chunk_count
                ),

            "chunk_size":
                int(
                    self.chunk_size
                ),

            "schema": {
                "available_columns":
                    columns,

                "required_columns":
                    list(
                        self.REQUIRED_COLUMNS
                    ),

                "missing_required_columns":
                    [],
            },

            "validation_error_counts":
                validation_error_counts,

            "variant_descriptor_statistics":
                dict(
                    sorted(
                        variant_descriptor_statistics
                        .items()
                    )
                ),

            "missing_values":
                missing_values,

            "mapping_statistics":
                mapping_statistics,

            "pvalue_statistics":
                pvalue_statistics,

            "duplicates": {
                "duplicate_excess_rows":
                    int(
                        duplicate_excess_rows
                    ),
            },

            "anomalies":
                list(
                    self.anomalies
                ),

            "validation_summary": {
                "severity_counts":
                    severity_counts,

                "total_anomaly_occurrences":
                    int(
                        self
                        .total_anomaly_occurrences
                    ),

                "stored_anomaly_samples":
                    int(
                        len(
                            self.anomalies
                        )
                    ),

                "anomaly_sample_limit":
                    int(
                        self
                        .anomaly_sample_limit
                    ),
            },
        }