"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/runners/moradi_de.py

Description:
    Dataset-level validation runner for Moradi differential-expression
    supplementary datasets S5-S10.

    Expected source datasets:

        S5  - differential exon expression, cis-related analysis
        S6  - differential exon expression, trans-related analysis
        S7  - differential intron expression, cis-related analysis
        S8  - differential intron expression, trans-related analysis
        S9  - differential isoform expression, cis-related analysis
        S10 - differential isoform expression, trans-related analysis

    Expected columns:

        Unnamed: 0
        baseMean
        log2FoldChange
        lfcSE
        stat
        pvalue
        padj

    Validation includes:

        - required schema
        - feature identifier format
        - baseMean validity
        - finite DE statistics
        - pvalue and padj ranges
        - missing-value statistics
        - exact duplicate statistics
        - source-level malformed feature detection
        - anomaly generation through AnomalyFactory

    Known source behavior:
        Intron identifiers resembling scientific notation, such as
        "INT1e+05", are preserved as source anomalies. They are not
        converted to canonical INT identifiers and are reported as
        WARNING rather than ERROR.

    Raw source files are never modified.

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
from pcatrfqtl.validation.validators.moradi_de import (
    MoradiDEValidator,
)


class MoradiDERunner:
    """Validate Moradi differential-expression datasets S5-S10."""

    DATASET_ID = "moradi_de"

    REQUIRED_COLUMNS = [
        "Unnamed: 0",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ]

    DATASETS = {
        "S5": {
            "filename":
                "Supplementary Table S5_DE_EP_cis.csv",
            "feature_type":
                "exon",
            "analysis_scope":
                "cis",
        },
        "S6": {
            "filename":
                "Supplementary Table S6_DE_EP_trans.csv",
            "feature_type":
                "exon",
            "analysis_scope":
                "trans",
        },
        "S7": {
            "filename":
                "Supplementary Table S7_DE_INT_cis.csv",
            "feature_type":
                "intron",
            "analysis_scope":
                "cis",
        },
        "S8": {
            "filename":
                "Supplementary Table S8_DE_INT_trans.csv",
            "feature_type":
                "intron",
            "analysis_scope":
                "trans",
        },
        "S9": {
            "filename":
                "Supplementary Table S9_DE_Isoform_cis.csv",
            "feature_type":
                "transcript",
            "analysis_scope":
                "cis",
        },
        "S10": {
            "filename":
                "Supplementary Table S10_DE_Isoform_trans_cancer_vs_control.csv",
            "feature_type":
                "transcript",
            "analysis_scope":
                "trans",
        },
    }

    ERROR_PATTERN = re.compile(
        r"^(?:\[[^\]]+\]\s*)?"
        r"(?P<field>.+?)"
        r"\[(?P<index>\d+)\]"
        r".*?:\s*"
        r"(?P<value>.+)$"
    )

    SCIENTIFIC_INTRON_PATTERN = re.compile(
        r"^INT\d+(?:\.\d+)?[eE][+-]?\d+$",
        re.IGNORECASE,
    )

    def __init__(
        self,
        source_directory: str | Path,
        *,
        anomaly_sample_limit: int = 1_000,
    ) -> None:
        """Initialize Moradi DE validation."""

        self.source_directory = Path(
            source_directory
        )

        self.validator = (
            MoradiDEValidator()
        )

        self.anomaly_sample_limit = int(
            anomaly_sample_limit
        )

        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies: list[
            dict[str, Any]
        ] = []

        self.total_anomaly_occurrences = 0

    # ==================================================================
    # Status helpers
    # ==================================================================

    @staticmethod
    def _status_from_severity_counts(
        severity_counts: dict[str, int],
    ) -> str:
        """Determine PASS, WARNING, or FAIL."""

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
    # Source anomaly helpers
    # ==================================================================

    @classmethod
    def _is_scientific_notation_like_intron(
        cls,
        value: Any,
    ) -> bool:
        """
        Detect source intron IDs resembling scientific notation.

        Example:
            INT1e+05

        These are preserved exactly and never automatically converted.
        """

        if value is None:
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

    # ==================================================================
    # Anomaly handling
    # ==================================================================

    def _add_anomaly(
        self,
        *,
        file_path: Path,
        row: int | None,
        column: str | None,
        value: Any,
        anomaly_type: str,
        category: str,
        severity: str,
        message: str,
        action: str,
    ) -> None:
        """Register one validation anomaly occurrence."""

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
                    file_path
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
        """Extract row index and source value from validator messages."""

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
        file_path: Path,
        dataframe: pd.DataFrame,
        errors: list[str],
        column: str,
        anomaly_type: str,
        category: str,
    ) -> None:
        """Convert field-validation errors to AnomalyFactory records."""

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

                # CSV line 1 contains the header.
                row = (
                    local_index
                    + 2
                )

                if (
                    0 <= local_index
                    < len(dataframe)
                    and column
                    in dataframe.columns
                ):
                    value = (
                        dataframe.iloc[
                            local_index
                        ][column]
                    )

            self._add_anomaly(
                file_path=file_path,
                row=row,
                column=column,
                value=value,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity="ERROR",
                message=error,
                action=(
                    "FLAG_FOR_REVIEW"
                ),
            )

    # ==================================================================
    # Statistics helpers
    # ==================================================================

    @staticmethod
    def _missing_value_statistics(
        dataframe: pd.DataFrame,
    ) -> dict[str, int]:
        """Return missing-value counts for every source column."""

        return {
            str(column):
                int(
                    dataframe[
                        column
                    ]
                    .isna()
                    .sum()
                )
            for column
            in dataframe.columns
        }

    @staticmethod
    def _duplicate_statistics(
        dataframe: pd.DataFrame,
    ) -> dict[str, int]:
        """
        Return exact duplicate statistics.

        duplicate_rows:
            Total rows participating in duplicate groups.

        duplicate_excess_rows:
            Rows exceeding the first occurrence of each duplicate group.
        """

        duplicate_mask = (
            dataframe.duplicated(
                keep=False
            )
        )

        duplicate_excess_mask = (
            dataframe.duplicated(
                keep="first"
            )
        )

        return {
            "duplicate_rows":
                int(
                    duplicate_mask.sum()
                ),
            "duplicate_excess_rows":
                int(
                    duplicate_excess_mask.sum()
                ),
        }

    @classmethod
    def _missing_required_columns(
        cls,
        dataframe: pd.DataFrame,
    ) -> list[str]:
        """Return required columns that are absent."""

        available = set(
            dataframe.columns
        )

        return [
            column
            for column
            in cls.REQUIRED_COLUMNS
            if column not in available
        ]

    # ==================================================================
    # Feature anomaly handling
    # ==================================================================

    def _classify_feature_errors(
        self,
        *,
        file_path: Path,
        dataframe: pd.DataFrame,
        feature_column: str,
        feature_type: str,
        feature_errors: list[str],
        severity_counts: dict[str, int],
    ) -> tuple[
        list[str],
        int,
    ]:
        """
        Separate genuine feature-validation errors from known source
        anomalies.

        Scientific-notation-like intron IDs are retained as WARNING.
        """

        if (
            feature_type
            != "intron"
        ):
            return (
                feature_errors,
                0,
            )

        retained_errors: list[str] = []

        source_feature_anomalies = 0

        for error in feature_errors:

            (
                local_index,
                _
            ) = (
                self
                ._extract_error_metadata(
                    error
                )
            )

            if (
                local_index is not None
                and 0 <= local_index
                < len(dataframe)
            ):
                value = (
                    dataframe.iloc[
                        local_index
                    ][feature_column]
                )

                if (
                    self
                    ._is_scientific_notation_like_intron(
                        value
                    )
                ):
                    source_feature_anomalies += 1

                    severity_counts[
                        "WARNING"
                    ] += 1

                    self._add_anomaly(
                        file_path=file_path,
                        row=(
                            local_index
                            + 2
                        ),
                        column=feature_column,
                        value=value,
                        anomaly_type=(
                            "MALFORMED_INTRON_FEATURE_ID"
                        ),
                        category="SOURCE",
                        severity="WARNING",
                        message=(
                            "Source intron identifier "
                            f"{value!r} resembles scientific "
                            "notation and does not match the "
                            "canonical INT<digits> format."
                        ),
                        action=(
                            "RETAIN_RAW_VALUE_AND_EXCLUDE_FROM_"
                            "CANONICAL_FEATURE_ID_ANALYSIS"
                        ),
                    )

                    continue

            retained_errors.append(
                error
            )

        return (
            retained_errors,
            source_feature_anomalies,
        )

    # ==================================================================
    # Single-file validation
    # ==================================================================

    def _validate_file(
        self,
        dataset_name: str,
        specification: dict[
            str,
            str,
        ],
    ) -> dict[str, Any]:
        """Validate one Moradi S5-S10 CSV file."""

        file_path = (
            self.source_directory
            / specification[
                "filename"
            ]
        )

        if not file_path.exists():

            self._add_anomaly(
                file_path=file_path,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_REQUIRED_FILE"
                ),
                category="SOURCE",
                severity="ERROR",
                message=(
                    f"Required Moradi DE file "
                    f"{file_path.name!r} "
                    "was not found."
                ),
                action=(
                    "STOP_DATASET_VALIDATION"
                ),
            )

            return {
                "dataset":
                    dataset_name,
                "file":
                    str(
                        file_path
                    ),
                "status":
                    "FAIL",
                "rows":
                    0,
                "columns":
                    0,
                "feature_type":
                    specification[
                        "feature_type"
                    ],
                "analysis_scope":
                    specification[
                        "analysis_scope"
                    ],
                "validation_error_counts":
                    {},
                "missing_values":
                    {},
                "duplicates":
                    {},
            }

        dataframe = pd.read_csv(
            file_path,
            low_memory=False,
        )

        dataframe.columns = [
            str(column).strip()
            for column
            in dataframe.columns
        ]

        missing_columns = (
            self
            ._missing_required_columns(
                dataframe
            )
        )

        severity_counts = {
            "WARNING": 0,
            "ERROR": 0,
        }

        if missing_columns:

            for column in missing_columns:

                severity_counts[
                    "ERROR"
                ] += 1

                self._add_anomaly(
                    file_path=file_path,
                    row=None,
                    column=column,
                    value=None,
                    anomaly_type=(
                        "MISSING_REQUIRED_COLUMN"
                    ),
                    category="SCHEMA",
                    severity="ERROR",
                    message=(
                        "Required Moradi DE column "
                        f"{column!r} was not found."
                    ),
                    action=(
                        "STOP_FILE_VALIDATION"
                    ),
                )

            return {
                "dataset":
                    dataset_name,
                "file":
                    str(
                        file_path
                    ),
                "status":
                    "FAIL",
                "rows":
                    int(
                        len(
                            dataframe
                        )
                    ),
                "columns":
                    int(
                        len(
                            dataframe.columns
                        )
                    ),
                "feature_type":
                    specification[
                        "feature_type"
                    ],
                "analysis_scope":
                    specification[
                        "analysis_scope"
                    ],
                "schema": {
                    "available_columns":
                        list(
                            dataframe.columns
                        ),
                    "required_columns":
                        list(
                            self.REQUIRED_COLUMNS
                        ),
                    "missing_required_columns":
                        missing_columns,
                },
                "validation_error_counts":
                    {},
                "missing_values":
                    (
                        self
                        ._missing_value_statistics(
                            dataframe
                        )
                    ),
                "duplicates":
                    (
                        self
                        ._duplicate_statistics(
                            dataframe
                        )
                    ),
                "validation_summary": {
                    "severity_counts":
                        severity_counts,
                },
            }

        feature_column = (
            "Unnamed: 0"
        )

        # --------------------------------------------------------------
        # Field validation
        # --------------------------------------------------------------

        feature_errors = (
            self.validator
            .validate_features(
                dataframe[
                    feature_column
                ].tolist(),
                specification[
                    "feature_type"
                ],
                feature_column,
                allow_missing=False,
            )
        )

        (
            feature_errors,
            source_feature_anomalies,
        ) = (
            self
            ._classify_feature_errors(
                file_path=file_path,
                dataframe=dataframe,
                feature_column=feature_column,
                feature_type=(
                    specification[
                        "feature_type"
                    ]
                ),
                feature_errors=(
                    feature_errors
                ),
                severity_counts=(
                    severity_counts
                ),
            )
        )

        base_mean_errors = (
            self.validator
            .validate_base_mean(
                dataframe[
                    "baseMean"
                ].tolist(),
                "baseMean",
                allow_missing=False,
            )
        )

        log2fc_errors = (
            self.validator
            .validate_numeric(
                dataframe[
                    "log2FoldChange"
                ].tolist(),
                "log2FoldChange",
                allow_missing=False,
            )
        )

        lfcse_errors = (
            self.validator
            .validate_numeric(
                dataframe[
                    "lfcSE"
                ].tolist(),
                "lfcSE",
                allow_missing=False,
            )
        )

        stat_errors = (
            self.validator
            .validate_numeric(
                dataframe[
                    "stat"
                ].tolist(),
                "stat",
                allow_missing=False,
            )
        )

        pvalue_errors = (
            self.validator
            .validate_probability(
                dataframe[
                    "pvalue"
                ].tolist(),
                "pvalue",
                allow_missing=False,
            )
        )

        padj_errors = (
            self.validator
            .validate_probability(
                dataframe[
                    "padj"
                ].tolist(),
                "padj",
                allow_missing=False,
            )
        )

        validation_error_counts = {
            "features":
                len(
                    feature_errors
                ),
            "source_feature_anomalies":
                int(
                    source_feature_anomalies
                ),
            "base_mean":
                len(
                    base_mean_errors
                ),
            "log2_fold_change":
                len(
                    log2fc_errors
                ),
            "lfc_se":
                len(
                    lfcse_errors
                ),
            "stat":
                len(
                    stat_errors
                ),
            "p_values":
                len(
                    pvalue_errors
                ),
            "adjusted_p_values":
                len(
                    padj_errors
                ),
        }

        error_groups = [
            (
                feature_errors,
                feature_column,
                "INVALID_FEATURE_ID",
                "FORMAT",
            ),
            (
                base_mean_errors,
                "baseMean",
                "INVALID_BASE_MEAN",
                "RANGE",
            ),
            (
                log2fc_errors,
                "log2FoldChange",
                "INVALID_LOG2_FOLD_CHANGE",
                "NUMERIC",
            ),
            (
                lfcse_errors,
                "lfcSE",
                "INVALID_LFC_SE",
                "NUMERIC",
            ),
            (
                stat_errors,
                "stat",
                "INVALID_STATISTIC",
                "NUMERIC",
            ),
            (
                pvalue_errors,
                "pvalue",
                "INVALID_P_VALUE",
                "RANGE",
            ),
            (
                padj_errors,
                "padj",
                "INVALID_ADJUSTED_P_VALUE",
                "RANGE",
            ),
        ]

        for (
            errors,
            column,
            anomaly_type,
            category,
        ) in error_groups:

            if not errors:
                continue

            severity_counts[
                "ERROR"
            ] += len(
                errors
            )

            self._add_validation_errors(
                file_path=file_path,
                dataframe=dataframe,
                errors=errors,
                column=column,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
            )

        # --------------------------------------------------------------
        # Duplicate statistics
        # --------------------------------------------------------------

        duplicates = (
            self
            ._duplicate_statistics(
                dataframe
            )
        )

        if (
            duplicates[
                "duplicate_excess_rows"
            ]
            > 0
        ):

            severity_counts[
                "WARNING"
            ] += 1

            self._add_anomaly(
                file_path=file_path,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "EXACT_DUPLICATE_DE_RECORD"
                ),
                category="DUPLICATE",
                severity="WARNING",
                message=(
                    f"{duplicates['duplicate_excess_rows']:,} "
                    "exact duplicate excess rows were "
                    f"detected in {dataset_name}."
                ),
                action=(
                    "RETAIN_RAW_AND_REVIEW_DURING_"
                    "NORMALIZATION"
                ),
            )

        status = (
            self
            ._status_from_severity_counts(
                severity_counts
            )
        )

        return {
            "dataset":
                dataset_name,

            "file":
                str(
                    file_path
                ),

            "status":
                status,

            "rows":
                int(
                    len(
                        dataframe
                    )
                ),

            "columns":
                int(
                    len(
                        dataframe.columns
                    )
                ),

            "feature_type":
                specification[
                    "feature_type"
                ],

            "analysis_scope":
                specification[
                    "analysis_scope"
                ],

            "schema": {
                "available_columns":
                    list(
                        dataframe.columns
                    ),

                "required_columns":
                    list(
                        self.REQUIRED_COLUMNS
                    ),

                "missing_required_columns":
                    [],
            },

            "validation_error_counts":
                validation_error_counts,

            "missing_values":
                (
                    self
                    ._missing_value_statistics(
                        dataframe
                    )
                ),

            "duplicates":
                duplicates,

            "validation_summary": {
                "severity_counts":
                    severity_counts,
            },
        }

    # ==================================================================
    # Full S5-S10 validation
    # ==================================================================

    def validate(
        self,
    ) -> dict[str, Any]:
        """Validate all Moradi DE datasets S5-S10."""

        if not self.source_directory.exists():
            raise FileNotFoundError(
                "Moradi source directory not found: "
                f"{self.source_directory}"
            )

        if not self.source_directory.is_dir():
            raise NotADirectoryError(
                "Moradi source path is not a directory: "
                f"{self.source_directory}"
            )

        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies = []

        self.total_anomaly_occurrences = 0

        files: dict[
            str,
            dict[str, Any],
        ] = {}

        overall_error_count = 0
        overall_warning_count = 0

        total_rows = 0

        for (
            dataset_name,
            specification,
        ) in (
            self.DATASETS
            .items()
        ):

            result = (
                self
                ._validate_file(
                    dataset_name,
                    specification,
                )
            )

            files[
                dataset_name
            ] = result

            total_rows += int(
                result.get(
                    "rows",
                    0,
                )
            )

            summary = result.get(
                "validation_summary",
                {},
            )

            severity = summary.get(
                "severity_counts",
                {},
            )

            overall_error_count += int(
                severity.get(
                    "ERROR",
                    0,
                )
            )

            overall_warning_count += int(
                severity.get(
                    "WARNING",
                    0,
                )
            )

            if (
                result.get(
                    "status"
                )
                == "FAIL"
                and not severity
            ):
                overall_error_count += 1

        severity_counts = {
            "WARNING":
                overall_warning_count,
            "ERROR":
                overall_error_count,
        }

        status = (
            self
            ._status_from_severity_counts(
                severity_counts
            )
        )

        return {
            "dataset_id":
                self.DATASET_ID,

            "source_directory":
                str(
                    self.source_directory
                ),

            "status":
                status,

            "files_validated":
                len(
                    files
                ),

            "total_rows":
                int(
                    total_rows
                ),

            "files":
                files,

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