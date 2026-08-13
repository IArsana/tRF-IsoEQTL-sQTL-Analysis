"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/runners/moradi_workbook.py

Description:
    Workbook-level validation runner for the Moradi supplementary
    regulatory-QTL datasets.

    The runner currently validates core Moradi datasets:

        Supplementary Table S1
            - cis intron-retention sQTL
            - cis cassette-exon sQTL
            - cis iso-eQTL

        Supplementary Table S2
            - GWAS-linked cis intron sQTL
            - GWAS-linked cis exon sQTL
            - GWAS-linked cis iso-eQTL

    Validation findings are converted into standardized anomalies using
    AnomalyFactory. Occurrence identifiers follow the globally sequential
    format:

        ANOM-000001
        ANOM-000002
        ...

    Status semantics:

        PASS
            No anomalies.

        WARNING
            Warning-level anomalies exist, but no error-level anomalies.

        FAIL
            At least one error-level anomaly exists.

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

from pcatrfqtl.validation.anomaly import AnomalyFactory
from pcatrfqtl.validation.validators.moradi import (
    MoradiQTLValidator,
)


class MoradiWorkbookValidator:
    """Validate core Moradi supplementary QTL workbooks."""

    DATASET_ID = "moradi"

    S1_FILE = "Supplementary Table S1.xlsx"
    S2_FILE = "Supplementary Table S2.xlsx"
    S3_FILE = "Supplementary Table S3.xlsx"
    S4_FILE = "Supplementary Table S4.xlsx"

    S1_SHEETS = {
        "Cis-intron-retention-sQTL-0.05": {
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Cis-cassette-exon-sQTL-0.05": {
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Cis-iso-eQTL-0.05": {
            "feature_column": "splicing_event",
            "feature_validator": "transcript",
        },
    }

    S2_SHEETS = {
        "Cis-intron-sQTL-GWAS-LD_0.5": {
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Cis-exon-sQTL-GWAS-LD_0.5": {
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Cis-iso-eQTl-GWAS-LD_0.5": {
            "feature_column": "mRNA_isoform",
            "feature_validator": "transcript",
        },
    }

    S3_SHEETS = {
        "Trans-intron-sQTL": {
            "snp_column": "SNP",
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Trans-exon-sQTl": {
            # Source schema uses "T" instead of "SNP".
            "snp_column": "T",
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Trans-iso-eQTL": {
            "snp_column": "SNP",
            "feature_column": "splicing_event",
            "feature_validator": "transcript",
        },
    }

    S4_SHEETS = {
        "Trans-intron-sQTL-LD_0.5": {
            "feature_column": "mRNA_isoform",
            "feature_validator": "transcript",
        },
        "Trans-exon-sQTL-LD_0.5": {
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
        "Trans-iso-eQTL-LD_0.5": {
            "feature_column": "splicing_event",
            "feature_validator": "splicing",
        },
    }

    S1_REQUIRED_COLUMNS = [
        "SNP",
        "R",
        "A",
        "SNP_pos",
        "splicing_event",
        "Stat",
        "p_value",
        "FDR",
        "beta",
        "MAF",
        "AvgCall",
        "Rsq",
    ]

    S2_COMMON_REQUIRED_COLUMNS = [
        "sQTL_SNP",
        "sQTL_SNP-pos",
        "tag_SNP",
        "tag_SNP_pos",
        "LD",
        "gwas_cancer",
    ]

    ERROR_PATTERN = re.compile(
        r"^(?:\[[^\]]+\]\s*)?"
        r"(?P<field>.+?)"
        r"\[(?P<index>\d+)\]"
        r".*?:\s*"
        r"(?P<value>.+?)"
        r"\.?$"
    )

    def __init__(
        self,
        moradi_dir: str | Path,
    ) -> None:
        """Initialize Moradi workbook validation."""

        self.moradi_dir = Path(
            moradi_dir
        )

        self.validator = (
            MoradiQTLValidator()
        )

        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies: list[
            dict[str, Any]
        ] = []

    # ==================================================================
    # Dataframe helpers
    # ==================================================================

    @staticmethod
    def _clean_dataframe(
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Remove completely empty rows/columns and clean headers."""

        cleaned = dataframe.dropna(
            axis=0,
            how="all",
        )

        cleaned = cleaned.dropna(
            axis=1,
            how="all",
        )

        cleaned = cleaned.copy()

        cleaned.columns = [
            str(column).strip()
            for column
            in cleaned.columns
        ]

        return cleaned.reset_index(
            drop=True
        )

    @staticmethod
    def _missing_value_statistics(
        dataframe: pd.DataFrame,
    ) -> dict[str, int]:
        """Return missing-value counts per column."""

        return {
            str(column): int(
                dataframe[column]
                .isna()
                .sum()
            )
            for column
            in dataframe.columns
        }

    @staticmethod
    def _unique_statistics(
        dataframe: pd.DataFrame,
        columns: list[str],
    ) -> dict[str, int]:
        """Return unique non-null counts for selected columns."""

        result: dict[str, int] = {}

        for column in columns:
            if column not in dataframe.columns:
                continue

            series = (
                dataframe[column]
                .dropna()
            )

            if (
                pd.api.types
                .is_object_dtype(
                    series.dtype
                )
                or
                pd.api.types
                .is_string_dtype(
                    series.dtype
                )
            ):
                series = (
                    series
                    .astype(str)
                    .str.strip()
                )

            result[column] = int(
                series.nunique()
            )

        return result

    @staticmethod
    def _duplicate_statistics(
        dataframe: pd.DataFrame,
    ) -> dict[str, int]:
        """Return exact-row duplicate statistics."""

        duplicate_rows = int(
            dataframe
            .duplicated(
                keep=False
            )
            .sum()
        )

        duplicate_excess_rows = int(
            dataframe
            .duplicated(
                keep="first"
            )
            .sum()
        )

        return {
            "duplicate_rows":
                duplicate_rows,
            "duplicate_excess_rows":
                duplicate_excess_rows,
        }

    @staticmethod
    def _missing_required_columns(
        dataframe: pd.DataFrame,
        required_columns: list[str],
    ) -> list[str]:
        """Return required columns absent from the dataframe."""

        available = set(
            dataframe.columns
        )

        return [
            column
            for column
            in required_columns
            if column not in available
        ]

    @staticmethod
    def _validate_required_columns(
        dataframe: pd.DataFrame,
        required_columns: list[str],
    ) -> list[str]:
        """
        Backward-compatible alias for required-column validation.

        Returns a list of required columns that are absent from the
        dataframe.
        """

        return (
            MoradiWorkbookValidator
            ._missing_required_columns(
                dataframe,
                required_columns,
            )
        )

    # ==================================================================
    # Status helpers
    # ==================================================================

    @staticmethod
    def _status_from_anomalies(
        anomalies: list[
            dict[str, Any]
        ],
    ) -> str:
        """Determine status from anomaly severities."""

        if not anomalies:
            return "PASS"

        severities = {
            str(
                anomaly.get(
                    "severity",
                    "",
                )
            ).upper()
            for anomaly
            in anomalies
        }

        if "ERROR" in severities:
            return "FAIL"

        return "WARNING"

    @staticmethod
    def _aggregate_status(
        statuses: list[str],
    ) -> str:
        """Aggregate child validation statuses."""

        normalized = {
            status.upper()
            for status
            in statuses
        }

        if "FAIL" in normalized:
            return "FAIL"

        if "WARNING" in normalized:
            return "WARNING"

        return "PASS"

    # ==================================================================
    # Anomaly helpers
    # ==================================================================

    def _add_anomaly(
        self,
        *,
        file: str | Path,
        sheet: str | None,
        row: int | None,
        column: str | None,
        value: Any,
        anomaly_type: str,
        category: str,
        severity: str,
        message: str,
        action: str,
    ) -> None:
        """Create and register one standardized anomaly."""

        anomaly = (
            self.anomaly_factory
            .create(
                type=anomaly_type,
                category=category,
                severity=severity,
                dataset=self.DATASET_ID,
                file=str(file),
                sheet=sheet,
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
        Extract dataframe index and source value from a validator message.

        Expected validator format resembles:

            [MORADI-EVENT-001]
            splicing_event[38278]
            invalid splicing-event identifier: 'INT1e+05'.

        Returns
        -------
        tuple
            (dataframe_index, value)
        """

        match = (
            cls.ERROR_PATTERN.search(
                error
            )
        )

        if not match:
            return None, None

        index = int(
            match.group(
                "index"
            )
        )

        raw_value = (
            match.group(
                "value"
            )
            .strip()
            .rstrip(".")
        )

        if (
            len(raw_value) >= 2
            and raw_value[0]
            == raw_value[-1]
            and raw_value[0]
            in {"'", '"'}
        ):
            raw_value = (
                raw_value[1:-1]
            )

        return index, raw_value

    def _add_validation_errors(
        self,
        *,
        errors: list[str],
        workbook: Path,
        sheet: str,
        column: str,
        anomaly_type: str,
        category: str,
        severity: str = "ERROR",
        action: str = "FLAG_FOR_REVIEW",
        source_dataframe: (
            pd.DataFrame | None
        ) = None,
    ) -> None:
        """
        Convert field-level validation errors into Anomaly objects.

        If validator messages contain the dataframe index, row/value
        metadata are recovered where possible.

        Excel row number is calculated assuming:
            header=0
            dataframe index 0 == Excel row 2
        """

        for error in errors:
            (
                dataframe_index,
                parsed_value,
            ) = (
                self
                ._extract_error_metadata(
                    error
                )
            )

            row: int | None = None
            value: Any = parsed_value

            if dataframe_index is not None:
                # header row = Excel row 1,
                # first data record = Excel row 2.
                row = (
                    dataframe_index
                    + 2
                )

                if (
                    source_dataframe
                    is not None
                    and column
                    in source_dataframe.columns
                    and 0
                    <= dataframe_index
                    < len(
                        source_dataframe
                    )
                ):
                    value = (
                        source_dataframe.iloc[
                            dataframe_index
                        ][column]
                    )

            self._add_anomaly(
                file=workbook,
                sheet=sheet,
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
    # Source anomaly policies
    # ==================================================================

    @staticmethod
    def _is_known_malformed_intron_event(
        value: Any,
    ) -> bool:
        """
        Detect malformed intron identifiers that appear to contain
        scientific notation.

        Example:
            INT1e+05

        Such values remain invalid identifiers but are treated as
        source-data warnings rather than structural validation failures.
        """

        if value is None:
            return False

        value_str = str(
            value
        ).strip()

        return bool(
            re.fullmatch(
                r"INT\d+(?:\.\d+)?[eE][+-]?\d+",
                value_str,
            )
        )

    # ==================================================================
    # Reading
    # ==================================================================

    @staticmethod
    def _read_sheet(
        workbook: Path,
        sheet_name: str,
    ) -> pd.DataFrame:
        """Read and clean one Moradi Excel worksheet."""

        dataframe = pd.read_excel(
            workbook,
            sheet_name=sheet_name,
            header=0,
        )

        return (
            MoradiWorkbookValidator
            ._clean_dataframe(
                dataframe
            )
        )

    # ==================================================================
    # S1 validation
    # ==================================================================

    def _validate_s1_sheet(
        self,
        dataframe: pd.DataFrame,
        sheet_name: str,
        workbook: Path,
    ) -> dict[str, Any]:
        """Validate one Supplementary Table S1 worksheet."""

        anomalies_before = len(
            self.anomalies
        )

        missing_columns = (
            self
            ._missing_required_columns(
                dataframe,
                self.S1_REQUIRED_COLUMNS,
            )
        )

        if missing_columns:
            self._add_anomaly(
                file=workbook,
                sheet=sheet_name,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_REQUIRED_COLUMN"
                ),
                category="MISSING",
                severity="ERROR",
                message=(
                    "Missing required columns: "
                    + ", ".join(
                        missing_columns
                    )
                ),
                action=(
                    "STOP_SHEET_VALIDATION"
                ),
            )

            sheet_anomalies = (
                self.anomalies[
                    anomalies_before:
                ]
            )

            return {
                "status": "FAIL",
                "rows": int(
                    len(dataframe)
                ),
                "columns": int(
                    len(
                        dataframe.columns
                    )
                ),
                "missing_required_columns":
                    missing_columns,
                "missing_values": (
                    self
                    ._missing_value_statistics(
                        dataframe
                    )
                ),
                "unique_values": {},
                "duplicates": (
                    self
                    ._duplicate_statistics(
                        dataframe
                    )
                ),
                "validation_error_counts":
                    {},
                "anomalies":
                    sheet_anomalies,
            }

        errors: dict[
            str,
            list[str],
        ] = {}

        errors["snp_ids"] = (
            self.validator
            .validate_snp_ids(
                dataframe[
                    "SNP"
                ].tolist(),
                "SNP",
            )
        )

        errors["coordinates"] = (
            self.validator
            .validate_coordinates(
                dataframe[
                    "SNP_pos"
                ].tolist(),
                "SNP_pos",
            )
        )

        errors[
            "reference_alleles"
        ] = (
            self.validator
            .validate_alleles(
                dataframe[
                    "R"
                ].tolist(),
                "R",
            )
        )

        errors[
            "alternate_alleles"
        ] = (
            self.validator
            .validate_alleles(
                dataframe[
                    "A"
                ].tolist(),
                "A",
            )
        )

        feature_validator = (
            self.S1_SHEETS[
                sheet_name
            ][
                "feature_validator"
            ]
        )

        if (
            feature_validator
            == "transcript"
        ):
            errors["features"] = (
                self.validator
                .validate_transcripts(
                    dataframe[
                        "splicing_event"
                    ].tolist(),
                    "splicing_event",
                )
            )
        else:
            errors["features"] = (
                self.validator
                .validate_splicing_events(
                    dataframe[
                        "splicing_event"
                    ].tolist(),
                    "splicing_event",
                )
            )

        errors["stat"] = (
            self.validator
            .validate_numeric(
                dataframe[
                    "Stat"
                ].tolist(),
                "Stat",
            )
        )

        errors["p_values"] = (
            self.validator
            .validate_p_values(
                dataframe[
                    "p_value"
                ].tolist(),
                "p_value",
            )
        )

        errors["fdr"] = (
            self.validator
            .validate_fdr(
                dataframe[
                    "FDR"
                ].tolist(),
                "FDR",
            )
        )

        errors["beta"] = (
            self.validator
            .validate_numeric(
                dataframe[
                    "beta"
                ].tolist(),
                "beta",
            )
        )

        errors["maf"] = (
            self.validator
            .validate_maf(
                dataframe[
                    "MAF"
                ].tolist(),
                "MAF",
            )
        )

        errors["call_rate"] = (
            self.validator
            .validate_call_rate(
                dataframe[
                    "AvgCall"
                ].tolist(),
                "AvgCall",
            )
        )

        errors["rsq"] = (
            self.validator
            .validate_rsq(
                dataframe[
                    "Rsq"
                ].tolist(),
                "Rsq",
            )
        )

        # --------------------------------------------------------------
        # Standard field errors
        # --------------------------------------------------------------

        standard_fields = [
            (
                "snp_ids",
                "SNP",
                "INVALID_SNP_ID",
                "FORMAT",
            ),
            (
                "coordinates",
                "SNP_pos",
                "INVALID_COORDINATE",
                "COORDINATE",
            ),
            (
                "reference_alleles",
                "R",
                "INVALID_REFERENCE_ALLELE",
                "ALLELE",
            ),
            (
                "alternate_alleles",
                "A",
                "INVALID_ALTERNATE_ALLELE",
                "ALLELE",
            ),
            (
                "stat",
                "Stat",
                "INVALID_STATISTIC",
                "RANGE",
            ),
            (
                "p_values",
                "p_value",
                "INVALID_P_VALUE",
                "RANGE",
            ),
            (
                "fdr",
                "FDR",
                "INVALID_FDR",
                "RANGE",
            ),
            (
                "beta",
                "beta",
                "INVALID_BETA",
                "RANGE",
            ),
            (
                "maf",
                "MAF",
                "INVALID_MAF",
                "RANGE",
            ),
            (
                "call_rate",
                "AvgCall",
                "INVALID_CALL_RATE",
                "RANGE",
            ),
            (
                "rsq",
                "Rsq",
                "INVALID_RSQ",
                "RANGE",
            ),
        ]

        for (
            error_key,
            column,
            anomaly_type,
            category,
        ) in standard_fields:
            self._add_validation_errors(
                errors=errors[
                    error_key
                ],
                workbook=workbook,
                sheet=sheet_name,
                column=column,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity="ERROR",
                source_dataframe=(
                    dataframe
                ),
            )

        # --------------------------------------------------------------
        # Feature anomalies
        #
        # Scientific-notation-like intron identifiers such as INT1e+05
        # are invalid IDs but are classified as source warnings.
        # Other feature errors remain ERROR-level findings.
        # --------------------------------------------------------------

        for error in errors[
            "features"
        ]:
            (
                dataframe_index,
                parsed_value,
            ) = (
                self
                ._extract_error_metadata(
                    error
                )
            )

            source_value: Any = (
                parsed_value
            )

            row: int | None = None

            if (
                dataframe_index
                is not None
            ):
                row = (
                    dataframe_index
                    + 2
                )

                if (
                    0
                    <= dataframe_index
                    < len(dataframe)
                ):
                    source_value = (
                        dataframe.iloc[
                            dataframe_index
                        ][
                            "splicing_event"
                        ]
                    )

            if (
                feature_validator
                == "splicing"
                and self
                ._is_known_malformed_intron_event(
                    source_value
                )
            ):
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=row,
                    column=(
                        "splicing_event"
                    ),
                    value=source_value,
                    anomaly_type=(
                        "MALFORMED_SPLICING_EVENT_ID"
                    ),
                    category="SOURCE",
                    severity="WARNING",
                    message=error,
                    action=(
                        "FLAG_FOR_REVIEW_AND_EXCLUDE_FROM_"
                        "EVENT_ID_DEPENDENT_ANALYSIS"
                    ),
                )

            else:
                anomaly_type = (
                    "INVALID_TRANSCRIPT_ID"
                    if feature_validator
                    == "transcript"
                    else
                    "INVALID_SPLICING_EVENT_ID"
                )

                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=row,
                    column=(
                        "splicing_event"
                    ),
                    value=source_value,
                    anomaly_type=(
                        anomaly_type
                    ),
                    category="FORMAT",
                    severity="ERROR",
                    message=error,
                    action=(
                        "FLAG_FOR_REVIEW"
                    ),
                )

        error_counts = {
            name: len(values)
            for name, values
            in errors.items()
        }

        sheet_anomalies = (
            self.anomalies[
                anomalies_before:
            ]
        )

        return {
            "status": (
                self
                ._status_from_anomalies(
                    sheet_anomalies
                )
            ),
            "rows": int(
                len(dataframe)
            ),
            "columns": int(
                len(
                    dataframe.columns
                )
            ),
            "missing_required_columns":
                [],
            "missing_values": (
                self
                ._missing_value_statistics(
                    dataframe
                )
            ),
            "unique_values": (
                self
                ._unique_statistics(
                    dataframe,
                    [
                        "SNP",
                        "SNP_pos",
                        "splicing_event",
                    ],
                )
            ),
            "duplicates": (
                self
                ._duplicate_statistics(
                    dataframe
                )
            ),
            "validation_error_counts":
                error_counts,
            "anomalies":
                sheet_anomalies,
        }

    # ==================================================================
    # S2 validation
    # ==================================================================

    def _validate_s2_sheet(
        self,
        dataframe: pd.DataFrame,
        sheet_name: str,
        workbook: Path,
    ) -> dict[str, Any]:
        """
        Validate one Supplementary Table S2 worksheet.

        Source semantic mapping
        -----------------------
        sQTL_SNP
            Genomic coordinate.

        sQTL_SNP-pos
            dbSNP rs identifier.

        tag_SNP
            GWAS genomic coordinate.

        tag_SNP_pos
            GWAS dbSNP identifier.

        The raw source headers are preserved.
        """

        anomalies_before = len(
            self.anomalies
        )

        feature_column = (
            self.S2_SHEETS[
                sheet_name
            ][
                "feature_column"
            ]
        )

        required_columns = [
            *self
            .S2_COMMON_REQUIRED_COLUMNS,
            feature_column,
        ]

        missing_columns = (
            self
            ._missing_required_columns(
                dataframe,
                required_columns,
            )
        )

        if missing_columns:
            self._add_anomaly(
                file=workbook,
                sheet=sheet_name,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_REQUIRED_COLUMN"
                ),
                category="MISSING",
                severity="ERROR",
                message=(
                    "Missing required columns: "
                    + ", ".join(
                        missing_columns
                    )
                ),
                action=(
                    "STOP_SHEET_VALIDATION"
                ),
            )

            sheet_anomalies = (
                self.anomalies[
                    anomalies_before:
                ]
            )

            return {
                "status": "FAIL",
                "rows": int(
                    len(dataframe)
                ),
                "columns": int(
                    len(
                        dataframe.columns
                    )
                ),
                "missing_required_columns":
                    missing_columns,
                "missing_values": (
                    self
                    ._missing_value_statistics(
                        dataframe
                    )
                ),
                "unique_values": {},
                "duplicates": (
                    self
                    ._duplicate_statistics(
                        dataframe
                    )
                ),
                "validation_error_counts":
                    {},
                "anomalies":
                    sheet_anomalies,
            }

        errors: dict[
            str,
            list[str],
        ] = {}

        errors[
            "qtl_coordinates"
        ] = (
            self.validator
            .validate_coordinates(
                dataframe[
                    "sQTL_SNP"
                ].tolist(),
                "sQTL_SNP",
            )
        )

        errors[
            "qtl_snp_ids"
        ] = (
            self.validator
            .validate_snp_ids(
                dataframe[
                    "sQTL_SNP-pos"
                ].tolist(),
                "sQTL_SNP-pos",
            )
        )

        errors[
            "gwas_coordinates"
        ] = (
            self.validator
            .validate_coordinates(
                dataframe[
                    "tag_SNP"
                ].tolist(),
                "tag_SNP",
            )
        )

        errors[
            "gwas_tag_snp_ids"
        ] = (
            self.validator
            .validate_composite_snp_ids(
                dataframe[
                    "tag_SNP_pos"
                ].tolist(),
                "tag_SNP_pos",
            )
        )

        feature_validator = (
            self.S2_SHEETS[
                sheet_name
            ][
                "feature_validator"
            ]
        )

        if (
            feature_validator
            == "transcript"
        ):
            errors["features"] = (
                self.validator
                .validate_transcripts(
                    dataframe[
                        feature_column
                    ].tolist(),
                    feature_column,
                )
            )
        else:
            errors["features"] = (
                self.validator
                .validate_splicing_events(
                    dataframe[
                        feature_column
                    ].tolist(),
                    feature_column,
                )
            )

        errors["ld"] = (
            self.validator
            .validate_ld(
                dataframe[
                    "LD"
                ].tolist(),
                "LD",
            )
        )

        field_definitions = [
            (
                "qtl_coordinates",
                "sQTL_SNP",
                "INVALID_QTL_COORDINATE",
                "COORDINATE",
            ),
            (
                "qtl_snp_ids",
                "sQTL_SNP-pos",
                "INVALID_QTL_SNP_ID",
                "FORMAT",
            ),
            (
                "gwas_coordinates",
                "tag_SNP",
                "INVALID_GWAS_COORDINATE",
                "COORDINATE",
            ),
            (
                "gwas_tag_snp_ids",
                "tag_SNP_pos",
                "INVALID_GWAS_SNP_ID",
                "FORMAT",
            ),
            (
                "ld",
                "LD",
                "INVALID_LD_R2",
                "RANGE",
            ),
        ]

        for (
            error_key,
            column,
            anomaly_type,
            category,
        ) in field_definitions:
            self._add_validation_errors(
                errors=errors[
                    error_key
                ],
                workbook=workbook,
                sheet=sheet_name,
                column=column,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity="ERROR",
                source_dataframe=(
                    dataframe
                ),
            )

        feature_anomaly_type = (
            "INVALID_TRANSCRIPT_ID"
            if feature_validator
            == "transcript"
            else
            "INVALID_SPLICING_EVENT_ID"
        )

        self._add_validation_errors(
            errors=errors[
                "features"
            ],
            workbook=workbook,
            sheet=sheet_name,
            column=feature_column,
            anomaly_type=(
                feature_anomaly_type
            ),
            category="FORMAT",
            severity="ERROR",
            source_dataframe=(
                dataframe
            ),
        )

        error_counts = {
            name: len(values)
            for name, values
            in errors.items()
        }

        sheet_anomalies = (
            self.anomalies[
                anomalies_before:
            ]
        )

        return {
            "status": (
                self
                ._status_from_anomalies(
                    sheet_anomalies
                )
            ),
            "rows": int(
                len(dataframe)
            ),
            "columns": int(
                len(
                    dataframe.columns
                )
            ),
            "missing_required_columns":
                [],
            "missing_values": (
                self
                ._missing_value_statistics(
                    dataframe
                )
            ),
            "unique_values": (
                self
                ._unique_statistics(
                    dataframe,
                    [
                        "sQTL_SNP",
                        "sQTL_SNP-pos",
                        feature_column,
                        "tag_SNP",
                        "tag_SNP_pos",
                        "gwas_cancer",
                    ],
                )
            ),
            "duplicates": (
                self
                ._duplicate_statistics(
                    dataframe
                )
            ),
            "validation_error_counts":
                error_counts,
            "anomalies":
                sheet_anomalies,
            "source_schema_notes": [
                (
                    "sQTL_SNP contains genomic coordinates "
                    "while sQTL_SNP-pos contains dbSNP identifiers."
                ),
                (
                    "tag_SNP contains genomic coordinates "
                    "while tag_SNP_pos contains dbSNP identifiers."
                ),
                (
                    "Raw source column names are preserved; "
                    "semantic renaming belongs to normalization."
                ),
            ],
        }

    # ==================================================================
    # s3 validation
    # ==================================================================

    def _validate_s3_sheet(
        self,
        dataframe: pd.DataFrame,
        sheet_name: str,
        workbook: Path,
    ) -> dict[str, Any]:
        """
        Validate one Moradi Supplementary Table S3 worksheet.

        S3 contains trans-QTL associations.

        Source-specific behavior:
            - Most worksheets use SNP as the dbSNP column.
            - Trans-exon-sQTl-0.05 uses T instead of SNP.
            - Allele columns support SNVs, indels, and comma-separated
            multi-allelic representations.
        """

        anomalies_before = len(
            self.anomalies
        )

        schema = self.S3_SHEETS[
            sheet_name
        ]

        snp_column = schema[
            "snp_column"
        ]

        feature_column = schema[
            "feature_column"
        ]

        feature_validator = schema[
            "feature_validator"
        ]

        required_columns = [
            snp_column,
            "R",
            "A",
            "SNP_pos",
            feature_column,
            "Stat",
            "p_value",
            "FDR",
            "beta",
            "MAF",
            "AvgCall",
            "Rsq",
        ]

        missing_columns = (
            self._missing_required_columns(
                dataframe,
                required_columns,
            )
        )

        if missing_columns:
            self._add_anomaly(
                file=workbook,
                sheet=sheet_name,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_REQUIRED_COLUMN"
                ),
                category="MISSING",
                severity="ERROR",
                message=(
                    "Missing required columns: "
                    + ", ".join(
                        missing_columns
                    )
                ),
                action=(
                    "STOP_SHEET_VALIDATION"
                ),
            )

            sheet_anomalies = (
                self.anomalies[
                    anomalies_before:
                ]
            )

            return {
                "status": "FAIL",
                "rows": int(
                    len(dataframe)
                ),
                "columns": int(
                    len(
                        dataframe.columns
                    )
                ),
                "missing_required_columns":
                    missing_columns,
                "missing_values":
                    self
                    ._missing_value_statistics(
                        dataframe
                    ),
                "unique_values": {},
                "duplicates":
                    self
                    ._duplicate_statistics(
                        dataframe
                    ),
                "validation_error_counts":
                    {},
                "anomalies":
                    sheet_anomalies,
            }

        errors: dict[
            str,
            list[str],
        ] = {}

        errors["snp_ids"] = (
            self.validator
            .validate_snp_ids(
                dataframe[
                    snp_column
                ].tolist(),
                snp_column,
            )
        )

        errors["coordinates"] = (
            self.validator
            .validate_coordinates(
                dataframe[
                    "SNP_pos"
                ].tolist(),
                "SNP_pos",
            )
        )

        errors[
            "reference_alleles"
        ] = (
            self.validator
            .validate_alleles(
                dataframe[
                    "R"
                ].tolist(),
                "R",
            )
        )

        errors[
            "alternate_alleles"
        ] = (
            self.validator
            .validate_alleles(
                dataframe[
                    "A"
                ].tolist(),
                "A",
            )
        )

        if (
            feature_validator
            == "transcript"
        ):
            errors["features"] = (
                self.validator
                .validate_transcripts(
                    dataframe[
                        feature_column
                    ].tolist(),
                    feature_column,
                )
            )
        else:
            errors["features"] = (
                self.validator
                .validate_splicing_events(
                    dataframe[
                        feature_column
                    ].tolist(),
                    feature_column,
                )
            )

        errors["stat"] = (
            self.validator
            .validate_numeric(
                dataframe[
                    "Stat"
                ].tolist(),
                "Stat",
            )
        )

        errors["p_values"] = (
            self.validator
            .validate_p_values(
                dataframe[
                    "p_value"
                ].tolist(),
                "p_value",
            )
        )

        errors["fdr"] = (
            self.validator
            .validate_fdr(
                dataframe[
                    "FDR"
                ].tolist(),
                "FDR",
            )
        )

        errors["beta"] = (
            self.validator
            .validate_numeric(
                dataframe[
                    "beta"
                ].tolist(),
                "beta",
            )
        )

        errors["maf"] = (
            self.validator
            .validate_maf(
                dataframe[
                    "MAF"
                ].tolist(),
                "MAF",
            )
        )

        errors["call_rate"] = (
            self.validator
            .validate_call_rate(
                dataframe[
                    "AvgCall"
                ].tolist(),
                "AvgCall",
            )
        )

        errors["rsq"] = (
            self.validator
            .validate_rsq(
                dataframe[
                    "Rsq"
                ].tolist(),
                "Rsq",
            )
        )

        field_definitions = [
            (
                "snp_ids",
                snp_column,
                "INVALID_SNP_ID",
                "FORMAT",
            ),
            (
                "coordinates",
                "SNP_pos",
                "INVALID_COORDINATE",
                "COORDINATE",
            ),
            (
                "reference_alleles",
                "R",
                "INVALID_REFERENCE_ALLELE",
                "ALLELE",
            ),
            (
                "alternate_alleles",
                "A",
                "INVALID_ALTERNATE_ALLELE",
                "ALLELE",
            ),
            (
                "stat",
                "Stat",
                "INVALID_STATISTIC",
                "RANGE",
            ),
            (
                "p_values",
                "p_value",
                "INVALID_P_VALUE",
                "RANGE",
            ),
            (
                "fdr",
                "FDR",
                "INVALID_FDR",
                "RANGE",
            ),
            (
                "beta",
                "beta",
                "INVALID_BETA",
                "RANGE",
            ),
            (
                "maf",
                "MAF",
                "INVALID_MAF",
                "RANGE",
            ),
            (
                "call_rate",
                "AvgCall",
                "INVALID_CALL_RATE",
                "RANGE",
            ),
            (
                "rsq",
                "Rsq",
                "INVALID_RSQ",
                "RANGE",
            ),
        ]

        for (
            error_key,
            column,
            anomaly_type,
            category,
        ) in field_definitions:

            self._add_validation_errors(
                errors=errors[
                    error_key
                ],
                workbook=workbook,
                sheet=sheet_name,
                column=column,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity="ERROR",
                source_dataframe=(
                    dataframe
                ),
            )

        feature_anomaly_type = (
            "INVALID_TRANSCRIPT_ID"
            if feature_validator
            == "transcript"
            else
            "INVALID_SPLICING_EVENT_ID"
        )

        self._add_validation_errors(
            errors=errors[
                "features"
            ],
            workbook=workbook,
            sheet=sheet_name,
            column=feature_column,
            anomaly_type=(
                feature_anomaly_type
            ),
            category="FORMAT",
            severity="ERROR",
            source_dataframe=(
                dataframe
            ),
        )

        error_counts = {
            name: len(values)
            for name, values
            in errors.items()
        }

        sheet_anomalies = (
            self.anomalies[
                anomalies_before:
            ]
        )

        source_schema_notes: list[
            str
        ] = []

        if snp_column == "T":
            source_schema_notes.append(
                (
                    "The source worksheet uses column 'T' "
                    "for dbSNP identifiers instead of 'SNP'. "
                    "The raw header is preserved."
                )
            )

        return {
            "status": (
                self
                ._status_from_anomalies(
                    sheet_anomalies
                )
            ),
            "rows": int(
                len(dataframe)
            ),
            "columns": int(
                len(
                    dataframe.columns
                )
            ),
            "missing_required_columns":
                [],
            "missing_values":
                self
                ._missing_value_statistics(
                    dataframe
                ),
            "unique_values":
                self
                ._unique_statistics(
                    dataframe,
                    [
                        snp_column,
                        "SNP_pos",
                        feature_column,
                    ],
                ),
            "duplicates":
                self
                ._duplicate_statistics(
                    dataframe
                ),
            "validation_error_counts":
                error_counts,
            "anomalies":
                sheet_anomalies,
            "source_schema_notes":
                source_schema_notes,
        }

    def _validate_s4_sheet(
    self,
    dataframe: pd.DataFrame,
    sheet_name: str,
    workbook: Path,
    ) -> dict[str, Any]:
        """
        Validate one Moradi Supplementary Table S4 worksheet.

        S4 contains GWAS-linked trans-QTL associations.

        Source semantics:
            sQTL_SNP      -> genomic coordinate
            sQTL_SNP-pos  -> dbSNP identifier
            tag_SNP       -> genomic coordinate
            tag_SNP_pos   -> single/composite dbSNP identifier

        Feature-column names are source-specific and are interpreted
        according to the worksheet schema rather than column name alone.
        """

        anomalies_before = len(
            self.anomalies
        )

        schema = self.S4_SHEETS[
            sheet_name
        ]

        feature_column = schema[
            "feature_column"
        ]

        feature_validator = schema[
            "feature_validator"
        ]

        required_columns = [
            "sQTL_SNP",
            "sQTL_SNP-pos",
            feature_column,
            "tag_SNP",
            "tag_SNP_pos",
            "LD",
            "gwas_cancer",
        ]

        missing_columns = (
            self._missing_required_columns(
                dataframe,
                required_columns,
            )
        )

        if missing_columns:
            self._add_anomaly(
                file=workbook,
                sheet=sheet_name,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_REQUIRED_COLUMN"
                ),
                category="MISSING",
                severity="ERROR",
                message=(
                    "Missing required columns: "
                    + ", ".join(
                        missing_columns
                    )
                ),
                action=(
                    "STOP_SHEET_VALIDATION"
                ),
            )

            sheet_anomalies = (
                self.anomalies[
                    anomalies_before:
                ]
            )

            return {
                "status": "FAIL",
                "rows": int(
                    len(dataframe)
                ),
                "columns": int(
                    len(
                        dataframe.columns
                    )
                ),
                "missing_required_columns":
                    missing_columns,
                "missing_values":
                    self
                    ._missing_value_statistics(
                        dataframe
                    ),
                "unique_values": {},
                "duplicates":
                    self
                    ._duplicate_statistics(
                        dataframe
                    ),
                "validation_error_counts":
                    {},
                "anomalies":
                    sheet_anomalies,
            }

        errors: dict[
            str,
            list[str],
        ] = {}

        errors[
            "qtl_coordinates"
        ] = (
            self.validator
            .validate_coordinates(
                dataframe[
                    "sQTL_SNP"
                ].tolist(),
                "sQTL_SNP",
            )
        )

        errors[
            "qtl_snp_ids"
        ] = (
            self.validator
            .validate_snp_ids(
                dataframe[
                    "sQTL_SNP-pos"
                ].tolist(),
                "sQTL_SNP-pos",
            )
        )

        errors[
            "gwas_coordinates"
        ] = (
            self.validator
            .validate_coordinates(
                dataframe[
                    "tag_SNP"
                ].tolist(),
                "tag_SNP",
            )
        )

        errors[
            "gwas_tag_snp_ids"
        ] = (
            self.validator
            .validate_composite_snp_ids(
                dataframe[
                    "tag_SNP_pos"
                ].tolist(),
                "tag_SNP_pos",
            )
        )

        if (
            feature_validator
            == "transcript"
        ):
            errors["features"] = (
                self.validator
                .validate_transcripts(
                    dataframe[
                        feature_column
                    ].tolist(),
                    feature_column,
                )
            )
        else:
            errors["features"] = (
                self.validator
                .validate_splicing_events(
                    dataframe[
                        feature_column
                    ].tolist(),
                    feature_column,
                )
            )

        errors["ld"] = (
            self.validator
            .validate_ld(
                dataframe[
                    "LD"
                ].tolist(),
                "LD",
            )
        )

        field_definitions = [
            (
                "qtl_coordinates",
                "sQTL_SNP",
                "INVALID_QTL_COORDINATE",
                "COORDINATE",
            ),
            (
                "qtl_snp_ids",
                "sQTL_SNP-pos",
                "INVALID_QTL_SNP_ID",
                "FORMAT",
            ),
            (
                "gwas_coordinates",
                "tag_SNP",
                "INVALID_GWAS_COORDINATE",
                "COORDINATE",
            ),
            (
                "gwas_tag_snp_ids",
                "tag_SNP_pos",
                "INVALID_GWAS_SNP_ID",
                "FORMAT",
            ),
            (
                "ld",
                "LD",
                "INVALID_LD_R2",
                "RANGE",
            ),
        ]

        for (
            error_key,
            column,
            anomaly_type,
            category,
        ) in field_definitions:

            self._add_validation_errors(
                errors=errors[
                    error_key
                ],
                workbook=workbook,
                sheet=sheet_name,
                column=column,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity="ERROR",
                source_dataframe=(
                    dataframe
                ),
            )

        feature_anomaly_type = (
            "INVALID_TRANSCRIPT_ID"
            if feature_validator
            == "transcript"
            else
            "INVALID_SPLICING_EVENT_ID"
        )

        self._add_validation_errors(
            errors=errors[
                "features"
            ],
            workbook=workbook,
            sheet=sheet_name,
            column=feature_column,
            anomaly_type=(
                feature_anomaly_type
            ),
            category="FORMAT",
            severity="ERROR",
            source_dataframe=(
                dataframe
            ),
        )

        error_counts = {
            name: len(values)
            for name, values
            in errors.items()
        }

        sheet_anomalies = (
            self.anomalies[
                anomalies_before:
            ]
        )

        return {
            "status": (
                self
                ._status_from_anomalies(
                    sheet_anomalies
                )
            ),
            "rows": int(
                len(dataframe)
            ),
            "columns": int(
                len(
                    dataframe.columns
                )
            ),
            "missing_required_columns":
                [],
            "missing_values":
                self
                ._missing_value_statistics(
                    dataframe
                ),
            "unique_values":
                self
                ._unique_statistics(
                    dataframe,
                    [
                        "sQTL_SNP",
                        "sQTL_SNP-pos",
                        feature_column,
                        "tag_SNP",
                        "tag_SNP_pos",
                        "gwas_cancer",
                    ],
                ),
            "duplicates":
                self
                ._duplicate_statistics(
                    dataframe
                ),
            "validation_error_counts":
                error_counts,
            "anomalies":
                sheet_anomalies,
            "source_schema_notes": [
                (
                    "sQTL_SNP contains genomic coordinates "
                    "while sQTL_SNP-pos contains dbSNP identifiers."
                ),
                (
                    "tag_SNP contains genomic coordinates "
                    "while tag_SNP_pos contains single or composite "
                    "dbSNP identifiers."
                ),
                (
                    "Feature-column semantics are interpreted "
                    "per worksheet because the source headers are "
                    "not fully consistent across S4."
                ),
            ],
        }

    # ==================================================================
    # Workbook validation
    # ==================================================================

    def validate_s1(
        self,
    ) -> dict[str, Any]:
        """Validate all core Supplementary Table S1 worksheets."""

        workbook = (
            self.moradi_dir
            / self.S1_FILE
        )

        if not workbook.exists():
            self._add_anomaly(
                file=workbook,
                sheet=None,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_WORKBOOK"
                ),
                category="MISSING",
                severity="ERROR",
                message=(
                    f"Workbook not found: "
                    f"{workbook}"
                ),
                action=(
                    "STOP_WORKBOOK_VALIDATION"
                ),
            )

            return {
                "status": "FAIL",
                "path": str(
                    workbook
                ),
                "sheets": {},
                "anomalies": (
                    self.anomalies[-1:]
                ),
            }

        workbook_anomalies_before = (
            len(
                self.anomalies
            )
        )

        excel_file = pd.ExcelFile(
            workbook
        )

        available_sheets = set(
            excel_file.sheet_names
        )

        sheets: dict[
            str,
            Any,
        ] = {}

        for sheet_name in (
            self.S1_SHEETS
        ):
            if (
                sheet_name
                not in available_sheets
            ):
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "MISSING_REQUIRED_SHEET"
                    ),
                    category="MISSING",
                    severity="ERROR",
                    message=(
                        "Required worksheet "
                        f"{sheet_name!r} "
                        "was not found."
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status":
                        "FAIL",
                    "rows":
                        0,
                    "columns":
                        0,
                    "anomalies":
                        self
                        .anomalies[-1:],
                }

                continue

            try:
                dataframe = (
                    self._read_sheet(
                        workbook,
                        sheet_name,
                    )
                )

                sheets[
                    sheet_name
                ] = (
                    self
                    ._validate_s1_sheet(
                        dataframe,
                        sheet_name,
                        workbook,
                    )
                )

            except Exception as exc:
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "VALIDATION_EXCEPTION"
                    ),
                    category="SYSTEM",
                    severity="ERROR",
                    message=str(
                        exc
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status":
                        "FAIL",
                    "rows":
                        0,
                    "columns":
                        0,
                    "anomalies":
                        self
                        .anomalies[-1:],
                }

        workbook_anomalies = (
            self.anomalies[
                workbook_anomalies_before:
            ]
        )

        return {
            "status": (
                self
                ._aggregate_status(
                    [
                        result.get(
                            "status",
                            "FAIL",
                        )
                        for result
                        in sheets.values()
                    ]
                )
            ),
            "path": str(
                workbook
            ),
            "sheets": sheets,
            "anomalies":
                workbook_anomalies,
        }

    def validate_s2(
        self,
    ) -> dict[str, Any]:
        """Validate all core Supplementary Table S2 worksheets."""

        workbook = (
            self.moradi_dir
            / self.S2_FILE
        )

        if not workbook.exists():
            self._add_anomaly(
                file=workbook,
                sheet=None,
                row=None,
                column=None,
                value=None,
                anomaly_type=(
                    "MISSING_WORKBOOK"
                ),
                category="MISSING",
                severity="ERROR",
                message=(
                    f"Workbook not found: "
                    f"{workbook}"
                ),
                action=(
                    "STOP_WORKBOOK_VALIDATION"
                ),
            )

            return {
                "status": "FAIL",
                "path": str(
                    workbook
                ),
                "sheets": {},
                "anomalies": (
                    self.anomalies[-1:]
                ),
            }

        workbook_anomalies_before = (
            len(
                self.anomalies
            )
        )

        excel_file = pd.ExcelFile(
            workbook
        )

        available_sheets = set(
            excel_file.sheet_names
        )

        sheets: dict[
            str,
            Any,
        ] = {}

        for sheet_name in (
            self.S2_SHEETS
        ):
            if (
                sheet_name
                not in available_sheets
            ):
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "MISSING_REQUIRED_SHEET"
                    ),
                    category="MISSING",
                    severity="ERROR",
                    message=(
                        "Required worksheet "
                        f"{sheet_name!r} "
                        "was not found."
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status":
                        "FAIL",
                    "rows":
                        0,
                    "columns":
                        0,
                    "anomalies":
                        self
                        .anomalies[-1:],
                }

                continue

            try:
                dataframe = (
                    self._read_sheet(
                        workbook,
                        sheet_name,
                    )
                )

                sheets[
                    sheet_name
                ] = (
                    self
                    ._validate_s2_sheet(
                        dataframe,
                        sheet_name,
                        workbook,
                    )
                )

            except Exception as exc:
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "VALIDATION_EXCEPTION"
                    ),
                    category="SYSTEM",
                    severity="ERROR",
                    message=str(
                        exc
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status":
                        "FAIL",
                    "rows":
                        0,
                    "columns":
                        0,
                    "anomalies":
                        self
                        .anomalies[-1:],
                }

        workbook_anomalies = (
            self.anomalies[
                workbook_anomalies_before:
            ]
        )

        return {
            "status": (
                self
                ._aggregate_status(
                    [
                        result.get(
                            "status",
                            "FAIL",
                        )
                        for result
                        in sheets.values()
                    ]
                )
            ),
            "path": str(
                workbook
            ),
            "sheets": sheets,
            "anomalies":
                workbook_anomalies,
        }

    def validate_s3(
        self,
    ) -> dict[str, Any]:
        """Validate all Supplementary Table S3 worksheets."""

        workbook = (
            self.moradi_dir
            / self.S3_FILE
        )

        if not workbook.exists():
            self._add_anomaly(
                file=workbook,
                sheet=None,
                row=None,
                column=None,
                value=None,
                anomaly_type="MISSING_WORKBOOK",
                category="MISSING",
                severity="ERROR",
                message=(
                    f"Workbook not found: {workbook}"
                ),
                action="STOP_WORKBOOK_VALIDATION",
            )

            return {
                "status": "FAIL",
                "path": str(workbook),
                "sheets": {},
                "anomalies":
                    self.anomalies[-1:],
            }

        anomalies_before = len(
            self.anomalies
        )

        excel_file = pd.ExcelFile(
            workbook
        )

        available_sheets = set(
            excel_file.sheet_names
        )

        sheets: dict[str, Any] = {}

        for sheet_name in self.S3_SHEETS:

            if (
                sheet_name
                not in available_sheets
            ):
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "MISSING_REQUIRED_SHEET"
                    ),
                    category="MISSING",
                    severity="ERROR",
                    message=(
                        f"Required worksheet "
                        f"{sheet_name!r} was not found."
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status": "FAIL",
                    "rows": 0,
                    "columns": 0,
                    "anomalies":
                        self.anomalies[-1:],
                }

                continue

            try:
                dataframe = (
                    self._read_sheet(
                        workbook,
                        sheet_name,
                    )
                )

                sheets[
                    sheet_name
                ] = (
                    self._validate_s3_sheet(
                        dataframe,
                        sheet_name,
                        workbook,
                    )
                )

            except Exception as exc:
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "VALIDATION_EXCEPTION"
                    ),
                    category="SYSTEM",
                    severity="ERROR",
                    message=str(exc),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status": "FAIL",
                    "rows": 0,
                    "columns": 0,
                    "anomalies":
                        self.anomalies[-1:],
                }

        return {
            "status":
                self._aggregate_status(
                    [
                        result.get(
                            "status",
                            "FAIL",
                        )
                        for result
                        in sheets.values()
                    ]
                ),
            "path": str(workbook),
            "sheets": sheets,
            "anomalies":
                self.anomalies[
                    anomalies_before:
                ],
        }


    def validate_s4(
        self,
    ) -> dict[str, Any]:
        """Validate all Supplementary Table S4 worksheets."""

        workbook = (
            self.moradi_dir
            / self.S4_FILE
        )

        if not workbook.exists():
            self._add_anomaly(
                file=workbook,
                sheet=None,
                row=None,
                column=None,
                value=None,
                anomaly_type="MISSING_WORKBOOK",
                category="MISSING",
                severity="ERROR",
                message=(
                    f"Workbook not found: {workbook}"
                ),
                action="STOP_WORKBOOK_VALIDATION",
            )

            return {
                "status": "FAIL",
                "path": str(workbook),
                "sheets": {},
                "anomalies":
                    self.anomalies[-1:],
            }

        anomalies_before = len(
            self.anomalies
        )

        excel_file = pd.ExcelFile(
            workbook
        )

        available_sheets = set(
            excel_file.sheet_names
        )

        sheets: dict[str, Any] = {}

        for sheet_name in self.S4_SHEETS:

            if (
                sheet_name
                not in available_sheets
            ):
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "MISSING_REQUIRED_SHEET"
                    ),
                    category="MISSING",
                    severity="ERROR",
                    message=(
                        f"Required worksheet "
                        f"{sheet_name!r} was not found."
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status": "FAIL",
                    "rows": 0,
                    "columns": 0,
                    "anomalies":
                        self.anomalies[-1:],
                }

                continue

            try:
                dataframe = (
                    self._read_sheet(
                        workbook,
                        sheet_name,
                    )
                )

                sheets[
                    sheet_name
                ] = (
                    self._validate_s4_sheet(
                        dataframe,
                        sheet_name,
                        workbook,
                    )
                )

            except Exception as exc:
                self._add_anomaly(
                    file=workbook,
                    sheet=sheet_name,
                    row=None,
                    column=None,
                    value=None,
                    anomaly_type=(
                        "VALIDATION_EXCEPTION"
                    ),
                    category="SYSTEM",
                    severity="ERROR",
                    message=str(exc),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                sheets[
                    sheet_name
                ] = {
                    "status": "FAIL",
                    "rows": 0,
                    "columns": 0,
                    "anomalies":
                        self.anomalies[-1:],
                }

        return {
            "status":
                self._aggregate_status(
                    [
                        result.get(
                            "status",
                            "FAIL",
                        )
                        for result
                        in sheets.values()
                    ]
                ),
            "path": str(workbook),
            "sheets": sheets,
            "anomalies":
                self.anomalies[
                    anomalies_before:
                ],
        }

    # ==================================================================
    # Summary
    # ==================================================================

    @staticmethod
    def _summarize_anomalies(
        anomalies: list[
            dict[str, Any]
        ],
    ) -> dict[str, Any]:
        """Create anomaly summary statistics."""

        by_category: dict[
            str,
            int,
        ] = {}

        by_severity: dict[
            str,
            int,
        ] = {}

        by_type: dict[
            str,
            int,
        ] = {}

        for anomaly in anomalies:
            category = str(
                anomaly.get(
                    "category",
                    "UNKNOWN",
                )
            )

            severity = str(
                anomaly.get(
                    "severity",
                    "UNKNOWN",
                )
            )

            anomaly_type = str(
                anomaly.get(
                    "type",
                    "UNKNOWN",
                )
            )

            by_category[
                category
            ] = (
                by_category.get(
                    category,
                    0,
                )
                + 1
            )

            by_severity[
                severity
            ] = (
                by_severity.get(
                    severity,
                    0,
                )
                + 1
            )

            by_type[
                anomaly_type
            ] = (
                by_type.get(
                    anomaly_type,
                    0,
                )
                + 1
            )

        return {
            "total_anomalies":
                len(anomalies),
            "by_category":
                by_category,
            "by_severity":
                by_severity,
            "by_type":
                by_type,
        }

    # ==================================================================
    # Main API
    # ==================================================================

    def validate_core(
        self,
    ) -> dict[str, Any]:
        """
        Validate core Moradi S1 and S2 datasets.

        A fresh anomaly sequence is created on every validation run.
        """

        # Reset state.
        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies = []

        s1 = self.validate_s1()
        s2 = self.validate_s2()

        status = (
            self
            ._aggregate_status(
                [
                    s1["status"],
                    s2["status"],
                ]
            )
        )

        return {
            "dataset_id":
                self.DATASET_ID,
            "status":
                status,
            "files": {
                "supplementary_s1":
                    s1,
                "supplementary_s2":
                    s2,
            },
            "anomalies":
                list(
                    self.anomalies
                ),
            "validation_summary":
                self
                ._summarize_anomalies(
                    self.anomalies
                ),
        }

    def validate_extended_qtl(
        self,
    ) -> dict[str, Any]:
        """
        Validate Moradi S1-S4 regulatory-QTL datasets.
        """

        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies = []

        s1 = self.validate_s1()
        s2 = self.validate_s2()
        s3 = self.validate_s3()
        s4 = self.validate_s4()

        status = (
            self._aggregate_status(
                [
                    s1["status"],
                    s2["status"],
                    s3["status"],
                    s4["status"],
                ]
            )
        )

        return {
            "dataset_id":
                self.DATASET_ID,
            "status":
                status,
            "files": {
                "supplementary_s1": s1,
                "supplementary_s2": s2,
                "supplementary_s3": s3,
                "supplementary_s4": s4,
            },
            "anomalies":
                list(
                    self.anomalies
                ),
            "validation_summary":
                self._summarize_anomalies(
                    self.anomalies
                ),
        }