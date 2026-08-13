"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/runners/trfqtl_workbook.py

Description:
    Validates the Cancer-tRFQTL supplementary workbook used by the
    pcatRFQTL research pipeline.

    The validator performs sheet-specific validation for:
        - SNP identifiers
        - Genomic coordinates
        - Alleles
        - tRF identifiers
        - P-values
        - FDR values
        - LD r2 values
        - GWAS effect alleles
        - Effect estimates

    Validation results include standardized data-quality anomalies
    with globally sequential identifiers in the form:

        ANOM-000001
        ANOM-000002
        ...

    Validation severity determines the final sheet and dataset status:

        PASS
            No validation anomalies.

        WARNING
            One or more warning-level anomalies exist, but no
            error-level anomalies are present.

        FAIL
            At least one error-level anomaly exists.

    Raw source files are never modified by this validator.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from pcatrfqtl.validation.anomaly import AnomalyFactory
from pcatrfqtl.validation.validators.trfqtl import TRFQTLValidator


class TRFQTLWorkbookValidator:
    """
    Validate the Cancer-tRFQTL supplementary workbook.

    The workbook contains supplementary tables S1-S11. Only sheets
    containing structured molecular or genetic association data are
    subjected to formal validation.

    Currently validated sheets
    --------------------------
    S1
        Survival-associated tRFQTLs.

    S2
        GWAS-associated tRFQTLs.

    S10
        Cancer-risk-associated tRFQTLs from GWAS meta-analysis.

    Parameters
    ----------
    workbook:
        Path to the Cancer-tRFQTL supplementary workbook.
    """

    DATASET_ID = "cancer_trfqtl"

    SUPPORTED_SHEETS = {
        "S1",
        "S2",
        "S10",
    }

    def __init__(
        self,
        workbook: str | Path,
    ) -> None:
        """Initialize the workbook validator."""

        self.workbook = Path(workbook)

        self.validator = TRFQTLValidator()

        # One factory per validation run.
        # This produces globally sequential occurrence IDs.
        self.anomaly_factory = AnomalyFactory()

        self.anomalies: list[dict[str, Any]] = []

    # ==================================================================
    # Reading and dataframe helpers
    # ==================================================================

    def _read_sheet(
        self,
        sheet_name: str,
    ) -> pd.DataFrame:
        """
        Read a supplementary worksheet.

        The Cancer-tRFQTL workbook contains a descriptive title in
        Excel row 1 and the real table header in Excel row 2.

        Therefore Pandas must read the worksheet with header=1.
        """

        return pd.read_excel(
            self.workbook,
            sheet_name=sheet_name,
            header=1,
        )

    @staticmethod
    def _clean_dataframe(
        df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Clean a dataframe without modifying raw source data.

        Operations
        ----------
        - Remove completely empty rows.
        - Remove completely empty columns.
        - Strip whitespace from column names.
        - Reset the dataframe index.

        No source values are modified.
        """

        cleaned = df.dropna(
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
            for column in cleaned.columns
        ]

        return cleaned.reset_index(
            drop=True
        )

    @staticmethod
    def _missing_value_statistics(
        df: pd.DataFrame,
    ) -> dict[str, int]:
        """
        Calculate missing-value counts for every dataframe column.
        """

        return {
            str(column): int(
                df[column]
                .isna()
                .sum()
            )
            for column in df.columns
        }

    @staticmethod
    def _unique_statistics(
        df: pd.DataFrame,
        columns: list[str],
    ) -> dict[str, int]:
        """
        Calculate unique non-null values for selected columns.

        String values are stripped before uniqueness is calculated.
        """

        result: dict[str, int] = {}

        for column in columns:
            if column not in df.columns:
                continue

            series = df[column].dropna()

            if (
                pd.api.types
                .is_object_dtype(series.dtype)
                or pd.api.types
                .is_string_dtype(series.dtype)
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
    def _unique_values(
        df: pd.DataFrame,
        columns: list[str],
    ) -> dict[str, int]:
        """
        Backward-compatible alias for `_unique_statistics`.

        Retained because earlier versions of the workbook runner used
        `_unique_values`.
        """

        return (
            TRFQTLWorkbookValidator
            ._unique_statistics(
                df,
                columns,
            )
        )

    @staticmethod
    def _duplicate_statistics(
        df: pd.DataFrame,
    ) -> dict[str, int]:
        """
        Calculate exact-row duplicate statistics.

        duplicate_rows
            Number of rows participating in duplicate groups.

        duplicate_groups
            Number of unique duplicated row patterns.
        """

        duplicated = df.duplicated(
            keep=False
        )

        duplicate_rows = int(
            duplicated.sum()
        )

        if duplicate_rows == 0:
            duplicate_groups = 0
        else:
            duplicate_groups = int(
                df.loc[duplicated]
                .drop_duplicates()
                .shape[0]
            )

        return {
            "duplicate_rows": duplicate_rows,
            "duplicate_groups": duplicate_groups,
        }

    # ==================================================================
    # Status helpers
    # ==================================================================

    @staticmethod
    def _status_from_anomalies(
        anomalies: list[dict[str, Any]],
    ) -> str:
        """
        Determine validation status from anomaly severity.

        Returns
        -------
        str
            PASS
                No anomalies.

            WARNING
                Warning-level anomalies exist but no errors.

            FAIL
                At least one ERROR-level anomaly exists.
        """

        if not anomalies:
            return "PASS"

        severities = {
            str(
                anomaly.get(
                    "severity",
                    "",
                )
            ).upper()
            for anomaly in anomalies
        }

        if "ERROR" in severities:
            return "FAIL"

        return "WARNING"

    @staticmethod
    def _aggregate_status(
        statuses: list[str],
    ) -> str:
        """
        Aggregate multiple sheet statuses into one dataset status.
        """

        normalized = {
            status.upper()
            for status in statuses
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
        sheet: str,
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
        Register a standardized data-quality anomaly.

        Parameters
        ----------
        sheet:
            Workbook sheet containing the anomaly.

        row:
            One-based worksheet row number if available.

        column:
            Column associated with the anomaly.

        value:
            Original source value, when available.

        anomaly_type:
            Stable semantic anomaly type.

        category:
            General anomaly category.

        severity:
            INFO, WARNING, or ERROR.

        message:
            Human-readable validation explanation.

        action:
            Recommended downstream handling.
        """

        anomaly = self.anomaly_factory.create(
            type=anomaly_type,
            category=category,
            severity=severity,
            dataset=self.DATASET_ID,
            file=str(self.workbook),
            sheet=sheet,
            row=row,
            column=column,
            value=value,
            message=message,
            action=action,
        )

        self.anomalies.append(
            anomaly.to_dict()
        )

    def _add_validation_errors(
        self,
        *,
        errors: list[str],
        sheet: str,
        column: str,
        anomaly_type: str,
        category: str,
        severity: str = "ERROR",
        action: str = "FLAG_FOR_REVIEW",
    ) -> None:
        """
        Convert field-level validation messages into anomalies.

        Field validators currently return human-readable validation
        messages. Therefore row/value metadata may not yet be available
        in structured form and are recorded as None.

        A future validator API may return structured field-validation
        results containing exact dataframe row and source value.
        """

        for error in errors:
            self._add_anomaly(
                sheet=sheet,
                row=None,
                column=column,
                value=None,
                anomaly_type=anomaly_type,
                category=category,
                severity=severity,
                message=error,
                action=action,
            )

    # ==================================================================
    # Required-column helper
    # ==================================================================

    @staticmethod
    def _missing_required_columns(
        df: pd.DataFrame,
        required_columns: list[str],
    ) -> list[str]:
        """Return required columns absent from the dataframe."""

        available = set(
            df.columns
        )

        return [
            column
            for column in required_columns
            if column not in available
        ]

    # ==================================================================
    # Supplementary Table S1
    # ==================================================================

    def _validate_s1(
        self,
        df: pd.DataFrame,
    ) -> dict[str, Any]:
        """
        Validate Supplementary Table S1.

        S1 contains survival-associated tRFQTLs.
        """

        required_columns = [
            "Cancer type",
            "SNP ID",
            "SNP position (hg19)",
            "Alleles",
            "tRF",
            "P-value (survival)",
        ]

        errors = {
            "snp_ids": 0,
            "coordinates": 0,
            "alleles": 0,
            "trf_ids": 0,
            "p_values": 0,
        }

        anomalies_before = len(
            self.anomalies
        )

        missing = (
            self._missing_required_columns(
                df,
                required_columns,
            )
        )

        if missing:
            self._add_anomaly(
                sheet="S1",
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
                    + ", ".join(missing)
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
                "rows": int(len(df)),
                "columns": int(
                    len(df.columns)
                ),
                "unique": {},
                "missing_values": (
                    self
                    ._missing_value_statistics(
                        df
                    )
                ),
                "duplicates": (
                    self
                    ._duplicate_statistics(
                        df
                    )
                ),
                "validation_error_counts": {
                    **errors,
                    "missing_columns": len(
                        missing
                    ),
                },
                "anomalies": sheet_anomalies,
            }

        snp_errors = (
            self.validator
            .validate_snp_ids(
                df["SNP ID"].tolist(),
                "SNP ID",
            )
        )

        coordinate_errors = (
            self.validator
            .validate_coordinates(
                df[
                    "SNP position (hg19)"
                ].tolist(),
                "SNP position (hg19)",
            )
        )

        allele_errors = (
            self.validator
            .validate_alleles(
                df["Alleles"].tolist()
            )
        )

        trf_errors = (
            self.validator
            .validate_trf_ids(
                df["tRF"].tolist()
            )
        )

        pvalue_errors = (
            self.validator
            .validate_p_values(
                df[
                    "P-value (survival)"
                ].tolist(),
                "P-value (survival)",
            )
        )

        errors["snp_ids"] = len(
            snp_errors
        )

        errors["coordinates"] = len(
            coordinate_errors
        )

        errors["alleles"] = len(
            allele_errors
        )

        errors["trf_ids"] = len(
            trf_errors
        )

        errors["p_values"] = len(
            pvalue_errors
        )

        self._add_validation_errors(
            errors=snp_errors,
            sheet="S1",
            column="SNP ID",
            anomaly_type=(
                "INVALID_SNP_ID"
            ),
            category="FORMAT",
        )

        self._add_validation_errors(
            errors=coordinate_errors,
            sheet="S1",
            column=(
                "SNP position (hg19)"
            ),
            anomaly_type=(
                "INVALID_COORDINATE"
            ),
            category="COORDINATE",
        )

        self._add_validation_errors(
            errors=allele_errors,
            sheet="S1",
            column="Alleles",
            anomaly_type=(
                "INVALID_ALLELE"
            ),
            category="ALLELE",
        )

        self._add_validation_errors(
            errors=trf_errors,
            sheet="S1",
            column="tRF",
            anomaly_type=(
                "INVALID_TRF_ID"
            ),
            category="FORMAT",
        )

        self._add_validation_errors(
            errors=pvalue_errors,
            sheet="S1",
            column=(
                "P-value (survival)"
            ),
            anomaly_type=(
                "INVALID_P_VALUE"
            ),
            category="RANGE",
        )

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
            "rows": int(len(df)),
            "columns": int(
                len(df.columns)
            ),
            "unique": (
                self._unique_statistics(
                    df,
                    [
                        "Cancer type",
                        "SNP ID",
                        "tRF",
                    ],
                )
            ),
            "missing_values": (
                self
                ._missing_value_statistics(
                    df
                )
            ),
            "duplicates": (
                self
                ._duplicate_statistics(
                    df
                )
            ),
            "validation_error_counts": (
                errors
            ),
            "anomalies": (
                sheet_anomalies
            ),
        }

    # ==================================================================
    # Supplementary Table S2
    # ==================================================================

    def _validate_s2(
        self,
        df: pd.DataFrame,
    ) -> dict[str, Any]:
        """
        Validate Supplementary Table S2.

        S2 contains GWAS-associated tRFQTLs.
        """

        required_columns = [
            "Cancer type",
            "SNP ID",
            "SNP position (hg19)",
            "Alleles",
            "tRF",
            "statistic",
            "P-value",
            "FDR",
            "GWAS tagSnp",
            "LD (r2)",
        ]

        errors = {
            "snp_ids": 0,
            "coordinates": 0,
            "alleles": 0,
            "trf_ids": 0,
            "p_values": 0,
            "fdr": 0,
            "ld_r2": 0,
        }

        anomalies_before = len(
            self.anomalies
        )

        missing = (
            self._missing_required_columns(
                df,
                required_columns,
            )
        )

        if missing:
            self._add_anomaly(
                sheet="S2",
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
                    + ", ".join(missing)
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
                "rows": int(len(df)),
                "columns": int(
                    len(df.columns)
                ),
                "unique": {},
                "missing_values": (
                    self
                    ._missing_value_statistics(
                        df
                    )
                ),
                "duplicates": (
                    self
                    ._duplicate_statistics(
                        df
                    )
                ),
                "validation_error_counts": {
                    **errors,
                    "missing_columns": len(
                        missing
                    ),
                },
                "anomalies": sheet_anomalies,
            }

        validators = [
            (
                "snp_ids",
                self.validator
                .validate_snp_ids(
                    df["SNP ID"].tolist(),
                    "SNP ID",
                ),
                "SNP ID",
                "INVALID_SNP_ID",
                "FORMAT",
                "ERROR",
            ),
            (
                "coordinates",
                self.validator
                .validate_coordinates(
                    df[
                        "SNP position (hg19)"
                    ].tolist(),
                    "SNP position (hg19)",
                ),
                "SNP position (hg19)",
                "INVALID_COORDINATE",
                "COORDINATE",
                "ERROR",
            ),
            (
                "alleles",
                self.validator
                .validate_alleles(
                    df["Alleles"].tolist()
                ),
                "Alleles",
                "INVALID_ALLELE",
                "ALLELE",
                "ERROR",
            ),
            (
                "trf_ids",
                self.validator
                .validate_trf_ids(
                    df["tRF"].tolist()
                ),
                "tRF",
                "INVALID_TRF_ID",
                "FORMAT",
                "ERROR",
            ),
            (
                "p_values",
                self.validator
                .validate_p_values(
                    df[
                        "P-value"
                    ].tolist(),
                    "P-value",
                ),
                "P-value",
                "INVALID_P_VALUE",
                "RANGE",
                "ERROR",
            ),
            (
                "fdr",
                self.validator
                .validate_fdr(
                    df["FDR"].tolist()
                ),
                "FDR",
                "INVALID_FDR",
                "RANGE",
                "ERROR",
            ),
            (
                "ld_r2",
                self.validator
                .validate_ld_r2(
                    df[
                        "LD (r2)"
                    ].tolist()
                ),
                "LD (r2)",
                "INVALID_LD_R2",
                "RANGE",
                "ERROR",
            ),
        ]

        for (
            key,
            validation_errors,
            column,
            anomaly_type,
            category,
            severity,
        ) in validators:

            errors[key] = len(
                validation_errors
            )

            self._add_validation_errors(
                errors=validation_errors,
                sheet="S2",
                column=column,
                anomaly_type=(
                    anomaly_type
                ),
                category=category,
                severity=severity,
            )

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
            "rows": int(len(df)),
            "columns": int(
                len(df.columns)
            ),
            "unique": (
                self._unique_statistics(
                    df,
                    [
                        "Cancer type",
                        "SNP ID",
                        "tRF",
                        "GWAS tagSnp",
                    ],
                )
            ),
            "missing_values": (
                self
                ._missing_value_statistics(
                    df
                )
            ),
            "duplicates": (
                self
                ._duplicate_statistics(
                    df
                )
            ),
            "validation_error_counts": (
                errors
            ),
            "anomalies": (
                sheet_anomalies
            ),
        }

    # ==================================================================
    # Supplementary Table S10
    # ==================================================================

    def _validate_s10(
        self,
        df: pd.DataFrame,
    ) -> dict[str, Any]:
        """
        Validate Supplementary Table S10.

        S10 contains tRFQTLs associated with cancer risk in GWAS
        meta-analyses.

        Invalid SNP identifiers in this source table are classified
        as WARNING-level source-data anomalies rather than fatal
        structural errors.

        This policy preserves the source value while allowing the
        remainder of the sheet to remain usable.
        """

        required_columns = [
            "Cancer type",
            "SNP ID",
            "Position (hg19)",
            "A1 (effect allele)",
            "Effect",
            "GWAS P-value",
        ]

        errors = {
            "snp_ids": 0,
            "coordinates": 0,
            "effect_alleles": 0,
            "effects": 0,
            "p_values": 0,
        }

        anomalies_before = len(
            self.anomalies
        )

        missing = (
            self._missing_required_columns(
                df,
                required_columns,
            )
        )

        if missing:
            self._add_anomaly(
                sheet="S10",
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
                    + ", ".join(missing)
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
                "rows": int(len(df)),
                "columns": int(
                    len(df.columns)
                ),
                "unique": {},
                "missing_values": (
                    self
                    ._missing_value_statistics(
                        df
                    )
                ),
                "duplicates": (
                    self
                    ._duplicate_statistics(
                        df
                    )
                ),
                "validation_error_counts": {
                    **errors,
                    "missing_columns": len(
                        missing
                    ),
                },
                "anomalies": sheet_anomalies,
            }

        snp_errors = (
            self.validator
            .validate_snp_ids(
                df["SNP ID"].tolist(),
                "SNP ID",
            )
        )

        coordinate_errors = (
            self.validator
            .validate_coordinates(
                df[
                    "Position (hg19)"
                ].tolist(),
                "Position (hg19)",
            )
        )

        allele_errors = (
            self.validator
            .validate_effect_alleles(
                df[
                    "A1 (effect allele)"
                ].tolist()
            )
        )

        effect_errors = (
            self.validator
            .validate_effects(
                df["Effect"].tolist()
            )
        )

        pvalue_errors = (
            self.validator
            .validate_p_values(
                df[
                    "GWAS P-value"
                ].tolist(),
                "GWAS P-value",
            )
        )

        errors["snp_ids"] = len(
            snp_errors
        )

        errors["coordinates"] = len(
            coordinate_errors
        )

        errors["effect_alleles"] = len(
            allele_errors
        )

        errors["effects"] = len(
            effect_errors
        )

        errors["p_values"] = len(
            pvalue_errors
        )

        # --------------------------------------------------------------
        # SNP-ID anomalies in S10
        #
        # The source workbook contains at least one known non-rsID
        # value ("6"). This is retained as a source-data anomaly.
        # --------------------------------------------------------------

        self._add_validation_errors(
            errors=snp_errors,
            sheet="S10",
            column="SNP ID",
            anomaly_type=(
                "INVALID_SNP_ID"
            ),
            category="SOURCE",
            severity="WARNING",
            action=(
                "FLAG_AND_EXCLUDE_FROM_RSID_DEPENDENT_ANALYSIS"
            ),
        )

        self._add_validation_errors(
            errors=coordinate_errors,
            sheet="S10",
            column="Position (hg19)",
            anomaly_type=(
                "INVALID_COORDINATE"
            ),
            category="COORDINATE",
            severity="ERROR",
        )

        self._add_validation_errors(
            errors=allele_errors,
            sheet="S10",
            column=(
                "A1 (effect allele)"
            ),
            anomaly_type=(
                "INVALID_EFFECT_ALLELE"
            ),
            category="ALLELE",
            severity="ERROR",
        )

        self._add_validation_errors(
            errors=effect_errors,
            sheet="S10",
            column="Effect",
            anomaly_type=(
                "INVALID_EFFECT"
            ),
            category="RANGE",
            severity="ERROR",
        )

        self._add_validation_errors(
            errors=pvalue_errors,
            sheet="S10",
            column="GWAS P-value",
            anomaly_type=(
                "INVALID_P_VALUE"
            ),
            category="RANGE",
            severity="ERROR",
        )

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
            "rows": int(len(df)),
            "columns": int(
                len(df.columns)
            ),
            "unique": (
                self._unique_statistics(
                    df,
                    [
                        "Cancer type",
                        "SNP ID",
                    ],
                )
            ),
            "missing_values": (
                self
                ._missing_value_statistics(
                    df
                )
            ),
            "duplicates": (
                self
                ._duplicate_statistics(
                    df
                )
            ),
            "validation_error_counts": (
                errors
            ),
            "anomalies": (
                sheet_anomalies
            ),
        }

    # ==================================================================
    # Validation summary
    # ==================================================================

    @staticmethod
    def _summarize_anomalies(
        anomalies: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """
        Create anomaly summary statistics.
        """

        by_category: dict[str, int] = {}
        by_severity: dict[str, int] = {}
        by_type: dict[str, int] = {}

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

            by_category[category] = (
                by_category.get(
                    category,
                    0,
                )
                + 1
            )

            by_severity[severity] = (
                by_severity.get(
                    severity,
                    0,
                )
                + 1
            )

            by_type[anomaly_type] = (
                by_type.get(
                    anomaly_type,
                    0,
                )
                + 1
            )

        return {
            "total_anomalies": len(
                anomalies
            ),
            "by_category": (
                by_category
            ),
            "by_severity": (
                by_severity
            ),
            "by_type": (
                by_type
            ),
        }

    # ==================================================================
    # Main validation
    # ==================================================================

    def validate_all(
        self,
    ) -> dict[str, Any]:
        """
        Validate all supported structured sheets.

        Returns
        -------
        dict
            Complete Cancer-tRFQTL workbook validation report.
        """

        if not self.workbook.exists():
            raise FileNotFoundError(
                "Workbook not found: "
                f"{self.workbook}"
            )

        if not self.workbook.is_file():
            raise ValueError(
                "Workbook path is not a file: "
                f"{self.workbook}"
            )

        # --------------------------------------------------------------
        # Reset state.
        #
        # Repeated calls must generate a fresh, deterministic sequence
        # starting again from ANOM-000001.
        # --------------------------------------------------------------

        self.anomaly_factory = (
            AnomalyFactory()
        )

        self.anomalies = []

        report: dict[str, Any] = {
            "dataset_id": (
                self.DATASET_ID
            ),
            "file": str(
                self.workbook
            ),
            "status": "PASS",
            "sheets": {},
            "anomalies": [],
            "validation_summary": {},
        }

        excel_file = pd.ExcelFile(
            self.workbook
        )

        available_sheets = set(
            excel_file.sheet_names
        )

        # --------------------------------------------------------------
        # Validate only explicitly supported structured sheets.
        # --------------------------------------------------------------

        for sheet_name in [
            "S1",
            "S2",
            "S10",
        ]:

            if (
                sheet_name
                not in available_sheets
            ):
                self._add_anomaly(
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
                        f"Required supported sheet "
                        f"{sheet_name!r} "
                        "was not found in the workbook."
                    ),
                    action=(
                        "STOP_SHEET_VALIDATION"
                    ),
                )

                report["sheets"][
                    sheet_name
                ] = {
                    "status": "FAIL",
                    "rows": 0,
                    "columns": 0,
                    "unique": {},
                    "missing_values": {},
                    "duplicates": {
                        "duplicate_rows": 0,
                        "duplicate_groups": 0,
                    },
                    "validation_error_counts": {
                        "missing_sheet": 1,
                    },
                    "anomalies": (
                        self.anomalies[-1:]
                    ),
                }

                continue

            anomalies_before = len(
                self.anomalies
            )

            try:
                df = self._read_sheet(
                    sheet_name
                )

                df = (
                    self
                    ._clean_dataframe(
                        df
                    )
                )

                if sheet_name == "S1":
                    result = (
                        self._validate_s1(
                            df
                        )
                    )

                elif sheet_name == "S2":
                    result = (
                        self._validate_s2(
                            df
                        )
                    )

                elif sheet_name == "S10":
                    result = (
                        self._validate_s10(
                            df
                        )
                    )

                else:
                    continue

                report["sheets"][
                    sheet_name
                ] = result

            except Exception as exc:
                self._add_anomaly(
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

                sheet_anomalies = (
                    self.anomalies[
                        anomalies_before:
                    ]
                )

                report["sheets"][
                    sheet_name
                ] = {
                    "status": "FAIL",
                    "rows": 0,
                    "columns": 0,
                    "unique": {},
                    "missing_values": {},
                    "duplicates": {
                        "duplicate_rows": 0,
                        "duplicate_groups": 0,
                    },
                    "validation_error_counts": {},
                    "anomalies": (
                        sheet_anomalies
                    ),
                }

        # --------------------------------------------------------------
        # Final anomaly collection
        # --------------------------------------------------------------

        report["anomalies"] = list(
            self.anomalies
        )

        report[
            "validation_summary"
        ] = self._summarize_anomalies(
            self.anomalies
        )

        # --------------------------------------------------------------
        # Determine final status from sheet statuses.
        # --------------------------------------------------------------

        sheet_statuses = [
            str(
                result.get(
                    "status",
                    "FAIL",
                )
            )
            for result
            in report["sheets"].values()
        ]

        report["status"] = (
            self._aggregate_status(
                sheet_statuses
            )
        )

        return report