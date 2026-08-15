"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/runners/standardize_trfqtl_all.py

Description:
    Backfill standardization runner for the complete Cancer-tRFQTL
    supplementary workbook.

    The workbook contains eleven supplementary tables (S1-S11).

    Existing core standardized tables:

        S2  - GWAS-associated tRFQTLs
        S10 - cancer-risk-associated variants

    are treated as locked and are not rewritten by this runner.

    Newly standardized tables:

        S1
        S3
        S4
        S5
        S6
        S7
        S8
        S9
        S11

    are generated using CancerTRFQTLSupplementaryStandardizer.

    The runner preserves source semantics and row cardinality and does
    not perform disease filtering, genome liftover, cross-source joins,
    or statistical filtering.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from pcatrfqtl.io.parquet import (
    read_parquet,
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)
from pcatrfqtl.standardization.trfqtl_supplementary import (
    CancerTRFQTLSupplementaryStandardizer,
    DATASET_ID,
    SHEET_REGISTRY,
)


logger = get_logger(
    __name__
)


@dataclass(frozen=True)
class TRFQTLAllInputs:
    """Input paths required for complete Cancer-tRFQTL standardization."""

    workbook: Path


class CancerTRFQTLAllStandardizationRunner:
    """
    Standardize and audit all Cancer-tRFQTL supplementary tables.

    S2 and S10 are expected to already exist from the locked core
    standardization pipeline.
    """

    EXPECTED_SHEETS = tuple(
        f"S{index}"
        for index in range(
            1,
            12,
        )
    )

    LOCKED_CORE_SHEETS = (
        "S2",
        "S10",
    )

    BACKFILL_SHEETS = (
        "S1",
        "S3",
        "S4",
        "S5",
        "S6",
        "S7",
        "S8",
        "S9",
        "S11",
    )

    SUMMARY_FILENAME = (
        "cancer_trfqtl_all_tables.json"
    )

    def __init__(
        self,
        *,
        inputs: TRFQTLAllInputs,
        output_directory: str | Path,
        qc_directory: str | Path | None = None,
    ) -> None:
        """Initialize complete Cancer-tRFQTL standardization runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

        if qc_directory is None:
            self.qc_directory = (
                self.output_directory
                .parent
                / "qc"
            )
        else:
            self.qc_directory = Path(
                qc_directory
            )

    # ------------------------------------------------------------------
    # Input checks
    # ------------------------------------------------------------------

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require an existing input file."""

        if not path.exists():
            raise FileNotFoundError(
                "Cancer-tRFQTL workbook not found: "
                f"{path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                "Cancer-tRFQTL workbook path is not a file: "
                f"{path}"
            )

    @classmethod
    def _inspect_workbook_sheets(
        cls,
        workbook: Path,
    ) -> list[str]:
        """Return workbook sheet names and verify S1-S11 presence."""

        with pd.ExcelFile(
            workbook
        ) as excel_file:
            sheet_names = list(
                excel_file.sheet_names
            )

        missing = [
            sheet
            for sheet
            in cls.EXPECTED_SHEETS
            if sheet
            not in sheet_names
        ]

        if missing:
            raise ValueError(
                "Cancer-tRFQTL workbook is missing expected sheets: "
                f"{missing}. Available sheets: {sheet_names}"
            )

        return sheet_names

    # ------------------------------------------------------------------
    # Workbook reading
    # ------------------------------------------------------------------

    @staticmethod
    def _read_sheet(
        workbook: Path,
        sheet: str,
    ) -> pd.DataFrame:
        """
        Read a Cancer-tRFQTL supplementary sheet.

        Header row 1 means the second Excel row is used as the actual
        dataframe header.
        """

        return pd.read_excel(
            workbook,
            sheet_name=sheet,
            header=1,
        )

    # ------------------------------------------------------------------
    # QC helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_nunique(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count unique non-missing values."""

        if (
            column
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                column
            ]
            .dropna()
            .nunique()
        )

    @staticmethod
    def _nonmissing_count(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count non-missing values."""

        if (
            column
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                column
            ]
            .notna()
            .sum()
        )

    @staticmethod
    def _numeric_usable_count(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count non-missing standardized numeric values."""

        if (
            column
            not in dataframe.columns
        ):
            return 0

        return int(
            pd.to_numeric(
                dataframe[
                    column
                ],
                errors="coerce",
            )
            .notna()
            .sum()
        )

    @classmethod
    def _table_qc(
        cls,
        dataframe: pd.DataFrame,
        *,
        sheet: str,
        input_rows: int,
        output_path: Path,
        locked_existing: bool,
    ) -> dict[str, Any]:
        """Build sheet-level standardization QC."""

        report: dict[
            str,
            Any,
        ] = {
            "sheet":
                sheet,

            "record_type":
                (
                    SHEET_REGISTRY[
                        sheet
                    ].record_type
                    if sheet
                    in SHEET_REGISTRY
                    else (
                        "gwas_associated_trfqtl"
                        if sheet == "S2"
                        else "cancer_risk_variant"
                    )
                ),

            "input_rows":
                int(
                    input_rows
                ),

            "output_rows":
                int(
                    len(
                        dataframe
                    )
                ),

            "cardinality_preserved":
                bool(
                    input_rows
                    == len(
                        dataframe
                    )
                ),

            "output":
                str(
                    output_path
                ),

            "columns":
                int(
                    len(
                        dataframe.columns
                    )
                ),

            "locked_existing":
                bool(
                    locked_existing
                ),
        }

        if (
            "disease_raw"
            in dataframe.columns
        ):
            report[
                "disease_nonmissing_rows"
            ] = cls._nonmissing_count(
                dataframe,
                "disease_raw",
            )

            report[
                "unique_disease_labels"
            ] = cls._safe_nunique(
                dataframe,
                "disease_raw",
            )

        if (
            "variant_raw"
            in dataframe.columns
        ):
            report[
                "variant_nonmissing_rows"
            ] = cls._nonmissing_count(
                dataframe,
                "variant_raw",
            )

            report[
                "unique_variants"
            ] = cls._safe_nunique(
                dataframe,
                "variant_raw",
            )

        if (
            "trf_raw"
            in dataframe.columns
        ):
            report[
                "trf_nonmissing_rows"
            ] = cls._nonmissing_count(
                dataframe,
                "trf_raw",
            )

            report[
                "unique_trfs"
            ] = cls._safe_nunique(
                dataframe,
                "trf_raw",
            )

        if (
            "gene_raw"
            in dataframe.columns
        ):
            report[
                "gene_nonmissing_rows"
            ] = cls._nonmissing_count(
                dataframe,
                "gene_raw",
            )

            report[
                "unique_genes"
            ] = cls._safe_nunique(
                dataframe,
                "gene_raw",
            )

        for column in (
            "p_value",
            "survival_p_value",
            "correlation_p_value",
            "log2_fold_change",
            "effect_value",
        ):
            if (
                column
                in dataframe.columns
            ):
                report[
                    f"{column}_usable_rows"
                ] = (
                    cls._numeric_usable_count(
                        dataframe,
                        column,
                    )
                )

        return report

    # ------------------------------------------------------------------
    # Locked S2/S10 audit
    # ------------------------------------------------------------------

    def _audit_locked_table(
        self,
        *,
        workbook: Path,
        sheet: str,
    ) -> dict[str, Any]:
        """
        Audit a previously standardized core table without rewriting it.
        """

        output_path = (
            self.output_directory
            / f"{sheet.lower()}.parquet"
        )

        if not output_path.exists():
            raise FileNotFoundError(
                f"Locked Cancer-tRFQTL {sheet} standardized output "
                f"was not found: {output_path}. "
                "Run the existing core S2/S10 standardization first."
            )

        raw = (
            CancerTRFQTLSupplementaryStandardizer
            .clean_source_dataframe(
                self._read_sheet(
                    workbook,
                    sheet,
                )
            )
        )

        standardized = read_parquet(
            output_path
        )

        if (
            len(
                raw
            )
            != len(
                standardized
            )
        ):
            raise RuntimeError(
                f"Locked Cancer-tRFQTL {sheet} cardinality mismatch: "
                f"workbook={len(raw)}, standardized={len(standardized)}"
            )

        return self._table_qc(
            standardized,
            sheet=sheet,
            input_rows=len(
                raw
            ),
            output_path=output_path,
            locked_existing=True,
        )

    # ------------------------------------------------------------------
    # Backfill sheet processing
    # ------------------------------------------------------------------

    def _standardize_backfill_sheet(
        self,
        *,
        workbook: Path,
        sheet: str,
    ) -> dict[str, Any]:
        """Standardize one non-core Cancer-tRFQTL sheet."""

        raw = self._read_sheet(
            workbook,
            sheet,
        )

        cleaned = (
            CancerTRFQTLSupplementaryStandardizer
            .clean_source_dataframe(
                raw
            )
        )

        input_rows = len(
            cleaned
        )

        standardized = (
            CancerTRFQTLSupplementaryStandardizer
            .standardize(
                cleaned,
                sheet=sheet,
            )
        )

        if (
            len(
                standardized
            )
            != input_rows
        ):
            raise RuntimeError(
                f"Cancer-tRFQTL {sheet} cardinality changed: "
                f"{input_rows} -> {len(standardized)}"
            )

        output_path = (
            self.output_directory
            / f"{sheet.lower()}.parquet"
        )

        write_parquet(
            standardized,
            output_path,
            index=False,
        )

        logger.info(
            "Cancer-tRFQTL %s standardized: %d rows -> %s",
            sheet,
            len(
                standardized
            ),
            output_path,
        )

        return self._table_qc(
            standardized,
            sheet=sheet,
            input_rows=input_rows,
            output_path=output_path,
            locked_existing=False,
        )

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute complete Cancer-tRFQTL supplementary standardization."""

        workbook = Path(
            self.inputs.workbook
        )

        self._require_file(
            workbook
        )

        workbook_sheets = (
            self._inspect_workbook_sheets(
                workbook
            )
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        logger.info(
            "Starting complete Cancer-tRFQTL supplementary "
            "standardization."
        )

        logger.info(
            "Workbook sheets detected: %d",
            len(
                workbook_sheets
            ),
        )

        tables: dict[
            str,
            Any,
        ] = {}

        # --------------------------------------------------------------
        # Process S1-S11 in natural order.
        # --------------------------------------------------------------

        for sheet in self.EXPECTED_SHEETS:

            if sheet in self.LOCKED_CORE_SHEETS:

                logger.info(
                    "Auditing locked Cancer-tRFQTL %s.",
                    sheet,
                )

                tables[
                    sheet
                ] = self._audit_locked_table(
                    workbook=workbook,
                    sheet=sheet,
                )

                continue

            logger.info(
                "Standardizing Cancer-tRFQTL %s.",
                sheet,
            )

            tables[
                sheet
            ] = self._standardize_backfill_sheet(
                workbook=workbook,
                sheet=sheet,
            )

        # --------------------------------------------------------------
        # Global invariants.
        # --------------------------------------------------------------

        all_expected_outputs_exist = all(
            (
                self.output_directory
                / f"{sheet.lower()}.parquet"
            ).exists()
            for sheet
            in self.EXPECTED_SHEETS
        )

        if not all_expected_outputs_exist:
            raise RuntimeError(
                "Complete Cancer-tRFQTL standardization finished but "
                "one or more expected S1-S11 Parquet outputs are missing."
            )

        all_cardinality_preserved = all(
            bool(
                table[
                    "cardinality_preserved"
                ]
            )
            for table
            in tables.values()
        )

        if not all_cardinality_preserved:
            raise RuntimeError(
                "Cancer-tRFQTL all-table standardization detected a "
                "cardinality mismatch."
            )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M2-backfill",

            "stage":
                "cancer_trfqtl_all_tables_standardized",

            "dataset":
                DATASET_ID,

            "input": {
                "workbook":
                    str(
                        workbook
                    ),

                "expected_sheet_count":
                    len(
                        self.EXPECTED_SHEETS
                    ),

                "workbook_sheet_count":
                    len(
                        workbook_sheets
                    ),

                "expected_sheets":
                    list(
                        self.EXPECTED_SHEETS
                    ),

                "workbook_sheets":
                    workbook_sheets,
            },

            "policy": {
                "raw_file_modified":
                    False,

                "source_values_corrected":
                    False,

                "empty_spreadsheet_rows_removed":
                    True,

                "source_columns_preserved":
                    True,

                "disease_filtering_performed":
                    False,

                "prad_filtering_performed":
                    False,

                "liftover_performed":
                    False,

                "cross_source_join_performed":
                    False,

                "statistical_filtering_performed":
                    False,

                "candidate_ranking_performed":
                    False,

                "s2_rewritten":
                    False,

                "s10_rewritten":
                    False,
            },

            "tables":
                tables,

            "qc": {
                "all_11_sheets_present":
                    all(
                        sheet
                        in workbook_sheets
                        for sheet
                        in self.EXPECTED_SHEETS
                    ),

                "all_11_outputs_exist":
                    all_expected_outputs_exist,

                "all_cardinality_preserved":
                    all_cardinality_preserved,

                "locked_core_tables":
                    list(
                        self.LOCKED_CORE_SHEETS
                    ),

                "newly_standardized_tables":
                    list(
                        self.BACKFILL_SHEETS
                    ),
            },
        }

        summary_path = (
            self.qc_directory
            / self.SUMMARY_FILENAME
        )

        report[
            "report_path"
        ] = str(
            summary_path
        )

        with summary_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "Cancer-tRFQTL complete standardization finished."
        )

        logger.info(
            "Cancer-tRFQTL standardized tables available: %d",
            len(
                tables
            ),
        )

        logger.info(
            "QC report: %s",
            summary_path,
        )

        return report