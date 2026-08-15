"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/runners/harmonize_trfqtl_all.py

Description:
    Backfill harmonization runner for all Cancer-tRFQTL supplementary
    tables.

    S2 and S10 are locked core harmonized tables and are audited only.
    S1, S3-S9, and S11 are harmonized from their standardized Parquet
    representations.

    The runner is schema-compatible with both the supplementary
    harmonization fields and the existing locked core harmonization
    fields.

    No liftover, disease filtering, deduplication, statistical
    filtering, cross-source joining, or candidate ranking is performed.

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

from pcatrfqtl.harmonization.trfqtl_supplementary import (
    CancerTRFQTLSupplementaryHarmonizer,
    DATASET_ID,
    SHEET_REGISTRY,
)
from pcatrfqtl.io.parquet import (
    read_parquet,
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


@dataclass(frozen=True)
class TRFQTLAllHarmonizationInputs:
    """Input directory for standardized Cancer-tRFQTL tables."""

    standardized_directory: Path


class CancerTRFQTLAllHarmonizationRunner:
    """Harmonize the complete Cancer-tRFQTL supplementary resource."""

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
        inputs: TRFQTLAllHarmonizationInputs,
        output_directory: str | Path,
        qc_directory: str | Path | None = None,
    ) -> None:
        """Initialize all-table Cancer-tRFQTL harmonization."""

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
    # Paths
    # ------------------------------------------------------------------

    def _standardized_path(
        self,
        sheet: str,
    ) -> Path:
        """Return standardized table path."""

        return (
            Path(
                self.inputs.standardized_directory
            )
            / f"{sheet.lower()}.parquet"
        )

    def _harmonized_path(
        self,
        sheet: str,
    ) -> Path:
        """Return harmonized table path."""

        return (
            self.output_directory
            / f"{sheet.lower()}.parquet"
        )

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require an existing file."""

        if not path.exists():
            raise FileNotFoundError(
                f"Required Cancer-tRFQTL file not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"Expected file but found non-file path: {path}"
            )

    # ------------------------------------------------------------------
    # Generic dataframe helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _first_existing_column(
        dataframe: pd.DataFrame,
        columns: tuple[str, ...],
    ) -> str | None:
        """
        Return the first schema-compatible column that exists.

        This allows the backfill QC runner to audit both newer
        supplementary harmonization outputs and previously locked core
        S2/S10 outputs without rewriting their schemas.
        """

        for column in columns:
            if column in dataframe.columns:
                return column

        return None

    @staticmethod
    def _count_true(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count explicit true values in a boolean-like column."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[
                column
            ]
            .fillna(False)
            .eq(True)
            .sum()
        )

    @staticmethod
    def _count_nonmissing(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count non-missing values."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[
                column
            ]
            .notna()
            .sum()
        )

    @staticmethod
    def _value_counts(
        dataframe: pd.DataFrame,
        column: str,
    ) -> dict[str, int]:
        """Return JSON-safe value counts."""

        if column not in dataframe.columns:
            return {}

        counts = (
            dataframe[
                column
            ]
            .dropna()
            .astype(str)
            .value_counts()
        )

        return {
            str(
                key
            ): int(
                value
            )
            for (
                key,
                value,
            ) in counts.items()
        }

    # ------------------------------------------------------------------
    # QC field resolution
    # ------------------------------------------------------------------

    @classmethod
    def _resolve_variant_rsid_usable_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """
        Resolve rsID usability column across harmonization schemas.
        """

        return cls._first_existing_column(
            dataframe,
            (
                "harm_variant_rsid_usable",
                "harm_rsid_usable",
                "variant_rsid_usable",
            ),
        )

    @classmethod
    def _resolve_coordinate_join_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """
        Resolve coordinate-join eligibility column across schemas.

        Locked core S2/S10 may use a field name different from the
        supplementary harmonizer.
        """

        return cls._first_existing_column(
            dataframe,
            (
                "harm_variant_coordinate_join_allowed",
                "harm_coordinate_join_allowed",
                "coordinate_join_allowed",
            ),
        )

    @classmethod
    def _resolve_variant_status_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """Resolve variant harmonization status column."""

        return cls._first_existing_column(
            dataframe,
            (
                "harm_variant_status",
                "variant_harmonization_status",
            ),
        )

    @classmethod
    def _resolve_feature_usable_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """Resolve feature-usability column."""

        return cls._first_existing_column(
            dataframe,
            (
                "harm_feature_usable",
                "feature_usable",
            ),
        )

    @classmethod
    def _resolve_feature_status_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """Resolve feature harmonization status column."""

        return cls._first_existing_column(
            dataframe,
            (
                "harm_feature_status",
                "feature_harmonization_status",
            ),
        )

    @classmethod
    def _resolve_feature_type_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """Resolve harmonized feature type column."""

        return cls._first_existing_column(
            dataframe,
            (
                "harm_feature_type",
                "feature_type",
            ),
        )

    @classmethod
    def _resolve_primary_disease_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """Resolve primary-disease boolean field."""

        return cls._first_existing_column(
            dataframe,
            (
                "harm_is_primary_disease",
                "is_primary_disease",
            ),
        )

    @classmethod
    def _resolve_disease_status_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str | None:
        """Resolve disease harmonization status field."""

        return cls._first_existing_column(
            dataframe,
            (
                "harm_disease_status",
                "disease_harmonization_status",
            ),
        )

    # ------------------------------------------------------------------
    # Sheet-level QC
    # ------------------------------------------------------------------

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
        """Create sheet-level harmonization QC."""

        report: dict[str, Any] = {
            "sheet":
                sheet,

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

        if sheet in SHEET_REGISTRY:
            report[
                "metadata_only"
            ] = bool(
                SHEET_REGISTRY[
                    sheet
                ].metadata_only
            )

        # --------------------------------------------------------------
        # Disease QC
        # --------------------------------------------------------------

        primary_disease_column = (
            cls._resolve_primary_disease_column(
                dataframe
            )
        )

        disease_status_column = (
            cls._resolve_disease_status_column(
                dataframe
            )
        )

        if primary_disease_column is not None:
            report[
                "primary_disease_rows"
            ] = cls._count_true(
                dataframe,
                primary_disease_column,
            )

        if disease_status_column is not None:
            report[
                "disease_status_counts"
            ] = cls._value_counts(
                dataframe,
                disease_status_column,
            )

        # --------------------------------------------------------------
        # Variant QC
        # --------------------------------------------------------------

        variant_rsid_column = (
            cls._first_existing_column(
                dataframe,
                (
                    "harm_variant_rsid",
                    "variant_rsid",
                ),
            )
        )

        if variant_rsid_column is not None:

            rsid_usable_column = (
                cls._resolve_variant_rsid_usable_column(
                    dataframe
                )
            )

            coordinate_join_column = (
                cls._resolve_coordinate_join_column(
                    dataframe
                )
            )

            variant_status_column = (
                cls._resolve_variant_status_column(
                    dataframe
                )
            )

            if rsid_usable_column is not None:
                report[
                    "variant_rsid_usable_rows"
                ] = cls._count_true(
                    dataframe,
                    rsid_usable_column,
                )
            else:
                report[
                    "variant_rsid_usable_rows"
                ] = None

                report[
                    "variant_rsid_qc_status"
                ] = "COLUMN_NOT_AVAILABLE"

            if coordinate_join_column is not None:
                report[
                    "coordinate_join_allowed_rows"
                ] = cls._count_true(
                    dataframe,
                    coordinate_join_column,
                )

                report[
                    "coordinate_join_qc_status"
                ] = "AVAILABLE"

                report[
                    "coordinate_join_source_column"
                ] = coordinate_join_column
            else:
                report[
                    "coordinate_join_allowed_rows"
                ] = None

                report[
                    "coordinate_join_qc_status"
                ] = "COLUMN_NOT_AVAILABLE"

            if variant_status_column is not None:
                report[
                    "variant_status_counts"
                ] = cls._value_counts(
                    dataframe,
                    variant_status_column,
                )

        # --------------------------------------------------------------
        # Feature QC
        # --------------------------------------------------------------

        feature_id_column = (
            cls._first_existing_column(
                dataframe,
                (
                    "harm_feature_id",
                    "feature_id",
                ),
            )
        )

        if feature_id_column is not None:

            feature_usable_column = (
                cls._resolve_feature_usable_column(
                    dataframe
                )
            )

            feature_status_column = (
                cls._resolve_feature_status_column(
                    dataframe
                )
            )

            feature_type_column = (
                cls._resolve_feature_type_column(
                    dataframe
                )
            )

            if feature_usable_column is not None:
                report[
                    "feature_usable_rows"
                ] = cls._count_true(
                    dataframe,
                    feature_usable_column,
                )

            if feature_status_column is not None:
                report[
                    "feature_status_counts"
                ] = cls._value_counts(
                    dataframe,
                    feature_status_column,
                )

            if feature_type_column is not None:
                report[
                    "feature_type_counts"
                ] = cls._value_counts(
                    dataframe,
                    feature_type_column,
                )

        return report

    # ------------------------------------------------------------------
    # Locked core tables
    # ------------------------------------------------------------------

    def _audit_locked_table(
        self,
        *,
        sheet: str,
    ) -> dict[str, Any]:
        """
        Audit S2/S10 without rewriting their harmonized outputs.

        Both standardized and harmonized tables must already exist.
        """

        standardized_path = (
            self._standardized_path(
                sheet
            )
        )

        harmonized_path = (
            self._harmonized_path(
                sheet
            )
        )

        self._require_file(
            standardized_path
        )

        self._require_file(
            harmonized_path
        )

        standardized = read_parquet(
            standardized_path
        )

        harmonized = read_parquet(
            harmonized_path
        )

        if (
            len(
                standardized
            )
            != len(
                harmonized
            )
        ):
            raise RuntimeError(
                f"Locked {sheet} harmonization cardinality mismatch: "
                f"{len(standardized)} -> {len(harmonized)}"
            )

        logger.info(
            "Locked Cancer-tRFQTL %s audit: %d rows.",
            sheet,
            len(
                harmonized
            ),
        )

        return self._table_qc(
            harmonized,
            sheet=sheet,
            input_rows=len(
                standardized
            ),
            output_path=harmonized_path,
            locked_existing=True,
        )

    # ------------------------------------------------------------------
    # Backfill tables
    # ------------------------------------------------------------------

    def _harmonize_backfill_sheet(
        self,
        *,
        sheet: str,
    ) -> dict[str, Any]:
        """Harmonize one supplementary table."""

        standardized_path = (
            self._standardized_path(
                sheet
            )
        )

        self._require_file(
            standardized_path
        )

        standardized = read_parquet(
            standardized_path
        )

        input_rows = len(
            standardized
        )

        harmonized = (
            CancerTRFQTLSupplementaryHarmonizer
            .harmonize(
                standardized,
                sheet=sheet,
            )
        )

        if (
            len(
                harmonized
            )
            != input_rows
        ):
            raise RuntimeError(
                f"Cancer-tRFQTL {sheet} harmonization changed "
                f"cardinality: {input_rows} -> {len(harmonized)}"
            )

        output_path = (
            self._harmonized_path(
                sheet
            )
        )

        write_parquet(
            harmonized,
            output_path,
            index=False,
        )

        logger.info(
            "Cancer-tRFQTL %s harmonized: %d rows -> %s",
            sheet,
            len(
                harmonized
            ),
            output_path,
        )

        return self._table_qc(
            harmonized,
            sheet=sheet,
            input_rows=input_rows,
            output_path=output_path,
            locked_existing=False,
        )

    # ------------------------------------------------------------------
    # Global QC
    # ------------------------------------------------------------------

    def _all_outputs_exist(
        self,
    ) -> bool:
        """Check whether harmonized S1-S11 outputs all exist."""

        return all(
            self._harmonized_path(
                sheet
            ).exists()
            for sheet
            in self.EXPECTED_SHEETS
        )

    @staticmethod
    def _all_cardinality_preserved(
        tables: dict[str, Any],
    ) -> bool:
        """Check cardinality across every audited/harmonized table."""

        return all(
            bool(
                table[
                    "cardinality_preserved"
                ]
            )
            for table
            in tables.values()
        )

    @staticmethod
    def _coordinate_qc_summary(
        tables: dict[str, Any],
    ) -> dict[str, Any]:
        """Summarize coordinate-join QC availability."""

        result: dict[
            str,
            Any,
        ] = {}

        for (
            sheet,
            table,
        ) in tables.items():

            if (
                "coordinate_join_allowed_rows"
                not in table
            ):
                continue

            result[
                sheet
            ] = {
                "rows":
                    table.get(
                        "coordinate_join_allowed_rows"
                    ),

                "status":
                    table.get(
                        "coordinate_join_qc_status"
                    ),

                "source_column":
                    table.get(
                        "coordinate_join_source_column"
                    ),
            }

        return result

    # ------------------------------------------------------------------
    # Runner
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute complete Cancer-tRFQTL supplementary harmonization."""

        standardized_directory = Path(
            self.inputs.standardized_directory
        )

        if not standardized_directory.exists():
            raise FileNotFoundError(
                "Cancer-tRFQTL standardized directory not found: "
                f"{standardized_directory}"
            )

        if not standardized_directory.is_dir():
            raise NotADirectoryError(
                "Cancer-tRFQTL standardized path is not a directory: "
                f"{standardized_directory}"
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
            "harmonization."
        )

        tables: dict[
            str,
            Any,
        ] = {}

        # --------------------------------------------------------------
        # Natural S1-S11 order
        # --------------------------------------------------------------

        for sheet in self.EXPECTED_SHEETS:

            if sheet in self.LOCKED_CORE_SHEETS:

                logger.info(
                    "Auditing locked harmonized Cancer-tRFQTL %s.",
                    sheet,
                )

                tables[
                    sheet
                ] = self._audit_locked_table(
                    sheet=sheet
                )

                continue

            logger.info(
                "Harmonizing Cancer-tRFQTL %s.",
                sheet,
            )

            tables[
                sheet
            ] = self._harmonize_backfill_sheet(
                sheet=sheet
            )

        # --------------------------------------------------------------
        # Global invariants
        # --------------------------------------------------------------

        all_outputs_exist = (
            self._all_outputs_exist()
        )

        all_cardinality_preserved = (
            self._all_cardinality_preserved(
                tables
            )
        )

        if not all_outputs_exist:
            raise RuntimeError(
                "One or more Cancer-tRFQTL harmonized outputs "
                "are missing."
            )

        if not all_cardinality_preserved:
            raise RuntimeError(
                "Cancer-tRFQTL harmonization cardinality invariant "
                "failed."
            )

        # --------------------------------------------------------------
        # Report
        # --------------------------------------------------------------

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M2-backfill",

            "stage":
                "cancer_trfqtl_all_tables_harmonized",

            "dataset":
                DATASET_ID,

            "inputs": {
                "standardized_directory":
                    str(
                        standardized_directory
                    ),
            },

            "output_directory":
                str(
                    self.output_directory
                ),

            "policy": {
                "liftover_performed":
                    False,

                "disease_filtering_performed":
                    False,

                "prad_filtering_performed":
                    False,

                "deduplication_performed":
                    False,

                "cross_source_join_performed":
                    False,

                "statistical_filtering_performed":
                    False,

                "candidate_ranking_performed":
                    False,

                "row_cardinality_preserved":
                    True,

                "metadata_tables_preserved":
                    True,

                "source_identifiers_preserved":
                    True,

                "s2_rewritten":
                    False,

                "s10_rewritten":
                    False,
            },

            "tables":
                tables,

            "qc": {
                "all_11_outputs_exist":
                    all_outputs_exist,

                "all_cardinality_preserved":
                    all_cardinality_preserved,

                "locked_core_tables":
                    list(
                        self.LOCKED_CORE_SHEETS
                    ),

                "newly_harmonized_tables":
                    list(
                        self.BACKFILL_SHEETS
                    ),

                "primary_disease_rows_by_sheet": {
                    sheet: int(
                        table.get(
                            "primary_disease_rows",
                            0,
                        )
                        or 0
                    )
                    for (
                        sheet,
                        table,
                    ) in tables.items()
                },

                "coordinate_join_by_sheet":
                    self._coordinate_qc_summary(
                        tables
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
            "Cancer-tRFQTL complete harmonization finished."
        )

        logger.info(
            "Harmonized tables available: %d",
            len(
                tables
            ),
        )

        logger.info(
            "QC report: %s",
            summary_path,
        )

        return report