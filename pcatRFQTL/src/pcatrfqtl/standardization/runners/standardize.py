"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/runners/standardize.py

Description:
    Dataset-level execution runner for M2 data standardization.

    Supported datasets:

        - GWAS Catalog association data
        - Cancer-tRFQTL S2 and S10
        - Moradi QTL S1-S4
        - Moradi differential-expression S5-S10

    Responsibilities:

        - read validated source datasets
        - execute source-specific standardizers
        - process GWAS Catalog in chunks
        - verify row cardinality
        - prepare serialization-safe Parquet representations
        - write standardized Parquet outputs
        - collect output-level statistics
        - generate a JSON standardization summary

    Scientific transformation rules are implemented separately in:

        pcatrfqtl.standardization.gwas
        pcatrfqtl.standardization.trfqtl
        pcatrfqtl.standardization.moradi_qtl
        pcatrfqtl.standardization.moradi_de

    Important:
        Parquet serialization may require physical-type coercion for
        mixed object columns originating from Excel files. Such
        coercion affects only the standardized storage representation.

        Raw source files are never modified.

    This runner does not:

        - perform genome-build liftover
        - perform cross-source harmonization
        - filter datasets to prostate cancer
        - remove duplicate source records
        - silently correct malformed source identifiers

    Dataset locations should ultimately be resolved from
    configs/datasets.yaml by the CLI or configuration layer.

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
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from pcatrfqtl.logging.logger import get_logger
from pcatrfqtl.standardization.gwas import GWASStandardizer
from pcatrfqtl.standardization.moradi_de import MoradiDEStandardizer
from pcatrfqtl.standardization.moradi_qtl import MoradiQTLStandardizer
from pcatrfqtl.standardization.trfqtl import TRFQTLStandardizer


logger = get_logger(__name__)


@dataclass(frozen=True)
class StandardizationInputs:
    """
    Source paths required for M2 dataset-level standardization.

    These paths are execution inputs only. Dataset configuration remains
    centralized in configs/datasets.yaml.
    """

    gwas_catalog: Path
    cancer_trfqtl_workbook: Path
    moradi_directory: Path


class StandardizationRunner:
    """Execute dataset-level standardization for all validated sources."""

    DEFAULT_GWAS_CHUNK_SIZE = 100_000

    # ==================================================================
    # Cancer-tRFQTL
    # ==================================================================

    TRFQTL_TABLES: dict[str, dict[str, str]] = {
        "S2": {
            "output": "s2.parquet",
        },
        "S10": {
            "output": "s10.parquet",
        },
    }

    # ==================================================================
    # Moradi S1
    # ==================================================================

    MORADI_S1: dict[str, dict[str, str]] = {
        "Cis-intron-retention-sQTL-0.05": {
            "output": "s1_cis_intron.parquet",
            "feature_type": "intron",
        },
        "Cis-cassette-exon-sQTL-0.05": {
            "output": "s1_cis_exon.parquet",
            "feature_type": "exon",
        },
        "Cis-iso-eQTL-0.05": {
            "output": "s1_cis_isoform.parquet",
            "feature_type": "transcript",
        },
    }

    # ==================================================================
    # Moradi S2
    # ==================================================================

    MORADI_S2: dict[str, dict[str, str]] = {
        "Cis-intron-sQTL-GWAS-LD_0.5": {
            "output": "s2_cis_intron.parquet",
            "feature_type": "intron",
            "feature_column": "splicing_event",
        },
        "Cis-exon-sQTL-GWAS-LD_0.5": {
            "output": "s2_cis_exon.parquet",
            "feature_type": "exon",
            "feature_column": "splicing_event",
        },
        "Cis-iso-eQTl-GWAS-LD_0.5": {
            "output": "s2_cis_isoform.parquet",
            "feature_type": "transcript",
            "feature_column": "mRNA_isoform",
        },
    }

    # ==================================================================
    # Moradi S3
    # ==================================================================

    MORADI_S3: dict[str, dict[str, str]] = {
        "Trans-intron-sQTL": {
            "output": "s3_trans_intron.parquet",
            "feature_type": "intron",
            "feature_column": "splicing_event",
            "variant_column": "SNP",
        },
        "Trans-exon-sQTl": {
            "output": "s3_trans_exon.parquet",
            "feature_type": "exon",
            "feature_column": "splicing_event",
            "variant_column": "T",
        },
        "Trans-iso-eQTL": {
            "output": "s3_trans_isoform.parquet",
            "feature_type": "transcript",
            "feature_column": "splicing_event",
            "variant_column": "SNP",
        },
    }

    # ==================================================================
    # Moradi S4
    # ==================================================================

    MORADI_S4: dict[str, dict[str, str]] = {
        "Trans-intron-sQTL-LD_0.5": {
            "output": "s4_trans_intron.parquet",
            "feature_type": "transcript",
            "feature_column": "mRNA_isoform",
        },
        "Trans-exon-sQTL-LD_0.5": {
            "output": "s4_trans_exon.parquet",
            "feature_type": "exon",
            "feature_column": "splicing_event",
        },
        "Trans-iso-eQTL-LD_0.5": {
            "output": "s4_trans_isoform.parquet",
            "feature_type": "intron",
            "feature_column": "splicing_event",
        },
    }

    # ==================================================================
    # Moradi DE S5-S10
    # ==================================================================

    MORADI_DE: dict[str, dict[str, str]] = {
        "S5": {
            "filename": "Supplementary Table S5_DE_EP_cis.csv",
            "output": "s5.parquet",
            "feature_type": "exon",
            "analysis_scope": "cis",
        },
        "S6": {
            "filename": "Supplementary Table S6_DE_EP_trans.csv",
            "output": "s6.parquet",
            "feature_type": "exon",
            "analysis_scope": "trans",
        },
        "S7": {
            "filename": "Supplementary Table S7_DE_INT_cis.csv",
            "output": "s7.parquet",
            "feature_type": "intron",
            "analysis_scope": "cis",
        },
        "S8": {
            "filename": "Supplementary Table S8_DE_INT_trans.csv",
            "output": "s8.parquet",
            "feature_type": "intron",
            "analysis_scope": "trans",
        },
        "S9": {
            "filename": "Supplementary Table S9_DE_Isoform_cis.csv",
            "output": "s9.parquet",
            "feature_type": "transcript",
            "analysis_scope": "cis",
        },
        "S10": {
            "filename": (
                "Supplementary Table "
                "S10_DE_Isoform_trans_cancer_vs_control.csv"
            ),
            "output": "s10.parquet",
            "feature_type": "transcript",
            "analysis_scope": "trans",
        },
    }

    def __init__(
        self,
        inputs: StandardizationInputs,
        output_directory: str | Path,
        *,
        gwas_chunk_size: int = DEFAULT_GWAS_CHUNK_SIZE,
    ) -> None:
        """Initialize the M2 standardization runner."""

        self.inputs = inputs
        self.output_directory = Path(output_directory)
        self.gwas_chunk_size = int(gwas_chunk_size)

        if self.gwas_chunk_size <= 0:
            raise ValueError(
                "GWAS chunk size must be greater than zero."
            )

    # ==================================================================
    # Generic filesystem helpers
    # ==================================================================

    @staticmethod
    def _ensure_file(
        path: Path,
    ) -> None:
        """Ensure that a required source file exists."""

        if not path.exists():
            raise FileNotFoundError(
                f"Required source file not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"Expected a source file but found: {path}"
            )

    @staticmethod
    def _ensure_directory(
        path: Path,
    ) -> None:
        """Ensure that a required source directory exists."""

        if not path.exists():
            raise FileNotFoundError(
                f"Required source directory not found: {path}"
            )

        if not path.is_dir():
            raise NotADirectoryError(
                f"Expected a directory but found: {path}"
            )

    @staticmethod
    def _prepare_directory(
        path: Path,
        *,
        clear: bool = False,
    ) -> None:
        """Create an output directory, optionally replacing it."""

        if clear and path.exists():
            shutil.rmtree(
                path
            )

        path.mkdir(
            parents=True,
            exist_ok=True,
        )

    # ==================================================================
    # Parquet serialization
    # ==================================================================

    @staticmethod
    def _prepare_parquet_dataframe(
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Prepare a DataFrame for stable PyArrow serialization.

        Excel and heterogeneous source files can produce pandas object
        columns containing mixed Python scalar types. For example:

            snp_id_raw:
                "rs123"
                6
                None

        PyArrow requires one compatible physical type per column.

        Mixed scalar object columns are therefore represented as pandas
        nullable strings in standardized Parquet outputs.

        List-like object columns are retained unchanged because PyArrow
        can serialize homogeneous list values as nested Arrow arrays.

        This method never modifies the input DataFrame in-place.
        """

        prepared = dataframe.copy()

        for column in prepared.columns:
            series = prepared[column]

            if series.dtype != object:
                continue

            non_missing = (
                series[
                    series.notna()
                ]
            )

            if non_missing.empty:
                prepared[column] = (
                    series.astype(
                        "string"
                    )
                )
                continue

            values = (
                non_missing.tolist()
            )

            contains_container = any(
                isinstance(
                    value,
                    (
                        list,
                        tuple,
                        dict,
                        set,
                    ),
                )
                for value in values
            )

            # Nested/list fields such as tag_rsids must remain nested.
            if contains_container:
                continue

            python_types = {
                type(value)
                for value in values
            }

            # Mixed scalar object columns are unsafe for PyArrow type
            # inference. Convert only the standardized storage
            # representation to nullable strings.
            if len(python_types) > 1:
                prepared[column] = (
                    series.astype(
                        "string"
                    )
                )

        return prepared

    @classmethod
    def _write_parquet(
        cls,
        dataframe: pd.DataFrame,
        path: Path,
    ) -> None:
        """Write one standardized DataFrame as Parquet."""

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        prepared = (
            cls
            ._prepare_parquet_dataframe(
                dataframe
            )
        )

        prepared.to_parquet(
            path,
            index=False,
            engine="pyarrow",
        )

    # ==================================================================
    # Summary helpers
    # ==================================================================

    @staticmethod
    def _status_counts(
        dataframe: pd.DataFrame,
    ) -> dict[str, int]:
        """Count row-level standardization status values."""

        column = "standardization_status"

        if column not in dataframe.columns:
            return {}

        counts = (
            dataframe[column]
            .fillna("UNKNOWN")
            .astype(str)
            .value_counts(
                dropna=False
            )
        )

        return {
            str(status): int(count)
            for status, count in counts.items()
        }

    @staticmethod
    def _boolean_count(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count explicit True values."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[column]
            .fillna(False)
            .eq(True)
            .sum()
        )

    @staticmethod
    def _verify_cardinality(
        *,
        dataset_name: str,
        input_rows: int,
        output_rows: int,
    ) -> None:
        """Ensure standardization never silently changes row count."""

        if input_rows != output_rows:
            raise RuntimeError(
                f"{dataset_name} standardization changed row "
                f"cardinality: {input_rows} -> {output_rows}"
            )

    @classmethod
    def _dataframe_summary(
        cls,
        *,
        dataframe: pd.DataFrame,
        output: Path,
        usable_columns: tuple[str, ...] = (),
    ) -> dict[str, Any]:
        """Generate common summary metadata for one output."""

        summary: dict[str, Any] = {
            "output": str(output),
            "rows": int(
                len(
                    dataframe
                )
            ),
            "columns": int(
                len(
                    dataframe.columns
                )
            ),
            "status_counts": (
                cls._status_counts(
                    dataframe
                )
            ),
        }

        for column in usable_columns:
            summary[
                f"{column}_rows"
            ] = cls._boolean_count(
                dataframe,
                column,
            )

        return summary

    # ==================================================================
    # GWAS Catalog
    # ==================================================================

    def standardize_gwas(
        self,
    ) -> dict[str, Any]:
        """
        Standardize the complete GWAS Catalog using chunked execution.

        Each source chunk is written as an independent Parquet part.
        """

        source = (
            self.inputs
            .gwas_catalog
        )

        self._ensure_file(
            source
        )

        output_directory = (
            self.output_directory
            / "gwas_catalog"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        logger.info(
            "Starting GWAS Catalog standardization: %s",
            source,
        )

        total_input_rows = 0
        total_output_rows = 0

        coordinate_usable_rows = 0
        variant_usable_rows = 0
        p_value_usable_rows = 0

        status_counts: dict[
            str,
            int,
        ] = {}

        parts = 0

        reader = pd.read_csv(
            source,
            sep="\t",
            chunksize=self.gwas_chunk_size,
            low_memory=False,
        )

        for (
            chunk_index,
            chunk,
        ) in enumerate(
            reader
        ):

            input_rows = len(
                chunk
            )

            standardized = (
                GWASStandardizer
                .standardize_dataframe(
                    chunk
                )
            )

            output_rows = len(
                standardized
            )

            self._verify_cardinality(
                dataset_name=(
                    f"GWAS chunk {chunk_index}"
                ),
                input_rows=input_rows,
                output_rows=output_rows,
            )

            output = (
                output_directory
                / f"part-{chunk_index:05d}.parquet"
            )

            self._write_parquet(
                standardized,
                output,
            )

            total_input_rows += (
                input_rows
            )

            total_output_rows += (
                output_rows
            )

            coordinate_usable_rows += (
                self._boolean_count(
                    standardized,
                    "coordinate_usable",
                )
            )

            variant_usable_rows += (
                self._boolean_count(
                    standardized,
                    "variant_usable",
                )
            )

            p_value_usable_rows += (
                self._boolean_count(
                    standardized,
                    "p_value_usable",
                )
            )

            current_status_counts = (
                self._status_counts(
                    standardized
                )
            )

            for (
                status,
                count,
            ) in current_status_counts.items():

                status_counts[
                    status
                ] = (
                    status_counts.get(
                        status,
                        0,
                    )
                    + count
                )

            parts += 1

            logger.info(
                "GWAS chunk %d standardized: %d rows",
                chunk_index + 1,
                output_rows,
            )

        self._verify_cardinality(
            dataset_name="GWAS Catalog",
            input_rows=total_input_rows,
            output_rows=total_output_rows,
        )

        logger.info(
            "GWAS Catalog standardization complete: "
            "%d rows across %d parts",
            total_output_rows,
            parts,
        )

        return {
            "dataset":
                "gwas_catalog",

            "source":
                str(
                    source
                ),

            "output":
                str(
                    output_directory
                ),

            "output_format":
                "parquet_dataset",

            "parts":
                parts,

            "input_rows":
                total_input_rows,

            "output_rows":
                total_output_rows,

            "cardinality_preserved":
                True,

            "status_counts":
                status_counts,

            "coordinate_usable_rows":
                coordinate_usable_rows,

            "variant_usable_rows":
                variant_usable_rows,

            "p_value_usable_rows":
                p_value_usable_rows,
        }

    # ==================================================================
    # Cancer-tRFQTL
    # ==================================================================

    def standardize_trfqtl(
        self,
    ) -> dict[str, Any]:
        """Standardize Cancer-tRFQTL S2 and S10."""

        source = (
            self.inputs
            .cancer_trfqtl_workbook
        )

        self._ensure_file(
            source
        )

        output_directory = (
            self.output_directory
            / "cancer_trfqtl"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        logger.info(
            "Starting Cancer-tRFQTL standardization: %s",
            source,
        )

        results: dict[
            str,
            Any,
        ] = {}

        # --------------------------------------------------------------
        # S2
        # --------------------------------------------------------------

        s2 = pd.read_excel(
            source,
            sheet_name="S2",
            header=1,
        )

        s2_standardized = (
            TRFQTLStandardizer
            .standardize_s2_dataframe(
                s2
            )
        )

        self._verify_cardinality(
            dataset_name="Cancer-tRFQTL S2",
            input_rows=len(
                s2
            ),
            output_rows=len(
                s2_standardized
            ),
        )

        s2_output = (
            output_directory
            / self.TRFQTL_TABLES[
                "S2"
            ][
                "output"
            ]
        )

        self._write_parquet(
            s2_standardized,
            s2_output,
        )

        results[
            "S2"
        ] = {
            "input_rows":
                int(
                    len(
                        s2
                    )
                ),

            "output_rows":
                int(
                    len(
                        s2_standardized
                    )
                ),

            **self._dataframe_summary(
                dataframe=s2_standardized,
                output=s2_output,
                usable_columns=(
                    "variant_usable",
                    "qtl_usable",
                ),
            ),
        }

        # --------------------------------------------------------------
        # S10
        # --------------------------------------------------------------

        s10 = pd.read_excel(
            source,
            sheet_name="S10",
            header=1,
        )

        s10_standardized = (
            TRFQTLStandardizer
            .standardize_s10_dataframe(
                s10
            )
        )

        self._verify_cardinality(
            dataset_name="Cancer-tRFQTL S10",
            input_rows=len(
                s10
            ),
            output_rows=len(
                s10_standardized
            ),
        )

        s10_output = (
            output_directory
            / self.TRFQTL_TABLES[
                "S10"
            ][
                "output"
            ]
        )

        self._write_parquet(
            s10_standardized,
            s10_output,
        )

        results[
            "S10"
        ] = {
            "input_rows":
                int(
                    len(
                        s10
                    )
                ),

            "output_rows":
                int(
                    len(
                        s10_standardized
                    )
                ),

            **self._dataframe_summary(
                dataframe=s10_standardized,
                output=s10_output,
                usable_columns=(
                    "variant_usable",
                    "effect_usable",
                ),
            ),
        }

        logger.info(
            "Cancer-tRFQTL standardization complete."
        )

        return {
            "dataset":
                "cancer_trfqtl",

            "source":
                str(
                    source
                ),

            "tables":
                results,
        }

    # ==================================================================
    # Moradi S1
    # ==================================================================

    def _standardize_moradi_s1(
        self,
        workbook: Path,
        output_directory: Path,
    ) -> dict[str, Any]:
        """Standardize all Moradi S1 sheets."""

        results: dict[
            str,
            Any,
        ] = {}

        for (
            sheet,
            specification,
        ) in self.MORADI_S1.items():

            logger.info(
                "Standardizing Moradi S1: %s",
                sheet,
            )

            dataframe = pd.read_excel(
                workbook,
                sheet_name=sheet,
            )

            standardized = (
                MoradiQTLStandardizer
                .standardize_s1_dataframe(
                    dataframe,
                    feature_type=(
                        specification[
                            "feature_type"
                        ]
                    ),
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Moradi S1 {sheet}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    standardized
                ),
            )

            output = (
                output_directory
                / specification[
                    "output"
                ]
            )

            self._write_parquet(
                standardized,
                output,
            )

            results[
                sheet
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            standardized
                        )
                    ),

                "feature_type":
                    specification[
                        "feature_type"
                    ],

                **self._dataframe_summary(
                    dataframe=standardized,
                    output=output,
                    usable_columns=(
                        "variant_usable",
                        "qtl_usable",
                    ),
                ),
            }

        return results

    # ==================================================================
    # Moradi S2 / S4
    # ==================================================================

    def _standardize_moradi_ld_workbook(
        self,
        workbook: Path,
        output_directory: Path,
        *,
        source_table: str,
        specifications: dict[
            str,
            dict[str, str],
        ],
        qtl_scope: str,
    ) -> dict[str, Any]:
        """Standardize Moradi S2 or S4 LD-linked QTL sheets."""

        results: dict[
            str,
            Any,
        ] = {}

        for (
            sheet,
            specification,
        ) in specifications.items():

            logger.info(
                "Standardizing Moradi %s: %s",
                source_table,
                sheet,
            )

            dataframe = pd.read_excel(
                workbook,
                sheet_name=sheet,
            )

            standardized = (
                MoradiQTLStandardizer
                .standardize_ld_dataframe(
                    dataframe,
                    source_table=source_table,
                    feature_type=(
                        specification[
                            "feature_type"
                        ]
                    ),
                    qtl_scope=qtl_scope,
                    feature_column=(
                        specification[
                            "feature_column"
                        ]
                    ),
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Moradi {source_table} {sheet}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    standardized
                ),
            )

            output = (
                output_directory
                / specification[
                    "output"
                ]
            )

            self._write_parquet(
                standardized,
                output,
            )

            results[
                sheet
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            standardized
                        )
                    ),

                "feature_type":
                    specification[
                        "feature_type"
                    ],

                "qtl_scope":
                    qtl_scope,

                **self._dataframe_summary(
                    dataframe=standardized,
                    output=output,
                    usable_columns=(
                        "qtl_variant_usable",
                        "tag_variant_usable",
                        "qtl_usable",
                    ),
                ),
            }

        return results

    # ==================================================================
    # Moradi S3
    # ==================================================================

    def _standardize_moradi_s3(
        self,
        workbook: Path,
        output_directory: Path,
    ) -> dict[str, Any]:
        """Standardize all Moradi S3 trans-QTL sheets."""

        results: dict[
            str,
            Any,
        ] = {}

        for (
            sheet,
            specification,
        ) in self.MORADI_S3.items():

            logger.info(
                "Standardizing Moradi S3: %s",
                sheet,
            )

            dataframe = pd.read_excel(
                workbook,
                sheet_name=sheet,
            )

            standardized = (
                MoradiQTLStandardizer
                .standardize_s3_dataframe(
                    dataframe,
                    feature_type=(
                        specification[
                            "feature_type"
                        ]
                    ),
                    feature_column=(
                        specification[
                            "feature_column"
                        ]
                    ),
                    variant_column=(
                        specification[
                            "variant_column"
                        ]
                    ),
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Moradi S3 {sheet}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    standardized
                ),
            )

            output = (
                output_directory
                / specification[
                    "output"
                ]
            )

            self._write_parquet(
                standardized,
                output,
            )

            results[
                sheet
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            standardized
                        )
                    ),

                "feature_type":
                    specification[
                        "feature_type"
                    ],

                "variant_column":
                    specification[
                        "variant_column"
                    ],

                **self._dataframe_summary(
                    dataframe=standardized,
                    output=output,
                    usable_columns=(
                        "variant_usable",
                        "qtl_usable",
                    ),
                ),
            }

        return results

    # ==================================================================
    # Moradi QTL S1-S4
    # ==================================================================

    def standardize_moradi_qtl(
        self,
    ) -> dict[str, Any]:
        """Standardize complete Moradi QTL datasets S1-S4."""

        source_directory = (
            self.inputs
            .moradi_directory
        )

        self._ensure_directory(
            source_directory
        )

        output_directory = (
            self.output_directory
            / "moradi_qtl"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        workbooks = {
            "S1": (
                source_directory
                / "Supplementary Table S1.xlsx"
            ),
            "S2": (
                source_directory
                / "Supplementary Table S2.xlsx"
            ),
            "S3": (
                source_directory
                / "Supplementary Table S3.xlsx"
            ),
            "S4": (
                source_directory
                / "Supplementary Table S4.xlsx"
            ),
        }

        for workbook in workbooks.values():
            self._ensure_file(
                workbook
            )

        logger.info(
            "Starting Moradi QTL standardization."
        )

        result = {
            "dataset":
                "moradi_qtl",

            "source_directory":
                str(
                    source_directory
                ),

            "genome_build":
                (
                    MoradiQTLStandardizer
                    .GENOME_BUILD
                ),

            "tables": {
                "S1": (
                    self
                    ._standardize_moradi_s1(
                        workbooks[
                            "S1"
                        ],
                        output_directory,
                    )
                ),

                "S2": (
                    self
                    ._standardize_moradi_ld_workbook(
                        workbooks[
                            "S2"
                        ],
                        output_directory,
                        source_table="S2",
                        specifications=(
                            self.MORADI_S2
                        ),
                        qtl_scope="cis",
                    )
                ),

                "S3": (
                    self
                    ._standardize_moradi_s3(
                        workbooks[
                            "S3"
                        ],
                        output_directory,
                    )
                ),

                "S4": (
                    self
                    ._standardize_moradi_ld_workbook(
                        workbooks[
                            "S4"
                        ],
                        output_directory,
                        source_table="S4",
                        specifications=(
                            self.MORADI_S4
                        ),
                        qtl_scope="trans",
                    )
                ),
            },
        }

        logger.info(
            "Moradi QTL standardization complete."
        )

        return result

    # ==================================================================
    # Moradi DE
    # ==================================================================

    def standardize_moradi_de(
        self,
    ) -> dict[str, Any]:
        """Standardize Moradi differential-expression S5-S10."""

        source_directory = (
            self.inputs
            .moradi_directory
        )

        self._ensure_directory(
            source_directory
        )

        output_directory = (
            self.output_directory
            / "moradi_de"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        logger.info(
            "Starting Moradi DE standardization."
        )

        results: dict[
            str,
            Any,
        ] = {}

        for (
            table,
            specification,
        ) in self.MORADI_DE.items():

            source = (
                source_directory
                / specification[
                    "filename"
                ]
            )

            self._ensure_file(
                source
            )

            logger.info(
                "Standardizing Moradi DE %s: %s",
                table,
                source.name,
            )

            dataframe = pd.read_csv(
                source,
                low_memory=False,
            )

            standardized = (
                MoradiDEStandardizer
                .standardize_dataframe(
                    dataframe,
                    source_table=table,
                    feature_type=(
                        specification[
                            "feature_type"
                        ]
                    ),
                    analysis_scope=(
                        specification[
                            "analysis_scope"
                        ]
                    ),
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Moradi DE {table}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    standardized
                ),
            )

            output = (
                output_directory
                / specification[
                    "output"
                ]
            )

            self._write_parquet(
                standardized,
                output,
            )

            results[
                table
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            standardized
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

                **self._dataframe_summary(
                    dataframe=standardized,
                    output=output,
                    usable_columns=(
                        "feature_usable",
                        "de_statistics_usable",
                        "de_record_usable",
                        "source_feature_anomaly",
                    ),
                ),
            }

        logger.info(
            "Moradi DE standardization complete."
        )

        return {
            "dataset":
                "moradi_de",

            "source_directory":
                str(
                    source_directory
                ),

            "tables":
                results,
        }

    # ==================================================================
    # Complete M2 execution
    # ==================================================================

    def run_all(
        self,
    ) -> dict[str, Any]:
        """Execute complete M2 dataset-level standardization."""

        self._prepare_directory(
            self.output_directory
        )

        logger.info(
            "Starting M2 — Data Standardized."
        )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M2",

            "stage":
                "data_standardized",

            "inputs": {
                key:
                    str(
                        value
                    )
                for (
                    key,
                    value,
                ) in asdict(
                    self.inputs
                ).items()
            },

            "output_directory":
                str(
                    self.output_directory
                ),

            "datasets":
                {},
        }

        report[
            "datasets"
        ][
            "gwas_catalog"
        ] = self.standardize_gwas()

        report[
            "datasets"
        ][
            "cancer_trfqtl"
        ] = self.standardize_trfqtl()

        report[
            "datasets"
        ][
            "moradi_qtl"
        ] = self.standardize_moradi_qtl()

        report[
            "datasets"
        ][
            "moradi_de"
        ] = self.standardize_moradi_de()

        qc_directory = (
            self.output_directory
            / "qc"
        )

        self._prepare_directory(
            qc_directory
        )

        report_path = (
            qc_directory
            / "standardization_summary.json"
        )

        report[
            "report_path"
        ] = str(
            report_path
        )

        with report_path.open(
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
            "M2 standardization complete."
        )

        logger.info(
            "Standardization report: %s",
            report_path,
        )

        return report