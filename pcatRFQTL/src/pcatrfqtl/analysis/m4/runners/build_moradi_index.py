"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/runners/build_moradi_index.py

Description:
    Runner for M4.1 Unified Moradi QTL Index construction.

    Twelve harmonized Moradi QTL tables are converted into a common
    schema and written as independent Parquet parts.

    Using multiple parts avoids loading the complete ~1.4 million-row
    regulatory resource into one large in-memory dataframe.

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

from pcatrfqtl.analysis.m4.moradi_index import (
    MORADI_INDEX_SPECS,
    MoradiQTLIndexBuilder,
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
class M41MoradiIndexInputs:
    """Input directory containing harmonized Moradi QTL tables."""

    harmonized_moradi_directory: Path


class M41MoradiIndexRunner:
    """Build M4.1 Unified Moradi QTL Index."""

    SUMMARY_FILENAME = (
        "m4_1_moradi_qtl_index_summary.json"
    )

    def __init__(
        self,
        *,
        inputs: M41MoradiIndexInputs,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:
        """Initialize M4.1 runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    @staticmethod
    def _count_true(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count True values."""

        return int(
            dataframe[
                column
            ]
            .fillna(False)
            .eq(True)
            .sum()
        )

    @staticmethod
    def _safe_nunique(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count unique nonmissing values."""

        return int(
            dataframe[
                column
            ]
            .dropna()
            .nunique()
        )

    def _input_path(
        self,
        table_id: str,
    ) -> Path:
        """Return harmonized Moradi input path."""

        return (
            Path(
                self.inputs.harmonized_moradi_directory
            )
            / f"{table_id}.parquet"
        )

    def _output_path(
        self,
        table_id: str,
    ) -> Path:
        """Return index part output path."""

        return (
            self.output_directory
            / f"{table_id}.parquet"
        )

    def run(
        self,
    ) -> dict[str, Any]:
        """Build all twelve unified Moradi index parts."""

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        logger.info(
            "Starting M4.1 Unified Moradi QTL Index."
        )

        tables: dict[
            str,
            Any,
        ] = {}

        total_input_rows = 0
        total_output_rows = 0

        total_qtl_ready = 0
        total_tag_ready = 0

        source_unique_qtl_rsids: set[str] = set()
        source_unique_tag_rsids: set[str] = set()

        for (
            table_id,
            spec,
        ) in MORADI_INDEX_SPECS.items():

            input_path = self._input_path(
                table_id
            )

            if not input_path.exists():
                raise FileNotFoundError(
                    f"M4.1 Moradi input not found: {input_path}"
                )

            logger.info(
                "Indexing Moradi table %s.",
                table_id,
            )

            source = read_parquet(
                input_path
            )

            indexed = (
                MoradiQTLIndexBuilder
                .build_table(
                    source,
                    spec=spec,
                )
            )

            output_path = self._output_path(
                table_id
            )

            write_parquet(
                indexed,
                output_path,
                index=False,
            )

            input_rows = len(
                source
            )

            output_rows = len(
                indexed
            )

            qtl_ready = self._count_true(
                indexed,
                "m4_direct_variant_match_ready",
            )

            tag_ready = self._count_true(
                indexed,
                "m4_tag_variant_match_ready",
            )

            unique_qtl_rsids = (
                indexed[
                    "qtl_rsid"
                ]
                .dropna()
                .astype(str)
                .unique()
            )

            unique_tag_rsids = (
                indexed[
                    "tag_rsid"
                ]
                .dropna()
                .astype(str)
                .unique()
            )

            source_unique_qtl_rsids.update(
                unique_qtl_rsids
            )

            source_unique_tag_rsids.update(
                unique_tag_rsids
            )

            total_input_rows += input_rows
            total_output_rows += output_rows

            total_qtl_ready += qtl_ready
            total_tag_ready += tag_ready

            tables[
                table_id
            ] = {
                "source_sheet":
                    spec.source_sheet,

                "regulatory_scope":
                    spec.regulatory_scope,

                "feature_class":
                    spec.feature_class,

                "has_tag_variant":
                    spec.has_tag_variant,

                "input_rows":
                    input_rows,

                "output_rows":
                    output_rows,

                "cardinality_preserved":
                    input_rows
                    == output_rows,

                "direct_variant_match_ready_rows":
                    qtl_ready,

                "tag_variant_match_ready_rows":
                    tag_ready,

                "unique_qtl_rsids":
                    self._safe_nunique(
                        indexed,
                        "qtl_rsid",
                    ),

                "unique_tag_rsids":
                    self._safe_nunique(
                        indexed,
                        "tag_rsid",
                    ),

                "unique_features":
                    self._safe_nunique(
                        indexed,
                        "feature_identity_key",
                    ),

                "output":
                    str(
                        output_path
                    ),
            }

            logger.info(
                "Moradi %s indexed: %d rows, "
                "%d direct-match ready, %d tag-match ready.",
                table_id,
                output_rows,
                qtl_ready,
                tag_ready,
            )

            del source
            del indexed

        all_cardinality_preserved = all(
            bool(
                table[
                    "cardinality_preserved"
                ]
            )
            for table
            in tables.values()
        )

        if (
            total_input_rows
            != total_output_rows
        ):
            raise RuntimeError(
                "M4.1 total cardinality mismatch: "
                f"{total_input_rows} -> {total_output_rows}"
            )

        if not all_cardinality_preserved:
            raise RuntimeError(
                "One or more M4.1 Moradi index parts "
                "changed source cardinality."
            )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M4.1",

            "stage":
                "unified_moradi_qtl_index",

            "input_directory":
                str(
                    self.inputs.harmonized_moradi_directory
                ),

            "output_directory":
                str(
                    self.output_directory
                ),

            "policy": {
                "source_rows_filtered":
                    False,

                "deduplication_performed":
                    False,

                "liftover_performed":
                    False,

                "ld_expansion_performed":
                    False,

                "colocalization_performed":
                    False,

                "candidate_ranking_performed":
                    False,

                "qtl_and_tag_variants_distinguished":
                    True,

                "cis_and_trans_distinguished":
                    True,

                "feature_classes_distinguished":
                    True,
            },

            "summary": {
                "source_tables":
                    len(
                        tables
                    ),

                "input_rows":
                    total_input_rows,

                "output_rows":
                    total_output_rows,

                "cardinality_preserved":
                    (
                        total_input_rows
                        == total_output_rows
                    ),

                "direct_variant_match_ready_rows":
                    total_qtl_ready,

                "tag_variant_match_ready_rows":
                    total_tag_ready,

                "unique_qtl_rsids":
                    len(
                        source_unique_qtl_rsids
                    ),

                "unique_tag_rsids":
                    len(
                        source_unique_tag_rsids
                    ),
            },

            "tables":
                tables,
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
            "M4.1 Unified Moradi QTL Index completed: "
            "%d rows across %d parts.",
            total_output_rows,
            len(
                tables
            ),
        )

        logger.info(
            "M4.1 QC report: %s",
            summary_path,
        )

        return report