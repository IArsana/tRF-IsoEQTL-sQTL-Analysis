"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/extract_gwas.py

Description:
    Dataset-level runner for M3.1 prostate-cancer GWAS extraction.

    M3.1 reads all harmonized GWAS Catalog partitions and selects
    association records whose harmonized disease context corresponds
    to prostate cancer.

    All selected prostate-cancer association rows are retained,
    including records that do not contain a usable canonical single
    rsID.

    Records are annotated for eligibility in the later M3.3 direct
    GWAS × tRF-QTL overlap analysis.

    M3.1 intentionally does not:

        - join GWAS records with tRF-QTL
        - deduplicate variants
        - collapse multiple GWAS studies
        - perform coordinate matching
        - perform genome liftover
        - impose a new significance threshold
        - rank candidate variants

    Parquet serialization is delegated to the shared pipeline I/O
    utility so heterogeneous preserved source columns remain compatible
    with PyArrow.

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

from pcatrfqtl.analysis.m3.gwas import (
    ProstateGWASSelector,
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
class M31GWASInputs:
    """Input paths required for M3.1."""

    harmonized_gwas_directory: Path


class M31GWASRunner:
    """
    Execute M3.1 prostate-cancer GWAS extraction.

    Input:
        harmonized GWAS Catalog partitions

    Output:
        data/processed/m3/prostate_gwas.parquet

    QC:
        data/processed/m3/qc/m3_1_gwas_summary.json
    """

    OUTPUT_FILENAME = (
        "prostate_gwas.parquet"
    )

    SUMMARY_FILENAME = (
        "m3_1_gwas_summary.json"
    )

    def __init__(
        self,
        inputs: M31GWASInputs,
        output_directory: str | Path,
    ) -> None:
        """Initialize the M3.1 runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

    # ==================================================================
    # Input validation
    # ==================================================================

    @staticmethod
    def _require_directory(
        path: Path,
    ) -> None:
        """Require a valid input directory."""

        if not path.exists():
            raise FileNotFoundError(
                "M3.1 input directory not found: "
                f"{path}"
            )

        if not path.is_dir():
            raise NotADirectoryError(
                "M3.1 expected directory but found: "
                f"{path}"
            )

    # ==================================================================
    # QC helpers
    # ==================================================================

    @staticmethod
    def _count_true(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count explicit True values."""

        if (
            column
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                column
            ]
            .fillna(
                False
            )
            .eq(
                True
            )
            .sum()
        )

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
    def _value_counts(
        dataframe: pd.DataFrame,
        column: str,
        *,
        top_n: int | None = None,
    ) -> dict[str, int]:
        """Return JSON-safe value counts."""

        if (
            column
            not in dataframe.columns
        ):
            return {}

        counts = (
            dataframe[
                column
            ]
            .dropna()
            .astype(
                str
            )
            .value_counts()
        )

        if (
            top_n
            is not None
        ):
            counts = (
                counts.head(
                    top_n
                )
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

    # ==================================================================
    # Execution
    # ==================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M3.1 prostate-cancer GWAS extraction."""

        source_directory = (
            self.inputs
            .harmonized_gwas_directory
        )

        self._require_directory(
            source_directory
        )

        parts = sorted(
            source_directory.glob(
                "part-*.parquet"
            )
        )

        if not parts:
            raise FileNotFoundError(
                "No harmonized GWAS partitions found in "
                f"{source_directory}"
            )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        qc_directory = (
            self.output_directory
            / "qc"
        )

        qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        logger.info(
            "Starting M3.1 prostate-cancer GWAS extraction."
        )

        logger.info(
            "GWAS partitions discovered: %d",
            len(
                parts
            ),
        )

        selected_parts: list[
            pd.DataFrame
        ] = []

        total_input_rows = 0
        total_selected_rows = 0

        # --------------------------------------------------------------
        # Process harmonized GWAS partitions
        # --------------------------------------------------------------

        for (
            index,
            path,
        ) in enumerate(
            parts,
            start=1,
        ):

            dataframe = (
                read_parquet(
                    path
                )
            )

            input_rows = len(
                dataframe
            )

            total_input_rows += (
                input_rows
            )

            selected = (
                ProstateGWASSelector
                .prepare(
                    dataframe
                )
            )

            selected[
                "m3_source_partition"
            ] = (
                path.name
            )

            selected_rows = len(
                selected
            )

            total_selected_rows += (
                selected_rows
            )

            selected_parts.append(
                selected
            )

            logger.info(
                "M3.1 part %d/%d: %d input -> %d prostate rows",
                index,
                len(
                    parts
                ),
                input_rows,
                selected_rows,
            )

        # --------------------------------------------------------------
        # Concatenate without deduplication
        # --------------------------------------------------------------

        if not selected_parts:
            raise RuntimeError(
                "M3.1 did not process any GWAS partitions."
            )

        prostate_gwas = (
            pd.concat(
                selected_parts,
                axis=0,
                ignore_index=True,
            )
        )

        if (
            len(
                prostate_gwas
            )
            != total_selected_rows
        ):
            raise RuntimeError(
                "M3.1 concatenation changed selected-row cardinality: "
                f"{total_selected_rows} -> "
                f"{len(prostate_gwas)}"
            )

        # --------------------------------------------------------------
        # Scientific invariants
        # --------------------------------------------------------------

        if not (
            prostate_gwas[
                "harm_is_primary_disease"
            ]
            .fillna(
                False
            )
            .eq(
                True
            )
            .all()
        ):
            raise RuntimeError(
                "M3.1 output contains records that are not marked as "
                "the primary disease."
            )

        if not (
            prostate_gwas[
                "harm_disease_id"
            ]
            .eq(
                ProstateGWASSelector
                .PRIMARY_DISEASE_ID
            )
            .all()
        ):
            raise RuntimeError(
                "M3.1 output contains a disease identity other than "
                f"{ProstateGWASSelector.PRIMARY_DISEASE_ID}."
            )

        # --------------------------------------------------------------
        # Output
        # --------------------------------------------------------------

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        write_parquet(
            prostate_gwas,
            output_path,
            index=False,
        )

        logger.info(
            "M3.1 Parquet output written: %s",
            output_path,
        )

        # --------------------------------------------------------------
        # QC statistics
        # --------------------------------------------------------------

        direct_overlap_eligible_rows = (
            self._count_true(
                prostate_gwas,
                "m3_direct_overlap_eligible",
            )
        )

        direct_overlap_ineligible_rows = (
            len(
                prostate_gwas
            )
            - direct_overlap_eligible_rows
        )

        unique_canonical_rsids = (
            self._safe_nunique(
                prostate_gwas,
                "harm_variant_rsid",
            )
        )

        unique_overlap_eligible_rsids = (
            self._safe_nunique(
                prostate_gwas.loc[
                    prostate_gwas[
                        "m3_direct_overlap_eligible"
                    ]
                    .fillna(
                        False
                    )
                    .eq(
                        True
                    )
                ],
                "harm_variant_rsid",
            )
        )

        duplicate_rsid_association_rows = (
            direct_overlap_eligible_rows
            - unique_overlap_eligible_rsids
        )

        # --------------------------------------------------------------
        # QC report
        # --------------------------------------------------------------

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M3.1",

            "stage":
                "prostate_cancer_gwas_subset",

            "input": {
                "directory":
                    str(
                        source_directory
                    ),

                "parts":
                    int(
                        len(
                            parts
                        )
                    ),

                "rows":
                    int(
                        total_input_rows
                    ),
            },

            "output": {
                "path":
                    str(
                        output_path
                    ),

                "rows":
                    int(
                        len(
                            prostate_gwas
                        )
                    ),

                "unique_canonical_rsids":
                    int(
                        unique_canonical_rsids
                    ),

                "direct_overlap_eligible_rows":
                    int(
                        direct_overlap_eligible_rows
                    ),

                "direct_overlap_ineligible_rows":
                    int(
                        direct_overlap_ineligible_rows
                    ),

                "unique_overlap_eligible_rsids":
                    int(
                        unique_overlap_eligible_rsids
                    ),

                "duplicate_rsid_association_rows":
                    int(
                        duplicate_rsid_association_rows
                    ),
            },

            "selection": {
                "primary_disease_id":
                    (
                        ProstateGWASSelector
                        .PRIMARY_DISEASE_ID
                    ),

                "disease_filtering_performed":
                    True,

                "rsid_filtering_performed":
                    False,

                "deduplication_performed":
                    False,

                "study_collapsing_performed":
                    False,

                "coordinate_matching_performed":
                    False,

                "liftover_performed":
                    False,

                "significance_filtering_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "qc": {
                "selected_rows_match_partition_sum":
                    (
                        len(
                            prostate_gwas
                        )
                        == total_selected_rows
                    ),

                "all_rows_primary_disease":
                    bool(
                        prostate_gwas[
                            "harm_is_primary_disease"
                        ]
                        .fillna(
                            False
                        )
                        .eq(
                            True
                        )
                        .all()
                    ),

                "all_rows_canonical_prostate_cancer":
                    bool(
                        prostate_gwas[
                            "harm_disease_id"
                        ]
                        .eq(
                            ProstateGWASSelector
                            .PRIMARY_DISEASE_ID
                        )
                        .all()
                    ),

                "disease_ids":
                    self._value_counts(
                        prostate_gwas,
                        "harm_disease_id",
                    ),

                "variant_status_counts":
                    self._value_counts(
                        prostate_gwas,
                        "harm_variant_status",
                    ),

                "variant_identity_method_counts":
                    self._value_counts(
                        prostate_gwas,
                        "harm_variant_identity_method",
                    ),

                "top_source_disease_labels":
                    self._value_counts(
                        prostate_gwas,
                        "harm_disease_source_value",
                        top_n=20,
                    ),
            },
        }

        summary_path = (
            qc_directory
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

        # --------------------------------------------------------------
        # Final log
        # --------------------------------------------------------------

        logger.info(
            "M3.1 complete: %d/%d GWAS association rows selected.",
            len(
                prostate_gwas
            ),
            total_input_rows,
        )

        logger.info(
            "M3.1 direct-overlap eligible rows: %d.",
            direct_overlap_eligible_rows,
        )

        logger.info(
            "M3.1 unique canonical rsIDs: %d.",
            unique_canonical_rsids,
        )

        logger.info(
            "M3.1 unique direct-overlap eligible rsIDs: %d.",
            unique_overlap_eligible_rsids,
        )

        logger.info(
            "M3.1 QC report: %s",
            summary_path,
        )

        return report