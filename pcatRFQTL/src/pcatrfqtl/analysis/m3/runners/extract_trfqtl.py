"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/extract_trfqtl.py

Description:
    Dataset-level runner for M3.2 prostate-cancer tRF-QTL extraction.

    M3.2 reads the harmonized Cancer-tRFQTL S2 dataset and selects
    prostate-cancer association records.

    All selected SNP-to-tRF association rows are retained without
    deduplication or candidate ranking.

    Output:
        data/processed/m3/prostate_trfqtl.parquet

    QC:
        data/processed/m3/qc/m3_2_trfqtl_summary.json

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

from pcatrfqtl.analysis.m3.trfqtl import (
    ProstateTRFQTLSelector,
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
class M32TRFQTLInputs:
    """Input paths required for M3.2."""

    harmonized_trfqtl_s2: Path


class M32TRFQTLRunner:
    """Execute M3.2 prostate-cancer tRF-QTL extraction."""

    OUTPUT_FILENAME = (
        "prostate_trfqtl.parquet"
    )

    SUMMARY_FILENAME = (
        "m3_2_trfqtl_summary.json"
    )

    def __init__(
        self,
        inputs: M32TRFQTLInputs,
        output_directory: str | Path,
    ) -> None:
        """Initialize the M3.2 runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require a valid input file."""

        if not path.exists():
            raise FileNotFoundError(
                "M3.2 input file not found: "
                f"{path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                "M3.2 expected file but found: "
                f"{path}"
            )

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
            .fillna(False)
            .eq(True)
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

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M3.2 prostate-cancer tRF-QTL extraction."""

        source = (
            self.inputs
            .harmonized_trfqtl_s2
        )

        self._require_file(
            source
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
            "Starting M3.2 prostate-cancer tRF-QTL extraction."
        )

        dataframe = (
            read_parquet(
                source
            )
        )

        input_rows = len(
            dataframe
        )

        selected = (
            ProstateTRFQTLSelector
            .prepare(
                dataframe
            )
        )

        selected[
            "m3_source_table"
        ] = "Cancer-tRFQTL_S2"

        # --------------------------------------------------------------
        # Scientific invariants
        # --------------------------------------------------------------

        if not (
            selected[
                "harm_is_primary_disease"
            ]
            .fillna(False)
            .eq(True)
            .all()
        ):
            raise RuntimeError(
                "M3.2 output contains records that are not marked "
                "as the primary disease."
            )

        if not (
            selected[
                "harm_disease_id"
            ]
            .eq(
                ProstateTRFQTLSelector
                .PRIMARY_DISEASE_ID
            )
            .all()
        ):
            raise RuntimeError(
                "M3.2 output contains a disease identity other than "
                f"{ProstateTRFQTLSelector.PRIMARY_DISEASE_ID}."
            )

        if not (
            selected[
                "harm_feature_type"
            ]
            .eq(
                "TRF"
            )
            .all()
        ):
            raise RuntimeError(
                "M3.2 output contains a biological feature that is "
                "not classified as TRF."
            )

        # --------------------------------------------------------------
        # Output
        # --------------------------------------------------------------

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        write_parquet(
            selected,
            output_path,
            index=False,
        )

        # --------------------------------------------------------------
        # QC
        # --------------------------------------------------------------

        eligible_rows = (
            self._count_true(
                selected,
                "m3_direct_overlap_eligible",
            )
        )

        unique_rsids = (
            self._safe_nunique(
                selected,
                "harm_variant_rsid",
            )
        )

        unique_trfs = (
            self._safe_nunique(
                selected,
                "harm_feature_id",
            )
        )

        unique_variant_trf_pairs = 0

        if (
            "harm_variant_rsid"
            in selected.columns
            and "harm_feature_identity_key"
            in selected.columns
        ):
            unique_variant_trf_pairs = int(
                selected[
                    [
                        "harm_variant_rsid",
                        "harm_feature_identity_key",
                    ]
                ]
                .dropna()
                .drop_duplicates()
                .shape[
                    0
                ]
            )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M3.2",

            "stage":
                "prostate_cancer_trfqtl_subset",

            "input": {
                "path":
                    str(
                        source
                    ),

                "rows":
                    int(
                        input_rows
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
                            selected
                        )
                    ),

                "unique_canonical_rsids":
                    unique_rsids,

                "unique_trfs":
                    unique_trfs,

                "unique_variant_trf_pairs":
                    unique_variant_trf_pairs,

                "direct_overlap_eligible_rows":
                    eligible_rows,

                "direct_overlap_ineligible_rows":
                    int(
                        len(
                            selected
                        )
                        - eligible_rows
                    ),
            },

            "selection": {
                "primary_disease_id":
                    (
                        ProstateTRFQTLSelector
                        .PRIMARY_DISEASE_ID
                    ),

                "source_table":
                    "Cancer-tRFQTL S2",

                "disease_filtering_performed":
                    True,

                "rsid_filtering_performed":
                    False,

                "feature_filtering_performed":
                    False,

                "deduplication_performed":
                    False,

                "coordinate_matching_performed":
                    False,

                "liftover_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "qc": {
                "all_rows_primary_disease":
                    bool(
                        selected[
                            "harm_is_primary_disease"
                        ]
                        .fillna(False)
                        .eq(True)
                        .all()
                    ),

                "all_rows_canonical_prostate_cancer":
                    bool(
                        selected[
                            "harm_disease_id"
                        ]
                        .eq(
                            ProstateTRFQTLSelector
                            .PRIMARY_DISEASE_ID
                        )
                        .all()
                    ),

                "all_features_are_trf":
                    bool(
                        selected[
                            "harm_feature_type"
                        ]
                        .eq(
                            "TRF"
                        )
                        .all()
                    ),

                "variant_status_counts":
                    self._value_counts(
                        selected,
                        "harm_variant_status",
                    ),

                "feature_status_counts":
                    self._value_counts(
                        selected,
                        "harm_feature_status",
                    ),

                "disease_status_counts":
                    self._value_counts(
                        selected,
                        "harm_disease_status",
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

        logger.info(
            "M3.2 complete: %d/%d tRF-QTL rows selected.",
            len(
                selected
            ),
            input_rows,
        )

        logger.info(
            "M3.2 unique canonical rsIDs: %d.",
            unique_rsids,
        )

        logger.info(
            "M3.2 unique tRFs: %d.",
            unique_trfs,
        )

        logger.info(
            "M3.2 direct-overlap eligible rows: %d.",
            eligible_rows,
        )

        logger.info(
            "M3.2 QC report: %s",
            summary_path,
        )

        return report