"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/direct_overlap.py

Description:
    Dataset-level runner for M3.3 direct GWAS × tRF-QTL overlap.

    Inputs:

        data/processed/m3/prostate_gwas.parquet
        data/processed/m3/prostate_trfqtl.parquet

    Outputs:

        data/processed/m3/direct_overlap.parquet

    QC:

        data/processed/m3/qc/m3_3_overlap_summary.json

    Direct matching is performed exclusively by exact canonical rsID.

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

from pcatrfqtl.analysis.m3.overlap import (
    DirectOverlapAnalyzer,
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
class M33OverlapInputs:
    """Input paths required for M3.3."""

    prostate_gwas: Path
    prostate_trfqtl: Path


class M33OverlapRunner:
    """Execute exact canonical-rsID GWAS × tRF-QTL overlap."""

    OUTPUT_FILENAME = (
        "direct_overlap.parquet"
    )

    SUMMARY_FILENAME = (
        "m3_3_overlap_summary.json"
    )

    def __init__(
        self,
        inputs: M33OverlapInputs,
        output_directory: str | Path,
    ) -> None:
        """Initialize M3.3 runner."""

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
                f"M3.3 input file not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"M3.3 expected file but found: {path}"
            )

    @staticmethod
    def _safe_nunique(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count unique non-missing values."""

        if (
            dataframe.empty
            or column
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
    def _count_eligible(
        dataframe: pd.DataFrame,
    ) -> int:
        """Count direct-overlap eligible rows."""

        if (
            "m3_direct_overlap_eligible"
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                "m3_direct_overlap_eligible"
            ]
            .fillna(False)
            .eq(True)
            .sum()
        )

    @staticmethod
    def _value_counts(
        dataframe: pd.DataFrame,
        column: str,
    ) -> dict[str, int]:
        """Return JSON-safe value counts."""

        if (
            dataframe.empty
            or column
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
            str(key): int(value)
            for (
                key,
                value,
            ) in counts.items()
        }

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M3.3."""

        self._require_file(
            self.inputs.prostate_gwas
        )

        self._require_file(
            self.inputs.prostate_trfqtl
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
            "Starting M3.3 direct GWAS × tRF-QTL overlap."
        )

        gwas = (
            read_parquet(
                self.inputs.prostate_gwas
            )
        )

        trfqtl = (
            read_parquet(
                self.inputs.prostate_trfqtl
            )
        )

        gwas_eligible = (
            DirectOverlapAnalyzer
            .eligible_gwas(
                gwas
            )
        )

        trfqtl_eligible = (
            DirectOverlapAnalyzer
            .eligible_trfqtl(
                trfqtl
            )
        )

        shared_rsids = (
            DirectOverlapAnalyzer
            .shared_rsids(
                gwas,
                trfqtl,
            )
        )

        overlap = (
            DirectOverlapAnalyzer
            .build_overlap_table(
                gwas,
                trfqtl,
            )
        )

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        # Ensure an empty result still has a stable Parquet output.
        if overlap.empty:
            logger.warning(
                "M3.3 found zero direct canonical-rsID overlaps."
            )

            overlap = pd.DataFrame(
                {
                    "harm_variant_rsid":
                        pd.Series(
                            dtype="string"
                        ),
                    "m3_overlap_method":
                        pd.Series(
                            dtype="string"
                        ),
                    "m3_direct_variant_overlap":
                        pd.Series(
                            dtype="boolean"
                        ),
                }
            )

        write_parquet(
            overlap,
            output_path,
            index=False,
        )

        unique_shared_rsids = (
            len(
                shared_rsids
            )
        )

        unique_shared_trfs = (
            self._safe_nunique(
                overlap,
                "harm_feature_id_trfqtl",
            )
        )

        if (
            unique_shared_trfs
            == 0
            and "harm_feature_id"
            in overlap.columns
        ):
            unique_shared_trfs = (
                self._safe_nunique(
                    overlap,
                    "harm_feature_id",
                )
            )

        unique_variant_trf_pairs = 0

        feature_column = None

        if (
            "harm_feature_identity_key_trfqtl"
            in overlap.columns
        ):
            feature_column = (
                "harm_feature_identity_key_trfqtl"
            )

        elif (
            "harm_feature_identity_key"
            in overlap.columns
        ):
            feature_column = (
                "harm_feature_identity_key"
            )

        if (
            not overlap.empty
            and feature_column
            is not None
            and "harm_variant_rsid"
            in overlap.columns
        ):
            unique_variant_trf_pairs = int(
                overlap[
                    [
                        "harm_variant_rsid",
                        feature_column,
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
                "M3.3",

            "stage":
                "direct_gwas_trfqtl_overlap",

            "inputs": {
                "gwas": {
                    "path":
                        str(
                            self.inputs.prostate_gwas
                        ),

                    "rows":
                        int(
                            len(
                                gwas
                            )
                        ),

                    "eligible_rows":
                        self._count_eligible(
                            gwas
                        ),

                    "unique_eligible_rsids":
                        self._safe_nunique(
                            gwas_eligible,
                            "harm_variant_rsid",
                        ),
                },

                "trfqtl": {
                    "path":
                        str(
                            self.inputs.prostate_trfqtl
                        ),

                    "rows":
                        int(
                            len(
                                trfqtl
                            )
                        ),

                    "eligible_rows":
                        self._count_eligible(
                            trfqtl
                        ),

                    "unique_eligible_rsids":
                        self._safe_nunique(
                            trfqtl_eligible,
                            "harm_variant_rsid",
                        ),
                },
            },

            "output": {
                "path":
                    str(
                        output_path
                    ),

                "rows":
                    int(
                        len(
                            overlap
                        )
                    ),

                "unique_shared_rsids":
                    int(
                        unique_shared_rsids
                    ),

                "unique_shared_trfs":
                    int(
                        unique_shared_trfs
                    ),

                "unique_variant_trf_pairs":
                    int(
                        unique_variant_trf_pairs
                    ),

                "shared_rsids":
                    sorted(
                        shared_rsids
                    ),
            },

            "method": {
                "match_key":
                    "canonical_rsid",

                "exact_matching":
                    True,

                "coordinate_matching_performed":
                    False,

                "ld_proxy_matching_performed":
                    False,

                "liftover_performed":
                    False,

                "deduplication_performed":
                    False,

                "study_collapsing_performed":
                    False,

                "colocalization_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "qc": {
                "overlap_found":
                    bool(
                        unique_shared_rsids
                        > 0
                    ),

                "overlap_method_counts":
                    self._value_counts(
                        overlap,
                        "m3_overlap_method",
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
            "M3.3 complete: %d shared canonical rsIDs.",
            unique_shared_rsids,
        )

        logger.info(
            "M3.3 overlap rows: %d.",
            len(
                overlap
            ),
        )

        logger.info(
            "M3.3 unique shared tRFs: %d.",
            unique_shared_trfs,
        )

        logger.info(
            "M3.3 QC report: %s",
            summary_path,
        )

        return report