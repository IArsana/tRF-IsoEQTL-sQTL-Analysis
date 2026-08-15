"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/summarize_statistics.py

Description:
    Dataset-level runner for M3.5 descriptive statistical summary.

    Inputs:

        data/processed/m3/prostate_gwas.parquet
        data/processed/m3/prostate_trfqtl.parquet
        data/processed/m3/candidate_snptrf.parquet

    Outputs:

        data/processed/m3/m3_candidate_statistics.parquet

    QC:

        data/processed/m3/qc/m3_5_statistics_summary.json

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

from pcatrfqtl.analysis.m3.statistics import (
    M3StatisticsAnalyzer,
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
class M35StatisticsInputs:
    """Input paths required for M3.5."""

    prostate_gwas: Path
    prostate_trfqtl: Path
    candidate_snptrf: Path


class M35StatisticsRunner:
    """Execute M3.5 descriptive statistics."""

    OUTPUT_FILENAME = (
        "m3_candidate_statistics.parquet"
    )

    SUMMARY_FILENAME = (
        "m3_5_statistics_summary.json"
    )

    def __init__(
        self,
        inputs: M35StatisticsInputs,
        output_directory: str | Path,
    ) -> None:

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:

        if not path.exists():
            raise FileNotFoundError(
                f"M3.5 input file not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"M3.5 expected file but found: {path}"
            )

    def run(
        self,
    ) -> dict[str, Any]:

        for path in (
            self.inputs.prostate_gwas,
            self.inputs.prostate_trfqtl,
            self.inputs.candidate_snptrf,
        ):
            self._require_file(
                path
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
            "Starting M3.5 descriptive statistical summary."
        )

        gwas = read_parquet(
            self.inputs.prostate_gwas
        )

        trfqtl = read_parquet(
            self.inputs.prostate_trfqtl
        )

        candidates = read_parquet(
            self.inputs.candidate_snptrf
        )

        statistics = (
            M3StatisticsAnalyzer
            .summarize(
                gwas=gwas,
                trfqtl=trfqtl,
                candidates=candidates,
            )
        )

        candidate_statistics = (
            M3StatisticsAnalyzer
            .build_candidate_statistics_table(
                candidates
            )
        )

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        write_parquet(
            candidate_statistics,
            output_path,
            index=False,
        )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M3.5",

            "stage":
                "descriptive_statistical_summary",

            "statistics":
                statistics,

            "interpretation_policy": {
                "descriptive_statistics_only":
                    True,

                "candidate_ranking_performed":
                    False,

                "causal_inference_performed":
                    False,

                "ld_analysis_performed":
                    False,

                "colocalization_performed":
                    False,

                "zero_overlap_considered_valid_result":
                    True,
            },

            "output": {
                "candidate_statistics_path":
                    str(
                        output_path
                    ),

                "candidate_statistics_rows":
                    int(
                        len(
                            candidate_statistics
                        )
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
            "M3.5 complete."
        )

        logger.info(
            "Unique GWAS rsIDs: %d",
            statistics[
                "gwas"
            ][
                "unique_eligible_rsids"
            ],
        )

        logger.info(
            "Unique tRF-QTL rsIDs: %d",
            statistics[
                "trfqtl"
            ][
                "unique_eligible_rsids"
            ],
        )

        logger.info(
            "Direct shared rsIDs: %d",
            statistics[
                "direct_overlap"
            ][
                "unique_shared_rsids"
            ],
        )

        logger.info(
            "M3.5 QC report: %s",
            summary_path,
        )

        return report