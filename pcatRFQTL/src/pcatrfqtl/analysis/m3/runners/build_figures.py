"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/build_figures.py

Description:
    Dataset-level runner for M3.6 descriptive figure generation.

    M3.6 consumes the M3.5 summary and produces publication-oriented
    descriptive figures representing dataset reduction and direct
    canonical-rsID intersection.

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

from pcatrfqtl.analysis.m3.figures import (
    M3FigureBuilder,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


@dataclass(frozen=True)
class M36FigureInputs:
    """Input paths required for M3.6."""

    statistics_summary: Path


class M36FigureRunner:
    """Execute M3.6 figure generation."""

    SUMMARY_FILENAME = (
        "m3_6_figures_summary.json"
    )

    def __init__(
        self,
        inputs: M36FigureInputs,
        figure_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.inputs = inputs

        self.figure_directory = Path(
            figure_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    def run(
        self,
    ) -> dict[str, Any]:

        if not (
            self.inputs
            .statistics_summary
            .exists()
        ):
            raise FileNotFoundError(
                "M3.6 statistics summary not found: "
                f"{self.inputs.statistics_summary}"
            )

        with (
            self.inputs
            .statistics_summary
            .open(
                "r",
                encoding="utf-8",
            )
        ) as handle:

            report = json.load(
                handle
            )

        statistics = report[
            "statistics"
        ]

        logger.info(
            "Starting M3.6 figure generation."
        )

        flow_paths = (
            M3FigureBuilder
            .dataset_flow(
                statistics,
                self.figure_directory,
            )
        )

        intersection_paths = (
            M3FigureBuilder
            .variant_intersection(
                statistics,
                self.figure_directory,
            )
        )

        all_paths = (
            flow_paths
            + intersection_paths
        )

        if not all(
            path.exists()
            for path in all_paths
        ):
            raise RuntimeError(
                "M3.6 failed to generate one or more figures."
            )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        summary: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M3.6",

            "stage":
                "figures_and_qc",

            "figures": {
                "dataset_flow":
                    [
                        str(
                            path
                        )
                        for path
                        in flow_paths
                    ],

                "variant_intersection":
                    [
                        str(
                            path
                        )
                        for path
                        in intersection_paths
                    ],
            },

            "qc": {
                "expected_figure_files":
                    4,

                "generated_figure_files":
                    len(
                        all_paths
                    ),

                "all_figures_exist":
                    all(
                        path.exists()
                        for path
                        in all_paths
                    ),

                "direct_overlap_count":
                    statistics[
                        "direct_overlap"
                    ][
                        "unique_shared_rsids"
                    ],

                "ld_evidence_displayed":
                    False,

                "colocalization_displayed":
                    False,

                "causal_inference_displayed":
                    False,
            },
        }

        summary_path = (
            self.qc_directory
            / self.SUMMARY_FILENAME
        )

        summary[
            "report_path"
        ] = str(
            summary_path
        )

        with summary_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                summary,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M3.6 complete: %d figure files generated.",
            len(
                all_paths
            ),
        )

        return summary