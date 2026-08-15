"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/resolve_coloc_data_gap.py

Description:
    Runner for M5.6 colocalization data-gap resolution.

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
from pathlib import Path
from typing import Any

import yaml

from pcatrfqtl.analysis.m5.coloc_data_gap import (
    resolve_coloc_data_gap,
)
from pcatrfqtl.io.parquet import (
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


class M56ColocDataGapRunner:
    """Execute M5.6 colocalization data-gap resolution."""

    def __init__(
        self,
        *,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.config_path = Path(
            config_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    @property
    def source_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "coloc_data_gap_sources.parquet"
        )

    @property
    def resolution_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "coloc_data_gap_resolution.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_coloc_resolution.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_6_coloc_data_gap.json"
        )

    def run(
        self,
    ) -> dict[str, Any]:
        """Run M5.6."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M5.6 config not found: {self.config_path}"
            )

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            config = yaml.safe_load(
                handle
            )

        result = resolve_coloc_data_gap(
            config
        )

        sources = result.sources

        resolution = result.resolution

        candidates = result.candidates

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        write_parquet(
            sources,
            self.source_output,
            index=False,
        )

        write_parquet(
            resolution,
            self.resolution_output,
            index=False,
        )

        write_parquet(
            candidates,
            self.candidate_output,
            index=False,
        )

        formal_ready = int(
            candidates[
                "formal_coloc_ready"
            ]
            .fillna(
                False
            )
            .sum()
        )

        resolution_state = str(
            resolution.iloc[
                0
            ][
                "resolution_state"
            ]
        )

        report = {
            "milestone":
                "M5.6",

            "stage":
                "colocalization_data_gap_resolution",

            "policy": {
                "significant_only_qtl_allowed_for_coloc":
                    False,

                "missing_qtl_rows_interpreted_as_null":
                    False,

                "missing_qtl_effects_imputed":
                    False,

                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "source_options_assessed":
                    int(
                        len(
                            sources
                        )
                    ),

                "candidate_leads_assessed":
                    int(
                        len(
                            candidates
                        )
                    ),

                "moradi_resolution_state":
                    resolution_state,

                "formal_coloc_ready_candidates":
                    formal_ready,

                "formal_coloc_blocked_candidates":
                    int(
                        len(
                            candidates
                        )
                        -
                        formal_ready
                    ),
            },

            "resolution":
                resolution.to_dict(
                    orient="records"
                ),

            "source_status_counts":
                {
                    str(
                        key
                    ):
                        int(
                            value
                        )
                    for key, value
                    in sources[
                        "role"
                    ]
                    .value_counts()
                    .items()
                },

            "candidate_status_counts":
                {
                    str(
                        key
                    ):
                        int(
                            value
                        )
                    for key, value
                    in candidates[
                        "current_resolution"
                    ]
                    .value_counts()
                    .items()
                },

            "outputs": {
                "sources":
                    str(
                        self.source_output
                    ),

                "resolution":
                    str(
                        self.resolution_output
                    ),

                "candidate_resolution":
                    str(
                        self.candidate_output
                    ),
            },
        }

        with self.qc_output.open(
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
            "M5.6 complete."
        )

        logger.info(
            "Moradi resolution: %s.",
            resolution_state,
        )

        logger.info(
            "Formal coloc-ready candidates: %d/%d.",
            formal_ready,
            len(
                candidates
            ),
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report