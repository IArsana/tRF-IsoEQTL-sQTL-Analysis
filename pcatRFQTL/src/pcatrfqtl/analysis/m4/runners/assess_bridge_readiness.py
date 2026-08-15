"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/runners/assess_bridge_readiness.py

Description:
    Runner for M4.4 Regulatory Bridge Readiness assessment.

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

from pcatrfqtl.analysis.m4.bridge_readiness import (
    M44RegulatoryBridgeReadinessBuilder,
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
class M44BridgeReadinessInputs:
    """M4.4 input resources."""

    evidence_summary: Path


class M44BridgeReadinessRunner:
    """Run M4.4 Regulatory Bridge Readiness."""

    OUTPUT_FILENAME = (
        "m4_regulatory_bridge_readiness.parquet"
    )

    SUMMARY_FILENAME = (
        "m4_4_regulatory_bridge_readiness.json"
    )

    def __init__(
        self,
        *,
        inputs: M44BridgeReadinessInputs,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:
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
        return int(
            dataframe[
                column
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
        values = (
            dataframe[
                column
            ]
            .dropna()
            .astype(str)
            .value_counts()
        )

        return {
            str(key): int(value)
            for key, value
            in values.items()
        }

    def run(
        self,
    ) -> dict[str, Any]:

        if not self.inputs.evidence_summary.exists():
            raise FileNotFoundError(
                "M4.4 evidence summary not found: "
                f"{self.inputs.evidence_summary}"
            )

        evidence = read_parquet(
            self.inputs.evidence_summary
        )

        logger.info(
            "Starting M4.4 Regulatory Bridge Readiness."
        )

        result = (
            M44RegulatoryBridgeReadinessBuilder
            .build(
                evidence
            )
        )

        if (
            len(result)
            != len(evidence)
        ):
            raise RuntimeError(
                "M4.4 changed candidate cardinality: "
                f"{len(evidence)} -> {len(result)}"
            )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        write_parquet(
            result,
            output_path,
            index=False,
        )

        report: dict[str, Any] = {
            "milestone":
                "M4.4",

            "stage":
                "regulatory_bridge_readiness",

            "output":
                str(
                    output_path
                ),

            "policy": {
                "workflow_routing_only":
                    True,
                "candidate_ranking_performed":
                    False,
                "ld_lookup_performed":
                    False,
                "locus_mapping_performed":
                    False,
                "colocalization_performed":
                    False,
                "causal_inference_performed":
                    False,
                "candidate_cardinality_preserved":
                    True,
            },

            "summary": {
                "candidate_rows":
                    len(result),

                "bridge_state_counts":
                    self._value_counts(
                        result,
                        "regulatory_bridge_state",
                    ),

                "requires_ld_lookup":
                    self._count_true(
                        result,
                        "m5_requires_ld_lookup",
                    ),

                "requires_locus_mapping":
                    self._count_true(
                        result,
                        "m5_requires_locus_mapping",
                    ),

                "requires_colocalization_assessment":
                    self._count_true(
                        result,
                        "m5_requires_colocalization_assessment",
                    ),

                "search_gwas_neighbors":
                    self._count_true(
                        result,
                        "m5_search_gwas_neighbors",
                    ),

                "search_qtl_neighbors":
                    self._count_true(
                        result,
                        "m5_search_qtl_neighbors",
                    ),

                "search_cis_regulatory_bridge":
                    self._count_true(
                        result,
                        "m5_search_cis_regulatory_bridge",
                    ),

                "search_trans_regulatory_bridge":
                    self._count_true(
                        result,
                        "m5_search_trans_regulatory_bridge",
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
            "M4.4 completed: %d candidates assessed.",
            len(result),
        )

        return report