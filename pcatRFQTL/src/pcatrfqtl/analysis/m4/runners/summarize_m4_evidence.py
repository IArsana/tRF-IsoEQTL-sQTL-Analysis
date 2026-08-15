"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/runners/summarize_m4_evidence.py

Description:
    Runner for M4.3 Multi-omic Exact-Evidence Summary.

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

from pcatrfqtl.analysis.m4.evidence_summary import (
    M43ExactEvidenceSummaryBuilder,
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
class M43EvidenceSummaryInputs:
    """Inputs for M4.3."""

    prostate_trfqtl: Path

    moradi_matches: Path

    m3_candidate_scaffold: Path | None = None


class M43EvidenceSummaryRunner:
    """Run M4.3 candidate-level exact evidence summarization."""

    OUTPUT_FILENAME = (
        "m4_candidate_evidence_summary.parquet"
    )

    SUMMARY_FILENAME = (
        "m4_3_exact_evidence_summary.json"
    )

    def __init__(
        self,
        *,
        inputs: M43EvidenceSummaryInputs,
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
    def _require_file(
        path: Path,
    ) -> None:
        if not path.exists():
            raise FileNotFoundError(
                f"M4.3 required file not found: {path}"
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

        self._require_file(
            self.inputs.prostate_trfqtl
        )

        self._require_file(
            self.inputs.moradi_matches
        )

        trfqtl = read_parquet(
            self.inputs.prostate_trfqtl
        )

        moradi_matches = read_parquet(
            self.inputs.moradi_matches
        )

        m3_scaffold = None

        if (
            self.inputs.m3_candidate_scaffold
            is not None
            and self.inputs.m3_candidate_scaffold.exists()
        ):
            m3_scaffold = read_parquet(
                self.inputs.m3_candidate_scaffold
            )

        logger.info(
            "Starting M4.3 Multi-omic Exact-Evidence Summary."
        )

        result = (
            M43ExactEvidenceSummaryBuilder
            .build(
                trfqtl,
                moradi_matches,
                m3_candidate_scaffold=(
                    m3_scaffold
                ),
            )
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
                "M4.3",

            "stage":
                "multiomic_exact_evidence_summary",

            "output":
                str(
                    output_path
                ),

            "policy": {
                "candidate_level_summary":
                    True,
                "coordinate_matching_performed":
                    False,
                "liftover_performed":
                    False,
                "ld_inference_performed":
                    False,
                "colocalization_performed":
                    False,
                "candidate_ranking_performed":
                    False,
            },

            "summary": {
                "candidate_rows":
                    len(
                        result
                    ),

                "unique_trfqtl_rsids":
                    int(
                        result[
                            "trfqtl_rsid"
                        ]
                        .nunique()
                    ),

                "direct_gwas_match_candidates":
                    self._count_true(
                        result,
                        "direct_gwas_match",
                    ),

                "moradi_direct_qtl_match_candidates":
                    self._count_true(
                        result,
                        "moradi_direct_qtl_match",
                    ),

                "moradi_source_tag_match_candidates":
                    self._count_true(
                        result,
                        "moradi_source_tag_match",
                    ),

                "cis_exact_support_candidates":
                    self._count_true(
                        result,
                        "cis_exact_support",
                    ),

                "trans_exact_support_candidates":
                    self._count_true(
                        result,
                        "trans_exact_support",
                    ),

                "exact_multiomic_support_candidates":
                    self._count_true(
                        result,
                        "exact_multiomic_support",
                    ),

                "m5_ld_required_candidates":
                    self._count_true(
                        result,
                        "m5_ld_required",
                    ),

                "status_counts":
                    self._value_counts(
                        result,
                        "m4_exact_evidence_status",
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
            "M4.3 completed: %d candidate rows; "
            "%d exact multi-omic bridges.",
            len(result),
            self._count_true(
                result,
                "exact_multiomic_support",
            ),
        )

        return report