"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/integrate_candidate_genes.py

Description:
    Runner for M6.3 Candidate → Gene Integration.

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

import pandas as pd
import yaml

from pcatrfqtl.analysis.m6.candidate_gene_integration import (
    integrate_candidate_genes,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(
    __name__
)


# ============================================================================
# JSON utilities
# ============================================================================


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert pandas missing values to standards-compliant JSON null."""

    if dataframe.empty:
        return []

    cleaned = dataframe.astype(
        object
    ).where(
        pd.notna(
            dataframe
        ),
        None,
    )

    return cleaned.to_dict(
        orient="records"
    )


def _json_safe_summary(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert summary row to JSON-native values."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        try:

            if pd.isna(
                value
            ):

                output[
                    str(key)
                ] = None

                continue

        except (
            TypeError,
            ValueError,
        ):

            pass

        if hasattr(
            value,
            "item",
        ):

            value = value.item()

        output[
            str(key)
        ] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M63CandidateGeneIntegrationRunner:
    """Execute M6.3 candidate-gene integration."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        )

        self.config_path = Path(
            config_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ======================================================================
    # Outputs
    # ======================================================================

    @property
    def evidence_matrix_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_gene_evidence_matrix.parquet"
        )

    @property
    def final_integration_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "final_candidate_gene_integration.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_gene_integration_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m6_3_candidate_gene_integration.json"
        )

    # ======================================================================
    # Loading
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            config = yaml.safe_load(
                handle
            )

        if not isinstance(
            config,
            dict,
        ):

            raise ValueError(
                "M6.3 config must contain a YAML mapping."
            )

        return config

    def _load_upstream_qc(
        self,
        *,
        config: dict[str, Any],
        key: str,
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        properties = config[
            "upstream_qc"
        ][
            key
        ]

        path = (
            self.project_root
            /
            str(
                properties[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M6.3 upstream QC missing: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        if payload.get(
            "milestone"
        ) != properties[
            "expected_milestone"
        ]:

            raise RuntimeError(
                f"Unexpected milestone for upstream {key}."
            )

        if payload.get(
            "stage"
        ) != properties[
            "expected_stage"
        ]:

            raise RuntimeError(
                f"Unexpected stage for upstream {key}."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        (
            m6_1_qc,
            m6_1_path,
        ) = self._load_upstream_qc(
            config=config,
            key="candidate_prioritization",
        )

        (
            m6_2a_qc,
            m6_2a_path,
        ) = self._load_upstream_qc(
            config=config,
            key="regulatory_feature_mapping",
        )

        (
            m6_2b_qc,
            m6_2b_path,
        ) = self._load_upstream_qc(
            config=config,
            key="regulatory_gene_annotation",
        )

        (
            m6_2c_qc,
            m6_2c_path,
        ) = self._load_upstream_qc(
            config=config,
            key="event_reconstruction",
        )

        logger.info(
            "Loaded all locked M6.3 upstream QC artifacts."
        )

        result = integrate_candidate_genes(
            m6_1_qc=m6_1_qc,
            m6_2a_qc=m6_2a_qc,
            m6_2b_qc=m6_2b_qc,
            m6_2c_qc=m6_2c_qc,
            config=config,
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        write_parquet(
            result.evidence_matrix,
            self.evidence_matrix_output,
            index=False,
        )

        write_parquet(
            result.final_integration,
            self.final_integration_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        summary = _json_safe_summary(
            result.summary.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M6.3",

            "stage":
                "candidate_to_gene_integration",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m6_1":
                    str(
                        m6_1_path
                    ),

                "m6_2a":
                    str(
                        m6_2a_path
                    ),

                "m6_2b":
                    str(
                        m6_2b_path
                    ),

                "m6_2c_1":
                    str(
                        m6_2c_path
                    ),
            },

            "summary":
                summary,

            "candidate_gene_evidence":
                _records_with_json_nulls(
                    result.evidence_matrix
                ),

            "final_candidate_gene_integration":
                _records_with_json_nulls(
                    result.final_integration
                ),

            "outputs": {
                "candidate_gene_evidence":
                    str(
                        self.evidence_matrix_output
                    ),

                "final_candidate_gene_integration":
                    str(
                        self.final_integration_output
                    ),

                "summary":
                    str(
                        self.summary_output
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
                allow_nan=False,
            )

        logger.info(
            "M6.3 complete."
        )

        logger.info(
            "Candidates=%d | resolved gene=%d | "
            "feature-supported unresolved=%d | "
            "no defensible gene assignment=%d.",
            summary[
                "candidates_assessed"
            ],
            summary[
                "candidates_with_resolved_gene"
            ],
            summary[
                "candidates_with_supported_feature_unresolved_gene"
            ],
            summary[
                "candidates_without_defensible_gene_assignment"
            ],
        )

        for row in report[
            "final_candidate_gene_integration"
        ]:

            logger.info(
                "%s | feature=%s | gene=%s | status=%s.",
                row[
                    "lead_rsid"
                ],
                row[
                    "regulatory_feature_id"
                ],
                row[
                    "integrated_gene"
                ],
                row[
                    "gene_resolution_status"
                ],
            )

        logger.info(
            "Next stage: %s.",
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report