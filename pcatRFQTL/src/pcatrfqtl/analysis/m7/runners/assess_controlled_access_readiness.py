"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/assess_controlled_access_readiness.py

Description:
    Runner for M7.3C Controlled-Access Quantification Readiness.

    The GDC authentication token itself is never written to logs,
    parquet artifacts, or QC JSON.

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
import os
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pcatrfqtl.analysis.m7.controlled_access_readiness import (
    assess_controlled_access_readiness,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(
    __name__
)


# ============================================================================
# Helpers
# ============================================================================


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert dataframe records into JSON-safe values."""

    if dataframe.empty:
        return []

    cleaned = dataframe.astype(
        object
    ).where(
        pd.notna(dataframe),
        None,
    )

    records = cleaned.to_dict(
        orient="records"
    )

    for record in records:
        for key, value in list(
            record.items()
        ):
            if hasattr(
                value,
                "item",
            ):
                record[key] = value.item()

    return records


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert pandas row into JSON-native values."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:
            output[str(key)] = None
            continue

        try:
            if pd.isna(value):
                output[str(key)] = None
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

        output[str(key)] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M73CControlledAccessReadinessRunner:
    """Execute M7.3C controlled-access readiness audit."""

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

    @property
    def resource_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            /
            "controlled_resource_readiness.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            /
            "candidate_controlled_access_readiness.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            /
            "controlled_access_readiness_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            /
            "m7_3c_controlled_access_readiness.json"
        )

    def _load_config(
        self,
    ) -> dict[str, Any]:

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:
            config = yaml.safe_load(handle)

        if not isinstance(config, dict):
            raise RuntimeError(
                "M7.3C config must be a YAML mapping."
            )

        return config

    def _load_upstream_qc(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        upstream = config[
            "upstream_qc"
        ][
            "quantification_feasibility"
        ]

        path = (
            self.project_root
            /
            str(
                upstream[
                    "relative_path"
                ]
            )
        )

        if not path.exists():
            raise FileNotFoundError(
                f"M7.3B QC artifact missing: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:
            payload = json.load(handle)

        if (
            payload.get("milestone")
            !=
            upstream["expected_milestone"]
        ):
            raise RuntimeError(
                "Unexpected upstream milestone for M7.3C."
            )

        if (
            payload.get("stage")
            !=
            upstream["expected_stage"]
        ):
            raise RuntimeError(
                "Unexpected upstream stage for M7.3C."
            )

        observed_status = (
            payload
            .get(
                "summary",
                {},
            )
            .get(
                "overall_feasibility_status"
            )
        )

        expected_status = upstream.get(
            "expected_overall_status"
        )

        if (
            expected_status is not None
            and
            observed_status != expected_status
        ):
            raise RuntimeError(
                "Unexpected M7.3B feasibility status: "
                f"{observed_status!r}."
            )

        return payload, path

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        (
            m7_3b_qc,
            m7_3b_path,
        ) = self._load_upstream_qc(
            config=config,
        )

        result = assess_controlled_access_readiness(
            m7_3b_qc=m7_3b_qc,
            config=config,
            environment=dict(
                os.environ
            ),
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
            result.resource_readiness,
            self.resource_output,
            index=False,
        )

        write_parquet(
            result.candidate_readiness,
            self.candidate_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        summary = _json_safe_row(
            result.summary.iloc[0]
        )

        # SECURITY:
        # token path itself is deliberately not written to QC.
        report = {
            "milestone":
                "M7.3C",

            "stage":
                "controlled_access_quantification_readiness",

            "policy":
                config[
                    "scientific_policy"
                ],

            "security": {
                "token_environment_variable":
                    config[
                        "authentication"
                    ][
                        "token_path_environment_variable"
                    ],

                "token_contents_logged":
                    False,

                "token_contents_written_to_qc":
                    False,

                "token_hash_written_to_qc":
                    False,
            },

            "upstream_validation": {
                "m7_3b_loaded":
                    True,

                "m7_3b_path":
                    str(
                        m7_3b_path
                    ),

                "m7_3b_status":
                    (
                        m7_3b_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "overall_feasibility_status"
                        )
                    ),
            },

            "summary":
                summary,

            "candidate_controlled_access_readiness":
                _records_with_json_nulls(
                    result.candidate_readiness
                ),

            "controlled_resource_summary": {
                "resource_count":
                    int(
                        len(
                            result.resource_readiness
                        )
                    ),

                "ready_resource_count":
                    int(
                        result.resource_readiness[
                            "resource_ready_for_controlled_download"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                        if not result.resource_readiness.empty
                        else 0
                    ),
            },

            "outputs": {
                "resource_readiness":
                    str(
                        self.resource_output
                    ),

                "candidate_readiness":
                    str(
                        self.candidate_output
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
            "M7.3C complete."
        )

        logger.info(
            "Controlled BAMs=%d | resource-ready=%d | "
            "sequence-ready candidates=%d | "
            "access-configured candidates=%d | download-ready=%d.",
            summary[
                "controlled_bam_files_assessed"
            ],
            summary[
                "controlled_bam_files_resource_ready"
            ],
            summary[
                "candidate_sequences_ready"
            ],
            summary[
                "candidate_access_configurations_ready"
            ],
            summary[
                "candidates_ready_for_controlled_download"
            ],
        )

        logger.info(
            "Token configured=%s | exists=%s | readable=%s.",
            summary[
                "token_path_configured"
            ],
            summary[
                "token_file_exists"
            ],
            summary[
                "token_file_readable"
            ],
        )

        logger.info(
            "Overall readiness=%s | next=%s.",
            summary[
                "overall_readiness_status"
            ],
            summary[
                "next_stage"
            ],
        )

        return report