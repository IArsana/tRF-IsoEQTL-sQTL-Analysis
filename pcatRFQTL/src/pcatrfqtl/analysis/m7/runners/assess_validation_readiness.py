"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/assess_validation_readiness.py

Description:
    Runner for M7.1 TCGA-PRAD Clinical Validation Readiness Audit.

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

from pcatrfqtl.analysis.m7.validation_readiness import (
    assess_validation_readiness,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(
    __name__
)


# ============================================================================
# JSON helpers
# ============================================================================


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert pandas missing values to JSON null."""

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
    """Convert pandas scalar values to JSON-native values."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:

            output[
                str(key)
            ] = None

            continue

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


class M71ValidationReadinessRunner:
    """Execute M7.1 TCGA-PRAD validation readiness audit."""

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
    def resource_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "validation_resource_readiness.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_validation_readiness.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "validation_readiness_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m7_1_validation_readiness.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M7.1 configuration not found: "
                f"{self.config_path}"
            )

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
                "M7.1 configuration must contain a YAML mapping."
            )

        required_sections = {
            "target_cohort",
            "upstream_qc",
            "gdc_capabilities",
            "trf_validation_policy",
            "variant_validation_policy",
            "regulatory_event_policy",
            "statuses",
            "scientific_policy",
            "next_stage",
        }

        missing = (
            required_sections
            -
            set(
                config
            )
        )

        if missing:

            raise ValueError(
                "M7.1 configuration missing required sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Upstream validation
    # ======================================================================

    def _load_m6_4b_qc(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        properties = config[
            "upstream_qc"
        ][
            "pathway_readiness"
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
                "M6.4B QC not found: "
                f"{path}"
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
                "Unexpected M6.4B milestone."
            )

        if payload.get(
            "stage"
        ) != properties[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected M6.4B stage."
            )

        expected_status = properties.get(
            "expected_pathway_status"
        )

        if expected_status is not None:

            observed_status = (
                payload
                .get(
                    "summary",
                    {},
                )
                .get(
                    "pathway_readiness_status"
                )
            )

            if observed_status != expected_status:

                raise RuntimeError(
                    "Unexpected M6.4B pathway status: "
                    f"{observed_status}; "
                    f"expected {expected_status}."
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
        """Execute complete M7.1 workflow."""

        config = self._load_config()

        logger.info(
            "Loaded M7.1 validation-readiness configuration."
        )

        (
            m6_4b_qc,
            m6_4b_path,
        ) = self._load_m6_4b_qc(
            config=config,
        )

        logger.info(
            "Loaded locked M6.4B pathway-readiness QC."
        )

        result = assess_validation_readiness(
            m6_4b_qc=m6_4b_qc,
            config=config,
        )

        if result.candidate_readiness.empty:

            raise RuntimeError(
                "M7.1 candidate readiness is unexpectedly empty."
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

        summary = _json_safe_summary(
            result.summary.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M7.1",

            "stage":
                "tcga_prad_clinical_validation_readiness",

            "target_cohort":
                config[
                    "target_cohort"
                ],

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m6_4b_loaded":
                    True,

                "m6_4b_path":
                    str(
                        m6_4b_path
                    ),
            },

            "summary":
                summary,

            "validation_resource_readiness":
                _records_with_json_nulls(
                    result.resource_readiness
                ),

            "candidate_validation_readiness":
                _records_with_json_nulls(
                    result.candidate_readiness
                ),

            "outputs": {
                "validation_resource_readiness":
                    str(
                        self.resource_output
                    ),

                "candidate_validation_readiness":
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
            "M7.1 complete."
        )

        logger.info(
            "Candidates=%d | tRF identities=%d | "
            "tRF quantification ready=%d | "
            "variant validation ready=%d.",
            summary[
                "candidates_assessed"
            ],
            summary[
                "candidates_with_trf_identity"
            ],
            summary[
                "candidates_ready_for_direct_trf_quantification"
            ],
            summary[
                "candidates_ready_for_direct_variant_validation"
            ],
        )

        logger.info(
            "Overall readiness: %s.",
            summary[
                "overall_readiness_status"
            ],
        )

        for row in report[
            "candidate_validation_readiness"
        ]:

            logger.info(
                "%s | tRF=%s | trf_ready=%s | "
                "variant_ready=%s | regulatory_event=%s.",
                row[
                    "lead_rsid"
                ],
                row[
                    "trf_ids"
                ],
                row[
                    "direct_trf_quantification_ready"
                ],
                row[
                    "direct_variant_validation_ready"
                ],
                row[
                    "regulatory_event_validation_status"
                ],
            )

        logger.info(
            "Clinical next stage: %s.",
            summary[
                "clinical_next_stage"
            ],
        )

        logger.info(
            "tRF next stage: %s.",
            summary[
                "trf_next_stage"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report