"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/assess_event_reconstruction.py

Description:
    Runner for M6.2C.1 FASE/Event Reconstruction Feasibility Audit.

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

from pcatrfqtl.analysis.m6.event_reconstruction_readiness import (
    assess_event_reconstruction_readiness,
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
    """Convert pandas missing values to JSON-safe null."""

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
    """Convert summary Series to native JSON-safe values."""

    output: dict[
        str,
        Any,
    ] = {}

    for key, value in row.to_dict().items():

        if value is None:

            output[
                str(key)
            ] = None

            continue

        try:

            if pd.isna(value):

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


class M62C1EventReconstructionRunner:
    """Execute M6.2C.1 reconstruction readiness audit."""

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
    def components_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "reconstruction_components.parquet"
        )

    @property
    def readiness_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "feature_reconstruction_readiness.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "reconstruction_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m6_2c_1_event_reconstruction.json"
        )

    # ======================================================================
    # Config
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M6.2C.1 config not found: "
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
                "M6.2C.1 configuration must be a YAML mapping."
            )

        required_sections = {
            "upstream_qc",
            "target_policy",
            "reconstruction_components",
            "local_search",
            "component_states",
            "classification",
            "decision_policy",
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
                "M6.2C.1 config missing sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Upstream
    # ======================================================================

    def _load_m6_2b(
        self,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        properties = config[
            "upstream_qc"
        ][
            "regulatory_gene_annotation"
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
                "M6.2B QC not found: "
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
                "Unexpected M6.2B milestone."
            )

        if payload.get(
            "stage"
        ) != properties[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected M6.2B stage."
            )

        unresolved = [
            row
            for row
            in payload.get(
                "feature_resolution",
                []
            )
            if bool(
                row.get(
                    "continue_to_m6_2c",
                    False,
                )
            )
        ]

        if not unresolved:

            raise RuntimeError(
                "M6.2B contains no feature requiring M6.2C."
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
            m6_2b_qc,
            m6_2b_path,
        ) = self._load_m6_2b(
            config
        )

        logger.info(
            "Loaded M6.2B regulatory-gene annotation QC."
        )

        result = assess_event_reconstruction_readiness(
            project_root=self.project_root,
            m6_2b_qc=m6_2b_qc,
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
            result.components,
            self.components_output,
            index=False,
        )

        write_parquet(
            result.feature_readiness,
            self.readiness_output,
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
                "M6.2C.1",

            "stage":
                "fase_event_reconstruction_feasibility_audit",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m6_2b_loaded":
                    True,

                "m6_2b_path":
                    str(
                        m6_2b_path
                    ),
            },

            "summary":
                summary,

            "reconstruction_components":
                _records_with_json_nulls(
                    result.components
                ),

            "feature_readiness":
                _records_with_json_nulls(
                    result.feature_readiness
                ),

            "outputs": {
                "components":
                    str(
                        self.components_output
                    ),

                "feature_readiness":
                    str(
                        self.readiness_output
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
            "M6.2C.1 complete."
        )

        logger.info(
            "Reconstruction classification: %s.",
            summary[
                "reconstruction_classification"
            ],
        )

        logger.info(
            "Required components available: %d/%d.",
            summary[
                "required_components_available"
            ],
            summary[
                "required_components"
            ],
        )

        logger.info(
            "Exact reconstruction allowed: %s.",
            summary[
                "exact_reconstruction_allowed"
            ],
        )

        logger.info(
            "Source-compatible reconstruction allowed: %s.",
            summary[
                "source_compatible_reconstruction_allowed"
            ],
        )

        for row in report[
            "feature_readiness"
        ]:

            logger.info(
                "%s | reconstruction=%s | "
                "gene assignment allowed=%s.",
                row[
                    "regulatory_feature_id"
                ],
                row[
                    "reconstruction_classification"
                ],
                row[
                    "parent_gene_assignment_from_reconstruction_allowed"
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