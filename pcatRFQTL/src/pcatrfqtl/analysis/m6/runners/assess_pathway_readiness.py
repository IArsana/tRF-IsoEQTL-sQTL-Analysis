"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/assess_pathway_readiness.py

Description:
    Runner for M6.4B Pathway Readiness Assessment.

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

from pcatrfqtl.analysis.m6.pathway_readiness import (
    assess_pathway_readiness,
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
    """Convert pandas missing values to JSON-compatible null."""

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
    """Convert summary Series to JSON-native values."""

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


class M64BPathwayReadinessRunner:
    """Execute M6.4B Pathway Readiness Assessment."""

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
    def gene_set_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "defensible_gene_set.parquet"
        )

    @property
    def candidate_readiness_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_pathway_readiness.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "pathway_readiness_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m6_4b_pathway_readiness.json"
        )

    # ======================================================================
    # Config
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M6.4B config not found: "
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
                "M6.4B configuration must contain a YAML mapping."
            )

        return config

    # ======================================================================
    # Upstream
    # ======================================================================

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
                f"M6.4B upstream QC missing: {path}"
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
            m6_3_qc,
            m6_3_path,
        ) = self._load_upstream_qc(
            config=config,
            key="candidate_gene_integration",
        )

        (
            m6_4a_qc,
            m6_4a_path,
        ) = self._load_upstream_qc(
            config=config,
            key="trf_biological_annotation",
        )

        logger.info(
            "Loaded locked M6.3 and M6.4A QC artifacts."
        )

        result = assess_pathway_readiness(
            m6_3_qc=m6_3_qc,
            m6_4a_qc=m6_4a_qc,
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
            result.defensible_gene_set,
            self.gene_set_output,
            index=False,
        )

        write_parquet(
            result.candidate_readiness,
            self.candidate_readiness_output,
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
                "M6.4B",

            "stage":
                "pathway_readiness_assessment",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m6_3_loaded":
                    True,

                "m6_3_path":
                    str(
                        m6_3_path
                    ),

                "m6_4a_loaded":
                    True,

                "m6_4a_path":
                    str(
                        m6_4a_path
                    ),
            },

            "summary":
                summary,

            "defensible_gene_set":
                _records_with_json_nulls(
                    result.defensible_gene_set
                ),

            "candidate_pathway_readiness":
                _records_with_json_nulls(
                    result.candidate_readiness
                ),

            "outputs": {
                "defensible_gene_set":
                    str(
                        self.gene_set_output
                    ),

                "candidate_pathway_readiness":
                    str(
                        self.candidate_readiness_output
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
            "M6.4B complete."
        )

        logger.info(
            "Candidates=%d | resolved genes=%d | eligible genes=%d | "
            "unique genes=%d.",
            summary[
                "candidates_assessed"
            ],
            summary[
                "resolved_candidate_genes"
            ],
            summary[
                "eligible_candidate_genes"
            ],
            summary[
                "unique_defensible_genes"
            ],
        )

        logger.info(
            "Pathway enrichment ready=%s | status=%s.",
            summary[
                "pathway_enrichment_ready"
            ],
            summary[
                "pathway_readiness_status"
            ],
        )

        for row in report[
            "candidate_pathway_readiness"
        ]:

            logger.info(
                "%s | gene=%s | eligible=%s | exclusion=%s.",
                row[
                    "lead_rsid"
                ],
                row[
                    "integrated_gene"
                ],
                row[
                    "eligible_for_pathway_gene_set"
                ],
                row[
                    "pathway_exclusion_reason"
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