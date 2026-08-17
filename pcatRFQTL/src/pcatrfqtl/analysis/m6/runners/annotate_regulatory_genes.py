"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/annotate_regulatory_genes.py

Description:
    Runner for M6.2B Regulatory Feature → Source Gene Annotation.

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

from pcatrfqtl.analysis.m6.regulatory_gene_annotation import (
    resolve_source_gene_annotations,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(
    __name__
)


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


class M62BRegulatoryGeneAnnotationRunner:
    """Execute M6.2B."""

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
    def annotation_candidates_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "source_annotation_candidates.parquet"
        )

    @property
    def feature_resolution_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "regulatory_feature_gene_resolution.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "regulatory_gene_annotation_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m6_2b_regulatory_gene_annotation.json"
        )

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
                "M6.2B configuration must be a YAML mapping."
            )

        return config

    def _load_m6_2a(
        self,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        properties = config[
            "upstream_qc"
        ][
            "regulatory_feature_mapping"
        ]

        path = (
            self.project_root
            /
            properties[
                "relative_path"
            ]
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M6.2A QC not found: {path}"
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
                "Unexpected M6.2A upstream milestone."
            )

        if payload.get(
            "stage"
        ) != properties[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected M6.2A upstream stage."
            )

        return payload, path

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        (
            m6_2a_qc,
            m6_2a_path,
        ) = self._load_m6_2a(
            config
        )

        result = resolve_source_gene_annotations(
            project_root=self.project_root,
            m6_2a_qc=m6_2a_qc,
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
            result.annotation_candidates,
            self.annotation_candidates_output,
            index=False,
        )

        write_parquet(
            result.feature_resolution,
            self.feature_resolution_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        summary = result.summary.iloc[
            0
        ]

        report = {
            "milestone":
                "M6.2B",

            "stage":
                "regulatory_feature_to_source_gene_annotation",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m6_2a_loaded":
                    True,

                "m6_2a_path":
                    str(
                        m6_2a_path
                    ),
            },

            "summary": {
                key:
                    (
                        value.item()
                        if hasattr(
                            value,
                            "item"
                        )
                        else value
                    )
                for key, value
                in summary.to_dict().items()
            },

            "annotation_candidates":
                _records_with_json_nulls(
                    result.annotation_candidates
                ),

            "feature_resolution":
                _records_with_json_nulls(
                    result.feature_resolution
                ),

            "outputs": {
                "annotation_candidates":
                    str(
                        self.annotation_candidates_output
                    ),

                "feature_resolution":
                    str(
                        self.feature_resolution_output
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
            "M6.2B complete."
        )

        logger.info(
            "Features assessed=%d | resolved=%d | unresolved=%d.",
            report[
                "summary"
            ][
                "features_assessed"
            ],
            report[
                "summary"
            ][
                "features_source_gene_resolved"
            ],
            report[
                "summary"
            ][
                "features_source_gene_unresolved"
            ],
        )

        for row in report[
            "feature_resolution"
        ]:

            logger.info(
                "%s | gene=%s | status=%s.",
                row[
                    "regulatory_feature_id"
                ],
                row[
                    "parent_gene"
                ],
                row[
                    "gene_annotation_status"
                ],
            )

        logger.info(
            "Next stage: %s.",
            report[
                "summary"
            ][
                "next_stage"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report