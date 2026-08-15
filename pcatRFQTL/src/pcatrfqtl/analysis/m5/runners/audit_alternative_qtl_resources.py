"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/audit_alternative_qtl_resources.py

Description:
    Runner for M5.6C alternative PRAD QTL resource audit.

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

from pcatrfqtl.analysis.m5.alternative_qtl_audit import (
    audit_alternative_qtl_resources,
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


class M56CAlternativeQTLAuditRunner:
    """Execute M5.6C alternative-QTL resource audit."""

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
    def resources_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "alternative_qtl_resources.parquet"
        )

    @property
    def requirements_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "formal_coloc_requirements.parquet"
        )

    @property
    def candidates_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_alternative_qtl_resolution.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_6c_alternative_qtl_audit.json"
        )

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load M5.6C YAML."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M5.6C config not found: {self.config_path}"
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
                "M5.6C configuration must contain a YAML mapping."
            )

        return config

    def run(
        self,
    ) -> dict[str, Any]:
        """Run M5.6C."""

        config = self._load_config()

        result = audit_alternative_qtl_resources(
            config
        )

        resources = result.resources

        requirements = result.requirements

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
            resources,
            self.resources_output,
            index=False,
        )

        write_parquet(
            requirements,
            self.requirements_output,
            index=False,
        )

        write_parquet(
            candidates,
            self.candidates_output,
            index=False,
        )

        abf_ready_resources = int(
            resources[
                "abf_formal_coloc_ready"
            ]
            .fillna(
                False
            )
            .sum()
        )

        susie_ready_resources = int(
            resources[
                "susie_formal_coloc_ready"
            ]
            .fillna(
                False
            )
            .sum()
        )

        inspection_resources = int(
            resources[
                "requires_dataset_inspection"
            ]
            .fillna(
                False
            )
            .sum()
        )

        formal_ready_candidates = int(
            candidates[
                "formal_coloc_ready"
            ]
            .fillna(
                False
            )
            .sum()
        )

        report = {
            "milestone":
                "M5.6C",

            "stage":
                "alternative_prad_qtl_resource_audit",

            "policy": {
                "documentation_only_can_establish_coloc_readiness":
                    False,

                "direct_dataset_inspection_required":
                    True,

                "significant_only_qtl_allowed":
                    False,

                "missing_qtl_rows_interpreted_as_null":
                    False,

                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "alternative_resources_assessed":
                    int(
                        len(
                            resources
                        )
                    ),

                "abf_ready_resources":
                    abf_ready_resources,

                "susie_ready_resources":
                    susie_ready_resources,

                "resources_requiring_dataset_inspection":
                    inspection_resources,

                "candidate_leads_assessed":
                    int(
                        len(
                            candidates
                        )
                    ),

                "formal_coloc_ready_candidates":
                    formal_ready_candidates,

                "formal_coloc_blocked_or_pending_candidates":
                    int(
                        len(
                            candidates
                        )
                        -
                        formal_ready_candidates
                    ),
            },

            "resource_classification_counts":
                {
                    str(
                        key
                    ):
                        int(
                            value
                        )
                    for key, value
                    in resources[
                        "classification"
                    ]
                    .value_counts(
                        dropna=False
                    )
                    .items()
                },

            "candidate_resolution_counts":
                {
                    str(
                        key
                    ):
                        int(
                            value
                        )
                    for key, value
                    in candidates[
                        "resolution_status"
                    ]
                    .value_counts(
                        dropna=False
                    )
                    .items()
                },

            "outputs": {
                "resources":
                    str(
                        self.resources_output
                    ),

                "requirements":
                    str(
                        self.requirements_output
                    ),

                "candidate_resolution":
                    str(
                        self.candidates_output
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
            "M5.6C complete."
        )

        logger.info(
            "Alternative resources: %d | ABF-ready=%d | "
            "SuSiE-ready=%d | dataset inspection required=%d.",
            len(
                resources
            ),
            abf_ready_resources,
            susie_ready_resources,
            inspection_resources,
        )

        logger.info(
            "Formal coloc-ready candidates: %d/%d.",
            formal_ready_candidates,
            len(
                candidates
            ),
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report