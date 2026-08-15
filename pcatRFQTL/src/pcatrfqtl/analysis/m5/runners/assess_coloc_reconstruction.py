"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/assess_coloc_reconstruction.py

Description:
    Runner for M5.6B colocalization reconstruction feasibility audit.

    Outputs:
        reconstruction_components.parquet
        reconstruction_readiness.parquet
        candidate_reconstruction_readiness.parquet
        m5_6b_coloc_reconstruction.json

    M5.6B is an audit only. It does not download TCGA data, reconstruct
    expression phenotypes, execute MatrixEQTL, or perform colocalization.

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

from pcatrfqtl.analysis.m5.coloc_reconstruction_readiness import (
    assess_coloc_reconstruction_readiness,
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


# ============================================================================
# Runner
# ============================================================================


class M56BColocReconstructionRunner:
    """Execute M5.6B reconstruction feasibility assessment."""

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
            / "reconstruction_readiness.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_reconstruction_readiness.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_6b_coloc_reconstruction.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M5.6B YAML configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M5.6B config not found: "
                f"{self.config_path}"
            )

        if not self.config_path.is_file():

            raise RuntimeError(
                "M5.6B config path is not a file: "
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
                "M5.6B config must contain a YAML mapping."
            )

        required_sections = {
            "milestone",
            "candidate_leads",
            "target_analysis",
            "components",
            "readiness_policy",
            "scientific_policy",
        }

        missing = (
            required_sections
            - set(
                config
            )
        )

        if missing:

            raise ValueError(
                "M5.6B config missing required sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Run M5.6B."""

        config = self._load_config()

        result = (
            assess_coloc_reconstruction_readiness(
                config
            )
        )

        components = result.components

        readiness = result.summary

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
            components,
            self.components_output,
            index=False,
        )

        write_parquet(
            readiness,
            self.readiness_output,
            index=False,
        )

        write_parquet(
            candidates,
            self.candidate_output,
            index=False,
        )

        summary_row = readiness.iloc[
            0
        ]

        readiness_state = str(
            summary_row[
                "reconstruction_readiness"
            ]
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

        component_status_counts = {
            str(
                key
            ):
                int(
                    value
                )
            for key, value
            in components[
                "status"
            ]
            .value_counts(
                dropna=False
            )
            .items()
        }

        report = {
            "milestone":
                "M5.6B",

            "stage":
                "colocalization_reconstruction_feasibility_audit",

            "policy": {
                "reconstruction_executed":
                    False,

                "matrixeqtl_executed":
                    False,

                "missing_original_parameters_inferred":
                    False,

                "unknown_covariates_fabricated":
                    False,

                "approximate_reconstruction_treated_as_exact":
                    False,

                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "components_assessed":
                    int(
                        len(
                            components
                        )
                    ),

                "required_components":
                    int(
                        summary_row[
                            "required_component_count"
                        ]
                    ),

                "fully_available_components":
                    int(
                        summary_row[
                            "fully_available_component_count"
                        ]
                    ),

                "nonblocking_components":
                    int(
                        summary_row[
                            "nonblocking_component_count"
                        ]
                    ),

                "not_identified_components":
                    int(
                        summary_row[
                            "not_identified_component_count"
                        ]
                    ),

                "exact_availability_fraction":
                    float(
                        summary_row[
                            "exact_availability_fraction"
                        ]
                    ),

                "approximate_nonblocking_fraction":
                    float(
                        summary_row[
                            "approximate_nonblocking_fraction"
                        ]
                    ),

                "reconstruction_readiness":
                    readiness_state,

                "reconstruction_can_proceed":
                    bool(
                        summary_row[
                            "reconstruction_can_proceed"
                        ]
                    ),

                "exact_replication_claim_allowed":
                    bool(
                        summary_row[
                            "exact_replication_claim_allowed"
                        ]
                    ),

                "candidate_leads_assessed":
                    int(
                        len(
                            candidates
                        )
                    ),

                "formal_coloc_ready_candidates":
                    formal_ready_candidates,

                "formal_coloc_blocked_candidates":
                    int(
                        len(
                            candidates
                        )
                        -
                        formal_ready_candidates
                    ),
            },

            "component_status_counts":
                component_status_counts,

            "blocking_components":
                summary_row[
                    "blocking_components"
                ],

            "partial_components":
                summary_row[
                    "partial_components"
                ],

            "derivable_components":
                summary_row[
                    "derivable_components"
                ],

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
                        "reconstruction_readiness"
                    ]
                    .value_counts(
                        dropna=False
                    )
                    .items()
                },

            "outputs": {
                "components":
                    str(
                        self.components_output
                    ),

                "readiness":
                    str(
                        self.readiness_output
                    ),

                "candidate_readiness":
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
            "M5.6B complete."
        )

        logger.info(
            "Reconstruction readiness: %s.",
            readiness_state,
        )

        logger.info(
            "Required components: %d | "
            "available=%d | nonblocking=%d | not identified=%d.",
            report[
                "summary"
            ][
                "required_components"
            ],
            report[
                "summary"
            ][
                "fully_available_components"
            ],
            report[
                "summary"
            ][
                "nonblocking_components"
            ],
            report[
                "summary"
            ][
                "not_identified_components"
            ],
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