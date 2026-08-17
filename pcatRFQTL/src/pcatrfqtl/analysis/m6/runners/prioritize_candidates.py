"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/prioritize_candidates.py

Description:
    Runner for M6.1 final candidate prioritization.

    The runner consumes the locked M5.6D final colocalization QC artifact,
    combines it with the deterministic M6.1 candidate evidence definition,
    and generates the final candidate interpretation priority.

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

from pcatrfqtl.analysis.m6.candidate_prioritization import (
    run_candidate_prioritization,
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


class M61CandidatePrioritizationRunner:
    """Execute M6.1 final candidate prioritization."""

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
    def evidence_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_evidence_matrix.parquet"
        )

    @property
    def priority_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "final_candidate_priority.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_priority_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m6_1_candidate_prioritization.json"
        )

    # ======================================================================
    # Config
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M6.1 configuration not found: "
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
                "M6.1 configuration must be a YAML mapping."
            )

        required = {
            "candidate_leads",
            "upstream_qc",
            "priority_hierarchy",
            "priority_classes",
            "interpretation_policy",
            "scientific_policy",
        }

        missing = (
            required
            -
            set(
                config
            )
        )

        if missing:

            raise ValueError(
                "M6.1 configuration missing sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Upstream M5.6D
    # ======================================================================

    def _load_final_coloc_qc(
        self,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        properties = config[
            "upstream_qc"
        ][
            "final_coloc_resolution"
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
                "Locked M5.6D QC not found: "
                f"{path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        if not isinstance(
            payload,
            dict,
        ):

            raise ValueError(
                "M5.6D QC payload must be a JSON object."
            )

        final_status = (
            payload.get(
                "summary",
                {},
            ).get(
                "final_coloc_status"
            )
        )

        if final_status != (
            "FORMAL_COLOC_NOT_JUSTIFIED_DATA_LIMITATION"
        ):

            raise RuntimeError(
                "Unexpected M5.6D final status: "
                f"{final_status}"
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
            final_coloc_qc,
            final_coloc_path,
        ) = self._load_final_coloc_qc(
            config
        )

        result = run_candidate_prioritization(
            config=config,
            final_coloc_qc=final_coloc_qc,
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
            self.evidence_output,
            index=False,
        )

        write_parquet(
            result.priorities,
            self.priority_output,
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
                "M6.1",

            "stage":
                "final_candidate_prioritization",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m5_6d_loaded":
                    True,

                "m5_6d_path":
                    str(
                        final_coloc_path
                    ),

                "m5_6d_final_status":
                    final_coloc_qc[
                        "summary"
                    ][
                        "final_coloc_status"
                    ],
            },

            "summary": {
                "candidates_assessed":
                    int(
                        summary[
                            "candidates_assessed"
                        ]
                    ),

                "candidates_retained":
                    int(
                        summary[
                            "candidates_retained"
                        ]
                    ),

                "highest_priority_candidate":
                    str(
                        summary[
                            "highest_priority_candidate"
                        ]
                    ),

                "highest_priority_class":
                    str(
                        summary[
                            "highest_priority_class"
                        ]
                    ),

                "formal_coloc_ready_candidates":
                    int(
                        summary[
                            "formal_coloc_ready_candidates"
                        ]
                    ),

                "causal_claim_allowed_candidates":
                    int(
                        summary[
                            "causal_claim_allowed_candidates"
                        ]
                    ),

                "ranking_method":
                    str(
                        summary[
                            "ranking_method"
                        ]
                    ),

                "additive_score_used":
                    bool(
                        summary[
                            "additive_score_used"
                        ]
                    ),

                "candidate_exclusion_performed":
                    bool(
                        summary[
                            "candidate_exclusion_performed"
                        ]
                    ),

                "next_stage":
                    str(
                        summary[
                            "next_stage"
                        ]
                    ),
            },

            "candidate_priority":
                result.priorities.to_dict(
                    orient="records"
                ),

            "outputs": {
                "evidence_matrix":
                    str(
                        self.evidence_output
                    ),

                "priority":
                    str(
                        self.priority_output
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
            )

        logger.info(
            "M6.1 complete."
        )

        logger.info(
            "Candidates assessed=%d | retained=%d.",
            report[
                "summary"
            ][
                "candidates_assessed"
            ],
            report[
                "summary"
            ][
                "candidates_retained"
            ],
        )

        logger.info(
            "Highest-priority candidate=%s (%s).",
            report[
                "summary"
            ][
                "highest_priority_candidate"
            ],
            report[
                "summary"
            ][
                "highest_priority_class"
            ],
        )

        for row in report[
            "candidate_priority"
        ]:

            logger.info(
                "Rank %d: %s | %s | %s.",
                row[
                    "priority_rank"
                ],
                row[
                    "lead_rsid"
                ],
                row[
                    "priority_class"
                ],
                row[
                    "disease_evidence_class"
                ],
            )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report