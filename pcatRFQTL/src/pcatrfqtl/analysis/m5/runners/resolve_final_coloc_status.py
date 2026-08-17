"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/resolve_final_coloc_status.py

Description:
    Runner for M5.6D final colocalization resolution.

    The runner integrates locked QC outputs from M5.5, M5.6A, M5.6B,
    and M5.6C.1 and produces the final colocalization decision.

    M5.6D does not re-run upstream analyses.

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

from pcatrfqtl.analysis.m5.final_coloc_resolution import (
    resolve_final_coloc_status,
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


class M56DFinalColocResolutionRunner:
    """Execute M5.6D final colocalization resolution."""

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
            / "final_coloc_evidence_chain.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "final_candidate_coloc_resolution.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "final_coloc_resolution_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_6d_final_coloc_resolution.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load M5.6D configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M5.6D configuration not found: "
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
                "M5.6D config must contain a YAML mapping."
            )

        required_sections = {
            "candidate_leads",
            "upstream_qc",
            "resolution_policy",
            "evidence_preservation",
            "scientific_policy",
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
                "M5.6D config missing required sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # JSON loading
    # ======================================================================

    @staticmethod
    def _load_json(
        path: Path,
    ) -> dict[str, Any]:
        """Load one upstream QC JSON."""

        if not path.exists():

            raise FileNotFoundError(
                f"Required upstream QC not found: {path}"
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
                f"QC payload is not a mapping: {path}"
            )

        return payload

    def _load_upstream_qc(
        self,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, dict[str, Any]],
        dict[str, str],
    ]:
        """Load all configured upstream locked QC artifacts."""

        payloads: dict[
            str,
            dict[str, Any],
        ] = {}

        paths: dict[
            str,
            str,
        ] = {}

        for stage_id, properties in config[
            "upstream_qc"
        ].items():

            relative_path = Path(
                str(
                    properties[
                        "relative_path"
                    ]
                )
            )

            absolute_path = (
                self.project_root
                /
                relative_path
            )

            required = bool(
                properties.get(
                    "required",
                    True,
                )
            )

            if (
                not absolute_path.exists()
                and
                not required
            ):

                continue

            payloads[
                str(
                    stage_id
                )
            ] = self._load_json(
                absolute_path
            )

            paths[
                str(
                    stage_id
                )
            ] = str(
                absolute_path
            )

            logger.info(
                "Loaded upstream QC %s: %s",
                stage_id,
                absolute_path,
            )

        return (
            payloads,
            paths,
        )

    # ======================================================================
    # Locked-state validation
    # ======================================================================

    @staticmethod
    def _validate_expected_states(
        *,
        config: dict[str, Any],
        payloads: dict[str, dict[str, Any]],
    ) -> list[str]:
        """
        Validate key upstream states.

        Returns human-readable validation messages.
        """

        messages: list[
            str
        ] = []

        # M5.5
        m55 = payloads[
            "m5_5"
        ]

        formal_ready_pairs = int(
            m55[
                "summary"
            ][
                "formal_coloc_ready_pairs"
            ]
        )

        if formal_ready_pairs != 0:

            raise RuntimeError(
                "M5.5 locked state changed: formal coloc-ready pairs "
                f"is {formal_ready_pairs}, expected 0."
            )

        messages.append(
            "M5.5 formal coloc-ready pairs = 0."
        )

        # M5.6A
        m56a = payloads[
            "m5_6a"
        ]

        moradi_state = str(
            m56a[
                "summary"
            ][
                "moradi_resolution_state"
            ]
        )

        expected_moradi = set(
            config[
                "expected_states"
            ][
                "m5_6a"
            ][
                "moradi_resolution_state"
            ]
        )

        if moradi_state not in expected_moradi:

            raise RuntimeError(
                "M5.6A locked state changed: "
                f"{moradi_state}"
            )

        messages.append(
            f"M5.6A Moradi resolution = {moradi_state}."
        )

        # M5.6B
        m56b = payloads[
            "m5_6b"
        ]

        reconstruction_state = str(
            m56b[
                "summary"
            ][
                "reconstruction_readiness"
            ]
        )

        expected_reconstruction = set(
            config[
                "expected_states"
            ][
                "m5_6b"
            ][
                "reconstruction_readiness"
            ]
        )

        if reconstruction_state not in expected_reconstruction:

            raise RuntimeError(
                "M5.6B locked state changed: "
                f"{reconstruction_state}"
            )

        messages.append(
            f"M5.6B reconstruction = {reconstruction_state}."
        )

        # M5.6C.1
        m56c1 = payloads[
            "m5_6c_1"
        ]

        alternative_state = str(
            m56c1[
                "summary"
            ][
                "classification"
            ]
        )

        expected_alternative = set(
            config[
                "expected_states"
            ][
                "m5_6c_1"
            ][
                "classification"
            ]
        )

        if alternative_state not in expected_alternative:

            raise RuntimeError(
                "M5.6C.1 locked state changed: "
                f"{alternative_state}"
            )

        messages.append(
            f"M5.6C.1 alternative resource = {alternative_state}."
        )

        return messages

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M5.6D."""

        config = self._load_config()

        (
            payloads,
            upstream_paths,
        ) = self._load_upstream_qc(
            config
        )

        required_stage_ids = {
            "m5_5",
            "m5_6a",
            "m5_6b",
            "m5_6c_1",
        }

        missing_stage_ids = (
            required_stage_ids
            -
            set(
                payloads
            )
        )

        if missing_stage_ids:

            raise RuntimeError(
                "Missing required upstream M5.6D QC payloads: "
                f"{sorted(missing_stage_ids)}"
            )

        validation_messages = (
            self._validate_expected_states(
                config=config,
                payloads=payloads,
            )
        )

        result = resolve_final_coloc_status(
            config=config,
            m5_5=payloads[
                "m5_5"
            ],
            m5_6a=payloads[
                "m5_6a"
            ],
            m5_6b=payloads[
                "m5_6b"
            ],
            m5_6c_1=payloads[
                "m5_6c_1"
            ],
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
            result.evidence_chain,
            self.evidence_output,
            index=False,
        )

        write_parquet(
            result.candidates,
            self.candidate_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        summary_row = (
            result.summary.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M5.6D",

            "stage":
                "final_colocalization_resolution",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "all_required_qc_loaded":
                    True,

                "locked_states_validated":
                    True,

                "messages":
                    validation_messages,
            },

            "summary": {
                "candidate_leads_assessed":
                    int(
                        summary_row[
                            "candidate_leads_assessed"
                        ]
                    ),

                "formal_coloc_routes_assessed":
                    int(
                        summary_row[
                            "formal_coloc_routes_assessed"
                        ]
                    ),

                "formal_coloc_routes_available":
                    int(
                        summary_row[
                            "formal_coloc_routes_available"
                        ]
                    ),

                "valid_dense_qtl_route_available":
                    bool(
                        summary_row[
                            "valid_dense_qtl_route_available"
                        ]
                    ),

                "formal_coloc_ready_candidates":
                    int(
                        summary_row[
                            "formal_coloc_ready_candidates"
                        ]
                    ),

                "formal_coloc_blocked_candidates":
                    int(
                        summary_row[
                            "formal_coloc_blocked_candidates"
                        ]
                    ),

                "final_coloc_status":
                    str(
                        summary_row[
                            "final_coloc_status"
                        ]
                    ),
            },

            "evidence_chain":
                result.evidence_chain.to_dict(
                    orient="records"
                ),

            "candidate_status":
                result.candidates.to_dict(
                    orient="records"
                ),

            "upstream_qc":
                upstream_paths,

            "outputs": {
                "evidence_chain":
                    str(
                        self.evidence_output
                    ),

                "candidate_resolution":
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
            )

        logger.info(
            "M5.6D complete."
        )

        logger.info(
            "Final colocalization status: %s.",
            report[
                "summary"
            ][
                "final_coloc_status"
            ],
        )

        logger.info(
            "Formal-coloc routes available: %d/%d.",
            report[
                "summary"
            ][
                "formal_coloc_routes_available"
            ],
            report[
                "summary"
            ][
                "formal_coloc_routes_assessed"
            ],
        )

        logger.info(
            "Formal coloc-ready candidates: %d/%d.",
            report[
                "summary"
            ][
                "formal_coloc_ready_candidates"
            ],
            report[
                "summary"
            ][
                "candidate_leads_assessed"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report