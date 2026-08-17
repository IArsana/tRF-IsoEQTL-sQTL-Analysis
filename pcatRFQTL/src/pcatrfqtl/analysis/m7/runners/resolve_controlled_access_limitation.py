"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/resolve_controlled_access_limitation.py

Description:
    Runner for M7.3C Controlled-Access Limitation Resolution.

    This stage closes the direct TCGA tRF quantification branch when
    controlled-access aligned miRNA-Seq BAM resources exist but the required
    authorization is not available at analysis time.

    This limitation is documented as a data-access constraint and must not be
    interpreted as:
        - absence of the underlying TCGA sequencing resource;
        - negative tRF expression;
        - failed candidate validation;
        - evidence against the candidate.

    This runner does NOT:
        - download controlled BAM files;
        - use a GDC authentication token;
        - perform read counting;
        - quantify tRF expression;
        - substitute processed miRNA expression for tRF expression;
        - perform clinical association;
        - perform causal inference.

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

from pcatrfqtl.analysis.m7.controlled_access_limitation import (
    resolve_controlled_access_limitation,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


# ============================================================================
# JSON helpers
# ============================================================================


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """
    Convert a DataFrame into JSON-safe records.

    Pandas missing values are converted to JSON null.
    """

    if dataframe.empty:
        return []

    cleaned = dataframe.astype(object).where(
        pd.notna(dataframe),
        None,
    )

    records = cleaned.to_dict(
        orient="records"
    )

    for record in records:
        for key, value in list(record.items()):
            if hasattr(value, "item"):
                record[key] = value.item()

    return records


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """
    Convert one pandas Series into JSON-native scalar values.
    """

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

        if hasattr(value, "item"):
            value = value.item()

        output[str(key)] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M73CControlledAccessLimitationRunner:
    """
    Execute M7.3C Controlled-Access Limitation Resolution.

    M7.3C consumes the locked M7.3B quantification-feasibility artifact
    and records that direct TCGA tRF quantification is deferred because
    controlled-access authorization is not available at analysis time.
    """

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:
        self.project_root = Path(project_root)
        self.config_path = Path(config_path)
        self.output_directory = Path(output_directory)
        self.qc_directory = Path(qc_directory)

    # ======================================================================
    # Output paths
    # ======================================================================

    @property
    def candidate_output(
        self,
    ) -> Path:
        """Candidate-level M7.3C output."""

        return (
            self.output_directory
            / "candidate_controlled_access_limitation.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:
        """One-row M7.3C summary output."""

        return (
            self.output_directory
            / "controlled_access_limitation_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:
        """M7.3C QC JSON output."""

        return (
            self.qc_directory
            / "m7_3c_controlled_access_limitation.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M7.3C YAML configuration."""

        if not self.config_path.exists():
            raise FileNotFoundError(
                f"M7.3C config not found: {self.config_path}"
            )

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:
            config = yaml.safe_load(handle)

        if not isinstance(config, dict):
            raise RuntimeError(
                "M7.3C configuration must contain a YAML mapping."
            )

        required_sections = {
            "milestone",
            "stage",
            "upstream_qc",
            "access_state",
            "interpretation",
            "tcga_expression_policy",
            "clinical_association_policy",
            "statuses",
            "scientific_policy",
            "next_stage",
        }

        missing_sections = (
            required_sections
            - set(config.keys())
        )

        if missing_sections:
            raise RuntimeError(
                "M7.3C configuration missing required sections: "
                f"{sorted(missing_sections)}"
            )

        if config.get("milestone") != "M7.3C":
            raise RuntimeError(
                "M7.3C configuration has unexpected milestone: "
                f"{config.get('milestone')!r}"
            )

        if (
            config.get("stage")
            != "controlled_access_limitation_resolution"
        ):
            raise RuntimeError(
                "M7.3C configuration has unexpected stage: "
                f"{config.get('stage')!r}"
            )

        return config

    # ======================================================================
    # Upstream artifact
    # ======================================================================

    def _load_upstream_qc(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """
        Load and validate the locked M7.3B QC artifact.
        """

        upstream = config[
            "upstream_qc"
        ][
            "quantification_feasibility"
        ]

        path = (
            self.project_root
            / str(
                upstream[
                    "relative_path"
                ]
            )
        )

        if not path.exists():
            raise FileNotFoundError(
                f"M7.3B QC artifact not found: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:
            payload = json.load(handle)

        if not isinstance(payload, dict):
            raise RuntimeError(
                "M7.3B QC artifact must contain a JSON object."
            )

        observed_milestone = payload.get(
            "milestone"
        )

        expected_milestone = upstream[
            "expected_milestone"
        ]

        if observed_milestone != expected_milestone:
            raise RuntimeError(
                "Unexpected upstream milestone: "
                f"{observed_milestone!r}; "
                f"expected {expected_milestone!r}."
            )

        observed_stage = payload.get(
            "stage"
        )

        expected_stage = upstream[
            "expected_stage"
        ]

        if observed_stage != expected_stage:
            raise RuntimeError(
                "Unexpected upstream stage: "
                f"{observed_stage!r}; "
                f"expected {expected_stage!r}."
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
            and observed_status != expected_status
        ):
            raise RuntimeError(
                "Unexpected M7.3B feasibility status: "
                f"{observed_status!r}; "
                f"expected {expected_status!r}."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Declared access state
    # ======================================================================

    def _validate_declared_access_state(
        self,
        *,
        config: dict[str, Any],
    ) -> None:
        """
        Validate that M7.3C represents a documented access limitation.

        This stage is valid only when controlled resources exist but the
        researcher does not hold the authorization required to access them.
        """

        access_state = config[
            "access_state"
        ]

        controlled_resource_identified = bool(
            access_state.get(
                "controlled_resource_identified",
                False,
            )
        )

        authorized_access_available = bool(
            access_state.get(
                "authorized_access_available_at_analysis_time",
                False,
            )
        )

        dbgap_authorization_available = bool(
            access_state.get(
                "dbgap_project_authorization_available",
                False,
            )
        )

        limitation_confirmed = bool(
            access_state.get(
                "access_limitation_confirmed_by_researcher",
                False,
            )
        )

        if not controlled_resource_identified:
            raise RuntimeError(
                "M7.3C cannot resolve a controlled-access limitation "
                "when controlled_resource_identified=False."
            )

        if authorized_access_available:
            raise RuntimeError(
                "M7.3C limitation route is invalid because authorized "
                "controlled access is declared available."
            )

        if dbgap_authorization_available:
            raise RuntimeError(
                "M7.3C limitation route is invalid because dbGaP project "
                "authorization is declared available."
            )

        if not limitation_confirmed:
            raise RuntimeError(
                "M7.3C requires access_limitation_confirmed_by_researcher=True."
            )

    # ======================================================================
    # Upstream consistency
    # ======================================================================

    def _validate_upstream_resource_state(
        self,
        *,
        m7_3b_qc: dict[str, Any],
    ) -> dict[str, int]:
        """
        Verify that locked M7.3B actually contains a controlled-access route.
        """

        summary = m7_3b_qc.get(
            "summary",
            {},
        )

        if not isinstance(summary, dict):
            raise RuntimeError(
                "M7.3B QC summary must be a JSON object."
            )

        aligned_bam_count = int(
            summary.get(
                "aligned_mirna_bam_files",
                0,
            )
        )

        aligned_bam_case_count = int(
            summary.get(
                "aligned_mirna_bam_cases",
                0,
            )
        )

        controlled_candidate_count = int(
            summary.get(
                "candidates_requiring_controlled_access",
                0,
            )
        )

        directly_quantifiable_count = int(
            summary.get(
                "candidates_directly_quantifiable_now",
                0,
            )
        )

        if aligned_bam_count <= 0:
            raise RuntimeError(
                "M7.3C expected aligned miRNA BAM resources from M7.3B, "
                "but aligned_mirna_bam_files <= 0."
            )

        if controlled_candidate_count <= 0:
            raise RuntimeError(
                "M7.3C expected at least one candidate requiring "
                "controlled access according to M7.3B."
            )

        if directly_quantifiable_count != 0:
            raise RuntimeError(
                "M7.3C limitation route expected zero candidates directly "
                "quantifiable without controlled access."
            )

        return {
            "aligned_bam_count":
                aligned_bam_count,

            "aligned_bam_case_count":
                aligned_bam_case_count,

            "controlled_candidate_count":
                controlled_candidate_count,

            "directly_quantifiable_count":
                directly_quantifiable_count,
        }

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """
        Execute M7.3C Controlled-Access Limitation Resolution.
        """

        # ------------------------------------------------------------------
        # Configuration
        # ------------------------------------------------------------------

        config = self._load_config()

        self._validate_declared_access_state(
            config=config,
        )

        # ------------------------------------------------------------------
        # Locked upstream M7.3B
        # ------------------------------------------------------------------

        (
            m7_3b_qc,
            m7_3b_path,
        ) = self._load_upstream_qc(
            config=config,
        )

        logger.info(
            "Loaded locked M7.3B quantification-feasibility artifact."
        )

        upstream_state = self._validate_upstream_resource_state(
            m7_3b_qc=m7_3b_qc,
        )

        # ------------------------------------------------------------------
        # Core resolution
        # ------------------------------------------------------------------

        result = resolve_controlled_access_limitation(
            m7_3b_qc=m7_3b_qc,
            config=config,
        )

        if result.summary.empty:
            raise RuntimeError(
                "M7.3C core returned an empty summary."
            )

        if result.candidate_resolution.empty:
            raise RuntimeError(
                "M7.3C core returned no candidate-resolution rows."
            )

        # ------------------------------------------------------------------
        # Output directories
        # ------------------------------------------------------------------

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ------------------------------------------------------------------
        # Parquet outputs
        # ------------------------------------------------------------------

        write_parquet(
            result.candidate_resolution,
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

        # ------------------------------------------------------------------
        # Final QC report
        # ------------------------------------------------------------------

        report = {
            "milestone":
                "M7.3C",

            "stage":
                "controlled_access_limitation_resolution",

            "access_state": {
                "controlled_resource_identified":
                    True,

                "authorized_access_available_at_analysis_time":
                    False,

                "gdc_token_available":
                    bool(
                        config[
                            "access_state"
                        ].get(
                            "gdc_token_available",
                            False,
                        )
                    ),

                "dbgap_project_authorization_available":
                    False,

                "access_limitation_confirmed_by_researcher":
                    True,
            },

            "interpretation":
                config[
                    "interpretation"
                ],

            "tcga_expression_policy":
                config[
                    "tcga_expression_policy"
                ],

            "clinical_association_policy":
                config[
                    "clinical_association_policy"
                ],

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m7_3b_loaded":
                    True,

                "m7_3b_path":
                    str(
                        m7_3b_path
                    ),

                "m7_3b_expected_status":
                    (
                        config[
                            "upstream_qc"
                        ][
                            "quantification_feasibility"
                        ].get(
                            "expected_overall_status"
                        )
                    ),

                "m7_3b_observed_status":
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

                "aligned_mirna_bam_files":
                    upstream_state[
                        "aligned_bam_count"
                    ],

                "aligned_mirna_bam_cases":
                    upstream_state[
                        "aligned_bam_case_count"
                    ],

                "candidates_requiring_controlled_access":
                    upstream_state[
                        "controlled_candidate_count"
                    ],

                "candidates_directly_quantifiable_now":
                    upstream_state[
                        "directly_quantifiable_count"
                    ],
            },

            "summary":
                summary,

            "candidate_controlled_access_limitation":
                _records_with_json_nulls(
                    result.candidate_resolution
                ),

            "outputs": {
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

        # ------------------------------------------------------------------
        # QC JSON
        # ------------------------------------------------------------------

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

        # ------------------------------------------------------------------
        # Logging
        # ------------------------------------------------------------------

        logger.info(
            "M7.3C complete."
        )

        logger.info(
            "Candidates=%d | access-limited=%d | sequence-limited=%d | "
            "TCGA quantification deferred=%d.",
            summary[
                "candidate_trfs_assessed"
            ],
            summary[
                "candidates_access_limited"
            ],
            summary[
                "candidates_sequence_limited"
            ],
            summary[
                "candidates_tcga_quantification_deferred"
            ],
        )

        for row in report[
            "candidate_controlled_access_limitation"
        ]:

            logger.info(
                "%s | %s | sequence_ready=%s | "
                "resource_available=%s | access_limited=%s | "
                "sequence_limited=%s | deferred=%s | status=%s.",
                row[
                    "lead_rsid"
                ],
                row[
                    "trf_id"
                ],
                row[
                    "sequence_ready_for_quantification"
                ],
                row[
                    "resource_technically_available"
                ],
                row[
                    "access_limited"
                ],
                row[
                    "sequence_limited"
                ],
                row[
                    "tcga_direct_quantification_deferred"
                ],
                row[
                    "final_status"
                ],
            )

        logger.info(
            "Overall resolution=%s.",
            summary[
                "overall_resolution_status"
            ],
        )

        logger.info(
            "External validation next=%s.",
            summary[
                "next_external_validation_stage"
            ],
        )

        logger.info(
            "Variant validation next=%s.",
            summary[
                "next_variant_validation_stage"
            ],
        )

        logger.info(
            "Clinical route=%s.",
            summary[
                "clinical_route"
            ],
        )

        logger.info(
            "QC output=%s.",
            self.qc_output,
        )

        return report