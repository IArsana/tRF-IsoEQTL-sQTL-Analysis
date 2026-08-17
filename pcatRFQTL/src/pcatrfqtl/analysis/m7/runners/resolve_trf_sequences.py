"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/resolve_trf_sequences.py

Description:
    Runner for M7.3A Candidate tRF Sequence Resolution.

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

from pcatrfqtl.analysis.m7.trf_sequence_resolution import (
    resolve_trf_sequences,
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
    """Convert DataFrame into JSON-safe records."""

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
                record[
                    key
                ] = value.item()

    return records


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert pandas row into JSON-native values."""

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


class M73ATrfSequenceResolutionRunner:
    """Execute M7.3A Candidate tRF Sequence Resolution."""

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
            / "trf_sequence_evidence.parquet"
        )

    @property
    def resolution_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_trf_sequence_resolution.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "trf_sequence_resolution_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m7_3a_trf_sequence_resolution.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.3A config not found: {self.config_path}"
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
                "M7.3A configuration must contain a YAML mapping."
            )

        required = {
            "upstream_qc",
            "sequence_evidence",
            "sequence_validation",
            "resolution_policy",
            "quantification_policy",
            "statuses",
            "scientific_policy",
            "next_stage",
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
                "M7.3A config missing required sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Upstream QC
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
                f"M7.3A upstream QC missing: {path}"
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
                f"Unexpected milestone for M7.3A upstream {key}."
            )

        if payload.get(
            "stage"
        ) != properties[
            "expected_stage"
        ]:

            raise RuntimeError(
                f"Unexpected stage for M7.3A upstream {key}."
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
        """Execute M7.3A."""

        config = self._load_config()

        (
            m7_1_qc,
            m7_1_path,
        ) = self._load_upstream_qc(
            config=config,
            key="validation_readiness",
        )

        (
            m7_2b_qc,
            m7_2b_path,
        ) = self._load_upstream_qc(
            config=config,
            key="endpoint_resolution",
        )

        clinical_route_continues = bool(
            m7_2b_qc
            .get(
                "summary",
                {},
            )
            .get(
                "clinical_validation_route_continues",
                False,
            )
        )

        if not clinical_route_continues:

            raise RuntimeError(
                "M7.2B does not permit continued clinical validation."
            )

        logger.info(
            "Loaded locked M7.1 and M7.2B QC artifacts."
        )

        result = resolve_trf_sequences(
            m7_1_qc=m7_1_qc,
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
            result.sequence_evidence,
            self.evidence_output,
            index=False,
        )

        write_parquet(
            result.candidate_resolution,
            self.resolution_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        summary = _json_safe_row(
            result.summary.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M7.3A",

            "stage":
                "candidate_trf_sequence_resolution",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m7_1_loaded":
                    True,

                "m7_1_path":
                    str(
                        m7_1_path
                    ),

                "m7_2b_loaded":
                    True,

                "m7_2b_path":
                    str(
                        m7_2b_path
                    ),

                "clinical_validation_route_continues":
                    clinical_route_continues,
            },

            "summary":
                summary,

            "trf_sequence_evidence":
                _records_with_json_nulls(
                    result.sequence_evidence
                ),

            "candidate_trf_sequence_resolution":
                _records_with_json_nulls(
                    result.candidate_resolution
                ),

            "outputs": {
                "trf_sequence_evidence":
                    str(
                        self.evidence_output
                    ),

                "candidate_trf_sequence_resolution":
                    str(
                        self.resolution_output
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
            "M7.3A complete."
        )

        logger.info(
            "Candidate tRFs=%d | identifiers confirmed=%d | "
            "sequences resolved=%d | quantification audit ready=%d.",
            summary[
                "candidate_trfs_assessed"
            ],
            summary[
                "candidate_identifiers_source_confirmed"
            ],
            summary[
                "candidate_sequences_resolved"
            ],
            summary[
                "candidates_ready_for_quantification_audit"
            ],
        )

        for row in report[
            "candidate_trf_sequence_resolution"
        ]:

            logger.info(
                "%s | %s | sequence=%s | resolved=%s | "
                "quantification_ready=%s | status=%s.",
                row[
                    "lead_rsid"
                ],
                row[
                    "trf_id"
                ],
                row[
                    "sequence_dna"
                ],
                row[
                    "sequence_resolved"
                ],
                row[
                    "candidate_quantification_ready"
                ],
                row[
                    "sequence_resolution_status"
                ],
            )

        logger.info(
            "Overall resolution=%s.",
            summary[
                "overall_resolution_status"
            ],
        )

        logger.info(
            "Next stage=%s.",
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report