"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/resolve_clinical_endpoints.py

Description:
    Runner for M7.2B Clinical Endpoint Readiness Resolution.

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
import pyarrow.parquet as pq
import yaml

from pcatrfqtl.analysis.m7.endpoint_resolution import (
    resolve_clinical_endpoints,
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
    """Convert pandas missing values to JSON-safe null values."""

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
    """Convert pandas row into JSON-native scalars."""

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
# Safe Parquet reader
# ============================================================================


def _read_parquet_safe(
    path: Path,
) -> pd.DataFrame:
    """
    Read Parquet while ignoring pandas metadata.

    This avoids nested/list metadata reconstruction issues in M5+ artifacts.
    """

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    return table.to_pandas(
        ignore_metadata=True
    )


# ============================================================================
# Runner
# ============================================================================


class M72BEndpointResolutionRunner:
    """Execute M7.2B Clinical Endpoint Readiness Resolution."""

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
    # Paths
    # ======================================================================

    @property
    def time_to_event_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "time_to_event_resolution.parquet"
        )

    @property
    def clinicopathologic_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "clinicopathologic_endpoint_readiness.parquet"
        )

    @property
    def final_resolution_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "final_clinical_endpoint_resolution.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "clinical_endpoint_resolution_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m7_2b_endpoint_resolution.json"
        )

    # ======================================================================
    # Config
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.2B config not found: {self.config_path}"
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
                "M7.2B configuration must contain a YAML mapping."
            )

        required = {
            "upstream_qc",
            "time_to_event_policy",
            "time_to_event_endpoints",
            "clinicopathologic_policy",
            "clinicopathologic_variables",
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
                "M7.2B configuration missing required sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Upstream
    # ======================================================================

    def _load_m7_2_qc(
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
            "clinical_cohort"
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
                f"M7.2 QC not found: {path}"
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
                "Unexpected upstream milestone for M7.2B."
            )

        if payload.get(
            "stage"
        ) != properties[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected upstream stage for M7.2B."
            )

        observed_project = (
            payload
            .get(
                "summary",
                {},
            )
            .get(
                "project_id"
            )
        )

        if observed_project != properties[
            "expected_project"
        ]:

            raise RuntimeError(
                "Unexpected project in M7.2 QC: "
                f"{observed_project}"
            )

        return (
            payload,
            path,
        )

    def _load_clinical_cohort(
        self,
        *,
        m7_2_qc: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        Path,
    ]:

        cohort_path_raw = (
            m7_2_qc
            .get(
                "outputs",
                {},
            )
            .get(
                "cohort"
            )
        )

        if not cohort_path_raw:

            raise RuntimeError(
                "M7.2 QC does not contain outputs.cohort."
            )

        cohort_path = Path(
            str(
                cohort_path_raw
            )
        )

        if not cohort_path.exists():

            raise FileNotFoundError(
                f"M7.2 clinical cohort not found: {cohort_path}"
            )

        cohort = _read_parquet_safe(
            cohort_path
        )

        if cohort.empty:

            raise RuntimeError(
                "M7.2 clinical cohort is empty."
            )

        if (
            "case_id"
            not in cohort.columns
        ):

            raise RuntimeError(
                "M7.2 clinical cohort missing case_id."
            )

        if cohort[
            "case_id"
        ].duplicated().any():

            raise RuntimeError(
                "M7.2 clinical cohort contains duplicate cases."
            )

        return (
            cohort,
            cohort_path,
        )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.2B."""

        config = self._load_config()

        (
            m7_2_qc,
            m7_2_qc_path,
        ) = self._load_m7_2_qc(
            config=config,
        )

        (
            clinical_cohort,
            clinical_cohort_path,
        ) = self._load_clinical_cohort(
            m7_2_qc=m7_2_qc,
        )

        logger.info(
            "Loaded locked M7.2 clinical cohort: %d cases.",
            len(
                clinical_cohort
            ),
        )

        result = resolve_clinical_endpoints(
            m7_2_qc=m7_2_qc,
            clinical_cohort=clinical_cohort,
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
            result.time_to_event_resolution,
            self.time_to_event_output,
            index=False,
        )

        write_parquet(
            result.clinicopathologic_readiness,
            self.clinicopathologic_output,
            index=False,
        )

        write_parquet(
            result.final_endpoint_resolution,
            self.final_resolution_output,
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

        final_resolution = _json_safe_row(
            result.final_endpoint_resolution.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M7.2B",

            "stage":
                "clinical_endpoint_readiness_resolution",

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m7_2_loaded":
                    True,

                "m7_2_qc_path":
                    str(
                        m7_2_qc_path
                    ),

                "clinical_cohort_loaded":
                    True,

                "clinical_cohort_path":
                    str(
                        clinical_cohort_path
                    ),

                "clinical_case_count":
                    int(
                        len(
                            clinical_cohort
                        )
                    ),
            },

            "summary":
                summary,

            "final_endpoint_resolution":
                final_resolution,

            "time_to_event_resolution":
                _records_with_json_nulls(
                    result.time_to_event_resolution
                ),

            "clinicopathologic_endpoint_readiness":
                _records_with_json_nulls(
                    result.clinicopathologic_readiness
                ),

            "outputs": {
                "time_to_event_resolution":
                    str(
                        self.time_to_event_output
                    ),

                "clinicopathologic_endpoint_readiness":
                    str(
                        self.clinicopathologic_output
                    ),

                "final_endpoint_resolution":
                    str(
                        self.final_resolution_output
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
            "M7.2B complete."
        )

        logger.info(
            "Time-to-event ready=%d/%d | clinicopathologic ready=%d/%d.",
            summary[
                "time_to_event_endpoints_ready"
            ],
            summary[
                "time_to_event_endpoints_assessed"
            ],
            summary[
                "clinicopathologic_endpoints_ready"
            ],
            summary[
                "clinicopathologic_endpoints_assessed"
            ],
        )

        logger.info(
            "Primary time-to-event endpoint available=%s.",
            summary[
                "primary_time_to_event_endpoint_available"
            ],
        )

        logger.info(
            "Ready clinicopathologic endpoints=%s.",
            summary[
                "ready_clinicopathologic_endpoints"
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