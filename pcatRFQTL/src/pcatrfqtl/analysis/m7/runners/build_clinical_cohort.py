"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/build_clinical_cohort.py

Description:
    Runner for M7.2 TCGA-PRAD Clinical Cohort Construction.

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
import requests
import yaml

from pcatrfqtl.analysis.m7.clinical_cohort import (
    construct_clinical_cohort,
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
    """Convert pandas missing values to JSON null."""

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


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert pandas Series into JSON-safe scalars."""

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


class M72ClinicalCohortRunner:
    """Execute M7.2 TCGA-PRAD clinical cohort construction."""

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
    def cohort_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "tcga_prad_clinical_cohort.parquet"
        )

    @property
    def completeness_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "clinical_variable_completeness.parquet"
        )

    @property
    def endpoint_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "clinical_endpoint_readiness.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "clinical_cohort_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m7_2_clinical_cohort.json"
        )

    # ======================================================================
    # Config
    # ======================================================================

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
                "M7.2 configuration must contain a YAML mapping."
            )

        return config

    # ======================================================================
    # M7.1 validation
    # ======================================================================

    def _load_m7_1_qc(
        self,
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        path = (
            self.project_root
            / "data"
            / "processed"
            / "m7"
            / "qc"
            / "m7_1_validation_readiness.json"
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M7.1 QC not found: {path}"
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
        ) != "M7.1":

            raise RuntimeError(
                "Unexpected upstream milestone for M7.2."
            )

        if not bool(
            payload
            .get(
                "summary",
                {},
            )
            .get(
                "clinical_metadata_route_ready",
                False,
            )
        ):

            raise RuntimeError(
                "M7.1 did not authorize the clinical metadata route."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Raw acquisition
    # ======================================================================

    def _raw_output_path(
        self,
        *,
        config: dict[str, Any],
    ) -> Path:

        raw = config[
            "raw_data"
        ]

        return (
            self.project_root
            /
            str(
                raw[
                    "directory"
                ]
            )
            /
            str(
                raw[
                    "filename"
                ]
            )
        )

    def _query_gdc(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, Any]:
        """Retrieve TCGA-PRAD cases from GDC Cases endpoint."""

        api = config[
            "gdc_api"
        ]

        filter_config = api[
            "filter"
        ]

        filters = {
            "op":
                str(
                    filter_config[
                        "op"
                    ]
                ),

            "content": {
                "field":
                    str(
                        filter_config[
                            "field"
                        ]
                    ),

                "value":
                    [
                        str(
                            filter_config[
                                "value"
                            ]
                        )
                    ],
            },
        }

        params = {
            "filters":
                json.dumps(
                    filters
                ),

            "fields":
                ",".join(
                    api[
                        "fields"
                    ]
                ),

            "expand":
                ",".join(
                    api[
                        "expand"
                    ]
                ),

            "format":
                "JSON",

            "size":
                int(
                    api[
                        "page_size"
                    ]
                ),
        }

        response = requests.get(
            str(
                api[
                    "cases_endpoint"
                ]
            ),
            params=params,
            timeout=int(
                api[
                    "timeout_seconds"
                ]
            ),
            verify=bool(
                api[
                    "verify_ssl"
                ]
            ),
        )

        response.raise_for_status()

        payload = response.json()

        if not isinstance(
            payload,
            dict,
        ):

            raise RuntimeError(
                "Unexpected GDC response type."
            )

        return payload

    def _load_or_acquire_raw(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
        bool,
    ]:
        """Use immutable raw GDC response if already present."""

        path = self._raw_output_path(
            config=config
        )

        if path.exists():

            with path.open(
                "r",
                encoding="utf-8",
            ) as handle:

                payload = json.load(
                    handle
                )

            logger.info(
                "Using existing immutable M7.2 raw GDC response."
            )

            return (
                payload,
                path,
                False,
            )

        payload = self._query_gdc(
            config=config
        )

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with path.open(
            "x",
            encoding="utf-8",
        ) as handle:

            json.dump(
                payload,
                handle,
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            )

        logger.info(
            "Acquired and preserved immutable raw GDC clinical response."
        )

        return (
            payload,
            path,
            True,
        )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.2."""

        config = self._load_config()

        (
            _m7_1_qc,
            m7_1_path,
        ) = self._load_m7_1_qc()

        (
            raw_payload,
            raw_path,
            downloaded_now,
        ) = self._load_or_acquire_raw(
            config=config
        )

        hits = (
            raw_payload
            .get(
                "data",
                {},
            )
            .get(
                "hits",
                [],
            )
        )

        if not isinstance(
            hits,
            list,
        ):

            raise RuntimeError(
                "GDC cases response does not contain data.hits list."
            )

        logger.info(
            "Loaded %d TCGA-PRAD cases from GDC clinical metadata.",
            len(
                hits
            ),
        )

        result = construct_clinical_cohort(
            cases=hits,
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
            result.cohort,
            self.cohort_output,
            index=False,
        )

        write_parquet(
            result.completeness,
            self.completeness_output,
            index=False,
        )

        write_parquet(
            result.endpoint_readiness,
            self.endpoint_output,
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
                "M7.2",

            "stage":
                "tcga_prad_clinical_cohort_construction",

            "target_cohort":
                config[
                    "target_cohort"
                ],

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
            },

            "acquisition": {
                "raw_path":
                    str(
                        raw_path
                    ),

                "downloaded_in_this_run":
                    downloaded_now,

                "raw_case_count":
                    int(
                        len(
                            hits
                        )
                    ),

                "raw_immutable":
                    True,
            },

            "summary":
                summary,

            "endpoint_readiness":
                _records_with_json_nulls(
                    result.endpoint_readiness
                ),

            "clinical_variable_completeness":
                _records_with_json_nulls(
                    result.completeness
                ),

            "outputs": {
                "cohort":
                    str(
                        self.cohort_output
                    ),

                "completeness":
                    str(
                        self.completeness_output
                    ),

                "endpoint_readiness":
                    str(
                        self.endpoint_output
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
            "M7.2 complete."
        )

        logger.info(
            "Cases=%d | clinical usable=%d | OS usable=%d/%d events | "
            "recurrence usable=%d/%d events.",
            summary[
                "unique_cases"
            ],
            summary[
                "clinically_usable_cases"
            ],
            summary[
                "os_usable_cases"
            ],
            summary[
                "os_events"
            ],
            summary[
                "recurrence_usable_cases"
            ],
            summary[
                "recurrence_events"
            ],
        )

        for row in report[
            "endpoint_readiness"
        ]:

            logger.info(
                "%s | usable=%d | events=%d | ready=%s | status=%s.",
                row[
                    "endpoint_id"
                ],
                row[
                    "usable_cases"
                ],
                row[
                    "events"
                ],
                row[
                    "endpoint_ready"
                ],
                row[
                    "endpoint_status"
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