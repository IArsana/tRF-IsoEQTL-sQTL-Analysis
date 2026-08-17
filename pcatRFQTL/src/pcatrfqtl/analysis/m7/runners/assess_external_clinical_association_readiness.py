"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/assess_external_clinical_association_readiness.py

Description:
    Runner for M7.4D External Clinical Association Readiness.

    This runner validates the locked M7.4C.2C fragment interpretation and
    evaluates whether GSE80400 supports defensible inferential clinical
    association analyses.

    No clinical association model is fitted in this stage.

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

from pcatrfqtl.analysis.m7.external_clinical_association_readiness import (
    assess_external_clinical_association_readiness,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


def _safe_read_parquet(
    path: Path,
) -> pd.DataFrame:
    """Safely read M7 parquet artifacts."""

    if not path.exists():

        raise FileNotFoundError(
            f"Required parquet artifact not found: {path}"
        )

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    return table.to_pandas(
        ignore_metadata=True,
    )


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert pandas row to JSON-safe dict."""

    result: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:

            result[
                str(
                    key
                )
            ] = None

            continue

        try:

            if pd.isna(
                value
            ):

                result[
                    str(
                        key
                    )
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

        result[
            str(
                key
            )
        ] = value

    return result


def _records(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert dataframe to JSON-safe records."""

    return [
        _json_safe_row(
            row
        )
        for _, row in dataframe.iterrows()
    ]


class M74DExternalClinicalAssociationReadinessRunner:
    """Execute M7.4D."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        )

        self.config_path = Path(
            config_path
        )

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load M7.4D configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.4D config not found: {self.config_path}"
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

            raise RuntimeError(
                "M7.4D config root must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.4D":

            raise RuntimeError(
                "Unexpected M7.4D milestone."
            )

        if config.get(
            "stage"
        ) != "external_clinical_association_readiness":

            raise RuntimeError(
                "Unexpected M7.4D stage."
            )

        return config

    def _validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Validate locked M7.4C.2C QC."""

        upstream = (
            config[
                "upstream_qc"
            ][
                "fragment_interpretation"
            ]
        )

        path = (
            self.project_root
            /
            upstream[
                "relative_path"
            ]
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M7.4C.2C QC not found: {path}"
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
        ) != upstream[
            "expected_milestone"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.2C milestone."
            )

        if payload.get(
            "stage"
        ) != upstream[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.2C stage."
            )

        summary = payload.get(
            "summary",
            {},
        )

        if summary.get(
            "overall_status"
        ) != upstream[
            "expected_overall_status"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.2C overall status."
            )

        if summary.get(
            "evidence_class"
        ) != upstream[
            "expected_evidence_class"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.2C evidence class."
            )

        return (
            payload,
            path,
        )

    def _validate_molecular_artifact(
        self,
        *,
        config: dict[str, Any],
    ) -> Path:
        """Validate existence and shape of exact-sequence abundance artifact."""

        path = (
            self.project_root
            /
            config[
                "molecular_input"
            ][
                "run_counts_relative_path"
            ]
        )

        counts = _safe_read_parquet(
            path
        )

        expected_runs = int(
            config[
                "molecular_input"
            ][
                "expected_runs"
            ]
        )

        if len(
            counts
        ) != expected_runs:

            raise RuntimeError(
                "M7.4D molecular artifact has unexpected number "
                f"of runs: {len(counts)} != {expected_runs}."
            )

        abundance_variable = config[
            "molecular_input"
        ][
            "abundance_variable"
        ]

        if abundance_variable not in counts.columns:

            raise RuntimeError(
                "M7.4D molecular abundance variable missing: "
                f"{abundance_variable}"
            )

        return path

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.4D."""

        config = self._load_config()

        (
            upstream_qc,
            upstream_path,
        ) = self._validate_upstream(
            config=config,
        )

        molecular_path = self._validate_molecular_artifact(
            config=config,
        )

        result = assess_external_clinical_association_readiness(
            config=config,
        )

        output_config = config[
            "outputs"
        ]

        output_directory = (
            self.project_root
            /
            output_config[
                "processed_directory"
            ]
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        sample_path = (
            output_directory
            /
            output_config[
                "sample_inventory_filename"
            ]
        )

        endpoint_path = (
            output_directory
            /
            output_config[
                "endpoint_readiness_filename"
            ]
        )

        summary_path = (
            output_directory
            /
            output_config[
                "summary_filename"
            ]
        )

        write_parquet(
            result.sample_inventory,
            sample_path,
            index=False,
        )

        write_parquet(
            result.endpoint_readiness,
            endpoint_path,
            index=False,
        )

        write_parquet(
            result.summary,
            summary_path,
            index=False,
        )

        summary = _json_safe_row(
            result.summary.iloc[
                0
            ]
        )

        qc_path = (
            self.project_root
            /
            output_config[
                "qc_relative_path"
            ]
        )

        qc_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        report = {
            "milestone":
                "M7.4D",

            "stage":
                "external_clinical_association_readiness",

            "dataset":
                config[
                    "dataset"
                ],

            "upstream_validation": {
                "m7_4c2c_loaded":
                    True,

                "m7_4c2c_path":
                    str(
                        upstream_path
                    ),

                "m7_4c2c_evidence_class":
                    (
                        upstream_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "evidence_class"
                        )
                    ),
            },

            "molecular_input": {
                "path":
                    str(
                        molecular_path
                    ),

                "abundance_variable":
                    config[
                        "molecular_input"
                    ][
                        "abundance_variable"
                    ],

                "interpretation":
                    config[
                        "molecular_input"
                    ][
                        "abundance_interpretation"
                    ],
            },

            "run_mapping_policy":
                config[
                    "run_mapping_policy"
                ],

            "summary":
                summary,

            "sample_inventory":
                _records(
                    result.sample_inventory
                ),

            "endpoint_readiness":
                _records(
                    result.endpoint_readiness
                ),

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs": {
                "sample_inventory":
                    str(
                        sample_path
                    ),

                "endpoint_readiness":
                    str(
                        endpoint_path
                    ),

                "summary":
                    str(
                        summary_path
                    ),
            },
        }

        with qc_path.open(
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
            "M7.4D complete."
        )

        logger.info(
            "Samples=%d | endpoints=%d | ready=%d | not-ready=%d.",
            summary[
                "source_samples"
            ],
            summary[
                "endpoints_assessed"
            ],
            summary[
                "endpoints_association_ready"
            ],
            summary[
                "endpoints_not_ready"
            ],
        )

        for row in report[
            "endpoint_readiness"
        ]:

            logger.info(
                "%s | usable=%d | groups=%s | ready=%s | status=%s.",
                row[
                    "endpoint"
                ],
                row[
                    "usable_cases"
                ],
                row[
                    "group_counts"
                ],
                row[
                    "association_ready"
                ],
                row[
                    "readiness_status"
                ],
            )

        logger.info(
            "Run-sample mapping verified=%s.",
            summary[
                "run_sample_mapping_verified"
            ],
        )

        logger.info(
            "Clinical validation claimed=%s | association fitted=%s.",
            summary[
                "clinical_validation_claimed"
            ],
            summary[
                "association_model_fitted"
            ],
        )

        logger.info(
            "Overall=%s | next=%s.",
            summary[
                "overall_status"
            ],
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC=%s.",
            qc_path,
        )

        return report