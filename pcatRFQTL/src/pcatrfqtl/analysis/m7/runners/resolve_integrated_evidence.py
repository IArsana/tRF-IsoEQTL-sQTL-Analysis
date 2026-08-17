"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/resolve_integrated_evidence.py

Description:
    Runner for M7.6 Integrated Evidence Resolution.

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

from pcatrfqtl.analysis.m7.integrated_evidence_resolution import (
    assess_integrated_evidence_resolution,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


# ============================================================================
# Safe parquet reader
# ============================================================================


def _read_parquet_safe(
    path: Path,
) -> pd.DataFrame:
    """Read M7 parquet without pandas metadata reconstruction."""

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    return table.to_pandas(
        ignore_metadata=True,
    )


# ============================================================================
# JSON helpers
# ============================================================================


def _json_safe_value(
    value: Any,
) -> Any:

    if value is None:

        return None

    if isinstance(
        value,
        (
            str,
            int,
            float,
            bool,
            list,
            dict,
        ),
    ):

        return value

    if hasattr(
        value,
        "item",
    ):

        return value.item()

    return str(
        value
    )


def _records(
    frame: pd.DataFrame,
) -> list[dict[str, Any]]:

    rows: list[dict[str, Any]] = []

    for _, row in frame.iterrows():

        record: dict[str, Any] = {}

        for key, value in row.to_dict().items():

            try:

                if (
                    not isinstance(
                        value,
                        (
                            list,
                            dict,
                        ),
                    )
                    and
                    pd.isna(
                        value
                    )
                ):

                    value = None

            except (
                TypeError,
                ValueError,
            ):

                pass

            record[
                str(
                    key
                )
            ] = _json_safe_value(
                value
            )

        rows.append(
            record
        )

    return rows


# ============================================================================
# Runner
# ============================================================================


class M76IntegratedEvidenceResolutionRunner:
    """Execute M7.6."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        ).resolve()

        self.config_path = Path(
            config_path
        ).resolve()

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

            raise RuntimeError(
                "M7.6 config must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.6":

            raise RuntimeError(
                "Unexpected M7.6 milestone."
            )

        if config.get(
            "stage"
        ) != (
            "integrated_evidence_resolution"
        ):

            raise RuntimeError(
                "Unexpected M7.6 stage."
            )

        return config

    # ======================================================================
    # Upstream QC
    # ======================================================================

    def _validate_single_upstream(
        self,
        *,
        key: str,
        spec: dict[str, Any],
    ) -> dict[str, Any]:

        path = (
            self.project_root
            /
            spec[
                "relative_path"
            ]
        )

        if not path.exists():

            raise RuntimeError(
                f"Missing upstream QC for {key}: {path}"
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
        ) != spec[
            "expected_milestone"
        ]:

            raise RuntimeError(
                f"{key} milestone mismatch."
            )

        if payload.get(
            "stage"
        ) != spec[
            "expected_stage"
        ]:

            raise RuntimeError(
                f"{key} stage mismatch."
            )

        summary = payload.get(
            "summary"
        )

        if not isinstance(
            summary,
            dict,
        ):

            raise RuntimeError(
                f"{key} summary missing."
            )

        status_field = str(
            spec[
                "status_field"
            ]
        )

        observed_status = summary.get(
            status_field
        )

        if observed_status != spec[
            "expected_status"
        ]:

            raise RuntimeError(
                f"{key} status mismatch: "
                f"{observed_status!r} != "
                f"{spec['expected_status']!r}"
            )

        for field, expected in spec.get(
            "required_summary_state",
            {},
        ).items():

            observed = summary.get(
                field
            )

            if observed != expected:

                raise RuntimeError(
                    f"{key} locked-state mismatch: "
                    f"{field}: {observed!r} != {expected!r}"
                )

        return {
            "path":
                str(
                    path
                ),

            "status_field":
                status_field,

            "status":
                observed_status,
        }

    def _validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, Any]:

        validated: dict[str, Any] = {}

        for key, spec in config[
            "upstream_qc"
        ].items():

            validated[
                key
            ] = self._validate_single_upstream(
                key=key,
                spec=spec,
            )

        return validated

    # ======================================================================
    # Inputs
    # ======================================================================

    def _load_inputs(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        pd.DataFrame,
    ]:

        inputs = config[
            "inputs"
        ]

        m7_5b_resolution = _read_parquet_safe(
            self.project_root
            /
            inputs[
                "m7_5b_candidate_resolution"
            ]
        )

        m7_5c_resolution = _read_parquet_safe(
            self.project_root
            /
            inputs[
                "m7_5c_candidate_resolution"
            ]
        )

        expected_rsids = {
            str(
                value
            )
            for value in config[
                "validation"
            ][
                "expected_rsids"
            ]
        }

        observed_b = set(
            m7_5b_resolution[
                "rsid"
            ].astype(str)
        )

        observed_c = set(
            m7_5c_resolution[
                "rsid"
            ].astype(str)
        )

        if observed_b != expected_rsids:

            raise RuntimeError(
                "M7.6 M7.5B candidate set mismatch."
            )

        if observed_c != expected_rsids:

            raise RuntimeError(
                "M7.6 M7.5C candidate set mismatch."
            )

        return (
            m7_5b_resolution,
            m7_5c_resolution,
        )

    # ======================================================================
    # Persist
    # ======================================================================

    def _persist(
        self,
        *,
        result: Any,
        config: dict[str, Any],
    ) -> dict[str, str]:

        output = config[
            "outputs"
        ]

        directory = (
            self.project_root
            /
            output[
                "processed_directory"
            ]
        )

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        paths = {
            "candidate_evidence_matrix":
                directory
                /
                output[
                    "candidate_evidence_matrix_filename"
                ],

            "evidence_domain_matrix":
                directory
                /
                output[
                    "evidence_domain_matrix_filename"
                ],

            "candidate_limitation_matrix":
                directory
                /
                output[
                    "candidate_limitation_matrix_filename"
                ],

            "candidate_resolution":
                directory
                /
                output[
                    "candidate_resolution_filename"
                ],

            "summary":
                directory
                /
                output[
                    "summary_filename"
                ],
        }

        write_parquet(
            result.candidate_evidence_matrix,
            paths[
                "candidate_evidence_matrix"
            ],
            index=False,
        )

        write_parquet(
            result.evidence_domain_matrix,
            paths[
                "evidence_domain_matrix"
            ],
            index=False,
        )

        write_parquet(
            result.candidate_limitation_matrix,
            paths[
                "candidate_limitation_matrix"
            ],
            index=False,
        )

        write_parquet(
            result.candidate_resolution,
            paths[
                "candidate_resolution"
            ],
            index=False,
        )

        write_parquet(
            result.summary,
            paths[
                "summary"
            ],
            index=False,
        )

        return {
            key:
                str(
                    path
                )
            for key, path
            in paths.items()
        }

    # ======================================================================
    # Main
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        logger.info(
            "Starting M7.6 Integrated Evidence Resolution."
        )

        upstream = self._validate_upstream(
            config=config,
        )

        (
            m7_5b_resolution,
            m7_5c_resolution,
        ) = self._load_inputs(
            config=config,
        )

        result = assess_integrated_evidence_resolution(
            m7_5b_resolution=(
                m7_5b_resolution
            ),
            m7_5c_resolution=(
                m7_5c_resolution
            ),
            config=config,
        )

        output_paths = self._persist(
            result=result,
            config=config,
        )

        report = {
            "milestone":
                "M7.6",

            "stage":
                "integrated_evidence_resolution",

            "upstream_validation":
                upstream,

            "summary":
                _records(
                    result.summary
                )[0],

            "candidate_evidence_matrix":
                _records(
                    result.candidate_evidence_matrix
                ),

            "evidence_domain_matrix":
                _records(
                    result.evidence_domain_matrix
                ),

            "candidate_limitation_matrix":
                _records(
                    result.candidate_limitation_matrix
                ),

            "candidate_resolution":
                _records(
                    result.candidate_resolution
                ),

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs":
                output_paths,
        }

        qc_path = (
            self.project_root
            /
            config[
                "outputs"
            ][
                "qc_relative_path"
            ]
        )

        qc_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        qc_path.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            ),
            encoding="utf-8",
        )

        logger.info(
            "M7.6 complete | overall=%s.",
            report[
                "summary"
            ][
                "overall_status"
            ],
        )

        return report