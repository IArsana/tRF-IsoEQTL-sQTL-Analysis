"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/audit_study_independence_harmonization.py

Description:
    Runner for M7.5C Study Independence & Harmonization Audit.

    M7 inputs contain nested list columns, therefore this runner deliberately
    uses the M7-safe PyArrow parquet reader rather than pandas metadata
    reconstruction.

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

from pcatrfqtl.analysis.m7.study_independence_harmonization import (
    assess_study_independence_harmonization,
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
    """Read nested M7 parquet without pandas metadata reconstruction."""

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

    records: list[dict[str, Any]] = []

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

        records.append(
            record
        )

    return records


# ============================================================================
# Runner
# ============================================================================


class M75CStudyIndependenceHarmonizationRunner:
    """Execute M7.5C."""

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
                "M7.5C config root must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.5C":

            raise RuntimeError(
                "Unexpected M7.5C milestone."
            )

        if config.get(
            "stage"
        ) != (
            "study_independence_and_harmonization_audit"
        ):

            raise RuntimeError(
                "Unexpected M7.5C stage."
            )

        return config

    # ======================================================================
    # Upstream
    # ======================================================================

    def _validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, Any]:

        spec = config[
            "upstream_qc"
        ][
            "m7_5b"
        ]

        path = (
            self.project_root
            /
            spec[
                "relative_path"
            ]
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
                "M7.5B milestone mismatch."
            )

        if payload.get(
            "stage"
        ) != spec[
            "expected_stage"
        ]:

            raise RuntimeError(
                "M7.5B stage mismatch."
            )

        summary = payload.get(
            "summary"
        )

        if not isinstance(
            summary,
            dict,
        ):

            raise RuntimeError(
                "M7.5B summary missing."
            )

        status_field = str(
            spec[
                "status_field"
            ]
        )

        if summary.get(
            status_field
        ) != spec[
            "expected_status"
        ]:

            raise RuntimeError(
                "M7.5B overall status mismatch."
            )

        for field, expected in spec[
            "required_summary_state"
        ].items():

            observed = summary.get(
                field
            )

            if observed != expected:

                raise RuntimeError(
                    "M7.5B locked-state mismatch.\n"
                    f"Field: {field}\n"
                    f"Observed: {observed!r}\n"
                    f"Expected: {expected!r}"
                )

        return {
            "path":
                str(
                    path
                ),

            "status":
                summary[
                    status_field
                ],
        }

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
        pd.DataFrame,
    ]:

        inputs = config[
            "inputs"
        ]

        study_inventory = _read_parquet_safe(
            self.project_root
            /
            inputs[
                "study_inventory"
            ]
        )

        file_manifest = _read_parquet_safe(
            self.project_root
            /
            inputs[
                "file_manifest"
            ]
        )

        exact_hits = _read_parquet_safe(
            self.project_root
            /
            inputs[
                "exact_variant_hits"
            ]
        )

        if len(
            exact_hits
        ) != int(
            config[
                "validation"
            ][
                "expected_exact_hit_rows"
            ]
        ):

            raise RuntimeError(
                "M7.5C exact-hit count mismatch."
            )

        expected_rsids = {
            str(value)
            for value in config[
                "validation"
            ][
                "expected_rsids"
            ]
        }

        observed_rsids = set(
            exact_hits[
                "rs_id"
            ].astype(str)
        )

        if observed_rsids != expected_rsids:

            raise RuntimeError(
                "M7.5C exact-hit candidate set mismatch."
            )

        return (
            study_inventory,
            file_manifest,
            exact_hits,
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

        outputs = config[
            "outputs"
        ]

        directory = (
            self.project_root
            /
            outputs[
                "processed_directory"
            ]
        )

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        paths = {
            "publication_family_inventory":
                directory
                /
                outputs[
                    "publication_family_inventory_filename"
                ],

            "study_pair_overlap":
                directory
                /
                outputs[
                    "study_pair_overlap_filename"
                ],

            "candidate_study_harmonization":
                directory
                /
                outputs[
                    "candidate_study_harmonization_filename"
                ],

            "candidate_effect_harmonization":
                directory
                /
                outputs[
                    "candidate_effect_harmonization_filename"
                ],

            "publication_family_direction":
                directory
                /
                outputs[
                    "publication_family_direction_filename"
                ],

            "candidate_direction_resolution":
                directory
                /
                outputs[
                    "candidate_direction_resolution_filename"
                ],

            "candidate_resolution":
                directory
                /
                outputs[
                    "candidate_resolution_filename"
                ],

            "summary":
                directory
                /
                outputs[
                    "summary_filename"
                ],
        }

        write_parquet(
            result.publication_family_inventory,
            paths[
                "publication_family_inventory"
            ],
            index=False,
        )

        write_parquet(
            result.study_pair_overlap,
            paths[
                "study_pair_overlap"
            ],
            index=False,
        )

        write_parquet(
            result.candidate_study_harmonization,
            paths[
                "candidate_study_harmonization"
            ],
            index=False,
        )

        write_parquet(
            result.candidate_effect_harmonization,
            paths[
                "candidate_effect_harmonization"
            ],
            index=False,
        )

        write_parquet(
            result.publication_family_direction,
            paths[
                "publication_family_direction"
            ],
            index=False,
        )

        write_parquet(
            result.candidate_direction_resolution,
            paths[
                "candidate_direction_resolution"
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
            "Starting M7.5C Study Independence & Harmonization Audit."
        )

        upstream = self._validate_upstream(
            config=config,
        )

        (
            study_inventory,
            file_manifest,
            exact_hits,
        ) = self._load_inputs(
            config=config,
        )

        result = assess_study_independence_harmonization(
            study_inventory=study_inventory,
            file_manifest=file_manifest,
            exact_hits=exact_hits,
            config=config,
        )

        output_paths = self._persist(
            result=result,
            config=config,
        )

        report = {
            "milestone":
                "M7.5C",

            "stage":
                "study_independence_and_harmonization_audit",

            "upstream_validation":
                upstream,

            "summary":
                _records(
                    result.summary
                )[0],

            "publication_family_inventory":
                _records(
                    result.publication_family_inventory
                ),

            "study_pair_overlap":
                _records(
                    result.study_pair_overlap
                ),

            "candidate_effect_harmonization":
                _records(
                    result.candidate_effect_harmonization
                ),

            "publication_family_direction":
                _records(
                    result.publication_family_direction
                ),

            "candidate_direction_resolution":
                _records(
                    result.candidate_direction_resolution
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
            "M7.5C complete | overall=%s.",
            report[
                "summary"
            ][
                "overall_status"
            ],
        )

        return report