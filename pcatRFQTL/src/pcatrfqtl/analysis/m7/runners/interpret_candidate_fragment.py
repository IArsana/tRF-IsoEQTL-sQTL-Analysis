"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/interpret_candidate_fragment.py

Description:
    Runner for M7.4C.2C Candidate Fragment Interpretation.

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

from pcatrfqtl.analysis.m7.candidate_fragment_interpretation import (
    interpret_candidate_fragment,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


def _safe_read_parquet(
    path: Path,
) -> pd.DataFrame:
    """Read M7 parquet safely without pandas nested-column metadata issues."""

    if not path.exists():

        raise FileNotFoundError(
            f"Required M7.4C.2B artifact not found: {path}"
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
    """Convert one pandas row into JSON-safe scalars."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:

            output[str(key)] = None
            continue

        try:

            if pd.isna(
                value
            ):

                output[str(key)] = None
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

        output[str(key)] = value

    return output


def _records(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert dataframe into JSON-safe records."""

    if dataframe.empty:
        return []

    result: list[
        dict[str, Any]
    ] = []

    for _, row in dataframe.iterrows():

        result.append(
            _json_safe_row(
                row
            )
        )

    return result


class M74C2CCandidateFragmentInterpretationRunner:
    """Execute M7.4C.2C."""

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

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.4C.2C config not found: {self.config_path}"
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
                "M7.4C.2C YAML root must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.4C.2C":

            raise RuntimeError(
                "Unexpected M7.4C.2C milestone."
            )

        if config.get(
            "stage"
        ) != "candidate_fragment_interpretation":

            raise RuntimeError(
                "Unexpected M7.4C.2C stage."
            )

        return config

    def _validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> Path:

        upstream = (
            config[
                "upstream_qc"
            ][
                "positional_flanking_audit"
            ]
        )

        path = (
            self.project_root
            /
            str(
                upstream[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M7.4C.2B QC not found: {path}"
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
                "Unexpected upstream milestone."
            )

        if payload.get(
            "stage"
        ) != upstream[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected upstream stage."
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
                "Unexpected M7.4C.2B overall status."
            )

        if summary.get(
            "interpretation_status"
        ) != upstream[
            "expected_interpretation_status"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.2B interpretation status."
            )

        return path

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        upstream_path = self._validate_upstream(
            config=config,
        )

        run_summary_path = (
            self.project_root
            /
            config[
                "inputs"
            ][
                "run_summary"
            ]
        )

        aggregate_path = (
            self.project_root
            /
            config[
                "inputs"
            ][
                "motif_context_aggregate"
            ]
        )

        run_summary = _safe_read_parquet(
            run_summary_path
        )

        motif_context_aggregate = _safe_read_parquet(
            aggregate_path
        )

        logger.info(
            "Loaded M7.4C.2B artifacts | runs=%d | contexts=%d.",
            len(
                run_summary
            ),
            len(
                motif_context_aggregate
            ),
        )

        result = interpret_candidate_fragment(
            run_summary=run_summary,
            motif_context_aggregate=motif_context_aggregate,
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

        flank_path = (
            output_directory
            /
            output_config[
                "flank_architecture_filename"
            ]
        )

        interpretation_path = (
            output_directory
            /
            output_config[
                "interpretation_filename"
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
            result.flank_architecture,
            flank_path,
            index=False,
        )

        write_parquet(
            result.interpretation,
            interpretation_path,
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
                "M7.4C.2C",

            "stage":
                "candidate_fragment_interpretation",

            "upstream_qc":
                str(
                    upstream_path
                ),

            "candidate":
                config[
                    "candidate"
                ],

            "dataset":
                config[
                    "dataset"
                ],

            "summary":
                summary,

            "interpretation":
                _records(
                    result.interpretation
                ),

            "top_flank_architecture":
                _records(
                    result
                    .flank_architecture
                    .head(
                        int(
                            config[
                                "flank_analysis"
                            ][
                                "dominant_context_top_n"
                            ]
                        )
                    )
                ),

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs": {
                "flank_architecture":
                    str(
                        flank_path
                    ),

                "interpretation":
                    str(
                        interpretation_path
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
            "M7.4C.2C complete."
        )

        logger.info(
            "Evidence=%s | positive reads=%d | positive runs=%d.",
            summary[
                "evidence_class"
            ],
            summary[
                "positive_reads"
            ],
            summary[
                "positive_runs"
            ],
        )

        logger.info(
            "Full-length=%d | 5prime=%d | 3prime=%d | internal=%d.",
            summary[
                "full_length_24nt_reads"
            ],
            summary[
                "five_prime_boundary_reads"
            ],
            summary[
                "three_prime_boundary_reads"
            ],
            summary[
                "internal_motif_reads"
            ],
        )

        logger.info(
            "5prime fraction=%.6f | dominant context fraction=%.6f | "
            "dominant context runs=%d.",
            summary[
                "five_prime_boundary_fraction"
            ],
            summary[
                "dominant_context_fraction"
            ],
            summary[
                "dominant_context_runs"
            ],
        )

        logger.info(
            "Direct 24nt support=%s | longer fragment context=%s.",
            summary[
                "direct_24nt_fragment_supported"
            ],
            summary[
                "longer_fragment_context_present"
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