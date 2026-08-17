"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/annotate_trf_biology.py

Description:
    Runner for M6.4A Candidate/tRF Biological Annotation.

    This runner:
        - loads locked M6.3 candidate-gene integration QC;
        - loads the retained M3 candidate/tRF source;
        - executes M6.4A candidate/tRF biological annotation;
        - writes Parquet outputs;
        - writes a standards-compliant QC JSON report.

    Important safeguards:
        - raw/source artifacts are not modified;
        - tRF identifiers are not decoded;
        - biological tRF class is not inferred;
        - parent tRNA is not inferred;
        - gene identity is not inferred from SNP coordinates;
        - nearest-gene mapping is not performed;
        - formal colocalization is not performed;
        - causal inference is not performed.

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

from pcatrfqtl.analysis.m6.trf_biological_annotation import (
    annotate_candidate_trf_biology,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(
    __name__
)


# ============================================================================
# Parquet helpers
# ============================================================================


def _read_parquet_safely(
    path: Path,
) -> pd.DataFrame:
    """
    Read Parquet without pandas metadata reconstruction.

    This follows the project-safe Parquet reading policy for artifacts that
    may contain metadata or nested-compatible structures.
    """

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


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """
    Convert DataFrame records to standards-compliant JSON.

    Pandas/NumPy missing values are converted to Python None so that the
    resulting JSON contains null rather than NaN.
    """

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


def _json_safe_summary(
    row: pd.Series,
) -> dict[str, Any]:
    """
    Convert one summary row to JSON-native scalar values.
    """

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


class M64ATrfBiologicalAnnotationRunner:
    """Execute M6.4A Candidate/tRF Biological Annotation."""

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
    # Output paths
    # ======================================================================

    @property
    def evidence_output(
        self,
    ) -> Path:
        """Return candidate/tRF evidence output path."""

        return (
            self.output_directory
            / "candidate_trf_evidence.parquet"
        )

    @property
    def annotation_output(
        self,
    ) -> Path:
        """Return candidate biological annotation output path."""

        return (
            self.output_directory
            / "candidate_trf_biological_annotation.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:
        """Return M6.4A summary output path."""

        return (
            self.output_directory
            / "trf_biological_annotation_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:
        """Return M6.4A QC output path."""

        return (
            self.qc_directory
            / "m6_4a_trf_biological_annotation.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M6.4A configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M6.4A configuration not found: "
                f"{self.config_path}"
            )

        if not self.config_path.is_file():

            raise RuntimeError(
                "M6.4A configuration path is not a file: "
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
                "M6.4A configuration must contain a YAML mapping."
            )

        required_sections = {
            "upstream_qc",
            "candidate_trf_source",
            "coordinate_policy",
            "trf_annotation_policy",
            "rna_processing_context",
            "annotation_status",
            "scientific_policy",
            "next_stage",
        }

        missing_sections = (
            required_sections
            -
            set(
                config
            )
        )

        if missing_sections:

            raise ValueError(
                "M6.4A configuration missing required sections: "
                f"{sorted(missing_sections)}"
            )

        return config

    # ======================================================================
    # Upstream M6.3
    # ======================================================================

    def _load_m6_3_qc(
        self,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Load and validate locked M6.3 QC."""

        properties = config[
            "upstream_qc"
        ][
            "candidate_gene_integration"
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
                "M6.3 QC not found: "
                f"{path}"
            )

        if not path.is_file():

            raise RuntimeError(
                "M6.3 QC path is not a file: "
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
                "M6.3 QC must contain a JSON object."
            )

        expected_milestone = str(
            properties[
                "expected_milestone"
            ]
        )

        expected_stage = str(
            properties[
                "expected_stage"
            ]
        )

        if str(
            payload.get(
                "milestone"
            )
        ) != expected_milestone:

            raise RuntimeError(
                "Unexpected upstream M6.3 milestone: "
                f"{payload.get('milestone')}; "
                f"expected {expected_milestone}."
            )

        if str(
            payload.get(
                "stage"
            )
        ) != expected_stage:

            raise RuntimeError(
                "Unexpected upstream M6.3 stage: "
                f"{payload.get('stage')}; "
                f"expected {expected_stage}."
            )

        candidate_rows = payload.get(
            "final_candidate_gene_integration"
        )

        if not isinstance(
            candidate_rows,
            list,
        ):

            raise RuntimeError(
                "M6.3 QC does not contain "
                "final_candidate_gene_integration records."
            )

        if not candidate_rows:

            raise RuntimeError(
                "M6.3 final candidate integration is empty."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Candidate/tRF source
    # ======================================================================

    def _load_candidate_trf_source(
        self,
        config: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        Path,
    ]:
        """Load M3 candidate/tRF source artifact."""

        properties = config[
            "candidate_trf_source"
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
                "M6.4A candidate/tRF source not found: "
                f"{path}"
            )

        if not path.is_file():

            raise RuntimeError(
                "M6.4A candidate/tRF source path is not a file: "
                f"{path}"
            )

        dataframe = _read_parquet_safely(
            path
        )

        if dataframe.empty:

            raise RuntimeError(
                "M6.4A candidate/tRF source is empty."
            )

        expected_rows = properties.get(
            "expected_rows"
        )

        if expected_rows is not None:

            if len(
                dataframe
            ) != int(
                expected_rows
            ):

                logger.warning(
                    "Candidate/tRF source has %d rows; config expected %d.",
                    len(
                        dataframe
                    ),
                    int(
                        expected_rows
                    ),
                )

        return (
            dataframe,
            path,
        )

    # ======================================================================
    # Report
    # ======================================================================

    def _build_report(
        self,
        *,
        config: dict[str, Any],
        m6_3_path: Path,
        candidate_trf_path: Path,
        candidate_trf_source: pd.DataFrame,
        result: Any,
    ) -> dict[str, Any]:
        """Build standards-compliant M6.4A QC report."""

        if result.summary.empty:

            raise RuntimeError(
                "M6.4A summary is unexpectedly empty."
            )

        summary = _json_safe_summary(
            result.summary.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M6.4A",

            "stage":
                "candidate_trf_biological_annotation",

            "policy":
                dict(
                    config[
                        "scientific_policy"
                    ]
                ),

            "upstream_validation": {
                "m6_3_loaded":
                    True,

                "m6_3_path":
                    str(
                        m6_3_path
                    ),

                "candidate_trf_source_loaded":
                    True,

                "candidate_trf_source_path":
                    str(
                        candidate_trf_path
                    ),

                "candidate_trf_source_rows":
                    int(
                        len(
                            candidate_trf_source
                        )
                    ),

                "candidate_trf_source_columns":
                    int(
                        len(
                            candidate_trf_source.columns
                        )
                    ),
            },

            "summary":
                summary,

            "candidate_trf_evidence":
                _records_with_json_nulls(
                    result.candidate_trf_evidence
                ),

            "candidate_trf_biological_annotation":
                _records_with_json_nulls(
                    result.candidate_biological_annotation
                ),

            "outputs": {
                "candidate_trf_evidence":
                    str(
                        self.evidence_output
                    ),

                "candidate_trf_biological_annotation":
                    str(
                        self.annotation_output
                    ),

                "summary":
                    str(
                        self.summary_output
                    ),
            },
        }

        return report

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute complete M6.4A workflow."""

        # ------------------------------------------------------------------
        # Configuration
        # ------------------------------------------------------------------

        config = self._load_config()

        logger.info(
            "Loaded M6.4A configuration."
        )

        # ------------------------------------------------------------------
        # Locked upstream M6.3
        # ------------------------------------------------------------------

        (
            m6_3_qc,
            m6_3_path,
        ) = self._load_m6_3_qc(
            config
        )

        logger.info(
            "Loaded locked M6.3 candidate-gene integration QC."
        )

        # ------------------------------------------------------------------
        # Candidate/tRF source
        # ------------------------------------------------------------------

        (
            candidate_trf_source,
            candidate_trf_path,
        ) = self._load_candidate_trf_source(
            config
        )

        logger.info(
            "Loaded candidate/tRF source: %d rows × %d columns.",
            len(
                candidate_trf_source
            ),
            len(
                candidate_trf_source.columns
            ),
        )

        # ------------------------------------------------------------------
        # Core analysis
        # ------------------------------------------------------------------

        result = annotate_candidate_trf_biology(
            m6_3_qc=m6_3_qc,
            candidate_trf_source=candidate_trf_source,
            config=config,
        )

        if result.candidate_biological_annotation.empty:

            raise RuntimeError(
                "M6.4A candidate biological annotation is unexpectedly empty."
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
        # Analytical outputs
        # ------------------------------------------------------------------

        write_parquet(
            result.candidate_trf_evidence,
            self.evidence_output,
            index=False,
        )

        write_parquet(
            result.candidate_biological_annotation,
            self.annotation_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        # ------------------------------------------------------------------
        # QC
        # ------------------------------------------------------------------

        report = self._build_report(
            config=config,
            m6_3_path=m6_3_path,
            candidate_trf_path=candidate_trf_path,
            candidate_trf_source=candidate_trf_source,
            result=result,
        )

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

        summary = report[
            "summary"
        ]

        logger.info(
            "M6.4A complete."
        )

        logger.info(
            "Candidates=%d | with tRF identity=%d | "
            "regulatory context=%d | resolved gene=%d.",
            summary[
                "candidates_assessed"
            ],
            summary[
                "candidates_with_trf_identity"
            ],
            summary[
                "candidates_with_regulatory_feature_context"
            ],
            summary[
                "candidates_with_resolved_gene"
            ],
        )

        for row in report[
            "candidate_trf_biological_annotation"
        ]:

            logger.info(
                "%s | tRFs=%d | regulatory_feature=%s | "
                "RNA_context=%s | annotation_status=%s.",
                row[
                    "lead_rsid"
                ],
                row[
                    "trf_count"
                ],
                row[
                    "regulatory_feature_id"
                ],
                row[
                    "rna_processing_context"
                ],
                row[
                    "annotation_status"
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