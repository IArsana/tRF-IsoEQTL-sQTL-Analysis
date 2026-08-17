"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/assess_trf_quantification_feasibility.py

Description:
    Runner for M7.3B tRF Quantification Feasibility Audit.

    This runner:
        - validates locked M7.3A and M7.2B upstream artifacts;
        - queries GDC TCGA-PRAD miRNA-Seq file metadata;
        - preserves query scope separately from returned nested metadata;
        - audits aligned BAM accessibility;
        - writes candidate/resource feasibility artifacts.

    This runner does NOT:
        - download controlled BAM files;
        - perform read counting;
        - quantify tRF expression;
        - perform clinical association.

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

from pcatrfqtl.analysis.m7.trf_quantification_feasibility import (
    assess_trf_quantification_feasibility,
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
                str(
                    key
                )
            ] = None

            continue

        try:

            if pd.isna(
                value
            ):

                output[
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

        output[
            str(
                key
            )
        ] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M73BTrfQuantificationFeasibilityRunner:
    """Execute M7.3B tRF Quantification Feasibility Audit."""

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
    def resource_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            /
            "gdc_mirna_resource_inventory.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            /
            "candidate_trf_quantification_feasibility.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            /
            "trf_quantification_feasibility_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            /
            "m7_3b_trf_quantification_feasibility.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M7.3B configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.3B config not found: {self.config_path}"
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
                "M7.3B configuration must contain a YAML mapping."
            )

        required = {
            "upstream_qc",
            "target_resource",
            "gdc_api",
            "query_scope_policy",
            "resource_classes",
            "quantification_policy",
            "access_policy",
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
                "M7.3B configuration missing required sections: "
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
        """Load and validate one upstream QC artifact."""

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
                f"M7.3B upstream QC missing: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        observed_milestone = payload.get(
            "milestone"
        )

        expected_milestone = properties[
            "expected_milestone"
        ]

        if observed_milestone != expected_milestone:

            raise RuntimeError(
                f"Unexpected milestone for upstream {key}: "
                f"{observed_milestone!r}; "
                f"expected {expected_milestone!r}."
            )

        observed_stage = payload.get(
            "stage"
        )

        expected_stage = properties[
            "expected_stage"
        ]

        if observed_stage != expected_stage:

            raise RuntimeError(
                f"Unexpected stage for upstream {key}: "
                f"{observed_stage!r}; "
                f"expected {expected_stage!r}."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # GDC query construction
    # ======================================================================

    def _build_gdc_filters(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, Any]:
        """Build deterministic GDC file query filters."""

        target = config[
            "target_resource"
        ]

        return {
            "op":
                "and",

            "content": [
                {
                    "op":
                        "in",

                    "content": {
                        "field":
                            "cases.project.project_id",

                        "value": [
                            str(
                                target[
                                    "project_id"
                                ]
                            )
                        ],
                    },
                },
                {
                    "op":
                        "in",

                    "content": {
                        "field":
                            "files.experimental_strategy",

                        "value": [
                            str(
                                target[
                                    "experimental_strategy"
                                ]
                            )
                        ],
                    },
                },
            ],
        }

    # ======================================================================
    # GDC query
    # ======================================================================

    def _query_gdc_files(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        list[dict[str, Any]],
        dict[str, Any],
    ]:
        """
        Query TCGA-PRAD miRNA-Seq file metadata.

        M7.3B performs metadata retrieval only.
        Controlled files are not downloaded.

        Returns:
            (file_hits, query_metadata)
        """

        api = config[
            "gdc_api"
        ]

        target = config[
            "target_resource"
        ]

        filters = self._build_gdc_filters(
            config=config,
        )

        page_size = int(
            api[
                "page_size"
            ]
        )

        request_payload = {
            "filters":
                filters,

            "fields":
                ",".join(
                    str(
                        field
                    )
                    for field in api[
                        "fields"
                    ]
                ),

            "expand":
                ",".join(
                    str(
                        field
                    )
                    for field in api[
                        "expand"
                    ]
                ),

            "format":
                "JSON",

            "from":
                0,

            "size":
                page_size,
        }

        response = requests.post(
            str(
                api[
                    "files_endpoint"
                ]
            ),
            headers={
                "Content-Type":
                    "application/json",
            },
            json=request_payload,
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
                "GDC files response is not a JSON object."
            )

        data = payload.get(
            "data",
            {},
        )

        if not isinstance(
            data,
            dict,
        ):

            raise RuntimeError(
                "GDC files response does not contain a valid data object."
            )

        hits = data.get(
            "hits",
            [],
        )

        if not isinstance(
            hits,
            list,
        ):

            raise RuntimeError(
                "GDC files response does not contain data.hits list."
            )

        pagination = data.get(
            "pagination",
            {},
        )

        if not isinstance(
            pagination,
            dict,
        ):
            pagination = {}

        total_raw = pagination.get(
            "total"
        )

        total = (
            int(
                total_raw
            )
            if total_raw is not None
            else len(
                hits
            )
        )

        # ------------------------------------------------------------------
        # Do not silently truncate inventory.
        # ------------------------------------------------------------------

        if total > len(
            hits
        ):

            raise RuntimeError(
                "GDC returned more matching miRNA-Seq files than were "
                "retrieved in the configured page. "
                f"Retrieved={len(hits)}, total={total}, "
                f"page_size={page_size}. "
                "Increase page_size or implement pagination before locking "
                "M7.3B."
            )

        # ------------------------------------------------------------------
        # Optional API warnings
        # ------------------------------------------------------------------

        warnings = payload.get(
            "warnings",
            []
        )

        if warnings is None:
            warnings = []

        if not isinstance(
            warnings,
            list,
        ):
            warnings = [
                str(
                    warnings
                )
            ]

        query_metadata = {
            "repository":
                str(
                    target[
                        "repository"
                    ]
                ),

            "endpoint":
                str(
                    api[
                        "files_endpoint"
                    ]
                ),

            "project_filter":
                str(
                    target[
                        "project_id"
                    ]
                ),

            "experimental_strategy_filter":
                str(
                    target[
                        "experimental_strategy"
                    ]
                ),

            "project_query_scope_applied":
                True,

            "experimental_strategy_query_scope_applied":
                True,

            "metadata_only":
                True,

            "page_size":
                page_size,

            "retrieved_file_count":
                int(
                    len(
                        hits
                    )
                ),

            "reported_total_file_count":
                total,

            "all_matching_files_retrieved":
                bool(
                    len(
                        hits
                    )
                    ==
                    total
                ),

            "api_warning_count":
                int(
                    len(
                        warnings
                    )
                ),

            "api_warnings":
                warnings,

            "controlled_files_downloaded":
                False,

            "read_counting_performed":
                False,
        }

        return (
            hits,
            query_metadata,
        )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.3B."""

        config = self._load_config()

        # ------------------------------------------------------------------
        # Locked upstream M7.3A
        # ------------------------------------------------------------------

        (
            m7_3a_qc,
            m7_3a_path,
        ) = self._load_upstream_qc(
            config=config,
            key="trf_sequence_resolution",
        )

        # ------------------------------------------------------------------
        # Locked upstream M7.2B
        # ------------------------------------------------------------------

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

        # ------------------------------------------------------------------
        # M7.3A expected resolution state
        # ------------------------------------------------------------------

        expected_status = (
            config[
                "upstream_qc"
            ][
                "trf_sequence_resolution"
            ]
            .get(
                "expected_overall_status"
            )
        )

        if expected_status is not None:

            observed_status = (
                m7_3a_qc
                .get(
                    "summary",
                    {},
                )
                .get(
                    "overall_resolution_status"
                )
            )

            if observed_status != expected_status:

                raise RuntimeError(
                    "Unexpected M7.3A sequence-resolution status: "
                    f"{observed_status!r}; "
                    f"expected {expected_status!r}."
                )

        logger.info(
            "Loaded locked M7.3A and M7.2B QC artifacts."
        )

        # ------------------------------------------------------------------
        # Query GDC metadata
        # ------------------------------------------------------------------

        (
            gdc_files,
            query_metadata,
        ) = self._query_gdc_files(
            config=config,
        )

        logger.info(
            "GDC miRNA-Seq resource inventory returned %d files.",
            len(
                gdc_files
            ),
        )

        # ------------------------------------------------------------------
        # Core assessment
        # ------------------------------------------------------------------

        result = assess_trf_quantification_feasibility(
            m7_3a_qc=m7_3a_qc,
            gdc_files=gdc_files,
            config=config,
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
        # Parquet artifacts
        # ------------------------------------------------------------------

        write_parquet(
            result.resource_inventory,
            self.resource_output,
            index=False,
        )

        write_parquet(
            result.candidate_feasibility,
            self.candidate_output,
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

        # ------------------------------------------------------------------
        # QC report
        # ------------------------------------------------------------------

        report = {
            "milestone":
                "M7.3B",

            "stage":
                "trf_quantification_feasibility",

            "target_resource":
                config[
                    "target_resource"
                ],

            "resource_query":
                query_metadata,

            "policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m7_3a_loaded":
                    True,

                "m7_3a_path":
                    str(
                        m7_3a_path
                    ),

                "m7_3a_expected_status":
                    expected_status,

                "m7_3a_observed_status":
                    (
                        m7_3a_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "overall_resolution_status"
                        )
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

            "candidate_trf_quantification_feasibility":
                _records_with_json_nulls(
                    result.candidate_feasibility
                ),

            "gdc_mirna_resource_inventory":
                _records_with_json_nulls(
                    result.resource_inventory
                ),

            "outputs": {
                "resource_inventory":
                    str(
                        self.resource_output
                    ),

                "candidate_feasibility":
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
                allow_nan=False,
            )

        # ------------------------------------------------------------------
        # Logging
        # ------------------------------------------------------------------

        logger.info(
            "M7.3B complete."
        )

        logger.info(
            "Candidates=%d | sequence-ready=%d | "
            "aligned BAMs=%d | BAM cases=%d | "
            "direct ready=%d | controlled route=%d.",
            summary[
                "candidate_trfs_assessed"
            ],
            summary[
                "candidate_sequences_resolved"
            ],
            summary[
                "aligned_mirna_bam_files"
            ],
            summary[
                "aligned_mirna_bam_cases"
            ],
            summary[
                "candidates_directly_quantifiable_now"
            ],
            summary[
                "candidates_requiring_controlled_access"
            ],
        )

        logger.info(
            "Project provenance: query scoped=%s | "
            "files with nested project metadata=%d | "
            "files explicitly showing target project=%d.",
            summary[
                "project_query_scope_applied"
            ],
            summary[
                "files_with_project_metadata"
            ],
            summary[
                "files_with_target_project_observed"
            ],
        )

        for row in report[
            "candidate_trf_quantification_feasibility"
        ]:

            logger.info(
                "%s | %s | sequence_ready=%s | "
                "aligned_bams=%d | aligned_cases=%d | "
                "controlled_required=%s | direct_ready=%s | status=%s.",
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
                    "aligned_mirna_bam_file_count"
                ],
                row[
                    "aligned_mirna_bam_case_count"
                ],
                row[
                    "controlled_access_required"
                ],
                row[
                    "direct_quantification_ready_now"
                ],
                row[
                    "candidate_quantification_feasibility_status"
                ],
            )

        logger.info(
            "Overall feasibility=%s.",
            summary[
                "overall_feasibility_status"
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