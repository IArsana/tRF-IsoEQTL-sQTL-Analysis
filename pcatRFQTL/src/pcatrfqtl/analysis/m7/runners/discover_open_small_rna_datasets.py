"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/discover_open_small_rna_datasets.py

Description:
    Runner for M7.4A Open/Public Small-RNA Dataset Discovery.

    This runner:
        - validates the locked M7.3C controlled-access limitation artifact;
        - loads dataset metadata exclusively from configs/datasets.yaml;
        - extracts the top-level "datasets" registry;
        - selects only M7 external-validation datasets explicitly listed
          in the M7.4A configuration;
        - performs discovery-level prioritization through the M7.4A core;
        - writes reproducible Parquet and QC artifacts.

    Important dataset-registry rule:
        configs/datasets.yaml is the single source-of-truth for dataset
        metadata. M7.4A configuration contains dataset keys only and does
        not duplicate dataset metadata.

    M7.4A is discovery only.

    This runner does NOT:
        - download FASTQ/SRA sequencing files;
        - inspect sequencing reads;
        - search for the candidate tRF sequence;
        - quantify tRF expression;
        - establish technical eligibility;
        - perform differential expression;
        - perform clinical association;
        - make causal claims.

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

from pcatrfqtl.analysis.m7.open_small_rna_discovery import (
    discover_open_small_rna_datasets,
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

        if hasattr(
            value,
            "item",
        ):
            value = value.item()

        output[str(key)] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M74AOpenSmallRnaDiscoveryRunner:
    """
    Execute M7.4A Open/Public Small-RNA Dataset Discovery.

    Dataset metadata are loaded exclusively from the top-level "datasets"
    mapping in configs/datasets.yaml.
    """

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
    def discovery_output(
        self,
    ) -> Path:
        """Dataset-level discovery artifact."""

        return (
            self.output_directory
            /
            "open_small_rna_dataset_discovery.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:
        """One-row M7.4A summary."""

        return (
            self.output_directory
            /
            "open_small_rna_dataset_discovery_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:
        """M7.4A QC JSON artifact."""

        return (
            self.qc_directory
            /
            "m7_4a_open_small_rna_dataset_discovery.json"
        )

    # ======================================================================
    # Generic YAML loader
    # ======================================================================

    def _load_yaml(
        self,
        path: Path,
    ) -> dict[str, Any]:
        """
        Load a YAML document and require a mapping root.
        """

        if not path.exists():

            raise FileNotFoundError(
                f"Required YAML file not found: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = yaml.safe_load(
                handle
            )

        if not isinstance(
            payload,
            dict,
        ):

            raise RuntimeError(
                f"YAML root must be a mapping: {path}"
            )

        return payload

    # ======================================================================
    # M7.4A configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """
        Load and validate M7.4A configuration.
        """

        config = self._load_yaml(
            self.config_path
        )

        required_sections = {
            "milestone",
            "stage",
            "upstream_qc",
            "dataset_registry",
            "discovery_scope",
            "inclusion_criteria",
            "prioritization",
            "statuses",
            "scientific_policy",
            "next_stage",
        }

        missing_sections = (
            required_sections
            -
            set(
                config.keys()
            )
        )

        if missing_sections:

            raise RuntimeError(
                "M7.4A configuration missing required sections: "
                f"{sorted(missing_sections)}"
            )

        if config.get(
            "milestone"
        ) != "M7.4A":

            raise RuntimeError(
                "Unexpected M7.4A config milestone: "
                f"{config.get('milestone')!r}"
            )

        if (
            config.get(
                "stage"
            )
            !=
            "open_small_rna_dataset_discovery"
        ):

            raise RuntimeError(
                "Unexpected M7.4A config stage: "
                f"{config.get('stage')!r}"
            )

        registry_config = config[
            "dataset_registry"
        ]

        if not isinstance(
            registry_config,
            dict,
        ):

            raise RuntimeError(
                "M7.4A dataset_registry configuration must be a mapping."
            )

        relative_path = registry_config.get(
            "relative_path"
        )

        if not isinstance(
            relative_path,
            str,
        ) or not relative_path.strip():

            raise RuntimeError(
                "M7.4A dataset_registry.relative_path must be defined."
            )

        dataset_keys = registry_config.get(
            "dataset_keys"
        )

        if not isinstance(
            dataset_keys,
            list,
        ) or not dataset_keys:

            raise RuntimeError(
                "M7.4A dataset_registry.dataset_keys must contain "
                "at least one dataset key."
            )

        if any(
            not isinstance(
                key,
                str,
            )
            or not key.strip()
            for key in dataset_keys
        ):

            raise RuntimeError(
                "Every M7.4A dataset registry key must be a non-empty string."
            )

        if len(
            dataset_keys
        ) != len(
            set(
                dataset_keys
            )
        ):

            raise RuntimeError(
                "M7.4A dataset_registry.dataset_keys contains duplicates."
            )

        return config

    # ======================================================================
    # Dataset registry
    # ======================================================================

    def _load_dataset_registry(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """
        Load configs/datasets.yaml and extract its top-level "datasets" map.

        Expected structure:

            datasets:
              gwas_catalog:
                ...
              cancer_trfqtl:
                ...
              moradi:
                ...
              m7_gse80400:
                ...
              ...

            integration:
              ...

        Only the "datasets" mapping is supplied to the M7.4A core.
        """

        registry_path = (
            self.project_root
            /
            str(
                config[
                    "dataset_registry"
                ][
                    "relative_path"
                ]
            )
        )

        payload = self._load_yaml(
            registry_path
        )

        dataset_registry = payload.get(
            "datasets"
        )

        if not isinstance(
            dataset_registry,
            dict,
        ):

            raise RuntimeError(
                "configs/datasets.yaml must contain a top-level "
                "'datasets' mapping."
            )

        if not dataset_registry:

            raise RuntimeError(
                "configs/datasets.yaml contains an empty "
                "top-level 'datasets' mapping."
            )

        requested_keys = list(
            config[
                "dataset_registry"
            ][
                "dataset_keys"
            ]
        )

        missing_keys = [
            key
            for key in requested_keys
            if key not in dataset_registry
        ]

        if missing_keys:

            raise RuntimeError(
                "M7.4A requested dataset keys are missing from "
                "configs/datasets.yaml: "
                f"{missing_keys}"
            )

        for dataset_key in requested_keys:

            dataset = dataset_registry[
                dataset_key
            ]

            if not isinstance(
                dataset,
                dict,
            ):

                raise RuntimeError(
                    f"Dataset registry entry {dataset_key!r} "
                    "must be a mapping."
                )

        return (
            dataset_registry,
            registry_path,
        )

    # ======================================================================
    # Locked M7.3C upstream
    # ======================================================================

    def _load_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """
        Load and validate locked M7.3C controlled-access limitation QC.
        """

        upstream = (
            config[
                "upstream_qc"
            ][
                "controlled_access_limitation"
            ]
        )

        upstream_path = (
            self.project_root
            /
            str(
                upstream[
                    "relative_path"
                ]
            )
        )

        if not upstream_path.exists():

            raise FileNotFoundError(
                f"M7.3C QC artifact not found: {upstream_path}"
            )

        with upstream_path.open(
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

            raise RuntimeError(
                "M7.3C QC artifact must contain a JSON object."
            )

        observed_milestone = payload.get(
            "milestone"
        )

        expected_milestone = upstream[
            "expected_milestone"
        ]

        if observed_milestone != expected_milestone:

            raise RuntimeError(
                "Unexpected upstream milestone for M7.4A: "
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
                "Unexpected upstream stage for M7.4A: "
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
                "overall_resolution_status"
            )
        )

        expected_status = upstream.get(
            "expected_overall_status"
        )

        if (
            expected_status is not None
            and
            observed_status != expected_status
        ):

            raise RuntimeError(
                "Unexpected M7.3C resolution status for M7.4A: "
                f"{observed_status!r}; "
                f"expected {expected_status!r}."
            )

        return (
            payload,
            upstream_path,
        )

    # ======================================================================
    # Scientific-policy validation
    # ======================================================================

    def _validate_discovery_policy(
        self,
        *,
        config: dict[str, Any],
    ) -> None:
        """
        Ensure M7.4A remains a discovery-only stage.
        """

        policy = config[
            "scientific_policy"
        ]

        required_false_flags = [
            "raw_fastq_download_performed",
            "candidate_sequence_search_performed",
            "trf_quantification_performed",
            "differential_expression_performed",
            "clinical_association_performed",
            "technical_eligibility_claimed",
            "candidate_detected_claimed",
            "absent_candidate_claimed",
            "processed_mirna_used_as_trf_substitute",
            "sequence_inference_performed",
            "causal_inference_performed",
        ]

        if not bool(
            policy.get(
                "dataset_discovery_only",
                False,
            )
        ):

            raise RuntimeError(
                "M7.4A scientific policy must declare "
                "dataset_discovery_only=True."
            )

        invalid_flags = [
            flag
            for flag in required_false_flags
            if bool(
                policy.get(
                    flag,
                    False,
                )
            )
        ]

        if invalid_flags:

            raise RuntimeError(
                "M7.4A discovery-only policy violated. "
                "The following operations/claims must remain False: "
                f"{invalid_flags}"
            )

    # ======================================================================
    # Result validation
    # ======================================================================

    def _validate_result(
        self,
        *,
        discovery: pd.DataFrame,
        summary: pd.DataFrame,
        config: dict[str, Any],
    ) -> None:
        """
        Perform defensive validation before persisting M7.4A artifacts.
        """

        if discovery.empty:

            raise RuntimeError(
                "M7.4A core returned an empty dataset-discovery table."
            )

        if summary.empty:

            raise RuntimeError(
                "M7.4A core returned an empty summary table."
            )

        if len(
            summary
        ) != 1:

            raise RuntimeError(
                "M7.4A summary must contain exactly one row."
            )

        requested_keys = set(
            config[
                "dataset_registry"
            ][
                "dataset_keys"
            ]
        )

        if "dataset_key" not in discovery.columns:

            raise RuntimeError(
                "M7.4A discovery output is missing dataset_key."
            )

        observed_keys = set(
            discovery[
                "dataset_key"
            ]
            .dropna()
            .astype(
                str
            )
            .tolist()
        )

        if observed_keys != requested_keys:

            raise RuntimeError(
                "M7.4A dataset-key mismatch. "
                f"Requested={sorted(requested_keys)}, "
                f"observed={sorted(observed_keys)}."
            )

        if discovery[
            "dataset_key"
        ].duplicated().any():

            duplicated = (
                discovery.loc[
                    discovery[
                        "dataset_key"
                    ].duplicated(
                        keep=False
                    ),
                    "dataset_key",
                ]
                .astype(
                    str
                )
                .unique()
                .tolist()
            )

            raise RuntimeError(
                "M7.4A discovery output contains duplicate dataset keys: "
                f"{duplicated}"
            )

        # ------------------------------------------------------------------
        # Safeguards
        # ------------------------------------------------------------------

        safeguard_columns = [
            "technical_eligibility_established",
            "fastq_download_performed",
            "candidate_sequence_search_performed",
            "trf_quantification_performed",
            "clinical_association_performed",
        ]

        for column in safeguard_columns:

            if column not in discovery.columns:

                raise RuntimeError(
                    "M7.4A discovery output missing safeguard column: "
                    f"{column}"
                )

            if (
                discovery[
                    column
                ]
                .fillna(
                    False
                )
                .astype(
                    bool
                )
                .any()
            ):

                raise RuntimeError(
                    "M7.4A discovery-only safeguard violated: "
                    f"{column}=True observed."
                )

        summary_row = summary.iloc[
            0
        ]

        summary_false_flags = [
            "technical_eligibility_established",
            "candidate_sequence_search_performed",
            "trf_quantification_performed",
            "clinical_association_performed",
        ]

        for field in summary_false_flags:

            if field not in summary.columns:

                raise RuntimeError(
                    f"M7.4A summary missing safeguard field: {field}"
                )

            if bool(
                summary_row[
                    field
                ]
            ):

                raise RuntimeError(
                    "M7.4A summary violates discovery-only safeguard: "
                    f"{field}=True."
                )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """
        Execute M7.4A Open/Public Small-RNA Dataset Discovery.
        """

        # ------------------------------------------------------------------
        # Config
        # ------------------------------------------------------------------

        config = self._load_config()

        self._validate_discovery_policy(
            config=config,
        )

        # ------------------------------------------------------------------
        # Locked M7.3C
        # ------------------------------------------------------------------

        (
            m7_3c_qc,
            m7_3c_path,
        ) = self._load_upstream(
            config=config,
        )

        logger.info(
            "Loaded locked M7.3C controlled-access limitation artifact."
        )

        # ------------------------------------------------------------------
        # Dataset registry
        # ------------------------------------------------------------------

        (
            dataset_registry,
            dataset_registry_path,
        ) = self._load_dataset_registry(
            config=config,
        )

        requested_dataset_keys = list(
            config[
                "dataset_registry"
            ][
                "dataset_keys"
            ]
        )

        logger.info(
            "Loaded dataset registry from %s.",
            dataset_registry_path,
        )

        logger.info(
            "M7.4A requested %d external-validation dataset entries.",
            len(
                requested_dataset_keys
            ),
        )

        # ------------------------------------------------------------------
        # Core discovery
        # ------------------------------------------------------------------

        result = discover_open_small_rna_datasets(
            dataset_registry=dataset_registry,
            config=config,
        )

        # ------------------------------------------------------------------
        # Defensive validation
        # ------------------------------------------------------------------

        self._validate_result(
            discovery=result.dataset_discovery,
            summary=result.summary,
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
            result.dataset_discovery,
            self.discovery_output,
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
                "M7.4A",

            "stage":
                "open_small_rna_dataset_discovery",

            "discovery_scope":
                config[
                    "discovery_scope"
                ],

            "inclusion_criteria":
                config[
                    "inclusion_criteria"
                ],

            "prioritization_policy":
                config[
                    "prioritization"
                ],

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m7_3c_loaded":
                    True,

                "m7_3c_path":
                    str(
                        m7_3c_path
                    ),

                "m7_3c_expected_status":
                    (
                        config[
                            "upstream_qc"
                        ][
                            "controlled_access_limitation"
                        ].get(
                            "expected_overall_status"
                        )
                    ),

                "m7_3c_observed_status":
                    (
                        m7_3c_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "overall_resolution_status"
                        )
                    ),

                "tcga_direct_quantification_deferred":
                    bool(
                        m7_3c_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "candidates_tcga_quantification_deferred",
                            0,
                        )
                        >
                        0
                    ),
            },

            "dataset_registry": {
                "source_of_truth":
                    True,

                "relative_path":
                    config[
                        "dataset_registry"
                    ][
                        "relative_path"
                    ],

                "resolved_path":
                    str(
                        dataset_registry_path
                    ),

                "registry_root":
                    "datasets",

                "requested_dataset_keys":
                    requested_dataset_keys,

                "requested_dataset_count":
                    int(
                        len(
                            requested_dataset_keys
                        )
                    ),
            },

            "summary":
                summary,

            "dataset_discovery":
                _records_with_json_nulls(
                    result.dataset_discovery
                ),

            "outputs": {
                "dataset_discovery":
                    str(
                        self.discovery_output
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
            "M7.4A complete."
        )

        logger.info(
            "Datasets=%d | eligible=%d | high=%d | "
            "intermediate=%d | exploratory=%d.",
            summary[
                "datasets_assessed"
            ],
            summary[
                "datasets_discovery_eligible"
            ],
            summary[
                "high_priority_datasets"
            ],
            summary[
                "intermediate_priority_datasets"
            ],
            summary[
                "exploratory_datasets"
            ],
        )

        logger.info(
            "Raw-sequence routes=%d | explicit tRNA contexts=%d | "
            "recurrence-metadata datasets=%d.",
            summary[
                "datasets_with_raw_sequence_route"
            ],
            summary[
                "datasets_with_trna_context"
            ],
            summary[
                "datasets_with_recurrence_metadata"
            ],
        )

        for row in report[
            "dataset_discovery"
        ]:

            logger.info(
                "%s | n=%s | raw=%s | tRNA_context=%s | "
                "comparator=%s | clinical=%s | recurrence=%s | "
                "eligible=%s | priority=%s | status=%s.",
                row[
                    "dataset_id"
                ],
                row[
                    "sample_count_series"
                ],
                row[
                    "raw_sequence_route"
                ],
                row[
                    "trna_derived_rna_context_reported"
                ],
                row[
                    "tumor_comparator_available"
                ],
                row[
                    "clinicopathologic_metadata_reported"
                ],
                row[
                    "recurrence_metadata_reported"
                ],
                row[
                    "discovery_eligible"
                ],
                row[
                    "priority_class"
                ],
                row[
                    "discovery_status"
                ],
            )

        logger.info(
            "Leading dataset for technical audit=%s.",
            summary[
                "leading_dataset_for_technical_audit"
            ],
        )

        logger.info(
            "Technical eligibility established=%s | "
            "candidate sequence searched=%s | "
            "tRF quantification=%s | clinical association=%s.",
            summary[
                "technical_eligibility_established"
            ],
            summary[
                "candidate_sequence_search_performed"
            ],
            summary[
                "trf_quantification_performed"
            ],
            summary[
                "clinical_association_performed"
            ],
        )

        logger.info(
            "Overall discovery status=%s.",
            summary[
                "overall_discovery_status"
            ],
        )

        logger.info(
            "Next stage=%s.",
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC output=%s.",
            self.qc_output,
        )

        return report