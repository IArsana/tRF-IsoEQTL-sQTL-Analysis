"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/assess_external_small_rna_technical_eligibility.py

Description:
    Runner for M7.4B External Small-RNA Technical Eligibility.

    This runner:
        - validates the locked M7.4A discovery artifact;
        - loads the M7.4B technical-audit configuration;
        - evaluates external small-RNA datasets for suitability for
          exact-sequence auditing of the prioritized candidate tRF;
        - writes dataset-level and summary Parquet artifacts;
        - writes a QC JSON artifact.

    M7.4B is a metadata-level technical audit only.

    This runner does NOT:
        - download FASTQ or SRA reads;
        - inspect raw reads;
        - search for the candidate sequence;
        - count exact sequence matches;
        - quantify tRF abundance;
        - claim candidate presence or absence;
        - perform differential expression;
        - perform clinical association;
        - infer parent tRNA identity;
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

from pcatrfqtl.analysis.m7.external_small_rna_technical_eligibility import (
    assess_external_small_rna_technical_eligibility,
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
    """Convert a DataFrame into JSON-safe records."""

    if dataframe.empty:
        return []

    cleaned = dataframe.astype(object).where(
        pd.notna(dataframe),
        None,
    )

    records = cleaned.to_dict(
        orient="records"
    )

    for record in records:
        for key, value in list(record.items()):
            if hasattr(value, "item"):
                record[key] = value.item()

    return records


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert one pandas Series into JSON-safe scalar values."""

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

        if hasattr(value, "item"):
            value = value.item()

        output[str(key)] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M74BExternalSmallRnaTechnicalEligibilityRunner:
    """Execute M7.4B External Small-RNA Technical Eligibility."""

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
    def eligibility_output(
        self,
    ) -> Path:
        """Dataset-level technical eligibility artifact."""

        return (
            self.output_directory
            /
            "external_small_rna_technical_eligibility.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:
        """One-row M7.4B summary."""

        return (
            self.output_directory
            /
            "external_small_rna_technical_eligibility_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:
        """M7.4B QC JSON artifact."""

        return (
            self.qc_directory
            /
            "m7_4b_external_small_rna_technical_eligibility.json"
        )

    # ======================================================================
    # YAML loading
    # ======================================================================

    def _load_yaml(
        self,
        path: Path,
    ) -> dict[str, Any]:
        """Load YAML and require a mapping root."""

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
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M7.4B configuration."""

        config = self._load_yaml(
            self.config_path
        )

        required_sections = {
            "milestone",
            "stage",
            "upstream_qc",
            "dataset_registry",
            "candidate",
            "technical_criteria",
            "dataset_audit",
            "statuses",
            "scientific_policy",
            "next_stage",
        }

        missing = (
            required_sections
            -
            set(
                config.keys()
            )
        )

        if missing:

            raise RuntimeError(
                "M7.4B configuration missing required sections: "
                f"{sorted(missing)}"
            )

        if config.get(
            "milestone"
        ) != "M7.4B":

            raise RuntimeError(
                "Unexpected M7.4B milestone: "
                f"{config.get('milestone')!r}"
            )

        if (
            config.get(
                "stage"
            )
            !=
            "external_small_rna_technical_eligibility"
        ):

            raise RuntimeError(
                "Unexpected M7.4B stage: "
                f"{config.get('stage')!r}"
            )

        dataset_keys = (
            config[
                "dataset_registry"
            ].get(
                "dataset_keys"
            )
        )

        if not isinstance(
            dataset_keys,
            list,
        ) or not dataset_keys:

            raise RuntimeError(
                "M7.4B dataset_registry.dataset_keys must contain "
                "at least one dataset key."
            )

        dataset_audit = config[
            "dataset_audit"
        ]

        if not isinstance(
            dataset_audit,
            dict,
        ):

            raise RuntimeError(
                "M7.4B dataset_audit must be a mapping."
            )

        missing_audit_keys = [
            key
            for key in dataset_keys
            if key not in dataset_audit
        ]

        if missing_audit_keys:

            raise RuntimeError(
                "M7.4B dataset audit missing entries for: "
                f"{missing_audit_keys}"
            )

        extra_audit_keys = [
            key
            for key in dataset_audit
            if key not in dataset_keys
        ]

        if extra_audit_keys:

            raise RuntimeError(
                "M7.4B dataset audit contains unexpected entries: "
                f"{extra_audit_keys}"
            )

        return config

    # ======================================================================
    # Locked upstream M7.4A
    # ======================================================================

    def _load_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Load and validate the locked M7.4A QC artifact."""

        upstream = (
            config[
                "upstream_qc"
            ][
                "open_small_rna_discovery"
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
                f"M7.4A QC artifact not found: {path}"
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

            raise RuntimeError(
                "M7.4A QC artifact must contain a JSON object."
            )

        observed_milestone = payload.get(
            "milestone"
        )

        expected_milestone = upstream[
            "expected_milestone"
        ]

        if observed_milestone != expected_milestone:

            raise RuntimeError(
                "Unexpected M7.4A milestone: "
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
                "Unexpected M7.4A stage: "
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
                "overall_discovery_status"
            )
        )

        expected_status = upstream[
            "expected_overall_status"
        ]

        if observed_status != expected_status:

            raise RuntimeError(
                "Unexpected M7.4A discovery status: "
                f"{observed_status!r}; "
                f"expected {expected_status!r}."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Dataset registry consistency
    # ======================================================================

    def _validate_dataset_registry(
        self,
        *,
        config: dict[str, Any],
        m7_4a_qc: dict[str, Any],
    ) -> Path:
        """
        Verify that M7.4B uses the same canonical dataset registry entries
        that were discovered in locked M7.4A.
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

        datasets = payload.get(
            "datasets"
        )

        if not isinstance(
            datasets,
            dict,
        ):

            raise RuntimeError(
                "configs/datasets.yaml must contain top-level 'datasets'."
            )

        requested_keys = list(
            config[
                "dataset_registry"
            ][
                "dataset_keys"
            ]
        )

        missing_registry = [
            key
            for key in requested_keys
            if key not in datasets
        ]

        if missing_registry:

            raise RuntimeError(
                "M7.4B dataset keys missing from configs/datasets.yaml: "
                f"{missing_registry}"
            )

        upstream_rows = m7_4a_qc.get(
            "dataset_discovery",
            [],
        )

        if not isinstance(
            upstream_rows,
            list,
        ):

            raise RuntimeError(
                "M7.4A dataset_discovery must be a list."
            )

        upstream_keys = {
            str(
                row.get(
                    "dataset_key"
                )
            )
            for row in upstream_rows
            if isinstance(
                row,
                dict,
            )
            and row.get(
                "dataset_key"
            )
        }

        missing_upstream = [
            key
            for key in requested_keys
            if key not in upstream_keys
        ]

        if missing_upstream:

            raise RuntimeError(
                "M7.4B requested datasets not present in locked M7.4A: "
                f"{missing_upstream}"
            )

        upstream_by_key = {
            str(
                row[
                    "dataset_key"
                ]
            ):
            row
            for row in upstream_rows
            if isinstance(
                row,
                dict,
            )
            and row.get(
                "dataset_key"
            )
        }

        not_discovery_eligible = [
            key
            for key in requested_keys
            if not bool(
                upstream_by_key[
                    key
                ].get(
                    "discovery_eligible",
                    False,
                )
            )
        ]

        if not_discovery_eligible:

            raise RuntimeError(
                "M7.4B cannot audit datasets excluded by M7.4A: "
                f"{not_discovery_eligible}"
            )

        return registry_path

    # ======================================================================
    # Scientific policy
    # ======================================================================

    def _validate_policy(
        self,
        *,
        config: dict[str, Any],
    ) -> None:
        """Ensure M7.4B remains a metadata-only technical audit."""

        policy = config[
            "scientific_policy"
        ]

        if not bool(
            policy.get(
                "metadata_audit_only",
                False,
            )
        ):

            raise RuntimeError(
                "M7.4B scientific_policy.metadata_audit_only must be True."
            )

        required_false = [
            "fastq_download_performed",
            "raw_read_inspection_performed",
            "candidate_sequence_search_performed",
            "exact_match_counting_performed",
            "approximate_match_counting_performed",
            "sequence_alignment_performed",
            "trf_quantification_performed",
            "candidate_detected_claimed",
            "candidate_absent_claimed",
            "processed_mirna_used_as_trf_substitute",
            "parent_trna_inference_performed",
            "differential_expression_performed",
            "clinical_association_performed",
            "causal_inference_performed",
        ]

        violations = [
            field
            for field in required_false
            if bool(
                policy.get(
                    field,
                    False,
                )
            )
        ]

        if violations:

            raise RuntimeError(
                "M7.4B metadata-audit safeguard violation: "
                f"{violations}"
            )

    # ======================================================================
    # Result validation
    # ======================================================================

    def _validate_result(
        self,
        *,
        eligibility: pd.DataFrame,
        summary: pd.DataFrame,
        config: dict[str, Any],
    ) -> None:
        """Validate core outputs before persistence."""

        if eligibility.empty:

            raise RuntimeError(
                "M7.4B core returned empty eligibility output."
            )

        if summary.empty:

            raise RuntimeError(
                "M7.4B core returned empty summary output."
            )

        if len(
            summary
        ) != 1:

            raise RuntimeError(
                "M7.4B summary must contain exactly one row."
            )

        expected_keys = set(
            config[
                "dataset_registry"
            ][
                "dataset_keys"
            ]
        )

        observed_keys = set(
            eligibility[
                "dataset_key"
            ]
            .astype(
                str
            )
            .tolist()
        )

        if observed_keys != expected_keys:

            raise RuntimeError(
                "M7.4B dataset-key mismatch. "
                f"Expected={sorted(expected_keys)}, "
                f"observed={sorted(observed_keys)}."
            )

        if eligibility[
            "dataset_key"
        ].duplicated().any():

            raise RuntimeError(
                "M7.4B eligibility output contains duplicate dataset keys."
            )

        forbidden_true_columns = [
            "fastq_download_performed",
            "candidate_sequence_search_performed",
            "candidate_detected",
            "candidate_absent",
            "exact_match_counting_performed",
            "trf_quantification_performed",
            "clinical_association_performed",
        ]

        for column in forbidden_true_columns:

            if (
                eligibility[
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
                    "M7.4B technical-audit safeguard violated: "
                    f"{column}=True."
                )

        summary_row = summary.iloc[
            0
        ]

        summary_false_fields = [
            "fastq_download_performed",
            "candidate_sequence_search_performed",
            "candidate_detection_claimed",
            "candidate_absence_claimed",
            "trf_quantification_performed",
            "clinical_association_performed",
        ]

        for field in summary_false_fields:

            if bool(
                summary_row[
                    field
                ]
            ):

                raise RuntimeError(
                    "M7.4B summary safeguard violated: "
                    f"{field}=True."
                )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.4B."""

        config = self._load_config()

        self._validate_policy(
            config=config,
        )

        (
            m7_4a_qc,
            m7_4a_path,
        ) = self._load_upstream(
            config=config,
        )

        logger.info(
            "Loaded locked M7.4A open small-RNA discovery artifact."
        )

        registry_path = self._validate_dataset_registry(
            config=config,
            m7_4a_qc=m7_4a_qc,
        )

        logger.info(
            "Validated canonical dataset registry: %s.",
            registry_path,
        )

        result = assess_external_small_rna_technical_eligibility(
            config=config,
        )

        self._validate_result(
            eligibility=result.dataset_eligibility,
            summary=result.summary,
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
            result.dataset_eligibility,
            self.eligibility_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        summary = _json_safe_row(
            result.summary.iloc[0]
        )

        report = {
            "milestone":
                "M7.4B",

            "stage":
                "external_small_rna_technical_eligibility",

            "candidate":
                config[
                    "candidate"
                ],

            "technical_criteria":
                config[
                    "technical_criteria"
                ],

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "upstream_validation": {
                "m7_4a_loaded":
                    True,

                "m7_4a_path":
                    str(
                        m7_4a_path
                    ),

                "m7_4a_expected_status":
                    (
                        config[
                            "upstream_qc"
                        ][
                            "open_small_rna_discovery"
                        ][
                            "expected_overall_status"
                        ]
                    ),

                "m7_4a_observed_status":
                    (
                        m7_4a_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "overall_discovery_status"
                        )
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
                        registry_path
                    ),

                "dataset_keys":
                    config[
                        "dataset_registry"
                    ][
                        "dataset_keys"
                    ],
            },

            "summary":
                summary,

            "dataset_technical_eligibility":
                _records_with_json_nulls(
                    result.dataset_eligibility
                ),

            "outputs": {
                "dataset_eligibility":
                    str(
                        self.eligibility_output
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
            "M7.4B complete."
        )

        logger.info(
            "Datasets=%d | technically eligible=%d | "
            "fully eligible=%d | eligible with caveats=%d | "
            "metadata incomplete=%d | not eligible=%d.",
            summary[
                "datasets_assessed"
            ],
            summary[
                "datasets_technically_eligible"
            ],
            summary[
                "datasets_fully_eligible"
            ],
            summary[
                "datasets_eligible_with_caveats"
            ],
            summary[
                "datasets_metadata_incomplete"
            ],
            summary[
                "datasets_not_eligible"
            ],
        )

        for row in report[
            "dataset_technical_eligibility"
        ]:

            logger.info(
                "%s | raw=%s | library=%s | length=%s | "
                "orientation=%s | adapter=%s | ffpe=%s | "
                "eligible=%s | class=%s | status=%s.",
                row[
                    "dataset_key"
                ],
                row[
                    "raw_sequence_public"
                ],
                row[
                    "small_rna_library_supported"
                ],
                row[
                    "candidate_length_compatible"
                ],
                row[
                    "read_orientation_resolvable"
                ],
                row[
                    "adapter_handling_resolvable"
                ],
                row[
                    "ffpe"
                ],
                row[
                    "technical_eligibility"
                ],
                row[
                    "resolution_class"
                ],
                row[
                    "technical_status"
                ],
            )

        logger.info(
            "Overall technical status=%s.",
            summary[
                "overall_technical_status"
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