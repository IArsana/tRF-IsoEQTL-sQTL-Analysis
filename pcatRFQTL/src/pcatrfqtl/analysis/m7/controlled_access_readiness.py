"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/controlled_access_readiness.py

Description:
    Core logic for M7.3C Controlled-Access Quantification Readiness.

    This stage determines whether candidate-specific tRF quantification can
    proceed operationally after M7.3B established that the required
    TCGA-PRAD miRNA-Seq BAM resources are controlled-access.

    Security rules:
        - token contents are never persisted;
        - token contents are never logged;
        - token hashes are not persisted;
        - only token-path configuration state is assessed.

    M7.3C does NOT:
        - download BAM files;
        - perform read counting;
        - quantify tRF expression;
        - perform clinical association;
        - infer authorization merely from resource existence.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class ControlledAccessReadinessResult:
    """Container for M7.3C output tables."""

    resource_readiness: pd.DataFrame
    candidate_readiness: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional scalar text."""

    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except (
        TypeError,
        ValueError,
    ):
        pass

    normalized = str(value).strip()

    return normalized or None


# ============================================================================
# Controlled resource extraction
# ============================================================================


def extract_controlled_resources(
    *,
    m7_3b_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract controlled aligned miRNA BAM resources from locked M7.3B."""

    records = m7_3b_qc.get(
        "gdc_mirna_resource_inventory",
        [],
    )

    if not isinstance(records, list):
        raise RuntimeError(
            "M7.3B gdc_mirna_resource_inventory must be a list."
        )

    rows: list[dict[str, Any]] = []

    for record in records:

        if not isinstance(record, dict):
            continue

        if record.get("resource_class") != "ALIGNED_MIRNA_READS":
            continue

        if not bool(record.get("controlled_access", False)):
            continue

        rows.append(
            {
                "file_id":
                    _normalize_text(
                        record.get("file_id")
                    ),

                "file_name":
                    _normalize_text(
                        record.get("file_name")
                    ),

                "access":
                    _normalize_text(
                        record.get("access")
                    ),

                "state":
                    _normalize_text(
                        record.get("state")
                    ),

                "case_ids":
                    record.get(
                        "case_ids",
                        [],
                    ),

                "case_submitter_ids":
                    record.get(
                        "case_submitter_ids",
                        [],
                    ),

                "sample_types":
                    record.get(
                        "sample_types",
                        [],
                    ),

                "project_ids":
                    record.get(
                        "project_ids",
                        [],
                    ),

                "target_project_observed":
                    bool(
                        record.get(
                            "target_project_observed",
                            False,
                        )
                    ),

                "target_project_query_scoped":
                    bool(
                        record.get(
                            "target_project_query_scoped",
                            False,
                        )
                    ),

                "accepted_for_direct_trf_quantification":
                    bool(
                        record.get(
                            "accepted_for_direct_trf_quantification",
                            False,
                        )
                    ),
            }
        )

    result = pd.DataFrame(rows)

    if result.empty:
        return result

    if result["file_id"].isna().any():
        raise RuntimeError(
            "Controlled BAM inventory contains missing GDC file IDs."
        )

    if result["file_id"].duplicated().any():
        raise RuntimeError(
            "Controlled BAM inventory contains duplicate GDC file IDs."
        )

    return result.reset_index(drop=True)


# ============================================================================
# Candidate extraction
# ============================================================================


def extract_candidates(
    *,
    m7_3b_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract M7.3B candidate feasibility state."""

    records = m7_3b_qc.get(
        "candidate_trf_quantification_feasibility",
        [],
    )

    if not isinstance(records, list):
        raise RuntimeError(
            "M7.3B candidate feasibility artifact must be a list."
        )

    rows: list[dict[str, Any]] = []

    for record in records:

        if not isinstance(record, dict):
            continue

        rows.append(
            {
                "lead_rsid":
                    _normalize_text(
                        record.get("lead_rsid")
                    ),

                "priority_rank":
                    record.get("priority_rank"),

                "priority_class":
                    _normalize_text(
                        record.get("priority_class")
                    ),

                "trf_id":
                    _normalize_text(
                        record.get("trf_id")
                    ),

                "sequence_dna":
                    _normalize_text(
                        record.get("sequence_dna")
                    ),

                "sequence_ready_for_quantification":
                    bool(
                        record.get(
                            "sequence_ready_for_quantification",
                            False,
                        )
                    ),

                "aligned_mirna_bam_file_count":
                    int(
                        record.get(
                            "aligned_mirna_bam_file_count",
                            0,
                        )
                    ),

                "aligned_mirna_bam_case_count":
                    int(
                        record.get(
                            "aligned_mirna_bam_case_count",
                            0,
                        )
                    ),

                "controlled_aligned_bam_count":
                    int(
                        record.get(
                            "controlled_aligned_bam_count",
                            0,
                        )
                    ),

                "controlled_access_required":
                    bool(
                        record.get(
                            "controlled_access_required",
                            False,
                        )
                    ),

                "candidate_retained":
                    bool(
                        record.get(
                            "candidate_retained",
                            True,
                        )
                    ),
            }
        )

    result = pd.DataFrame(rows)

    if result.empty:
        raise RuntimeError(
            "No M7.3B candidates available for M7.3C."
        )

    return result.sort_values(
        [
            "priority_rank",
            "lead_rsid",
            "trf_id",
        ],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)


# ============================================================================
# Local authentication state
# ============================================================================


def evaluate_token_configuration(
    *,
    token_environment_variable: str,
    environment: dict[str, str],
) -> dict[str, Any]:
    """
    Evaluate local token-file configuration without reading token contents.

    The environment variable must contain a path to the GDC token file.
    """

    token_path_raw = _normalize_text(
        environment.get(
            token_environment_variable
        )
    )

    if token_path_raw is None:

        return {
            "token_environment_variable":
                token_environment_variable,

            "token_path_configured":
                False,

            "token_file_exists":
                False,

            "token_file_is_file":
                False,

            "token_file_readable":
                False,

            "token_path":
                None,

            "token_contents_read":
                False,
        }

    token_path = Path(
        token_path_raw
    ).expanduser()

    token_file_exists = token_path.exists()
    token_file_is_file = token_path.is_file()

    token_file_readable = False

    if token_file_exists and token_file_is_file:
        try:
            with token_path.open(
                "rb"
            ) as handle:
                handle.read(0)

            token_file_readable = True

        except OSError:
            token_file_readable = False

    return {
        "token_environment_variable":
            token_environment_variable,

        "token_path_configured":
            True,

        "token_file_exists":
            bool(
                token_file_exists
            ),

        "token_file_is_file":
            bool(
                token_file_is_file
            ),

        "token_file_readable":
            bool(
                token_file_readable
            ),

        # Intentionally store only the configured path.
        # Token contents are never read.
        "token_path":
            str(
                token_path
            ),

        "token_contents_read":
            False,
    }


# ============================================================================
# Resource readiness
# ============================================================================


def build_resource_readiness(
    *,
    resources: pd.DataFrame,
) -> pd.DataFrame:
    """Build one row per controlled aligned BAM."""

    if resources.empty:
        return resources.copy()

    result = resources.copy()

    result[
        "released_resource"
    ] = (
        result[
            "state"
        ]
        .fillna("")
        .astype(str)
        .str.lower()
        .eq("released")
    )

    result[
        "project_scope_valid"
    ] = (
        result[
            "target_project_query_scoped"
        ].astype(bool)
        &
        result[
            "target_project_observed"
        ].astype(bool)
    )

    result[
        "resource_ready_for_controlled_download"
    ] = (
        result[
            "accepted_for_direct_trf_quantification"
        ].astype(bool)
        &
        result[
            "released_resource"
        ].astype(bool)
        &
        result[
            "project_scope_valid"
        ].astype(bool)
    )

    result[
        "downloaded_in_m7_3c"
    ] = False

    result[
        "read_counting_performed"
    ] = False

    return result.reset_index(drop=True)


# ============================================================================
# Candidate readiness
# ============================================================================


def build_candidate_readiness(
    *,
    candidates: pd.DataFrame,
    resources: pd.DataFrame,
    token_state: dict[str, Any],
    access_probe_enabled: bool,
    statuses: dict[str, Any],
) -> pd.DataFrame:
    """Resolve operational controlled-access readiness per candidate."""

    ready_resource_count = int(
        resources[
            "resource_ready_for_controlled_download"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
        if not resources.empty
        else 0
    )

    rows: list[dict[str, Any]] = []

    for record in candidates.to_dict(
        orient="records"
    ):

        sequence_ready = bool(
            record[
                "sequence_ready_for_quantification"
            ]
        )

        if not sequence_ready:

            status = statuses[
                "sequence_blocked"
            ]

            access_configuration_ready = False
            ready_for_download = False

        elif ready_resource_count == 0:

            status = statuses[
                "controlled_resource_absent"
            ]

            access_configuration_ready = False
            ready_for_download = False

        elif not token_state[
            "token_path_configured"
        ]:

            status = statuses[
                "token_not_configured"
            ]

            access_configuration_ready = False
            ready_for_download = False

        elif not token_state[
            "token_file_exists"
        ]:

            status = statuses[
                "token_file_missing"
            ]

            access_configuration_ready = False
            ready_for_download = False

        elif not token_state[
            "token_file_readable"
        ]:

            status = statuses[
                "token_file_unreadable"
            ]

            access_configuration_ready = False
            ready_for_download = False

        elif not access_probe_enabled:

            status = statuses[
                "authorization_not_tested"
            ]

            access_configuration_ready = True
            ready_for_download = False

        else:

            # M7.3C v1 intentionally does not interpret authorization as
            # confirmed without an explicit remote access probe result.
            status = statuses[
                "authorization_not_tested"
            ]

            access_configuration_ready = True
            ready_for_download = False

        rows.append(
            {
                **record,

                "controlled_resource_files_ready":
                    ready_resource_count,

                "token_path_configured":
                    bool(
                        token_state[
                            "token_path_configured"
                        ]
                    ),

                "token_file_exists":
                    bool(
                        token_state[
                            "token_file_exists"
                        ]
                    ),

                "token_file_readable":
                    bool(
                        token_state[
                            "token_file_readable"
                        ]
                    ),

                "access_configuration_ready":
                    access_configuration_ready,

                "remote_authorization_tested":
                    False,

                "remote_authorization_confirmed":
                    False,

                "ready_for_controlled_download":
                    ready_for_download,

                "controlled_access_readiness_status":
                    str(
                        status
                    ),

                "bam_download_performed":
                    False,

                "read_counting_performed":
                    False,

                "trf_quantification_performed":
                    False,
            }
        )

    return pd.DataFrame(rows)


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    candidate_readiness: pd.DataFrame,
    resource_readiness: pd.DataFrame,
    token_state: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build M7.3C one-row summary."""

    resource_ready_count = int(
        resource_readiness[
            "resource_ready_for_controlled_download"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
        if not resource_readiness.empty
        else 0
    )

    sequence_ready_count = int(
        candidate_readiness[
            "sequence_ready_for_quantification"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    configuration_ready_count = int(
        candidate_readiness[
            "access_configuration_ready"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    ready_for_download_count = int(
        candidate_readiness[
            "ready_for_controlled_download"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    if ready_for_download_count > 0:

        overall_status = config[
            "statuses"
        ][
            "overall"
        ][
            "ready_for_download"
        ]

        next_stage = config[
            "next_stage"
        ][
            "if_ready_for_download"
        ]

    elif (
        token_state[
            "token_path_configured"
        ]
        and
        token_state[
            "token_file_exists"
        ]
        and
        token_state[
            "token_file_readable"
        ]
        and
        configuration_ready_count
        >
        0
    ):

        overall_status = config[
            "statuses"
        ][
            "overall"
        ][
            "authorization_pending"
        ]

        next_stage = config[
            "next_stage"
        ][
            "if_authorization_pending"
        ]

    elif sequence_ready_count > 0 and resource_ready_count > 0:

        overall_status = config[
            "statuses"
        ][
            "overall"
        ][
            "access_configuration_required"
        ]

        next_stage = config[
            "next_stage"
        ][
            "if_access_not_configured"
        ]

    else:

        overall_status = config[
            "statuses"
        ][
            "overall"
        ][
            "blocked"
        ]

        next_stage = config[
            "next_stage"
        ][
            "if_blocked"
        ]

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.3C",

                "candidate_trfs_assessed":
                    int(
                        len(candidate_readiness)
                    ),

                "candidate_sequences_ready":
                    sequence_ready_count,

                "controlled_bam_files_assessed":
                    int(
                        len(resource_readiness)
                    ),

                "controlled_bam_files_resource_ready":
                    resource_ready_count,

                "token_path_configured":
                    bool(
                        token_state[
                            "token_path_configured"
                        ]
                    ),

                "token_file_exists":
                    bool(
                        token_state[
                            "token_file_exists"
                        ]
                    ),

                "token_file_readable":
                    bool(
                        token_state[
                            "token_file_readable"
                        ]
                    ),

                "candidate_access_configurations_ready":
                    configuration_ready_count,

                "candidates_ready_for_controlled_download":
                    ready_for_download_count,

                "remote_authorization_tested":
                    False,

                "bam_download_performed":
                    False,

                "read_counting_performed":
                    False,

                "trf_quantification_performed":
                    False,

                "overall_readiness_status":
                    str(
                        overall_status
                    ),

                "next_stage":
                    str(
                        next_stage
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def assess_controlled_access_readiness(
    *,
    m7_3b_qc: dict[str, Any],
    config: dict[str, Any],
    environment: dict[str, str],
) -> ControlledAccessReadinessResult:
    """Execute M7.3C Controlled-Access Quantification Readiness."""

    resources = extract_controlled_resources(
        m7_3b_qc=m7_3b_qc,
    )

    candidates = extract_candidates(
        m7_3b_qc=m7_3b_qc,
    )

    token_variable = str(
        config[
            "authentication"
        ][
            "token_path_environment_variable"
        ]
    )

    token_state = evaluate_token_configuration(
        token_environment_variable=token_variable,
        environment=environment,
    )

    resource_readiness = build_resource_readiness(
        resources=resources,
    )

    candidate_readiness = build_candidate_readiness(
        candidates=candidates,
        resources=resource_readiness,
        token_state=token_state,
        access_probe_enabled=bool(
            config[
                "access_probe"
            ][
                "enabled"
            ]
        ),
        statuses=config[
            "statuses"
        ][
            "candidate"
        ],
    )

    summary = build_summary(
        candidate_readiness=candidate_readiness,
        resource_readiness=resource_readiness,
        token_state=token_state,
        config=config,
    )

    return ControlledAccessReadinessResult(
        resource_readiness=resource_readiness,
        candidate_readiness=candidate_readiness,
        summary=summary,
    )