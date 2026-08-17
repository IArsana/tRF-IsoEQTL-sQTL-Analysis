"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/endpoint_resolution.py

Description:
    Core logic for M7.2B Clinical Endpoint Readiness Resolution.

    This stage formally resolves which clinical endpoints are defensible
    for downstream candidate association analysis after M7.2 cohort
    construction.

    M7.2B distinguishes:

        1. time-to-event endpoints;
        2. clinicopathologic endpoints.

    A time-to-event endpoint must have:
        - sufficient usable cases;
        - sufficient events;
        - sufficient completeness;
        - analyzable censoring when required.

    Clinicopathologic endpoints are assessed based on:
        - source availability;
        - non-missing fraction;
        - usable case count;
        - observed variation.

    M7.2B does NOT:
        - fit Cox models;
        - perform Kaplan-Meier analysis;
        - test tRF associations;
        - impute missing clinical values;
        - reconstruct unavailable stages or grades;
        - select endpoints based on statistical significance;
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

from dataclasses import dataclass
from typing import Any

import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class EndpointResolutionResult:
    """Container for M7.2B outputs."""

    time_to_event_resolution: pd.DataFrame
    clinicopathologic_readiness: pd.DataFrame
    final_endpoint_resolution: pd.DataFrame
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

        if pd.isna(
            value
        ):
            return None

    except (
        TypeError,
        ValueError,
    ):
        pass

    normalized = str(
        value
    ).strip()

    return normalized or None


def _to_int(
    value: Any,
    *,
    default: int = 0,
) -> int:
    """Convert optional scalar to integer."""

    if value is None:
        return default

    try:

        return int(
            value
        )

    except (
        TypeError,
        ValueError,
    ):

        return default


def _to_float(
    value: Any,
    *,
    default: float = 0.0,
) -> float:
    """Convert optional scalar to float."""

    if value is None:
        return default

    try:

        return float(
            value
        )

    except (
        TypeError,
        ValueError,
    ):

        return default


# ============================================================================
# Time-to-event resolution
# ============================================================================


def resolve_time_to_event_endpoints(
    *,
    m7_2_qc: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Resolve M7.2 time-to-event endpoints.

    Upstream readiness is preserved. M7.2B does not relax thresholds after
    observing endpoint availability.
    """

    records = m7_2_qc.get(
        "endpoint_readiness",
        [],
    )

    if not isinstance(
        records,
        list,
    ):

        raise RuntimeError(
            "M7.2 endpoint_readiness must be a list."
        )

    policy = config[
        "time_to_event_policy"
    ]

    endpoint_config = config[
        "time_to_event_endpoints"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        endpoint_id = _normalize_text(
            record.get(
                "endpoint_id"
            )
        )

        if endpoint_id is None:

            continue

        if endpoint_id not in endpoint_config:

            continue

        usable_cases = _to_int(
            record.get(
                "usable_cases"
            )
        )

        usable_fraction = _to_float(
            record.get(
                "usable_fraction"
            )
        )

        events = _to_int(
            record.get(
                "events"
            )
        )

        censored = _to_int(
            record.get(
                "analyzable_censored_cases"
            )
        )

        raw_positive_events = _to_int(
            record.get(
                "raw_positive_event_labels"
            )
        )

        events_without_usable_time = _to_int(
            record.get(
                "events_without_usable_time"
            )
        )

        upstream_ready = bool(
            record.get(
                "endpoint_ready",
                False,
            )
        )

        upstream_status = _normalize_text(
            record.get(
                "endpoint_status"
            )
        )

        has_sufficient_cases = (
            usable_cases
            >=
            int(
                policy[
                    "minimum_usable_cases"
                ]
            )
        )

        has_sufficient_fraction = (
            usable_fraction
            >=
            float(
                policy[
                    "minimum_usable_fraction"
                ]
            )
        )

        has_sufficient_events = (
            events
            >=
            int(
                policy[
                    "minimum_events"
                ]
            )
        )

        has_censoring = (
            censored
            >
            0
        )

        require_censoring = bool(
            policy[
                "require_censored_observations"
            ]
        )

        final_ready = bool(
            upstream_ready
            and
            has_sufficient_cases
            and
            has_sufficient_fraction
            and
            has_sufficient_events
            and
            (
                has_censoring
                or
                not require_censoring
            )
        )

        endpoint_rules = endpoint_config[
            endpoint_id
        ]

        # ------------------------------------------------------------------
        # Final status
        # ------------------------------------------------------------------

        if final_ready:

            resolution_status = (
                "TIME_TO_EVENT_ENDPOINT_READY"
            )

            downstream_use = (
                "ELIGIBLE_FOR_TIME_TO_EVENT_MODELING"
            )

        elif not has_sufficient_events:

            resolution_status = str(
                endpoint_rules[
                    "insufficient_events_status"
                ]
            )

            downstream_use = str(
                endpoint_rules[
                    "if_not_ready"
                ]
            )

        elif not has_sufficient_fraction:

            resolution_status = str(
                endpoint_rules[
                    "insufficient_completeness_status"
                ]
            )

            downstream_use = str(
                endpoint_rules[
                    "if_not_ready"
                ]
            )

        elif require_censoring and not has_censoring:

            resolution_status = str(
                endpoint_rules[
                    "no_censoring_status"
                ]
            )

            downstream_use = str(
                endpoint_rules[
                    "if_not_ready"
                ]
            )

        elif not has_sufficient_cases:

            resolution_status = (
                "TIME_TO_EVENT_NOT_SELECTED_INSUFFICIENT_USABLE_CASES"
            )

            downstream_use = str(
                endpoint_rules[
                    "if_not_ready"
                ]
            )

        else:

            resolution_status = (
                "TIME_TO_EVENT_NOT_SELECTED_UPSTREAM_NOT_READY"
            )

            downstream_use = str(
                endpoint_rules[
                    "if_not_ready"
                ]
            )

        rows.append(
            {
                "endpoint_id":
                    endpoint_id,

                "endpoint_class":
                    str(
                        endpoint_rules[
                            "endpoint_class"
                        ]
                    ),

                "usable_cases":
                    usable_cases,

                "usable_fraction":
                    usable_fraction,

                "analyzable_events":
                    events,

                "analyzable_censored_cases":
                    censored,

                "raw_positive_event_labels":
                    raw_positive_events,

                "events_without_usable_time":
                    events_without_usable_time,

                "has_sufficient_cases":
                    has_sufficient_cases,

                "has_sufficient_completeness":
                    has_sufficient_fraction,

                "has_sufficient_events":
                    has_sufficient_events,

                "has_analyzable_censoring":
                    has_censoring,

                "upstream_endpoint_ready":
                    upstream_ready,

                "upstream_endpoint_status":
                    upstream_status,

                "final_time_to_event_ready":
                    final_ready,

                "primary_endpoint_selected":
                    False,

                "resolution_status":
                    resolution_status,

                "downstream_use":
                    downstream_use,

                "threshold_relaxed":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No supported time-to-event endpoints found in M7.2 QC."
        )

    return result


# ============================================================================
# Clinicopathologic readiness
# ============================================================================


def resolve_clinicopathologic_endpoints(
    *,
    clinical_cohort: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Assess clinicopathologic variables for downstream association analysis.

    This stage only determines readiness. No candidate association is tested.
    """

    policy = config[
        "clinicopathologic_policy"
    ]

    variables = config[
        "clinicopathologic_variables"
    ]

    statuses = config[
        "statuses"
    ][
        "clinicopathologic"
    ]

    minimum_fraction = float(
        policy[
            "minimum_non_missing_fraction"
        ]
    )

    minimum_cases = int(
        policy[
            "minimum_non_missing_cases"
        ]
    )

    require_variation = bool(
        policy[
            "require_variation"
        ]
    )

    total_cases = int(
        len(
            clinical_cohort
        )
    )

    rows: list[
        dict[str, Any]
    ] = []

    for endpoint_id, definition in variables.items():

        source_column = str(
            definition[
                "source_column"
            ]
        )

        if source_column not in clinical_cohort.columns:

            rows.append(
                {
                    "endpoint_id":
                        str(
                            endpoint_id
                        ),

                    "source_column":
                        source_column,

                    "endpoint_class":
                        str(
                            definition[
                                "endpoint_class"
                            ]
                        ),

                    "total_cases":
                        total_cases,

                    "non_missing_cases":
                        0,

                    "non_missing_fraction":
                        0.0,

                    "unique_non_missing_values":
                        0,

                    "has_variation":
                        False,

                    "endpoint_ready":
                        False,

                    "endpoint_status":
                        str(
                            statuses[
                                "insufficient_completeness"
                            ]
                        ),

                    "downstream_role":
                        str(
                            definition[
                                "downstream_role"
                            ]
                        ),

                    "missing_imputed":
                        False,
                }
            )

            continue

        series = clinical_cohort[
            source_column
        ]

        non_missing = series.dropna()

        non_missing_cases = int(
            len(
                non_missing
            )
        )

        non_missing_fraction = (
            non_missing_cases
            /
            total_cases
            if total_cases
            else 0.0
        )

        normalized_values = (
            non_missing
            .astype(
                str
            )
            .str.strip()
        )

        normalized_values = normalized_values.loc[
            normalized_values.ne(
                ""
            )
        ]

        unique_values = int(
            normalized_values.nunique(
                dropna=True
            )
        )

        has_variation = bool(
            unique_values
            >=
            2
        )

        enough_cases = (
            non_missing_cases
            >=
            minimum_cases
        )

        enough_fraction = (
            non_missing_fraction
            >=
            minimum_fraction
        )

        variation_pass = (
            has_variation
            or
            not require_variation
        )

        ready = bool(
            enough_cases
            and
            enough_fraction
            and
            variation_pass
        )

        if ready:

            status = str(
                statuses[
                    "ready"
                ]
            )

        elif not enough_cases:

            status = str(
                statuses[
                    "insufficient_cases"
                ]
            )

        elif not enough_fraction:

            status = str(
                statuses[
                    "insufficient_completeness"
                ]
            )

        else:

            status = str(
                statuses[
                    "insufficient_variation"
                ]
            )

        rows.append(
            {
                "endpoint_id":
                    str(
                        endpoint_id
                    ),

                "source_column":
                    source_column,

                "endpoint_class":
                    str(
                        definition[
                            "endpoint_class"
                        ]
                    ),

                "total_cases":
                    total_cases,

                "non_missing_cases":
                    non_missing_cases,

                "non_missing_fraction":
                    non_missing_fraction,

                "unique_non_missing_values":
                    unique_values,

                "has_variation":
                    has_variation,

                "endpoint_ready":
                    ready,

                "endpoint_status":
                    status,

                "downstream_role":
                    str(
                        definition[
                            "downstream_role"
                        ]
                    ),

                "missing_imputed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Final endpoint resolution
# ============================================================================


def build_final_endpoint_resolution(
    *,
    time_to_event: pd.DataFrame,
    clinicopathologic: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Resolve the final M7 clinical-validation route.
    """

    primary_time_to_event_available = bool(
        time_to_event[
            "final_time_to_event_ready"
        ]
        .astype(
            bool
        )
        .any()
    )

    ready_clinicopathologic = clinicopathologic.loc[
        clinicopathologic[
            "endpoint_ready"
        ]
        .astype(
            bool
        )
    ].copy()

    ready_ids = sorted(
        ready_clinicopathologic[
            "endpoint_id"
        ]
        .astype(
            str
        )
        .tolist()
    )

    time_status = (
        config[
            "statuses"
        ][
            "time_to_event"
        ][
            "primary_endpoint_available"
        ]
        if primary_time_to_event_available
        else
        config[
            "statuses"
        ][
            "time_to_event"
        ][
            "no_primary_endpoint"
        ]
    )

    clinical_route_continues = bool(
        ready_ids
    )

    overall_status = (
        config[
            "statuses"
        ][
            "overall"
        ][
            "clinical_route_continues"
        ]
        if clinical_route_continues
        else
        config[
            "statuses"
        ][
            "overall"
        ][
            "clinical_route_blocked"
        ]
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.2B",

                "primary_time_to_event_endpoint_available":
                    primary_time_to_event_available,

                "primary_time_to_event_endpoint":
                    None,

                "time_to_event_resolution_status":
                    str(
                        time_status
                    ),

                "ready_clinicopathologic_endpoint_count":
                    int(
                        len(
                            ready_ids
                        )
                    ),

                "ready_clinicopathologic_endpoints":
                    ready_ids,

                "clinical_validation_route_continues":
                    clinical_route_continues,

                "overall_resolution_status":
                    str(
                        overall_status
                    ),

                "primary_endpoint_selected":
                    False,

                "survival_model_allowed":
                    primary_time_to_event_available,

                "clinicopathologic_association_allowed":
                    clinical_route_continues,

                "endpoint_selected_by_significance":
                    False,

                "threshold_relaxed":
                    False,

                "next_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "primary"
                        ]
                    ),
            }
        ]
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    time_to_event: pd.DataFrame,
    clinicopathologic: pd.DataFrame,
    final_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.2B summary."""

    final = final_resolution.iloc[
        0
    ]

    ready_clinical = clinicopathologic.loc[
        clinicopathologic[
            "endpoint_ready"
        ]
        .astype(
            bool
        )
    ]

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.2B",

                "time_to_event_endpoints_assessed":
                    int(
                        len(
                            time_to_event
                        )
                    ),

                "time_to_event_endpoints_ready":
                    int(
                        time_to_event[
                            "final_time_to_event_ready"
                        ]
                        .astype(
                            bool
                        )
                        .sum()
                    ),

                "primary_time_to_event_endpoint_available":
                    bool(
                        final[
                            "primary_time_to_event_endpoint_available"
                        ]
                    ),

                "clinicopathologic_endpoints_assessed":
                    int(
                        len(
                            clinicopathologic
                        )
                    ),

                "clinicopathologic_endpoints_ready":
                    int(
                        len(
                            ready_clinical
                        )
                    ),

                "ready_clinicopathologic_endpoints":
                    sorted(
                        ready_clinical[
                            "endpoint_id"
                        ]
                        .astype(
                            str
                        )
                        .tolist()
                    ),

                "clinical_validation_route_continues":
                    bool(
                        final[
                            "clinical_validation_route_continues"
                        ]
                    ),

                "primary_endpoint_selected":
                    False,

                "survival_model_fitted":
                    False,

                "clinicopathologic_model_fitted":
                    False,

                "threshold_relaxed":
                    False,

                "overall_resolution_status":
                    str(
                        final[
                            "overall_resolution_status"
                        ]
                    ),

                "next_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "primary"
                        ]
                    ),

                "downstream_if_trf_quantifiable":
                    str(
                        config[
                            "next_stage"
                        ][
                            "downstream_if_trf_quantifiable"
                        ]
                    ),
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def resolve_clinical_endpoints(
    *,
    m7_2_qc: dict[str, Any],
    clinical_cohort: pd.DataFrame,
    config: dict[str, Any],
) -> EndpointResolutionResult:
    """Execute M7.2B Clinical Endpoint Readiness Resolution."""

    time_to_event = resolve_time_to_event_endpoints(
        m7_2_qc=m7_2_qc,
        config=config,
    )

    clinicopathologic = resolve_clinicopathologic_endpoints(
        clinical_cohort=clinical_cohort,
        config=config,
    )

    final_resolution = build_final_endpoint_resolution(
        time_to_event=time_to_event,
        clinicopathologic=clinicopathologic,
        config=config,
    )

    summary = build_summary(
        time_to_event=time_to_event,
        clinicopathologic=clinicopathologic,
        final_resolution=final_resolution,
        config=config,
    )

    return EndpointResolutionResult(
        time_to_event_resolution=time_to_event,
        clinicopathologic_readiness=clinicopathologic,
        final_endpoint_resolution=final_resolution,
        summary=summary,
    )