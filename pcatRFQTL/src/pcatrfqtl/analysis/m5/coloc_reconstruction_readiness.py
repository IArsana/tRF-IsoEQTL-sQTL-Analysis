"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/coloc_reconstruction_readiness.py

Description:
    Core logic for M5.6B colocalization reconstruction feasibility audit.

    This module evaluates whether the unavailable full Moradi regulatory-QTL
    summary statistics could be reconstructed with sufficient fidelity from
    underlying source data.

    It distinguishes:

        EXACTLY_REPRODUCIBLE
            All components required for close reproduction are available.

        APPROXIMATE_RECONSTRUCTION_FEASIBLE
            The majority of required components are available or derivable,
            but exact reproduction cannot be claimed.

        THEORETICALLY_POSSIBLE_BUT_NOT_REPRODUCIBLE
            Reconstruction is conceptually possible, but one or more critical
            original-analysis components remain unidentified or incomplete.

        RECONSTRUCTION_NOT_FEASIBLE
            Required underlying inputs are unavailable to a degree that
            prevents a meaningful reconstruction.

    Scientific safeguards:
        - Availability of TCGA data alone does not imply reproducibility.
        - Unknown original covariates are not silently inferred.
        - Approximate reconstruction is not treated as exact replication.
        - This stage does not execute MatrixEQTL.
        - This stage does not perform colocalization.
        - This stage does not perform causal inference.

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
# Status constants
# ============================================================================


AVAILABLE = "AVAILABLE"

DERIVABLE = "DERIVABLE"

PARTIAL = "PARTIAL"

NOT_IDENTIFIED = "NOT_IDENTIFIED"

NOT_REQUIRED = "NOT_REQUIRED"


VALID_COMPONENT_STATUSES = {
    AVAILABLE,
    DERIVABLE,
    PARTIAL,
    NOT_IDENTIFIED,
    NOT_REQUIRED,
}


EXACTLY_REPRODUCIBLE = (
    "EXACTLY_REPRODUCIBLE"
)

APPROXIMATE_RECONSTRUCTION_FEASIBLE = (
    "APPROXIMATE_RECONSTRUCTION_FEASIBLE"
)

THEORETICALLY_POSSIBLE_BUT_NOT_REPRODUCIBLE = (
    "THEORETICALLY_POSSIBLE_BUT_NOT_REPRODUCIBLE"
)

RECONSTRUCTION_NOT_FEASIBLE = (
    "RECONSTRUCTION_NOT_FEASIBLE"
)


# ============================================================================
# Result models
# ============================================================================


@dataclass(frozen=True)
class ReconstructionReadinessResult:
    """Container for M5.6B reconstruction-readiness outputs."""

    components: pd.DataFrame

    summary: pd.DataFrame

    candidates: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _normalize_status(
    value: Any,
) -> str:
    """Normalize and validate one reconstruction-component status."""

    status = (
        str(value)
        .strip()
        .upper()
    )

    if status not in VALID_COMPONENT_STATUSES:

        raise ValueError(
            "Invalid reconstruction component status: "
            f"{value}. "
            f"Allowed={sorted(VALID_COMPONENT_STATUSES)}"
        )

    return status


def _flatten_notes(
    value: Any,
) -> str | None:
    """Convert optional YAML note collection into deterministic text."""

    if value is None:

        return None

    if isinstance(
        value,
        list,
    ):

        return " | ".join(
            str(item)
            for item in value
        )

    return str(
        value
    )


# ============================================================================
# Component audit
# ============================================================================


def build_reconstruction_component_table(
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one row per reconstruction component."""

    components = config.get(
        "components"
    )

    if not isinstance(
        components,
        dict,
    ):

        raise ValueError(
            "M5.6B config must contain a 'components' mapping."
        )

    rows: list[
        dict[str, Any]
    ] = []

    for component_name, properties in components.items():

        if not isinstance(
            properties,
            dict,
        ):

            raise ValueError(
                "Component configuration must be a mapping: "
                f"{component_name}"
            )

        required = bool(
            properties.get(
                "required",
                True,
            )
        )

        status = _normalize_status(
            properties.get(
                "status",
                NOT_IDENTIFIED,
            )
        )

        rows.append(
            {
                "component":
                    str(
                        component_name
                    ),

                "required":
                    required,

                "status":
                    status,

                "is_available":
                    bool(
                        status == AVAILABLE
                    ),

                "is_derivable":
                    bool(
                        status == DERIVABLE
                    ),

                "is_partial":
                    bool(
                        status == PARTIAL
                    ),

                "is_not_identified":
                    bool(
                        status == NOT_IDENTIFIED
                    ),

                "is_nonblocking_for_approximate_reconstruction":
                    bool(
                        status
                        in {
                            AVAILABLE,
                            DERIVABLE,
                            PARTIAL,
                            NOT_REQUIRED,
                        }
                    ),

                "notes":
                    _flatten_notes(
                        properties.get(
                            "notes"
                        )
                    ),
            }
        )

    dataframe = pd.DataFrame(
        rows
    )

    if dataframe.empty:

        raise ValueError(
            "M5.6B reconstruction component table is empty."
        )

    return dataframe


# ============================================================================
# Overall reconstruction assessment
# ============================================================================


def assess_reconstruction_feasibility(
    *,
    components: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Determine overall reconstruction feasibility.

    Exact reproduction is intentionally strict.

    Approximate reconstruction allows AVAILABLE, DERIVABLE, and PARTIAL
    components, but cannot support an exact replication claim.
    """

    required_columns = {
        "component",
        "required",
        "status",
        "is_available",
        "is_nonblocking_for_approximate_reconstruction",
        "is_not_identified",
    }

    missing = (
        required_columns
        - set(
            components.columns
        )
    )

    if missing:

        raise ValueError(
            "Reconstruction component table missing columns: "
            f"{sorted(missing)}"
        )

    policy = config[
        "readiness_policy"
    ]

    required = components.loc[
        components[
            "required"
        ].astype(
            bool
        )
    ].copy()

    required_count = int(
        len(
            required
        )
    )

    if required_count == 0:

        raise ValueError(
            "No required reconstruction components were defined."
        )

    available_count = int(
        required[
            "is_available"
        ]
        .astype(
            bool
        )
        .sum()
    )

    nonblocking_count = int(
        required[
            "is_nonblocking_for_approximate_reconstruction"
        ]
        .astype(
            bool
        )
        .sum()
    )

    missing_count = int(
        required[
            "is_not_identified"
        ]
        .astype(
            bool
        )
        .sum()
    )

    exact_fraction = float(
        available_count
        / required_count
    )

    nonblocking_fraction = float(
        nonblocking_count
        / required_count
    )

    minimum_nonblocking_fraction = float(
        policy.get(
            "minimum_fraction_nonblocking_for_approximate_reconstruction",
            0.80,
        )
    )

    exact_components = {
        str(value)
        for value
        in policy.get(
            "exact_reproduction_requires",
            [],
        )
    }

    required_status_lookup = {
        str(
            row[
                "component"
            ]
        ):
            str(
                row[
                    "status"
                ]
            )
        for _, row
        in required.iterrows()
    }

    exact_requirements_met = bool(
        exact_components
        and
        all(
            required_status_lookup.get(
                component
            )
            == AVAILABLE
            for component
            in exact_components
        )
    )

    approximate_threshold_met = bool(
        nonblocking_fraction
        >=
        minimum_nonblocking_fraction
    )

    if exact_requirements_met:

        readiness_state = (
            EXACTLY_REPRODUCIBLE
        )

        reconstruction_can_proceed = True

        exact_replication_claim_allowed = True

    elif (
        approximate_threshold_met
        and
        missing_count == 0
    ):

        readiness_state = (
            APPROXIMATE_RECONSTRUCTION_FEASIBLE
        )

        reconstruction_can_proceed = True

        exact_replication_claim_allowed = False

    elif approximate_threshold_met:

        readiness_state = (
            THEORETICALLY_POSSIBLE_BUT_NOT_REPRODUCIBLE
        )

        reconstruction_can_proceed = False

        exact_replication_claim_allowed = False

    else:

        readiness_state = (
            RECONSTRUCTION_NOT_FEASIBLE
        )

        reconstruction_can_proceed = False

        exact_replication_claim_allowed = False

    blocking_components = (
        required.loc[
            required[
                "is_not_identified"
            ].astype(
                bool
            ),
            "component",
        ]
        .astype(
            str
        )
        .tolist()
    )

    partial_components = (
        required.loc[
            required[
                "status"
            ].eq(
                PARTIAL
            ),
            "component",
        ]
        .astype(
            str
        )
        .tolist()
    )

    derivable_components = (
        required.loc[
            required[
                "status"
            ].eq(
                DERIVABLE
            ),
            "component",
        ]
        .astype(
            str
        )
        .tolist()
    )

    return pd.DataFrame(
        [
            {
                "resource":
                    "Moradi_2022",

                "required_component_count":
                    required_count,

                "fully_available_component_count":
                    available_count,

                "nonblocking_component_count":
                    nonblocking_count,

                "not_identified_component_count":
                    missing_count,

                "exact_availability_fraction":
                    exact_fraction,

                "approximate_nonblocking_fraction":
                    nonblocking_fraction,

                "minimum_approximate_threshold":
                    minimum_nonblocking_fraction,

                "exact_requirements_met":
                    exact_requirements_met,

                "approximate_threshold_met":
                    approximate_threshold_met,

                "blocking_components":
                    blocking_components,

                "partial_components":
                    partial_components,

                "derivable_components":
                    derivable_components,

                "reconstruction_readiness":
                    readiness_state,

                "reconstruction_can_proceed":
                    reconstruction_can_proceed,

                "exact_replication_claim_allowed":
                    exact_replication_claim_allowed,

                "formal_coloc_input_generated":
                    False,

                "formal_coloc_can_proceed_from_reconstruction":
                    False,
            }
        ]
    )


# ============================================================================
# Candidate propagation
# ============================================================================


def build_candidate_reconstruction_readiness(
    *,
    summary: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Propagate reconstruction status to the three candidate loci."""

    if summary.empty:

        raise ValueError(
            "Reconstruction summary is empty."
        )

    row = summary.iloc[
        0
    ]

    readiness = str(
        row[
            "reconstruction_readiness"
        ]
    )

    reconstruction_can_proceed = bool(
        row[
            "reconstruction_can_proceed"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for lead_rsid in config[
        "candidate_leads"
    ]:

        rows.append(
            {
                "lead_rsid":
                    str(
                        lead_rsid
                    ),

                "reconstruction_readiness":
                    readiness,

                "reconstruction_can_proceed":
                    reconstruction_can_proceed,

                "full_qtl_summary_statistics_generated":
                    False,

                "formal_coloc_ready":
                    False,

                "coloc_abf_ready":
                    False,

                "coloc_susie_ready":
                    False,

                "recommended_action":
                    (
                        "RECONSTRUCT_AND_REVALIDATE"
                        if reconstruction_can_proceed
                        else
                        "DO_NOT_USE_RECONSTRUCTION_FOR_FORMAL_COLOC"
                    ),

                "exact_replication_claim_allowed":
                    bool(
                        row[
                            "exact_replication_claim_allowed"
                        ]
                    ),

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Public API
# ============================================================================


def assess_coloc_reconstruction_readiness(
    config: dict[str, Any],
) -> ReconstructionReadinessResult:
    """Execute deterministic M5.6B reconstruction feasibility audit."""

    components = (
        build_reconstruction_component_table(
            config
        )
    )

    summary = (
        assess_reconstruction_feasibility(
            components=components,
            config=config,
        )
    )

    candidates = (
        build_candidate_reconstruction_readiness(
            summary=summary,
            config=config,
        )
    )

    return ReconstructionReadinessResult(
        components=components,
        summary=summary,
        candidates=candidates,
    )