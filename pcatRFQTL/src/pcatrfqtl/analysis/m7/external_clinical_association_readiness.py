"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/external_clinical_association_readiness.py

Description:
    Core logic for M7.4D External Clinical Association Readiness.

    This stage evaluates whether the external GSE80400 cohort has sufficient
    clinical structure for inferential association between exact candidate
    sequence abundance and source-supported clinical endpoints.

    The stage performs readiness assessment only.

    It does NOT:
        - fit association models;
        - calculate p-values;
        - estimate effect sizes;
        - perform survival analysis;
        - infer missing recurrence status;
        - infer run-to-GEO mappings from accession order;
        - claim mature 24-nt tRF abundance;
        - claim clinical validation.

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


@dataclass(frozen=True)
class ExternalClinicalAssociationReadinessResult:
    """Container for M7.4D outputs."""

    sample_inventory: pd.DataFrame
    endpoint_readiness: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Sample inventory
# ============================================================================


def build_sample_inventory(
    *,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build conservative source-supported clinical sample inventory."""

    source_samples = config[
        "source_samples"
    ]

    annotations = config[
        "clinical_annotations"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for sample in source_samples:

        geo_accession = str(
            sample[
                "geo_accession"
            ]
        )

        if geo_accession not in annotations:

            raise RuntimeError(
                "Clinical annotation missing for "
                f"{geo_accession}."
            )

        annotation = annotations[
            geo_accession
        ]

        rows.append(
            {
                "geo_accession":
                    geo_accession,

                "source_title":
                    str(
                        sample[
                            "source_title"
                        ]
                    ),

                "source_group_number":
                    int(
                        sample[
                            "source_group_number"
                        ]
                    ),

                "broad_class":
                    annotation.get(
                        "broad_class"
                    ),

                "prostate_cancer":
                    annotation.get(
                        "prostate_cancer"
                    ),

                "bph":
                    annotation.get(
                        "bph"
                    ),

                "gleason_score":
                    annotation.get(
                        "gleason_score"
                    ),

                "recurrence_explicit":
                    annotation.get(
                        "recurrence_explicit"
                    ),

                "cured_explicit":
                    annotation.get(
                        "cured_explicit"
                    ),

                "hormone_refractory_explicit":
                    annotation.get(
                        "hormone_refractory_explicit"
                    ),

                "lymph_node_pca_explicit":
                    annotation.get(
                        "lymph_node_pca_explicit"
                    ),

                "ffpe":
                    annotation.get(
                        "ffpe"
                    ),

                "derived_from_group":
                    annotation.get(
                        "derived_from_group"
                    ),

                "run_accession":
                    None,

                "run_sample_mapping_verified":
                    False,

                "candidate_cpm":
                    None,

                "candidate_abundance_linked_to_sample":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "M7.4D sample inventory is empty."
        )

    if result[
        "geo_accession"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicate GEO sample accessions detected."
        )

    return result.sort_values(
        "source_group_number",
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Endpoint helpers
# ============================================================================


def _endpoint_row(
    *,
    endpoint: str,
    endpoint_type: str,
    usable_cases: int,
    group_counts: dict[str, int],
    mapping_verified: bool,
    minimum_group_support_met: bool,
    metadata_complete: bool,
    heterogeneous_definition: bool,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Construct one endpoint-readiness row."""

    statuses = config[
        "statuses"
    ][
        "endpoint"
    ]

    if not mapping_verified:

        ready = False

        status = statuses[
            "mapping_unverified"
        ]

    elif not metadata_complete:

        ready = False

        status = statuses[
            "incomplete_metadata"
        ]

    elif heterogeneous_definition:

        ready = False

        status = statuses[
            "heterogeneous"
        ]

    elif not minimum_group_support_met:

        ready = False

        status = statuses[
            "insufficient_group_support"
        ]

    else:

        ready = True

        status = statuses[
            "ready"
        ]

    return {
        "endpoint":
            endpoint,

        "endpoint_type":
            endpoint_type,

        "usable_cases":
            int(
                usable_cases
            ),

        "group_counts":
            [
                f"{key}:{value}"
                for key, value in sorted(
                    group_counts.items()
                )
            ],

        "number_of_groups":
            int(
                len(
                    group_counts
                )
            ),

        "run_sample_mapping_verified":
            bool(
                mapping_verified
            ),

        "metadata_complete":
            bool(
                metadata_complete
            ),

        "minimum_group_support_met":
            bool(
                minimum_group_support_met
            ),

        "heterogeneous_definition":
            bool(
                heterogeneous_definition
            ),

        "association_ready":
            bool(
                ready
            ),

        "readiness_status":
            str(
                status
            ),

        "association_model_fitted":
            False,

        "p_value_calculated":
            False,

        "effect_size_estimated":
            False,
    }


# ============================================================================
# Endpoint readiness
# ============================================================================


def build_endpoint_readiness(
    *,
    sample_inventory: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Assess all predefined external clinical endpoints."""

    mapping_verified = bool(
        config[
            "run_mapping_policy"
        ][
            "run_to_geo_mapping_verified"
        ]
    )

    endpoints = config[
        "endpoints"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    # ----------------------------------------------------------------------
    # PCa vs reference
    #
    # Do not combine NAP and BPH automatically.
    # ----------------------------------------------------------------------

    pca_count = int(
        sample_inventory[
            "prostate_cancer"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    reference_count = int(
        (
            sample_inventory[
                "broad_class"
            ]
            ==
            "REFERENCE"
        ).sum()
    )

    tumor_reference_min = int(
        endpoints[
            "tumor_reference"
        ][
            "minimum_cases_per_group"
        ]
    )

    rows.append(
        _endpoint_row(
            endpoint="tumor_reference",
            endpoint_type="BINARY",
            usable_cases=(
                pca_count
                +
                reference_count
            ),
            group_counts={
                "PCA":
                    pca_count,

                "REFERENCE":
                    reference_count,
            },
            mapping_verified=mapping_verified,
            minimum_group_support_met=bool(
                pca_count
                >=
                tumor_reference_min
                and
                reference_count
                >=
                tumor_reference_min
            ),
            metadata_complete=True,
            heterogeneous_definition=False,
            config=config,
        )
    )

    # ----------------------------------------------------------------------
    # Explicit recurrence vs cured
    # ----------------------------------------------------------------------

    recurrence_known = sample_inventory.loc[
        sample_inventory[
            "recurrence_explicit"
        ].notna()
        &
        sample_inventory[
            "prostate_cancer"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    ].copy()

    recurrent_count = int(
        (
            recurrence_known[
                "recurrence_explicit"
            ]
            ==
            True
        ).sum()
    )

    non_recurrent_count = int(
        (
            recurrence_known[
                "recurrence_explicit"
            ]
            ==
            False
        ).sum()
    )

    recurrence_min = int(
        endpoints[
            "recurrence"
        ][
            "minimum_cases_per_group"
        ]
    )

    rows.append(
        _endpoint_row(
            endpoint="recurrence",
            endpoint_type="BINARY",
            usable_cases=len(
                recurrence_known
            ),
            group_counts={
                "RECURRENT":
                    recurrent_count,

                "CURED_EXPLICIT":
                    non_recurrent_count,
            },
            mapping_verified=mapping_verified,
            minimum_group_support_met=bool(
                recurrent_count
                >=
                recurrence_min
                and
                non_recurrent_count
                >=
                recurrence_min
            ),
            metadata_complete=bool(
                len(
                    recurrence_known
                )
                >
                0
            ),
            heterogeneous_definition=False,
            config=config,
        )
    )

    # ----------------------------------------------------------------------
    # Gleason
    # ----------------------------------------------------------------------

    gleason = sample_inventory.loc[
        sample_inventory[
            "gleason_score"
        ].notna()
        &
        ~sample_inventory[
            "ffpe"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    ].copy()

    gleason_counts = (
        gleason[
            "gleason_score"
        ]
        .astype(
            int
        )
        .value_counts()
        .sort_index()
        .to_dict()
    )

    minimum_gleason_total = int(
        endpoints[
            "gleason"
        ][
            "minimum_total_samples"
        ]
    )

    minimum_gleason_group = int(
        endpoints[
            "gleason"
        ][
            "minimum_samples_per_category"
        ]
    )

    minimum_gleason_categories = int(
        endpoints[
            "gleason"
        ][
            "require_at_least_categories"
        ]
    )

    gleason_group_support = bool(
        len(
            gleason
        )
        >=
        minimum_gleason_total
        and
        len(
            gleason_counts
        )
        >=
        minimum_gleason_categories
        and
        all(
            int(
                count
            )
            >=
            minimum_gleason_group
            for count in gleason_counts.values()
        )
    )

    rows.append(
        _endpoint_row(
            endpoint="gleason",
            endpoint_type="ORDINAL",
            usable_cases=len(
                gleason
            ),
            group_counts={
                f"GS{int(key)}":
                    int(
                        value
                    )
                for key, value in gleason_counts.items()
            },
            mapping_verified=mapping_verified,
            minimum_group_support_met=gleason_group_support,
            metadata_complete=bool(
                len(
                    gleason
                )
                >
                0
            ),
            heterogeneous_definition=False,
            config=config,
        )
    )

    # ----------------------------------------------------------------------
    # Hormone refractory
    #
    # Only an explicit positive case is source-labelled.
    # Missing labels are NOT interpreted as negative.
    # ----------------------------------------------------------------------

    hormone_positive = int(
        sample_inventory[
            "hormone_refractory_explicit"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    rows.append(
        _endpoint_row(
            endpoint="hormone_refractory",
            endpoint_type="BINARY",
            usable_cases=hormone_positive,
            group_counts={
                "EXPLICIT_POSITIVE":
                    hormone_positive,
            },
            mapping_verified=mapping_verified,
            minimum_group_support_met=False,
            metadata_complete=False,
            heterogeneous_definition=False,
            config=config,
        )
    )

    # ----------------------------------------------------------------------
    # LN-PCa
    # ----------------------------------------------------------------------

    ln_positive = int(
        sample_inventory[
            "lymph_node_pca_explicit"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    rows.append(
        _endpoint_row(
            endpoint="lymph_node_pca",
            endpoint_type="BINARY",
            usable_cases=ln_positive,
            group_counts={
                "EXPLICIT_LN_PCA":
                    ln_positive,
            },
            mapping_verified=mapping_verified,
            minimum_group_support_met=False,
            metadata_complete=False,
            heterogeneous_definition=False,
            config=config,
        )
    )

    # ----------------------------------------------------------------------
    # TMPRSS2-ERG
    #
    # One explicitly labelled molecular subgroup is insufficient to define
    # positive and negative groups.
    # ----------------------------------------------------------------------

    tmprss2_positive = int(
        sample_inventory[
            "source_title"
        ]
        .str.contains(
            "TMPRSS2-ERG",
            case=False,
            regex=False,
        )
        .sum()
    )

    rows.append(
        _endpoint_row(
            endpoint="tmprss2_erg",
            endpoint_type="BINARY",
            usable_cases=tmprss2_positive,
            group_counts={
                "EXPLICIT_POSITIVE":
                    tmprss2_positive,
            },
            mapping_verified=mapping_verified,
            minimum_group_support_met=False,
            metadata_complete=False,
            heterogeneous_definition=False,
            config=config,
        )
    )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "M7.4D endpoint-readiness output is empty."
        )

    return result


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    sample_inventory: pd.DataFrame,
    endpoint_readiness: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.4D summary."""

    endpoints_assessed = int(
        len(
            endpoint_readiness
        )
    )

    endpoints_ready = int(
        endpoint_readiness[
            "association_ready"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    endpoints_not_ready = int(
        endpoints_assessed
        -
        endpoints_ready
    )

    mapping_verified = bool(
        config[
            "run_mapping_policy"
        ][
            "run_to_geo_mapping_verified"
        ]
    )

    ffpe_samples = int(
        sample_inventory[
            "ffpe"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    derived_samples = int(
        sample_inventory[
            "derived_from_group"
        ].notna().sum()
    )

    recurrence_explicit_cases = int(
        sample_inventory[
            "recurrence_explicit"
        ].notna().sum()
    )

    if endpoints_ready > 0:

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "ready"
            ]
        )

        next_stage = (
            config[
                "next_stage"
            ][
                "if_ready"
            ]
        )

    else:

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "descriptive_only"
            ]
        )

        next_stage = (
            config[
                "next_stage"
            ][
                "if_descriptive_only"
            ]
        )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.4D",

                "dataset_key":
                    config[
                        "dataset"
                    ][
                        "dataset_key"
                    ],

                "geo_accession":
                    config[
                        "dataset"
                    ][
                        "geo_accession"
                    ],

                "source_samples":
                    int(
                        len(
                            sample_inventory
                        )
                    ),

                "ffpe_samples":
                    ffpe_samples,

                "derived_or_nonindependent_samples":
                    derived_samples,

                "recurrence_status_explicit_cases":
                    recurrence_explicit_cases,

                "run_sample_mapping_verified":
                    mapping_verified,

                "candidate_abundance_variable":
                    config[
                        "molecular_input"
                    ][
                        "abundance_variable"
                    ],

                "candidate_abundance_interpretation":
                    config[
                        "molecular_input"
                    ][
                        "abundance_interpretation"
                    ],

                "endpoints_assessed":
                    endpoints_assessed,

                "endpoints_association_ready":
                    endpoints_ready,

                "endpoints_not_ready":
                    endpoints_not_ready,

                "association_model_fitted":
                    False,

                "hypothesis_test_performed":
                    False,

                "p_value_calculated":
                    False,

                "effect_size_estimated":
                    False,

                "survival_analysis_performed":
                    False,

                "recurrence_model_performed":
                    False,

                "clinical_validation_claimed":
                    False,

                "mature_24nt_trf_abundance_claimed":
                    False,

                "external_molecular_evidence_retained":
                    True,

                "overall_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Validation
# ============================================================================


def _validate_result(
    *,
    sample_inventory: pd.DataFrame,
    endpoint_readiness: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """Validate M7.4D safeguards."""

    if len(
        summary
    ) != 1:

        raise RuntimeError(
            "M7.4D summary must contain one row."
        )

    if len(
        sample_inventory
    ) != 11:

        raise RuntimeError(
            "M7.4D expected 11 GSE80400 source samples."
        )

    forbidden_endpoint_true = [
        "association_model_fitted",
        "p_value_calculated",
        "effect_size_estimated",
    ]

    for column in forbidden_endpoint_true:

        if (
            endpoint_readiness[
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
                "M7.4D safeguard violation: "
                f"{column}=True."
            )

    row = summary.iloc[
        0
    ]

    forbidden_summary_true = [
        "association_model_fitted",
        "hypothesis_test_performed",
        "p_value_calculated",
        "effect_size_estimated",
        "survival_analysis_performed",
        "recurrence_model_performed",
        "clinical_validation_claimed",
        "mature_24nt_trf_abundance_claimed",
    ]

    for field in forbidden_summary_true:

        if bool(
            row[
                field
            ]
        ):

            raise RuntimeError(
                "M7.4D summary safeguard violation: "
                f"{field}=True."
            )


# ============================================================================
# Public API
# ============================================================================


def assess_external_clinical_association_readiness(
    *,
    config: dict[str, Any],
) -> ExternalClinicalAssociationReadinessResult:
    """Execute M7.4D readiness assessment."""

    sample_inventory = build_sample_inventory(
        config=config,
    )

    endpoint_readiness = build_endpoint_readiness(
        sample_inventory=sample_inventory,
        config=config,
    )

    summary = build_summary(
        sample_inventory=sample_inventory,
        endpoint_readiness=endpoint_readiness,
        config=config,
    )

    _validate_result(
        sample_inventory=sample_inventory,
        endpoint_readiness=endpoint_readiness,
        summary=summary,
    )

    return ExternalClinicalAssociationReadinessResult(
        sample_inventory=sample_inventory,
        endpoint_readiness=endpoint_readiness,
        summary=summary,
    )