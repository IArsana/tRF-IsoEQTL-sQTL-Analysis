"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/alternative_qtl_audit.py

Description:
    Core logic for M5.6C alternative PRAD QTL resource audit.

    The audit determines whether an alternative prostate regulatory-QTL
    resource can resolve the colocalization input gap identified in M5.5
    and M5.6A/B.

    Formal colocalization readiness requires more than the presence of
    significant QTL associations. The resource must provide sufficiently
    complete feature-specific regional association statistics.

    Resource classifications:
        FORMAL_COLOC_READY
        ABF_COLOC_READY_SUSIE_NOT_READY
        ALTERNATIVE_RESOURCE_REQUIRES_DATASET_INSPECTION
        SIGNIFICANT_ONLY_NOT_COLOC_READY
        INSUFFICIENT_METADATA
        NOT_PRAD_RELEVANT

    Scientific safeguards:
        - Searchable/downloadable does not imply dense summary statistics.
        - Significant QTL associations are not treated as the full test set.
        - NOT_VERIFIED is not interpreted as unavailable.
        - Missing variants are not interpreted as null associations.
        - No colocalization is performed.
        - No fine-mapping is performed.
        - No causal inference is performed.

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


CONFIRMED = "CONFIRMED"
LIKELY_AVAILABLE = "LIKELY_AVAILABLE"
PARTIAL = "PARTIAL"
NOT_VERIFIED = "NOT_VERIFIED"
UNAVAILABLE = "UNAVAILABLE"


FORMAL_COLOC_READY = (
    "FORMAL_COLOC_READY"
)

ABF_COLOC_READY_SUSIE_NOT_READY = (
    "ABF_COLOC_READY_SUSIE_NOT_READY"
)

ALTERNATIVE_RESOURCE_REQUIRES_DATASET_INSPECTION = (
    "ALTERNATIVE_RESOURCE_REQUIRES_DATASET_INSPECTION"
)

SIGNIFICANT_ONLY_NOT_COLOC_READY = (
    "SIGNIFICANT_ONLY_NOT_COLOC_READY"
)

INSUFFICIENT_METADATA = (
    "INSUFFICIENT_METADATA"
)

NOT_PRAD_RELEVANT = (
    "NOT_PRAD_RELEVANT"
)


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class AlternativeQTLAuditResult:
    """Container for M5.6C outputs."""

    resources: pd.DataFrame

    requirements: pd.DataFrame

    candidates: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _normalize_status(
    value: Any,
    *,
    allowed: set[str],
) -> str:
    """Normalize one audit status."""

    status = (
        str(value)
        .strip()
        .upper()
    )

    if status not in allowed:

        raise ValueError(
            f"Invalid alternative-QTL status: {value}. "
            f"Allowed={sorted(allowed)}"
        )

    return status


def _is_confirmed(
    status: str,
) -> bool:
    """Return whether a field is formally confirmed."""

    return (
        status
        ==
        CONFIRMED
    )


# ============================================================================
# Requirements table
# ============================================================================


def build_requirement_table(
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build formal coloc requirement table."""

    requirements = config[
        "formal_coloc_requirements"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for requirement, properties in requirements.items():

        if not isinstance(
            properties,
            dict,
        ):

            raise ValueError(
                "Formal coloc requirement must be a mapping: "
                f"{requirement}"
            )

        rows.append(
            {
                "requirement":
                    str(
                        requirement
                    ),

                "required":
                    bool(
                        properties.get(
                            "required",
                            False,
                        )
                    ),

                "required_for_susie_only":
                    bool(
                        properties.get(
                            "required_for_susie_only",
                            False,
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Resource audit
# ============================================================================


def audit_alternative_resources(
    config: dict[str, Any],
) -> pd.DataFrame:
    """Audit configured alternative regulatory-QTL resources."""

    policy = config[
        "audit_policy"
    ]

    allowed = {
        str(value)
        .strip()
        .upper()
        for value
        in policy[
            "allowed_statuses"
        ]
    }

    resources = config[
        "resources"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for resource_id, resource in resources.items():

        evidence = resource[
            "evidence_status"
        ]

        def status_of(
            field: str,
        ) -> str:

            if field not in evidence:

                return NOT_VERIFIED

            value = evidence[
                field
            ]

            if isinstance(
                value,
                dict,
            ):

                value = value.get(
                    "status",
                    NOT_VERIFIED,
                )

            return _normalize_status(
                value,
                allowed=allowed,
            )

        prad_status = status_of(
            "prad_available"
        )

        feature_status = status_of(
            "feature_specific_qtl"
        )

        variant_status = status_of(
            "canonical_variant_identity"
        )

        p_status = status_of(
            "p_value"
        )

        beta_status = status_of(
            "effect_estimate"
        )

        se_status = status_of(
            "standard_error"
        )

        af_status = status_of(
            "allele_frequency"
        )

        sample_size_status = status_of(
            "sample_size"
        )

        build_status = status_of(
            "genome_build"
        )

        tested_matrix_status = status_of(
            "full_tested_variant_feature_matrix"
        )

        nonsignificant_status = status_of(
            "nonsignificant_associations"
        )

        dense_status = status_of(
            "dense_locus_summary_statistics"
        )

        downloadable_status = status_of(
            "downloadable"
        )

        ld_status = status_of(
            "ld_matrix"
        )

        core_confirmed = bool(
            _is_confirmed(
                prad_status
            )
            and
            _is_confirmed(
                feature_status
            )
            and
            _is_confirmed(
                variant_status
            )
            and
            _is_confirmed(
                p_status
            )
            and
            _is_confirmed(
                beta_status
            )
            and
            _is_confirmed(
                se_status
            )
            and
            _is_confirmed(
                sample_size_status
            )
            and
            _is_confirmed(
                build_status
            )
            and
            _is_confirmed(
                tested_matrix_status
            )
            and
            _is_confirmed(
                nonsignificant_status
            )
            and
            _is_confirmed(
                dense_status
            )
            and
            _is_confirmed(
                downloadable_status
            )
        )

        susie_ready = bool(
            core_confirmed
            and
            _is_confirmed(
                ld_status
            )
        )

        if prad_status != CONFIRMED:

            classification = (
                NOT_PRAD_RELEVANT
            )

        elif (
            tested_matrix_status
            ==
            UNAVAILABLE
            or
            nonsignificant_status
            ==
            UNAVAILABLE
            or
            dense_status
            ==
            UNAVAILABLE
        ):

            classification = (
                SIGNIFICANT_ONLY_NOT_COLOC_READY
            )

        elif core_confirmed and susie_ready:

            classification = (
                FORMAL_COLOC_READY
            )

        elif core_confirmed:

            classification = (
                ABF_COLOC_READY_SUSIE_NOT_READY
            )

        elif (
            tested_matrix_status
            ==
            NOT_VERIFIED
            or
            nonsignificant_status
            ==
            NOT_VERIFIED
            or
            dense_status
            ==
            NOT_VERIFIED
        ):

            classification = (
                ALTERNATIVE_RESOURCE_REQUIRES_DATASET_INSPECTION
            )

        else:

            classification = (
                INSUFFICIENT_METADATA
            )

        rows.append(
            {
                "resource_id":
                    str(
                        resource_id
                    ),

                "display_name":
                    str(
                        resource[
                            "display_name"
                        ]
                    ),

                "prad_status":
                    prad_status,

                "feature_specific_qtl_status":
                    feature_status,

                "variant_identity_status":
                    variant_status,

                "p_value_status":
                    p_status,

                "effect_estimate_status":
                    beta_status,

                "standard_error_status":
                    se_status,

                "allele_frequency_status":
                    af_status,

                "sample_size_status":
                    sample_size_status,

                "genome_build_status":
                    build_status,

                "full_tested_variant_feature_matrix_status":
                    tested_matrix_status,

                "nonsignificant_associations_status":
                    nonsignificant_status,

                "dense_locus_summary_statistics_status":
                    dense_status,

                "downloadable_status":
                    downloadable_status,

                "ld_matrix_status":
                    ld_status,

                "abf_formal_coloc_ready":
                    bool(
                        core_confirmed
                    ),

                "susie_formal_coloc_ready":
                    susie_ready,

                "classification":
                    classification,

                "requires_dataset_inspection":
                    bool(
                        classification
                        ==
                        ALTERNATIVE_RESOURCE_REQUIRES_DATASET_INSPECTION
                    ),

                "can_resolve_moradi_gap_now":
                    bool(
                        core_confirmed
                    ),
            }
        )

    dataframe = pd.DataFrame(
        rows
    )

    if dataframe.empty:

        raise ValueError(
            "No alternative QTL resources configured."
        )

    return dataframe


# ============================================================================
# Candidate-level propagation
# ============================================================================


def build_candidate_alternative_qtl_status(
    *,
    resources: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Propagate alternative-resource audit to candidate loci."""

    if resources.empty:

        raise ValueError(
            "Alternative resource audit is empty."
        )

    any_abf_ready = bool(
        resources[
            "abf_formal_coloc_ready"
        ]
        .fillna(
            False
        )
        .sum()
        >
        0
    )

    any_susie_ready = bool(
        resources[
            "susie_formal_coloc_ready"
        ]
        .fillna(
            False
        )
        .sum()
        >
        0
    )

    inspection_required = bool(
        resources[
            "requires_dataset_inspection"
        ]
        .fillna(
            False
        )
        .sum()
        >
        0
    )

    rows: list[
        dict[str, Any]
    ] = []

    for lead in config[
        "candidate_leads"
    ]:

        if any_abf_ready:

            status = (
                "ALTERNATIVE_QTL_RESOLUTION_AVAILABLE"
            )

            action = (
                "HARMONIZE_ALTERNATIVE_QTL_AND_REASSESS_COLOC"
            )

        elif inspection_required:

            status = (
                "ALTERNATIVE_QTL_RESOLUTION_PENDING_DATASET_INSPECTION"
            )

            action = (
                "INSPECT_DOWNLOAD_SCHEMA_AND_ASSOCIATION_COVERAGE"
            )

        else:

            status = (
                "NO_ALTERNATIVE_FORMAL_COLOC_RESOURCE"
            )

            action = (
                "RETAIN_LOCUS_LEVEL_REGULATORY_EVIDENCE"
            )

        rows.append(
            {
                "lead_rsid":
                    str(
                        lead
                    ),

                "alternative_abf_resource_available":
                    any_abf_ready,

                "alternative_susie_resource_available":
                    any_susie_ready,

                "dataset_inspection_required":
                    inspection_required,

                "formal_coloc_ready":
                    any_abf_ready,

                "coloc_abf_ready":
                    any_abf_ready,

                "coloc_susie_ready":
                    any_susie_ready,

                "resolution_status":
                    status,

                "recommended_action":
                    action,

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


def audit_alternative_qtl_resources(
    config: dict[str, Any],
) -> AlternativeQTLAuditResult:
    """Execute M5.6C resource audit."""

    requirements = (
        build_requirement_table(
            config
        )
    )

    resources = (
        audit_alternative_resources(
            config
        )
    )

    candidates = (
        build_candidate_alternative_qtl_status(
            resources=resources,
            config=config,
        )
    )

    return AlternativeQTLAuditResult(
        resources=resources,
        requirements=requirements,
        candidates=candidates,
    )