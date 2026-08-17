"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/external_small_rna_technical_eligibility.py

Description:
    Core logic for M7.4B External Small-RNA Technical Eligibility.

    M7.4B determines whether public external small-RNA sequencing datasets
    are technically suitable for exact-sequence auditing of the prioritized
    tRF candidate.

    This stage evaluates metadata only.

    It does NOT:
        - download FASTQ/SRA reads;
        - search for the candidate sequence;
        - count candidate reads;
        - quantify tRF expression;
        - claim candidate presence or absence;
        - perform clinical association;
        - infer parent tRNA origin.

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
class ExternalSmallRnaTechnicalEligibilityResult:
    """Container for M7.4B output tables."""

    dataset_eligibility: pd.DataFrame
    summary: pd.DataFrame


def _as_bool(
    value: Any,
) -> bool:
    """Convert optional scalar to bool conservatively."""

    if value is None:
        return False

    if isinstance(
        value,
        bool,
    ):
        return value

    if isinstance(
        value,
        str,
    ):

        normalized = value.strip().lower()

        if normalized in {
            "true",
            "yes",
            "1",
        }:
            return True

        if normalized in {
            "false",
            "no",
            "0",
            "",
        }:
            return False

    return bool(value)


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional text."""

    if value is None:
        return None

    normalized = str(value).strip()

    return normalized or None


def _candidate_length_supported(
    *,
    audit: dict[str, Any],
    candidate_length_nt: int,
) -> bool | None:
    """
    Determine whether candidate length is technically compatible.

    Explicit source-supported compatibility takes precedence.
    """

    explicit = audit.get(
        "candidate_length_compatible"
    )

    if explicit is not None:
        return _as_bool(
            explicit
        )

    min_nt = audit.get(
        "documented_processed_read_min_nt"
    )

    max_nt = audit.get(
        "documented_processed_read_max_nt"
    )

    if min_nt is None or max_nt is None:
        return None

    return bool(
        int(min_nt)
        <=
        candidate_length_nt
        <=
        int(max_nt)
    )


def classify_dataset_technical_eligibility(
    *,
    dataset_key: str,
    audit: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Classify one dataset using source-supported technical metadata."""

    candidate_length_nt = int(
        config[
            "candidate"
        ][
            "sequence_length_nt"
        ]
    )

    criteria = config[
        "technical_criteria"
    ]

    statuses = config[
        "statuses"
    ][
        "dataset"
    ]

    raw_public = bool(
        _as_bool(
            audit.get(
                "raw_sequence_public"
            )
        )
        and
        (
            _normalize_text(
                audit.get(
                    "raw_sequence_repository"
                )
            )
            is not None
        )
    )

    library_strategy = (
        _normalize_text(
            audit.get(
                "library_strategy"
            )
        )
        or ""
    ).lower()

    small_rna_library = bool(
        "mirna-seq"
        in library_strategy
        or
        "ncrna-seq"
        in library_strategy
        or
        "small"
        in library_strategy
    )

    size_fractionated = (
        (
            _normalize_text(
                audit.get(
                    "library_selection"
                )
            )
            or ""
        ).lower()
        ==
        "size fractionation"
    )

    explicit_small_rna_protocol = (
        _normalize_text(
            audit.get(
                "small_rna_protocol"
            )
        )
        is not None
    )

    library_preparation_supported = bool(
        size_fractionated
        or
        explicit_small_rna_protocol
    )

    candidate_length_compatible = _candidate_length_supported(
        audit=audit,
        candidate_length_nt=candidate_length_nt,
    )

    read_orientation_resolvable = _as_bool(
        audit.get(
            "read_orientation_resolvable"
        )
    )

    adapter_information_available = _as_bool(
        audit.get(
            "adapter_information_available"
        )
    )

    # For exact-sequence auditing, an explicitly identified standard
    # small-RNA protocol may permit adapter resolution during M7.4C even
    # when the exact adapter sequence has not yet been recorded here.
    adapter_handling_resolvable = bool(
        adapter_information_available
        or
        explicit_small_rna_protocol
    )

    metadata_source_verified = _as_bool(
        audit.get(
            "metadata_source_verified"
        )
    )

    requirements: list[
        bool | None
    ] = []

    if _as_bool(
        criteria.get(
            "require_public_raw_sequence"
        )
    ):
        requirements.append(
            raw_public
        )

    if _as_bool(
        criteria.get(
            "require_small_rna_or_ncrna_library"
        )
    ):
        requirements.append(
            small_rna_library
        )

    if _as_bool(
        criteria.get(
            "require_size_fractionation_or_explicit_small_rna_protocol"
        )
    ):
        requirements.append(
            library_preparation_supported
        )

    if _as_bool(
        criteria.get(
            "require_candidate_length_compatible"
        )
    ):
        requirements.append(
            candidate_length_compatible
        )

    if _as_bool(
        criteria.get(
            "require_read_orientation_resolvable"
        )
    ):
        requirements.append(
            read_orientation_resolvable
        )

    if _as_bool(
        criteria.get(
            "require_adapter_handling_resolvable"
        )
    ):
        requirements.append(
            adapter_handling_resolvable
        )

    any_unresolved = any(
        value is None
        for value in requirements
    )

    any_failed = any(
        value is False
        for value in requirements
    )

    ffpe = _as_bool(
        audit.get(
            "ffpe"
        )
    )

    caveats = audit.get(
        "technical_caveats",
        [],
    )

    if not isinstance(
        caveats,
        list,
    ):
        caveats = []

    if not metadata_source_verified:

        eligible = False

        status = str(
            statuses[
                "incomplete"
            ]
        )

        resolution_class = (
            "METADATA_INCOMPLETE"
        )

    elif any_unresolved:

        eligible = False

        status = str(
            statuses[
                "incomplete"
            ]
        )

        resolution_class = (
            "METADATA_INCOMPLETE"
        )

    elif any_failed:

        eligible = False

        status = str(
            statuses[
                "not_eligible"
            ]
        )

        resolution_class = (
            "NOT_TECHNICALLY_ELIGIBLE"
        )

    elif ffpe or caveats:

        eligible = True

        status = str(
            statuses[
                "eligible_with_caveats"
            ]
        )

        resolution_class = (
            "TECHNICALLY_ELIGIBLE_WITH_CAVEATS"
        )

    else:

        eligible = True

        status = str(
            statuses[
                "eligible"
            ]
        )

        resolution_class = (
            "TECHNICALLY_ELIGIBLE"
        )

    return {
        "dataset_key":
            dataset_key,

        "raw_sequence_public":
            raw_public,

        "raw_sequence_repository":
            _normalize_text(
                audit.get(
                    "raw_sequence_repository"
                )
            ),

        "library_strategy":
            _normalize_text(
                audit.get(
                    "library_strategy"
                )
            ),

        "library_source":
            _normalize_text(
                audit.get(
                    "library_source"
                )
            ),

        "library_selection":
            _normalize_text(
                audit.get(
                    "library_selection"
                )
            ),

        "small_rna_protocol":
            _normalize_text(
                audit.get(
                    "small_rna_protocol"
                )
            ),

        "read_mode":
            _normalize_text(
                audit.get(
                    "read_mode"
                )
            ),

        "small_rna_library_supported":
            small_rna_library,

        "library_preparation_supported":
            library_preparation_supported,

        "candidate_length_nt":
            candidate_length_nt,

        "candidate_length_compatible":
            candidate_length_compatible,

        "documented_processed_read_min_nt":
            audit.get(
                "documented_processed_read_min_nt"
            ),

        "documented_processed_read_max_nt":
            audit.get(
                "documented_processed_read_max_nt"
            ),

        "adapter_information_available":
            adapter_information_available,

        "documented_adapter_sequence":
            _normalize_text(
                audit.get(
                    "documented_adapter_sequence"
                )
            ),

        "adapter_handling_resolvable":
            adapter_handling_resolvable,

        "read_orientation_resolvable":
            read_orientation_resolvable,

        "original_pipeline_trf_aware":
            _as_bool(
                audit.get(
                    "original_pipeline_trf_aware"
                )
            ),

        "ffpe":
            ffpe,

        "metadata_source_verified":
            metadata_source_verified,

        "technical_caveat_count":
            int(
                len(
                    caveats
                )
            ),

        "technical_caveats":
            [
                str(
                    value
                )
                for value in caveats
            ],

        "technical_eligibility":
            eligible,

        "resolution_class":
            resolution_class,

        "technical_status":
            status,

        # M7.4B safeguards
        "fastq_download_performed":
            False,

        "candidate_sequence_search_performed":
            False,

        "candidate_detected":
            False,

        "candidate_absent":
            False,

        "exact_match_counting_performed":
            False,

        "trf_quantification_performed":
            False,

        "clinical_association_performed":
            False,
    }


def build_dataset_technical_eligibility(
    *,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build dataset-level M7.4B technical eligibility table."""

    audits = config[
        "dataset_audit"
    ]

    rows = [
        classify_dataset_technical_eligibility(
            dataset_key=dataset_key,
            audit=audit,
            config=config,
        )
        for dataset_key, audit in audits.items()
    ]

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        raise RuntimeError(
            "M7.4B produced no dataset audit rows."
        )

    return result.sort_values(
        [
            "technical_eligibility",
            "resolution_class",
            "dataset_key",
        ],
        ascending=[
            False,
            True,
            True,
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


def build_summary(
    *,
    eligibility: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.4B summary."""

    eligible = eligibility.loc[
        eligibility[
            "technical_eligibility"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    ].copy()

    fully_eligible = int(
        (
            eligible[
                "resolution_class"
            ]
            ==
            "TECHNICALLY_ELIGIBLE"
        ).sum()
    )

    eligible_with_caveats = int(
        (
            eligible[
                "resolution_class"
            ]
            ==
            "TECHNICALLY_ELIGIBLE_WITH_CAVEATS"
        ).sum()
    )

    incomplete = int(
        (
            eligibility[
                "resolution_class"
            ]
            ==
            "METADATA_INCOMPLETE"
        ).sum()
    )

    not_eligible = int(
        (
            eligibility[
                "resolution_class"
            ]
            ==
            "NOT_TECHNICALLY_ELIGIBLE"
        ).sum()
    )

    if not eligible.empty:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "route_available"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_route_available"
            ]
        )

    elif incomplete > 0:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "only_incomplete"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_no_route"
            ]
        )

    else:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "no_route"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_no_route"
            ]
        )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.4B",

                "datasets_assessed":
                    int(
                        len(
                            eligibility
                        )
                    ),

                "datasets_technically_eligible":
                    int(
                        len(
                            eligible
                        )
                    ),

                "datasets_fully_eligible":
                    fully_eligible,

                "datasets_eligible_with_caveats":
                    eligible_with_caveats,

                "datasets_metadata_incomplete":
                    incomplete,

                "datasets_not_eligible":
                    not_eligible,

                "candidate_trf":
                    config[
                        "candidate"
                    ][
                        "trf_id"
                    ],

                "candidate_sequence":
                    config[
                        "candidate"
                    ][
                        "exact_sequence"
                    ],

                "candidate_length_nt":
                    int(
                        config[
                            "candidate"
                        ][
                            "sequence_length_nt"
                        ]
                    ),

                "fastq_download_performed":
                    False,

                "candidate_sequence_search_performed":
                    False,

                "candidate_detection_claimed":
                    False,

                "candidate_absence_claimed":
                    False,

                "trf_quantification_performed":
                    False,

                "clinical_association_performed":
                    False,

                "overall_technical_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


def assess_external_small_rna_technical_eligibility(
    *,
    config: dict[str, Any],
) -> ExternalSmallRnaTechnicalEligibilityResult:
    """Execute M7.4B."""

    eligibility = build_dataset_technical_eligibility(
        config=config,
    )

    summary = build_summary(
        eligibility=eligibility,
        config=config,
    )

    return ExternalSmallRnaTechnicalEligibilityResult(
        dataset_eligibility=eligibility,
        summary=summary,
    )