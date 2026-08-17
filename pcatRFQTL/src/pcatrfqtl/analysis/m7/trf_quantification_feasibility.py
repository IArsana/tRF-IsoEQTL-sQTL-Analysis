"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/trf_quantification_feasibility.py

Description:
    Core logic for M7.3B tRF Quantification Feasibility Audit.

    This stage determines whether source-supported candidate tRF sequences
    can, in principle, be quantified from TCGA-PRAD miRNA-Seq resources.

    Processed GDC miRNA and miRNA-isoform quantification files are not
    accepted as direct tRF measurements.

    Candidate-specific tRF quantification requires:
        - an exact source-supported candidate sequence;
        - successful sequence validation;
        - aligned miRNA-Seq reads;
        - a defensible case mapping;
        - appropriate access to sequencing files.

    Important provenance rule:
        API query scope and returned nested project metadata are represented
        separately.

        target_project_query_scoped=True:
            the GDC request itself was restricted to TCGA-PRAD.

        target_project_observed=True:
            TCGA-PRAD was explicitly present in returned nested
            cases.project.project_id metadata.

        Missing project metadata is NOT treated as project confirmation.

    M7.3B does NOT:
        - download controlled BAM files;
        - count reads;
        - realign sequencing reads;
        - quantify tRF expression;
        - infer unresolved tRF sequences;
        - substitute miRNAs for tRFs;
        - substitute another tRF candidate;
        - perform expression association;
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

from dataclasses import dataclass
from typing import Any

import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class TrfQuantificationFeasibilityResult:
    """Container for M7.3B outputs."""

    resource_inventory: pd.DataFrame
    candidate_feasibility: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Generic helpers
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


def _as_dict_list(
    value: Any,
) -> list[dict[str, Any]]:
    """Normalize arbitrary nested value to list of dictionaries."""

    if not isinstance(
        value,
        list,
    ):
        return []

    return [
        item
        for item in value
        if isinstance(
            item,
            dict,
        )
    ]


# ============================================================================
# Candidate extraction
# ============================================================================


def extract_candidate_sequences(
    *,
    m7_3a_qc: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract candidate sequence-resolution state from locked M7.3A.

    Candidate identity and sequence state are inherited directly from M7.3A.
    No sequence inference or candidate substitution occurs here.
    """

    records = m7_3a_qc.get(
        "candidate_trf_sequence_resolution",
        [],
    )

    if not isinstance(
        records,
        list,
    ):
        raise RuntimeError(
            "M7.3A candidate_trf_sequence_resolution must be a list."
        )

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        lead_rsid = _normalize_text(
            record.get(
                "lead_rsid"
            )
        )

        trf_id = _normalize_text(
            record.get(
                "trf_id"
            )
        )

        if lead_rsid is None:
            raise RuntimeError(
                "M7.3B encountered candidate without lead_rsid."
            )

        if trf_id is None:
            raise RuntimeError(
                f"M7.3B candidate {lead_rsid} has no tRF identifier."
            )

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "priority_rank":
                    record.get(
                        "priority_rank"
                    ),

                "priority_class":
                    _normalize_text(
                        record.get(
                            "priority_class"
                        )
                    ),

                "trf_id":
                    trf_id,

                "sequence_dna":
                    _normalize_text(
                        record.get(
                            "sequence_dna"
                        )
                    ),

                "sequence_resolved":
                    bool(
                        record.get(
                            "sequence_resolved",
                            False,
                        )
                    ),

                "sequence_validation_passed":
                    bool(
                        record.get(
                            "sequence_validation_passed",
                            False,
                        )
                    ),

                "candidate_quantification_ready_m7_3a":
                    bool(
                        record.get(
                            "candidate_quantification_ready",
                            False,
                        )
                    ),

                "sequence_resolution_status":
                    _normalize_text(
                        record.get(
                            "sequence_resolution_status"
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

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        raise RuntimeError(
            "No M7.3A candidates available for M7.3B."
        )

    if result[
        [
            "lead_rsid",
            "trf_id",
        ]
    ].duplicated().any():
        raise RuntimeError(
            "M7.3B candidate source contains duplicate lead_rsid/tRF pairs."
        )

    return result.sort_values(
        [
            "priority_rank",
            "lead_rsid",
            "trf_id",
        ],
        kind="stable",
        na_position="last",
    ).reset_index(
        drop=True
    )


# ============================================================================
# GDC resource inventory
# ============================================================================


def build_resource_inventory(
    *,
    files: list[dict[str, Any]],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Harmonize GDC TCGA-PRAD miRNA-Seq file metadata.

    One output row represents one GDC file.

    Project provenance semantics:
        target_project_query_scoped
            True because the upstream GDC API query was explicitly filtered
            to the configured target project.

        project_metadata_available
            True only when returned nested case.project.project_id metadata
            exists.

        target_project_observed
            True only when the configured target project is explicitly found
            in returned nested project metadata.

    Missing nested project metadata does NOT imply project confirmation.
    """

    target_project = str(
        config[
            "target_resource"
        ][
            "project_id"
        ]
    )

    target_strategy = str(
        config[
            "target_resource"
        ][
            "experimental_strategy"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for file_record in files:

        if not isinstance(
            file_record,
            dict,
        ):
            continue

        # ------------------------------------------------------------------
        # Nested case metadata
        # ------------------------------------------------------------------

        cases = _as_dict_list(
            file_record.get(
                "cases"
            )
        )

        case_ids: set[
            str
        ] = set()

        case_submitter_ids: set[
            str
        ] = set()

        sample_types: set[
            str
        ] = set()

        project_ids: set[
            str
        ] = set()

        for case in cases:

            # --------------------------------------------------------------
            # Case identity
            # --------------------------------------------------------------

            case_id = _normalize_text(
                case.get(
                    "case_id"
                )
            )

            if case_id is not None:
                case_ids.add(
                    case_id
                )

            case_submitter_id = _normalize_text(
                case.get(
                    "submitter_id"
                )
            )

            if case_submitter_id is not None:
                case_submitter_ids.add(
                    case_submitter_id
                )

            # --------------------------------------------------------------
            # Explicit returned project metadata
            # --------------------------------------------------------------

            project = case.get(
                "project"
            )

            if isinstance(
                project,
                dict,
            ):

                project_id = _normalize_text(
                    project.get(
                        "project_id"
                    )
                )

                if project_id is not None:
                    project_ids.add(
                        project_id
                    )

            # --------------------------------------------------------------
            # Nested sample metadata
            # --------------------------------------------------------------

            samples = _as_dict_list(
                case.get(
                    "samples"
                )
            )

            for sample in samples:

                sample_type = _normalize_text(
                    sample.get(
                        "sample_type"
                    )
                )

                if sample_type is not None:
                    sample_types.add(
                        sample_type
                    )

        # ------------------------------------------------------------------
        # File metadata
        # ------------------------------------------------------------------

        data_type = _normalize_text(
            file_record.get(
                "data_type"
            )
        )

        data_format = _normalize_text(
            file_record.get(
                "data_format"
            )
        )

        experimental_strategy = _normalize_text(
            file_record.get(
                "experimental_strategy"
            )
        )

        access = _normalize_text(
            file_record.get(
                "access"
            )
        )

        normalized_access = (
            access.lower()
            if access is not None
            else None
        )

        # ------------------------------------------------------------------
        # Resource classification
        # ------------------------------------------------------------------

        is_aligned_mirna_bam = bool(
            experimental_strategy
            ==
            target_strategy
            and
            data_type
            ==
            "Aligned Reads"
            and
            data_format
            ==
            "BAM"
        )

        is_processed_mirna = bool(
            experimental_strategy
            ==
            target_strategy
            and
            data_type
            ==
            "miRNA Expression Quantification"
        )

        is_processed_isoform = bool(
            experimental_strategy
            ==
            target_strategy
            and
            data_type
            ==
            "Isoform Expression Quantification"
        )

        if is_aligned_mirna_bam:

            resource_class = (
                "ALIGNED_MIRNA_READS"
            )

            accepted_for_direct_trf_quantification = True

        elif is_processed_mirna:

            resource_class = (
                "PROCESSED_MIRNA_EXPRESSION"
            )

            accepted_for_direct_trf_quantification = False

        elif is_processed_isoform:

            resource_class = (
                "PROCESSED_MIRNA_ISOFORM_EXPRESSION"
            )

            accepted_for_direct_trf_quantification = False

        else:

            resource_class = (
                "OTHER_MIRNA_RESOURCE"
            )

            accepted_for_direct_trf_quantification = False

        # ------------------------------------------------------------------
        # Access classification
        # ------------------------------------------------------------------

        controlled_access = bool(
            normalized_access
            ==
            "controlled"
        )

        open_access = bool(
            normalized_access
            ==
            "open"
        )

        # ------------------------------------------------------------------
        # Project provenance
        # ------------------------------------------------------------------

        project_metadata_available = bool(
            project_ids
        )

        target_project_observed = bool(
            target_project
            in project_ids
        )

        # Query scope is guaranteed by the upstream GDC filter.
        target_project_query_scoped = True

        # ------------------------------------------------------------------
        # Output row
        # ------------------------------------------------------------------

        rows.append(
            {
                "file_id":
                    _normalize_text(
                        file_record.get(
                            "file_id"
                        )
                    ),

                "file_name":
                    _normalize_text(
                        file_record.get(
                            "file_name"
                        )
                    ),

                "data_category":
                    _normalize_text(
                        file_record.get(
                            "data_category"
                        )
                    ),

                "data_type":
                    data_type,

                "data_format":
                    data_format,

                "experimental_strategy":
                    experimental_strategy,

                "access":
                    access,

                "file_size":
                    file_record.get(
                        "file_size"
                    ),

                "state":
                    _normalize_text(
                        file_record.get(
                            "state"
                        )
                    ),

                # ----------------------------------------------------------
                # Case/sample metadata
                # ----------------------------------------------------------

                "case_ids":
                    sorted(
                        case_ids
                    ),

                "case_submitter_ids":
                    sorted(
                        case_submitter_ids
                    ),

                "case_count":
                    int(
                        len(
                            case_ids
                        )
                    ),

                "sample_types":
                    sorted(
                        sample_types
                    ),

                # ----------------------------------------------------------
                # Project provenance
                # ----------------------------------------------------------

                "project_ids":
                    sorted(
                        project_ids
                    ),

                "project_metadata_available":
                    project_metadata_available,

                "target_project_observed":
                    target_project_observed,

                "target_project_query_scoped":
                    target_project_query_scoped,

                # ----------------------------------------------------------
                # Resource classification
                # ----------------------------------------------------------

                "resource_class":
                    resource_class,

                "accepted_for_direct_trf_quantification":
                    accepted_for_direct_trf_quantification,

                # ----------------------------------------------------------
                # Access
                # ----------------------------------------------------------

                "controlled_access":
                    controlled_access,

                "open_access":
                    open_access,

                # ----------------------------------------------------------
                # M7.3B safeguards
                # ----------------------------------------------------------

                "downloaded_in_m7_3b":
                    False,

                "read_counting_performed":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        return result

    # ----------------------------------------------------------------------
    # Defensive validation
    # ----------------------------------------------------------------------

    if (
        "file_id"
        in result.columns
        and
        result[
            "file_id"
        ].notna().any()
        and
        result[
            "file_id"
        ].dropna().duplicated().any()
    ):

        duplicated = (
            result.loc[
                result[
                    "file_id"
                ].duplicated(
                    keep=False
                ),
                "file_id",
            ]
            .dropna()
            .astype(
                str
            )
            .unique()
            .tolist()
        )

        raise RuntimeError(
            "M7.3B GDC resource inventory contains duplicate file IDs: "
            f"{duplicated[:10]}"
        )

    if not result[
        "target_project_query_scoped"
    ].all():

        raise RuntimeError(
            "M7.3B contains resource rows not marked as "
            "target-project query scoped."
        )

    invalid_project_state = (
        result[
            "target_project_observed"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        &
        ~result[
            "project_metadata_available"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    if invalid_project_state.any():

        raise RuntimeError(
            "M7.3B project provenance inconsistency: "
            "target_project_observed=True while "
            "project_metadata_available=False."
        )

    return result.reset_index(
        drop=True
    )


# ============================================================================
# Candidate feasibility
# ============================================================================


def build_candidate_feasibility(
    *,
    candidates: pd.DataFrame,
    resources: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Determine candidate-specific tRF quantification feasibility.

    Resource availability is evaluated independently from candidate sequence
    readiness.

    A candidate with an unresolved sequence remains blocked even when aligned
    miRNA-Seq BAM resources exist.
    """

    statuses = config[
        "statuses"
    ][
        "candidate"
    ]

    # ----------------------------------------------------------------------
    # Aligned miRNA BAM inventory
    # ----------------------------------------------------------------------

    if resources.empty:

        aligned = resources.copy()

    else:

        aligned = resources.loc[
            resources[
                "resource_class"
            ]
            ==
            "ALIGNED_MIRNA_READS"
        ].copy()

    aligned_file_count = int(
        len(
            aligned
        )
    )

    aligned_case_ids: set[
        str
    ] = set()

    if not aligned.empty:

        for values in aligned[
            "case_ids"
        ].tolist():

            if not isinstance(
                values,
                (
                    list,
                    tuple,
                    set,
                ),
            ):
                continue

            aligned_case_ids.update(
                str(
                    value
                )
                for value in values
            )

    controlled_aligned_count = int(
        aligned[
            "controlled_access"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
        if not aligned.empty
        else 0
    )

    open_aligned_count = int(
        aligned[
            "open_access"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
        if not aligned.empty
        else 0
    )

    rows: list[
        dict[str, Any]
    ] = []

    for record in candidates.to_dict(
        orient="records"
    ):

        sequence_resolved = bool(
            record.get(
                "sequence_resolved",
                False,
            )
        )

        sequence_validation_passed = bool(
            record.get(
                "sequence_validation_passed",
                False,
            )
        )

        upstream_quantification_ready = bool(
            record.get(
                "candidate_quantification_ready_m7_3a",
                False,
            )
        )

        sequence_ready = bool(
            sequence_resolved
            and
            sequence_validation_passed
            and
            upstream_quantification_ready
        )

        # ------------------------------------------------------------------
        # Sequence unresolved / unusable
        # ------------------------------------------------------------------

        if not sequence_ready:

            status = str(
                statuses[
                    "sequence_unresolved"
                ]
            )

            resource_available = bool(
                aligned_file_count
                >
                0
            )

            controlled_access_required = False
            direct_quantification_ready = False

        # ------------------------------------------------------------------
        # No aligned miRNA reads
        # ------------------------------------------------------------------

        elif aligned_file_count == 0:

            status = str(
                statuses[
                    "no_aligned_reads"
                ]
            )

            resource_available = False
            controlled_access_required = False
            direct_quantification_ready = False

        # ------------------------------------------------------------------
        # Open aligned BAM available
        # ------------------------------------------------------------------

        elif open_aligned_count > 0:

            status = str(
                statuses[
                    "open_access_ready"
                ]
            )

            resource_available = True
            controlled_access_required = False
            direct_quantification_ready = True

        # ------------------------------------------------------------------
        # Controlled-access aligned BAM route
        # ------------------------------------------------------------------

        elif controlled_aligned_count > 0:

            status = str(
                statuses[
                    "controlled_access_required"
                ]
            )

            resource_available = True
            controlled_access_required = True
            direct_quantification_ready = False

        # ------------------------------------------------------------------
        # Unresolved inventory state
        # ------------------------------------------------------------------

        else:

            status = str(
                statuses[
                    "resource_inventory_incomplete"
                ]
            )

            resource_available = bool(
                aligned_file_count
                >
                0
            )

            controlled_access_required = False
            direct_quantification_ready = False

        rows.append(
            {
                "lead_rsid":
                    record[
                        "lead_rsid"
                    ],

                "priority_rank":
                    record[
                        "priority_rank"
                    ],

                "priority_class":
                    record[
                        "priority_class"
                    ],

                "trf_id":
                    record[
                        "trf_id"
                    ],

                "sequence_dna":
                    record[
                        "sequence_dna"
                    ],

                "sequence_resolved":
                    sequence_resolved,

                "sequence_validation_passed":
                    sequence_validation_passed,

                "sequence_ready_for_quantification":
                    sequence_ready,

                # ----------------------------------------------------------
                # Resource state
                # ----------------------------------------------------------

                "aligned_mirna_bam_file_count":
                    aligned_file_count,

                "aligned_mirna_bam_case_count":
                    int(
                        len(
                            aligned_case_ids
                        )
                    ),

                "controlled_aligned_bam_count":
                    controlled_aligned_count,

                "open_aligned_bam_count":
                    open_aligned_count,

                "aligned_read_resource_available":
                    resource_available,

                # ----------------------------------------------------------
                # Access / execution state
                # ----------------------------------------------------------

                "controlled_access_required":
                    controlled_access_required,

                "direct_quantification_ready_now":
                    direct_quantification_ready,

                "processed_mirna_accepted_as_direct_measurement":
                    False,

                # ----------------------------------------------------------
                # Resolution
                # ----------------------------------------------------------

                "candidate_quantification_feasibility_status":
                    status,

                "candidate_retained":
                    bool(
                        record.get(
                            "candidate_retained",
                            True,
                        )
                    ),

                "quantification_performed":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        raise RuntimeError(
            "M7.3B produced no candidate feasibility records."
        )

    return result.sort_values(
        [
            "priority_rank",
            "lead_rsid",
            "trf_id",
        ],
        kind="stable",
        na_position="last",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    resources: pd.DataFrame,
    candidate_feasibility: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.3B summary."""

    if resources.empty:

        aligned = resources.copy()
        processed_mirna = resources.copy()
        processed_isoform = resources.copy()

    else:

        aligned = resources.loc[
            resources[
                "resource_class"
            ]
            ==
            "ALIGNED_MIRNA_READS"
        ]

        processed_mirna = resources.loc[
            resources[
                "resource_class"
            ]
            ==
            "PROCESSED_MIRNA_EXPRESSION"
        ]

        processed_isoform = resources.loc[
            resources[
                "resource_class"
            ]
            ==
            "PROCESSED_MIRNA_ISOFORM_EXPRESSION"
        ]

    # ----------------------------------------------------------------------
    # Candidate counts
    # ----------------------------------------------------------------------

    sequence_ready = int(
        candidate_feasibility[
            "sequence_ready_for_quantification"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    direct_ready = int(
        candidate_feasibility[
            "direct_quantification_ready_now"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    controlled_required = int(
        candidate_feasibility[
            "controlled_access_required"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    blocked_sequence = int(
        (
            ~candidate_feasibility[
                "sequence_ready_for_quantification"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
        ).sum()
    )

    # ----------------------------------------------------------------------
    # Resource/case counts
    # ----------------------------------------------------------------------

    aligned_case_ids: set[
        str
    ] = set()

    if not aligned.empty:

        for values in aligned[
            "case_ids"
        ].tolist():

            if not isinstance(
                values,
                (
                    list,
                    tuple,
                    set,
                ),
            ):
                continue

            aligned_case_ids.update(
                str(
                    value
                )
                for value in values
            )

    project_metadata_available_files = int(
        resources[
            "project_metadata_available"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
        if not resources.empty
        else 0
    )

    target_project_observed_files = int(
        resources[
            "target_project_observed"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
        if not resources.empty
        else 0
    )

    # ----------------------------------------------------------------------
    # Final state
    # ----------------------------------------------------------------------

    if direct_ready > 0:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "direct_ready"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_direct_ready"
            ]
        )

    elif controlled_required > 0:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "controlled_route"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_controlled_access_required"
            ]
        )

    else:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "blocked"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_all_blocked"
            ]
        )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.3B",

                # ----------------------------------------------------------
                # Candidates
                # ----------------------------------------------------------

                "candidate_trfs_assessed":
                    int(
                        len(
                            candidate_feasibility
                        )
                    ),

                "candidate_sequences_resolved":
                    sequence_ready,

                "candidates_blocked_by_sequence_resolution":
                    blocked_sequence,

                # ----------------------------------------------------------
                # Resource inventory
                # ----------------------------------------------------------

                "gdc_mirna_files_observed":
                    int(
                        len(
                            resources
                        )
                    ),

                "aligned_mirna_bam_files":
                    int(
                        len(
                            aligned
                        )
                    ),

                "aligned_mirna_bam_cases":
                    int(
                        len(
                            aligned_case_ids
                        )
                    ),

                "processed_mirna_expression_files":
                    int(
                        len(
                            processed_mirna
                        )
                    ),

                "processed_mirna_isoform_files":
                    int(
                        len(
                            processed_isoform
                        )
                    ),

                # ----------------------------------------------------------
                # Project provenance
                # ----------------------------------------------------------

                "project_query_scope_applied":
                    True,

                "files_with_project_metadata":
                    project_metadata_available_files,

                "files_with_target_project_observed":
                    target_project_observed_files,

                "missing_project_metadata_interpreted_as_confirmation":
                    False,

                # ----------------------------------------------------------
                # Quantification readiness
                # ----------------------------------------------------------

                "candidates_directly_quantifiable_now":
                    direct_ready,

                "candidates_requiring_controlled_access":
                    controlled_required,

                # ----------------------------------------------------------
                # Safeguards
                # ----------------------------------------------------------

                "processed_mirna_used_as_trf_measurement":
                    False,

                "bam_download_performed":
                    False,

                "read_counting_performed":
                    False,

                # ----------------------------------------------------------
                # Resolution
                # ----------------------------------------------------------

                "overall_feasibility_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def assess_trf_quantification_feasibility(
    *,
    m7_3a_qc: dict[str, Any],
    gdc_files: list[dict[str, Any]],
    config: dict[str, Any],
) -> TrfQuantificationFeasibilityResult:
    """Execute M7.3B tRF Quantification Feasibility Audit."""

    candidates = extract_candidate_sequences(
        m7_3a_qc=m7_3a_qc,
    )

    resources = build_resource_inventory(
        files=gdc_files,
        config=config,
    )

    candidate_feasibility = build_candidate_feasibility(
        candidates=candidates,
        resources=resources,
        config=config,
    )

    summary = build_summary(
        resources=resources,
        candidate_feasibility=candidate_feasibility,
        config=config,
    )

    return TrfQuantificationFeasibilityResult(
        resource_inventory=resources,
        candidate_feasibility=candidate_feasibility,
        summary=summary,
    )