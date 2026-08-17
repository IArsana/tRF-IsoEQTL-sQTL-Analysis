"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/controlled_access_limitation.py

Description:
    Core logic for M7.3C Controlled-Access Limitation Resolution.

    M7.3B established that TCGA-PRAD contains aligned miRNA-Seq BAM
    resources suitable in principle for candidate-specific tRF
    quantification, but those resources require controlled access.

    M7.3C records the current analysis-time state in which the required
    project authorization is not available.

    This is treated as an access limitation, not as:
        - absence of the sequencing resource;
        - negative tRF expression;
        - failed candidate validation;
        - evidence against the candidate.

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
class ControlledAccessLimitationResult:
    """Container for M7.3C outputs."""

    candidate_resolution: pd.DataFrame
    summary: pd.DataFrame


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

    value = str(value).strip()

    return value or None


def extract_candidate_state(
    *,
    m7_3b_qc: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve candidate-level outcome after controlled-access limitation."""

    records = m7_3b_qc.get(
        "candidate_trf_quantification_feasibility",
        [],
    )

    if not isinstance(records, list):
        raise RuntimeError(
            "M7.3B candidate feasibility records must be a list."
        )

    statuses = config[
        "statuses"
    ]

    rows: list[dict[str, Any]] = []

    for record in records:

        if not isinstance(record, dict):
            continue

        sequence_ready = bool(
            record.get(
                "sequence_ready_for_quantification",
                False,
            )
        )

        controlled_required = bool(
            record.get(
                "controlled_access_required",
                False,
            )
        )

        if not sequence_ready:

            final_status = str(
                statuses[
                    "sequence_unresolved"
                ]
            )

            tcga_direct_quantification_deferred = False
            sequence_limited = True
            access_limited = False

        elif controlled_required:

            final_status = str(
                statuses[
                    "candidate_access_blocked"
                ]
            )

            tcga_direct_quantification_deferred = True
            sequence_limited = False
            access_limited = True

        else:

            raise RuntimeError(
                "M7.3C encountered a sequence-ready candidate without "
                "the expected controlled-access requirement."
            )

        rows.append(
            {
                "lead_rsid":
                    _normalize_text(
                        record.get(
                            "lead_rsid"
                        )
                    ),

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
                    _normalize_text(
                        record.get(
                            "trf_id"
                        )
                    ),

                "sequence_dna":
                    _normalize_text(
                        record.get(
                            "sequence_dna"
                        )
                    ),

                "sequence_ready_for_quantification":
                    sequence_ready,

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

                "resource_technically_available":
                    bool(
                        record.get(
                            "aligned_read_resource_available",
                            False,
                        )
                    ),

                "sequence_limited":
                    sequence_limited,

                "access_limited":
                    access_limited,

                "authorized_access_available":
                    False,

                "tcga_direct_quantification_deferred":
                    tcga_direct_quantification_deferred,

                "tcga_quantification_performed":
                    False,

                "tcga_expression_negative":
                    False,

                "candidate_retained":
                    bool(
                        record.get(
                            "candidate_retained",
                            True,
                        )
                    ),

                "final_status":
                    final_status,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        raise RuntimeError(
            "M7.3C produced no candidate records."
        )

    return result.sort_values(
        [
            "priority_rank",
            "lead_rsid",
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


def build_summary(
    *,
    candidate_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.3C summary."""

    access_limited = int(
        candidate_resolution[
            "access_limited"
        ]
        .astype(bool)
        .sum()
    )

    sequence_limited = int(
        candidate_resolution[
            "sequence_limited"
        ]
        .astype(bool)
        .sum()
    )

    deferred = int(
        candidate_resolution[
            "tcga_direct_quantification_deferred"
        ]
        .astype(bool)
        .sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.3C",

                "candidate_trfs_assessed":
                    int(
                        len(
                            candidate_resolution
                        )
                    ),

                "candidates_access_limited":
                    access_limited,

                "candidates_sequence_limited":
                    sequence_limited,

                "candidates_tcga_quantification_deferred":
                    deferred,

                "controlled_access_authorization_available":
                    False,

                "tcga_direct_trf_quantification_performed":
                    False,

                "processed_mirna_substitution_performed":
                    False,

                "clinical_association_performed":
                    False,

                "clinical_association_deferred":
                    bool(
                        deferred
                        >
                        0
                    ),

                "access_limitation_interpreted_as_negative_evidence":
                    False,

                "overall_resolution_status":
                    str(
                        config[
                            "statuses"
                        ][
                            "overall"
                        ]
                    ),

                "next_external_validation_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "external_validation"
                        ]
                    ),

                "next_variant_validation_stage":
                    str(
                        config[
                            "next_stage"
                        ][
                            "variant_validation"
                        ]
                    ),

                "clinical_route":
                    str(
                        config[
                            "next_stage"
                        ][
                            "clinical_route"
                        ]
                    ),
            }
        ]
    )


def resolve_controlled_access_limitation(
    *,
    m7_3b_qc: dict[str, Any],
    config: dict[str, Any],
) -> ControlledAccessLimitationResult:
    """Execute M7.3C."""

    candidate_resolution = extract_candidate_state(
        m7_3b_qc=m7_3b_qc,
        config=config,
    )

    summary = build_summary(
        candidate_resolution=candidate_resolution,
        config=config,
    )

    return ControlledAccessLimitationResult(
        candidate_resolution=candidate_resolution,
        summary=summary,
    )