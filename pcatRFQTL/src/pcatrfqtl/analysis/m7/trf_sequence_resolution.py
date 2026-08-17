"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/trf_sequence_resolution.py

Description:
    Core logic for M7.3A Candidate tRF Sequence Resolution.

    This stage resolves candidate tRF nucleotide sequences only when an exact
    sequence is explicitly supported by curated source evidence.

    Important semantic rule:
        absence of an exact sequence does NOT mean that a sequence failed
        alphabet or length validation. Such validation fields remain
        not assessable until sequence data are actually available.

    M7.3A does NOT:
        - decode sequence from a tRF identifier;
        - infer sequence from identifier characters;
        - infer sequence from a related tRF;
        - infer sequence from a parent tRNA;
        - infer parent tRNA;
        - infer tRF biological class;
        - substitute similar tRFs;
        - quantify TCGA small-RNA reads;
        - perform expression analysis;
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
import re
from typing import Any

import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class TrfSequenceResolutionResult:
    """Container for M7.3A outputs."""

    sequence_evidence: pd.DataFrame
    candidate_resolution: pd.DataFrame
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


def _normalize_sequence(
    sequence: Any,
    *,
    normalize_u_to_t: bool,
) -> str | None:
    """
    Normalize optional nucleotide sequence.

    Sequence normalization is formatting only and does not infer any missing
    nucleotide information.
    """

    normalized = _normalize_text(
        sequence
    )

    if normalized is None:
        return None

    normalized = (
        normalized
        .upper()
        .replace(
            " ",
            "",
        )
        .replace(
            "-",
            "",
        )
    )

    if normalize_u_to_t:

        normalized = normalized.replace(
            "U",
            "T",
        )

    return normalized or None


def _extract_identifier_length(
    trf_id: str,
    *,
    pattern: str,
) -> int | None:
    """
    Extract length metadata encoded in the identifier prefix.

    This is used only as a validation checksum after a source-supported
    sequence is available.

    It is NOT used to reconstruct or infer nucleotide sequence.
    """

    match = re.search(
        pattern,
        trf_id,
    )

    if match is None:
        return None

    try:

        return int(
            match.group(
                1
            )
        )

    except (
        TypeError,
        ValueError,
    ):

        return None


# ============================================================================
# Candidate extraction
# ============================================================================


def extract_candidates(
    *,
    m7_1_qc: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract candidate/tRF identities from locked M7.1.

    Candidate and tRF identity are inherited from upstream artifacts.
    """

    records = m7_1_qc.get(
        "candidate_validation_readiness",
        [],
    )

    if not isinstance(
        records,
        list,
    ):

        raise RuntimeError(
            "M7.1 candidate_validation_readiness must be a list."
        )

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        if not bool(
            record.get(
                "candidate_retained_for_m7",
                True,
            )
        ):

            continue

        lead_rsid = _normalize_text(
            record.get(
                "lead_rsid"
            )
        )

        if lead_rsid is None:

            raise RuntimeError(
                "M7.3A encountered a retained candidate without lead_rsid."
            )

        trf_ids_raw = record.get(
            "trf_ids",
            [],
        )

        if not isinstance(
            trf_ids_raw,
            list,
        ):

            trf_ids_raw = []

        trf_ids = sorted(
            {
                normalized
                for normalized in (
                    _normalize_text(
                        value
                    )
                    for value in trf_ids_raw
                )
                if normalized is not None
            }
        )

        if not trf_ids:

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
                        None,
                }
            )

            continue

        for trf_id in trf_ids:

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
                }
            )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No retained M7 candidate tRF identities are available."
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
# Curated evidence
# ============================================================================


def build_sequence_evidence(
    *,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Build curated candidate tRF sequence evidence.

    Validation behavior:
        - sequence available:
            alphabet and length consistency are evaluated;
        - sequence unavailable:
            alphabet and length validation are None / not assessable.
    """

    records = config.get(
        "sequence_evidence",
        [],
    )

    if not isinstance(
        records,
        list,
    ):

        raise RuntimeError(
            "sequence_evidence must be a YAML list."
        )

    validation = config[
        "sequence_validation"
    ]

    alphabet = {
        str(
            value
        ).upper()
        for value in validation[
            "alphabet"
        ]
    }

    normalize_u_to_t = bool(
        validation[
            "normalize_rna_u_to_dna_t"
        ]
    )

    require_identifier_length_match = bool(
        validation[
            "require_identifier_length_match"
        ]
    )

    length_pattern = str(
        validation[
            "length_prefix_pattern"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for record in records:

        trf_id = _normalize_text(
            record.get(
                "trf_id"
            )
        )

        if trf_id is None:

            raise RuntimeError(
                "M7.3A sequence evidence contains a missing tRF ID."
            )

        sequence = _normalize_sequence(
            record.get(
                "sequence_dna"
            ),
            normalize_u_to_t=normalize_u_to_t,
        )

        reported_length_raw = record.get(
            "sequence_length"
        )

        reported_length = (
            int(
                reported_length_raw
            )
            if reported_length_raw is not None
            else None
        )

        observed_length = (
            len(
                sequence
            )
            if sequence is not None
            else None
        )

        identifier_expected_length = _extract_identifier_length(
            trf_id,
            pattern=length_pattern,
        )

        exact_sequence_reported = bool(
            record.get(
                "exact_sequence_reported",
                False,
            )
        )

        exact_identifier_observed = bool(
            record.get(
                "exact_identifier_observed",
                False,
            )
        )

        # ------------------------------------------------------------------
        # Semantic validation
        #
        # Missing sequence = not assessable, not False.
        # ------------------------------------------------------------------

        if sequence is None:

            sequence_alphabet_valid: bool | None = None
            reported_length_consistent: bool | None = None
            identifier_length_consistent: bool | None = None
            sequence_validation_passed = False

        else:

            sequence_alphabet_valid = bool(
                all(
                    nucleotide in alphabet
                    for nucleotide in sequence
                )
            )

            reported_length_consistent = bool(
                reported_length is None
                or
                reported_length
                ==
                observed_length
            )

            if require_identifier_length_match:

                identifier_length_consistent = bool(
                    identifier_expected_length is not None
                    and
                    identifier_expected_length
                    ==
                    observed_length
                )

            else:

                identifier_length_consistent = True

            sequence_validation_passed = bool(
                exact_identifier_observed
                and
                exact_sequence_reported
                and
                sequence_alphabet_valid is True
                and
                reported_length_consistent is True
                and
                identifier_length_consistent is True
            )

        # ------------------------------------------------------------------
        # Internal source-evidence consistency
        # ------------------------------------------------------------------

        if exact_sequence_reported and sequence is None:

            raise RuntimeError(
                f"{trf_id}: exact_sequence_reported=True but "
                "sequence_dna is missing."
            )

        if (
            not exact_sequence_reported
            and
            sequence is not None
        ):

            raise RuntimeError(
                f"{trf_id}: sequence_dna is present but "
                "exact_sequence_reported=False."
            )

        rows.append(
            {
                "trf_id":
                    trf_id,

                "exact_identifier_observed":
                    exact_identifier_observed,

                "exact_sequence_reported":
                    exact_sequence_reported,

                "sequence_dna":
                    sequence,

                "reported_sequence_length":
                    reported_length,

                "observed_sequence_length":
                    observed_length,

                "identifier_expected_length":
                    identifier_expected_length,

                "sequence_alphabet_valid":
                    sequence_alphabet_valid,

                "reported_length_consistent":
                    reported_length_consistent,

                "identifier_length_consistent":
                    identifier_length_consistent,

                "sequence_validation_passed":
                    sequence_validation_passed,

                "sequence_validation_assessable":
                    bool(
                        sequence is not None
                    ),

                "source_type":
                    _normalize_text(
                        record.get(
                            "source_type"
                        )
                    ),

                "source_title":
                    _normalize_text(
                        record.get(
                            "source_title"
                        )
                    ),

                "source_locator":
                    _normalize_text(
                        record.get(
                            "source_locator"
                        )
                    ),

                "source_database_context":
                    _normalize_text(
                        record.get(
                            "source_database_context"
                        )
                    ),

                "evidence_status":
                    _normalize_text(
                        record.get(
                            "evidence_status"
                        )
                    ),

                "identifier_decoded":
                    False,

                "sequence_inferred":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No curated tRF sequence evidence configured."
        )

    if result[
        "trf_id"
    ].duplicated().any():

        duplicated = (
            result.loc[
                result[
                    "trf_id"
                ].duplicated(
                    keep=False
                ),
                "trf_id",
            ]
            .astype(
                str
            )
            .unique()
            .tolist()
        )

        raise RuntimeError(
            "Multiple curated sequence evidence records exist "
            f"for the same tRF ID: {duplicated}"
        )

    return result.sort_values(
        "trf_id",
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Candidate sequence resolution
# ============================================================================


def resolve_candidate_sequences(
    *,
    candidates: pd.DataFrame,
    evidence: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Resolve source-supported exact sequences for candidate tRFs.
    """

    statuses = config[
        "statuses"
    ]

    result = candidates.merge(
        evidence,
        on="trf_id",
        how="left",
        validate="many_to_one",
    )

    sequence_resolved: list[
        bool
    ] = []

    candidate_quantification_ready: list[
        bool
    ] = []

    resolution_status: list[
        str
    ] = []

    for row in result.to_dict(
        orient="records"
    ):

        exact_identifier_observed = bool(
            row.get(
                "exact_identifier_observed",
                False,
            )
        )

        exact_sequence_reported = bool(
            row.get(
                "exact_sequence_reported",
                False,
            )
        )

        sequence_validation_passed = bool(
            row.get(
                "sequence_validation_passed",
                False,
            )
        )

        sequence = _normalize_text(
            row.get(
                "sequence_dna"
            )
        )

        # ------------------------------------------------------------------
        # Exact source-supported resolution
        # ------------------------------------------------------------------

        if (
            exact_identifier_observed
            and
            exact_sequence_reported
            and
            sequence is not None
            and
            sequence_validation_passed
        ):

            sequence_resolved.append(
                True
            )

            candidate_quantification_ready.append(
                True
            )

            resolution_status.append(
                str(
                    statuses[
                        "resolved"
                    ]
                )
            )

            continue

        # ------------------------------------------------------------------
        # Source sequence exists but validation failed
        # ------------------------------------------------------------------

        if (
            exact_sequence_reported
            and
            sequence is not None
            and
            not sequence_validation_passed
        ):

            sequence_resolved.append(
                False
            )

            candidate_quantification_ready.append(
                False
            )

            resolution_status.append(
                str(
                    statuses[
                        "invalid_sequence"
                    ]
                )
            )

            continue

        # ------------------------------------------------------------------
        # Identifier confirmed but sequence unavailable
        # ------------------------------------------------------------------

        if exact_identifier_observed:

            sequence_resolved.append(
                False
            )

            candidate_quantification_ready.append(
                False
            )

            resolution_status.append(
                str(
                    statuses[
                        "unresolved_identifier_confirmed"
                    ]
                )
            )

            continue

        # ------------------------------------------------------------------
        # Identifier itself not externally confirmed
        # ------------------------------------------------------------------

        sequence_resolved.append(
            False
        )

        candidate_quantification_ready.append(
            False
        )

        resolution_status.append(
            str(
                statuses[
                    "unresolved_identifier_unconfirmed"
                ]
            )
        )

    result[
        "sequence_resolved"
    ] = sequence_resolved

    result[
        "candidate_quantification_ready"
    ] = candidate_quantification_ready

    result[
        "sequence_resolution_status"
    ] = resolution_status

    result[
        "parent_trna_inferred"
    ] = False

    result[
        "trf_class_inferred"
    ] = False

    result[
        "sequence_substitution_performed"
    ] = False

    result[
        "candidate_retained"
    ] = True

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
    candidate_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.3A summary."""

    total = int(
        len(
            candidate_resolution
        )
    )

    identifier_confirmed = int(
        candidate_resolution[
            "exact_identifier_observed"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    resolved = int(
        candidate_resolution[
            "sequence_resolved"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    unresolved = (
        total
        -
        resolved
    )

    quantification_ready = int(
        candidate_resolution[
            "candidate_quantification_ready"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    validation_assessable = int(
        candidate_resolution[
            "sequence_validation_assessable"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    if total > 0 and resolved == total:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall_complete"
            ]
        )

    elif resolved > 0:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall_partial"
            ]
        )

    else:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall_none"
            ]
        )

    if resolved > 0:

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_any_sequence_resolved"
            ]
        )

    else:

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_none_resolved"
            ]
        )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.3A",

                "candidate_trfs_assessed":
                    total,

                "candidate_identifiers_source_confirmed":
                    identifier_confirmed,

                "candidate_sequences_with_assessable_validation":
                    validation_assessable,

                "candidate_sequences_resolved":
                    resolved,

                "candidate_sequences_unresolved":
                    unresolved,

                "candidates_ready_for_quantification_audit":
                    quantification_ready,

                "identifier_decoding_performed":
                    False,

                "sequence_inference_performed":
                    False,

                "sequence_substitution_performed":
                    False,

                "unresolved_sequence_interpreted_as_invalid":
                    False,

                "tcga_quantification_performed":
                    False,

                "overall_resolution_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def resolve_trf_sequences(
    *,
    m7_1_qc: dict[str, Any],
    config: dict[str, Any],
) -> TrfSequenceResolutionResult:
    """
    Execute M7.3A Candidate tRF Sequence Resolution.
    """

    candidates = extract_candidates(
        m7_1_qc=m7_1_qc,
    )

    evidence = build_sequence_evidence(
        config=config,
    )

    candidate_resolution = resolve_candidate_sequences(
        candidates=candidates,
        evidence=evidence,
        config=config,
    )

    summary = build_summary(
        candidate_resolution=candidate_resolution,
        config=config,
    )

    return TrfSequenceResolutionResult(
        sequence_evidence=evidence,
        candidate_resolution=candidate_resolution,
        summary=summary,
    )