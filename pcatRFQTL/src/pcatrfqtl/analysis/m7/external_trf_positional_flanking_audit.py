"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/external_trf_positional_flanking_audit.py

Description:
    Core analysis for M7.4C.2B External tRF Positional and Flanking Audit.

    The analysis examines all FASTQ read lengths before applying the
    historical 15-35 nt canonical quantification window.

    It determines whether the exact candidate sequence occurs as:
        - a complete 24-nt read;
        - a 5' boundary sequence;
        - a 3' boundary sequence;
        - an internal motif;
        - multiple occurrences in one read.

    Reads shorter than 15 nt and longer than 35 nt are retained for
    descriptive QC and exact-sequence auditing rather than silently removed.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import pandas as pd


@dataclass(frozen=True)
class PositionalFlankingAuditResult:
    """Container for M7.4C.2B outputs."""

    run_summary: pd.DataFrame
    length_distribution: pd.DataFrame
    motif_context: pd.DataFrame
    motif_context_aggregate: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# FASTQ
# ============================================================================


def iter_fastq_sequences(
    path: Path,
) -> Iterator[str]:
    """Yield validated sequence records from FASTQ."""

    if not path.exists():

        raise FileNotFoundError(
            f"FASTQ not found: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:

        record_number = 0

        while True:

            identifier = handle.readline()

            if identifier == "":
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            record_number += 1

            if (
                sequence == ""
                or
                plus == ""
                or
                quality == ""
            ):

                raise RuntimeError(
                    "Incomplete FASTQ record "
                    f"{record_number} in {path}"
                )

            if not identifier.startswith("@"):

                raise RuntimeError(
                    "Invalid FASTQ identifier at "
                    f"record {record_number}."
                )

            if not plus.startswith("+"):

                raise RuntimeError(
                    "Invalid FASTQ '+' record at "
                    f"record {record_number}."
                )

            sequence = sequence.strip().upper()
            quality = quality.strip()

            if len(sequence) != len(quality):

                raise RuntimeError(
                    "FASTQ sequence/quality length mismatch at "
                    f"record {record_number}."
                )

            yield sequence


# ============================================================================
# Length bins
# ============================================================================


def classify_length_bin(
    length_nt: int,
    *,
    bins: list[dict[str, Any]],
) -> str:
    """Assign one configured read-length bin."""

    for definition in bins:

        minimum = definition.get(
            "min_nt"
        )

        maximum = definition.get(
            "max_nt"
        )

        minimum_ok = (
            minimum is None
            or
            length_nt >= int(
                minimum
            )
        )

        maximum_ok = (
            maximum is None
            or
            length_nt <= int(
                maximum
            )
        )

        if minimum_ok and maximum_ok:

            return str(
                definition[
                    "name"
                ]
            )

    raise RuntimeError(
        f"No configured read-length bin for {length_nt} nt."
    )


# ============================================================================
# Motif positions
# ============================================================================


def find_exact_occurrences(
    sequence: str,
    motif: str,
) -> list[int]:
    """
    Find all exact motif positions.

    Overlapping occurrences are retained.
    """

    positions: list[int] = []

    start = 0

    while True:

        index = sequence.find(
            motif,
            start,
        )

        if index < 0:
            break

        positions.append(
            index
        )

        start = index + 1

    return positions


def classify_boundary(
    *,
    read_length: int,
    candidate_length: int,
    start_0based: int,
    occurrence_count: int,
    classes: dict[str, str],
) -> str:
    """Classify candidate location relative to read boundaries."""

    end_0based_exclusive = (
        start_0based
        +
        candidate_length
    )

    if (
        read_length
        ==
        candidate_length
        and
        start_0based
        ==
        0
    ):

        return classes[
            "full_length"
        ]

    if occurrence_count > 1:

        return classes[
            "multiple_occurrences"
        ]

    if start_0based == 0:

        return classes[
            "five_prime_boundary"
        ]

    if end_0based_exclusive == read_length:

        return classes[
            "three_prime_boundary"
        ]

    return classes[
        "internal"
    ]


# ============================================================================
# One run
# ============================================================================


def audit_run(
    *,
    run_accession: str,
    fastq_path: Path,
    config: dict[str, Any],
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    """Audit all reads in one run."""

    candidate = str(
        config[
            "candidate"
        ][
            "exact_sequence"
        ]
    ).upper()

    candidate_length = len(
        candidate
    )

    policy = config[
        "read_length_policy"
    ]

    canonical_min = int(
        policy[
            "canonical_min_nt"
        ]
    )

    canonical_max = int(
        policy[
            "canonical_max_nt"
        ]
    )

    length_bins = policy[
        "length_bins"
    ]

    maximum_flank = int(
        config[
            "sequence_audit"
        ][
            "maximum_flank_sequence_nt"
        ]
    )

    boundary_classes = config[
        "boundary_classes"
    ]

    # ----------------------------------------------------------------------
    # Counters
    # ----------------------------------------------------------------------

    total_reads = 0

    canonical_reads = 0
    reads_lt15 = 0
    reads_gt35 = 0

    motif_positive_reads = 0
    motif_occurrences = 0

    full_length_reads = 0
    five_prime_boundary_reads = 0
    three_prime_boundary_reads = 0
    internal_reads = 0
    multiple_occurrence_reads = 0

    motif_positive_lt15 = 0
    motif_positive_canonical = 0
    motif_positive_gt35 = 0

    raw_length_counts: Counter[int] = Counter()

    context_counter: Counter[
        tuple[Any, ...]
    ] = Counter()

    # ----------------------------------------------------------------------
    # Reads
    # ----------------------------------------------------------------------

    for sequence in iter_fastq_sequences(
        fastq_path
    ):

        total_reads += 1

        read_length = len(
            sequence
        )

        raw_length_counts[
            read_length
        ] += 1

        if read_length < canonical_min:

            reads_lt15 += 1

        elif read_length > canonical_max:

            reads_gt35 += 1

        else:

            canonical_reads += 1

        positions = find_exact_occurrences(
            sequence,
            candidate,
        )

        if not positions:
            continue

        motif_positive_reads += 1
        motif_occurrences += len(
            positions
        )

        if read_length < canonical_min:

            motif_positive_lt15 += 1

        elif read_length > canonical_max:

            motif_positive_gt35 += 1

        else:

            motif_positive_canonical += 1

        occurrence_count = len(
            positions
        )

        read_classes: set[str] = set()

        for start_0based in positions:

            end_exclusive = (
                start_0based
                +
                candidate_length
            )

            left_flank = sequence[
                :start_0based
            ]

            right_flank = sequence[
                end_exclusive:
            ]

            boundary_class = classify_boundary(
                read_length=read_length,
                candidate_length=candidate_length,
                start_0based=start_0based,
                occurrence_count=occurrence_count,
                classes=boundary_classes,
            )

            read_classes.add(
                boundary_class
            )

            left_display = (
                left_flank[
                    -maximum_flank:
                ]
                if left_flank
                else ""
            )

            right_display = (
                right_flank[
                    :maximum_flank
                ]
                if right_flank
                else ""
            )

            key = (
                read_length,
                start_0based,
                end_exclusive,
                boundary_class,
                len(
                    left_flank
                ),
                len(
                    right_flank
                ),
                left_display,
                right_display,
            )

            context_counter[
                key
            ] += 1

        # One read contributes once to read-level boundary metrics.

        if boundary_classes[
            "full_length"
        ] in read_classes:

            full_length_reads += 1

        elif boundary_classes[
            "multiple_occurrences"
        ] in read_classes:

            multiple_occurrence_reads += 1

        elif boundary_classes[
            "five_prime_boundary"
        ] in read_classes:

            five_prime_boundary_reads += 1

        elif boundary_classes[
            "three_prime_boundary"
        ] in read_classes:

            three_prime_boundary_reads += 1

        else:

            internal_reads += 1

    # ----------------------------------------------------------------------
    # Length table
    # ----------------------------------------------------------------------

    length_rows: list[
        dict[str, Any]
    ] = []

    for length_nt, count in sorted(
        raw_length_counts.items()
    ):

        length_rows.append(
            {
                "run_accession":
                    run_accession,

                "read_length_nt":
                    int(
                        length_nt
                    ),

                "length_bin":
                    classify_length_bin(
                        length_nt,
                        bins=length_bins,
                    ),

                "read_count":
                    int(
                        count
                    ),

                "read_fraction":
                    (
                        float(
                            count
                            /
                            total_reads
                        )
                        if total_reads
                        else 0.0
                    ),

                "within_canonical_15_35":
                    bool(
                        canonical_min
                        <=
                        length_nt
                        <=
                        canonical_max
                    ),
            }
        )

    # ----------------------------------------------------------------------
    # Motif context
    # ----------------------------------------------------------------------

    context_rows: list[
        dict[str, Any]
    ] = []

    for (
        read_length,
        start,
        end,
        boundary_class,
        left_length,
        right_length,
        left_sequence,
        right_sequence,
    ), count in sorted(
        context_counter.items()
    ):

        context_rows.append(
            {
                "run_accession":
                    run_accession,

                "read_length_nt":
                    int(
                        read_length
                    ),

                "candidate_start_0based":
                    int(
                        start
                    ),

                "candidate_end_0based_exclusive":
                    int(
                        end
                    ),

                "candidate_start_1based":
                    int(
                        start
                        +
                        1
                    ),

                "candidate_end_1based_inclusive":
                    int(
                        end
                    ),

                "boundary_class":
                    str(
                        boundary_class
                    ),

                "candidate_at_5prime_end":
                    bool(
                        start
                        ==
                        0
                    ),

                "candidate_at_3prime_end":
                    bool(
                        end
                        ==
                        read_length
                    ),

                "left_flank_length":
                    int(
                        left_length
                    ),

                "right_flank_length":
                    int(
                        right_length
                    ),

                "left_flank_sequence":
                    left_sequence,

                "right_flank_sequence":
                    right_sequence,

                "occurrence_count":
                    int(
                        count
                    ),
            }
        )

    run_summary = {
        "run_accession":
            run_accession,

        "total_reads":
            int(
                total_reads
            ),

        "reads_lt15":
            int(
                reads_lt15
            ),

        "reads_15_to_35":
            int(
                canonical_reads
            ),

        "reads_gt35":
            int(
                reads_gt35
            ),

        "fraction_lt15":
            (
                float(
                    reads_lt15
                    /
                    total_reads
                )
                if total_reads
                else 0.0
            ),

        "fraction_15_to_35":
            (
                float(
                    canonical_reads
                    /
                    total_reads
                )
                if total_reads
                else 0.0
            ),

        "fraction_gt35":
            (
                float(
                    reads_gt35
                    /
                    total_reads
                )
                if total_reads
                else 0.0
            ),

        "motif_positive_reads":
            int(
                motif_positive_reads
            ),

        "motif_occurrences":
            int(
                motif_occurrences
            ),

        "motif_positive_reads_lt15":
            int(
                motif_positive_lt15
            ),

        "motif_positive_reads_15_to_35":
            int(
                motif_positive_canonical
            ),

        "motif_positive_reads_gt35":
            int(
                motif_positive_gt35
            ),

        "full_length_24nt_reads":
            int(
                full_length_reads
            ),

        "five_prime_boundary_reads":
            int(
                five_prime_boundary_reads
            ),

        "three_prime_boundary_reads":
            int(
                three_prime_boundary_reads
            ),

        "internal_motif_reads":
            int(
                internal_reads
            ),

        "multiple_occurrence_reads":
            int(
                multiple_occurrence_reads
            ),

        "motif_cpm_all_reads":
            (
                float(
                    motif_positive_reads
                    /
                    total_reads
                    *
                    1_000_000
                )
                if total_reads
                else 0.0
            ),
    }

    return (
        run_summary,
        length_rows,
        context_rows,
    )


# ============================================================================
# Multi-run
# ============================================================================


def build_audit_tables(
    *,
    run_fastq_paths: dict[str, Path],
    config: dict[str, Any],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """Audit all configured runs."""

    run_rows: list[
        dict[str, Any]
    ] = []

    length_rows: list[
        dict[str, Any]
    ] = []

    context_rows: list[
        dict[str, Any]
    ] = []

    for run_accession, path in sorted(
        run_fastq_paths.items()
    ):

        (
            run_summary,
            lengths,
            contexts,
        ) = audit_run(
            run_accession=run_accession,
            fastq_path=path,
            config=config,
        )

        run_rows.append(
            run_summary
        )

        length_rows.extend(
            lengths
        )

        context_rows.extend(
            contexts
        )

    return (
        pd.DataFrame(
            run_rows
        ),
        pd.DataFrame(
            length_rows
        ),
        pd.DataFrame(
            context_rows
        ),
    )


# ============================================================================
# Context aggregation
# ============================================================================


def aggregate_contexts(
    motif_context: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate identical positional/flanking contexts across runs."""

    if motif_context.empty:

        return pd.DataFrame(
            columns=[
                "read_length_nt",
                "candidate_start_0based",
                "candidate_end_0based_exclusive",
                "boundary_class",
                "left_flank_length",
                "right_flank_length",
                "left_flank_sequence",
                "right_flank_sequence",
                "occurrence_count",
                "runs_observed",
            ]
        )

    grouping = [
        "read_length_nt",
        "candidate_start_0based",
        "candidate_end_0based_exclusive",
        "boundary_class",
        "left_flank_length",
        "right_flank_length",
        "left_flank_sequence",
        "right_flank_sequence",
    ]

    result = (
        motif_context
        .groupby(
            grouping,
            dropna=False,
            as_index=False,
        )
        .agg(
            occurrence_count=(
                "occurrence_count",
                "sum",
            ),
            runs_observed=(
                "run_accession",
                "nunique",
            ),
        )
    )

    return result.sort_values(
        [
            "occurrence_count",
            "runs_observed",
        ],
        ascending=[
            False,
            False,
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    run_summary: pd.DataFrame,
    motif_context: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build overall M7.4C.2B summary."""

    total_reads = int(
        run_summary[
            "total_reads"
        ].sum()
    )

    reads_lt15 = int(
        run_summary[
            "reads_lt15"
        ].sum()
    )

    reads_15_35 = int(
        run_summary[
            "reads_15_to_35"
        ].sum()
    )

    reads_gt35 = int(
        run_summary[
            "reads_gt35"
        ].sum()
    )

    motif_positive = int(
        run_summary[
            "motif_positive_reads"
        ].sum()
    )

    motif_occurrences = int(
        run_summary[
            "motif_occurrences"
        ].sum()
    )

    full_length = int(
        run_summary[
            "full_length_24nt_reads"
        ].sum()
    )

    five_prime = int(
        run_summary[
            "five_prime_boundary_reads"
        ].sum()
    )

    three_prime = int(
        run_summary[
            "three_prime_boundary_reads"
        ].sum()
    )

    internal = int(
        run_summary[
            "internal_motif_reads"
        ].sum()
    )

    multiple = int(
        run_summary[
            "multiple_occurrence_reads"
        ].sum()
    )

    motif_lt15 = int(
        run_summary[
            "motif_positive_reads_lt15"
        ].sum()
    )

    motif_15_35 = int(
        run_summary[
            "motif_positive_reads_15_to_35"
        ].sum()
    )

    motif_gt35 = int(
        run_summary[
            "motif_positive_reads_gt35"
        ].sum()
    )

    # ----------------------------------------------------------------------
    # Interpretation class
    # ----------------------------------------------------------------------

    if full_length > 0:

        interpretation_status = (
            config[
                "statuses"
            ][
                "direct_fragment_support"
            ]
        )

    elif (
        five_prime > 0
        or
        three_prime > 0
    ):

        interpretation_status = (
            config[
                "statuses"
            ][
                "boundary_compatible"
            ]
        )

    elif motif_positive > 0:

        interpretation_status = (
            config[
                "statuses"
            ][
                "motif_only"
            ]
        )

    else:

        interpretation_status = (
            config[
                "statuses"
            ][
                "not_observed"
            ]
        )

    unique_contexts = int(
        len(
            motif_context
        )
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.4C.2B",

                "dataset_key":
                    config[
                        "dataset"
                    ][
                        "dataset_key"
                    ],

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

                "runs_audited":
                    int(
                        len(
                            run_summary
                        )
                    ),

                "total_reads":
                    total_reads,

                "reads_lt15":
                    reads_lt15,

                "reads_15_to_35":
                    reads_15_35,

                "reads_gt35":
                    reads_gt35,

                "fraction_lt15":
                    (
                        float(
                            reads_lt15
                            /
                            total_reads
                        )
                        if total_reads
                        else 0.0
                    ),

                "fraction_15_to_35":
                    (
                        float(
                            reads_15_35
                            /
                            total_reads
                        )
                        if total_reads
                        else 0.0
                    ),

                "fraction_gt35":
                    (
                        float(
                            reads_gt35
                            /
                            total_reads
                        )
                        if total_reads
                        else 0.0
                    ),

                "motif_positive_reads_all_lengths":
                    motif_positive,

                "motif_occurrences_all_lengths":
                    motif_occurrences,

                "motif_positive_reads_lt15":
                    motif_lt15,

                "motif_positive_reads_15_to_35":
                    motif_15_35,

                "motif_positive_reads_gt35":
                    motif_gt35,

                "full_length_24nt_reads":
                    full_length,

                "five_prime_boundary_reads":
                    five_prime,

                "three_prime_boundary_reads":
                    three_prime,

                "internal_motif_reads":
                    internal,

                "multiple_occurrence_reads":
                    multiple,

                "unique_positional_flanking_contexts":
                    unique_contexts,

                "all_length_exact_search_performed":
                    True,

                "canonical_window_retained":
                    True,

                "artificial_read_extension_performed":
                    False,

                "artificial_read_truncation_performed":
                    False,

                "approximate_matching_performed":
                    False,

                "parent_trna_inference_performed":
                    False,

                "genomic_origin_inference_performed":
                    False,

                "clinical_association_performed":
                    False,

                "interpretation_status":
                    interpretation_status,

                "overall_status":
                    config[
                        "statuses"
                    ][
                        "overall_complete"
                    ],

                "next_stage":
                    config[
                        "next_stage"
                    ],
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def audit_external_trf_positional_flanking_context(
    *,
    run_fastq_paths: dict[str, Path],
    config: dict[str, Any],
) -> PositionalFlankingAuditResult:
    """Execute M7.4C.2B."""

    (
        run_summary,
        length_distribution,
        motif_context,
    ) = build_audit_tables(
        run_fastq_paths=run_fastq_paths,
        config=config,
    )

    context_aggregate = aggregate_contexts(
        motif_context
    )

    summary = build_summary(
        run_summary=run_summary,
        motif_context=motif_context,
        config=config,
    )

    return PositionalFlankingAuditResult(
        run_summary=run_summary,
        length_distribution=length_distribution,
        motif_context=motif_context,
        motif_context_aggregate=context_aggregate,
        summary=summary,
    )