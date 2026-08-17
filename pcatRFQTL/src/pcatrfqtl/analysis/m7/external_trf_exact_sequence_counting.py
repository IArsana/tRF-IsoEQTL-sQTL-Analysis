"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/external_trf_exact_sequence_counting.py

Description:
    Core logic for M7.4C.2 Exact Candidate tRF Sequence Counting.

    This module performs conservative exact-sequence counting from FASTQ
    reads for the source-supported candidate tRF sequence.

    Candidate-positive reads are separated into:
        - exact full-length candidate reads;
        - longer usable reads containing the exact candidate sequence.

    This allows exact sequence-level evidence to be retained without
    conflating full-length candidate reads with exact substring evidence.

    The module additionally records preprocessing QC including:
        - total reads;
        - adapter-trimmed reads;
        - empty reads after trimming;
        - reads outside the allowed small-RNA length range;
        - usable reads;
        - usable-read fraction;
        - adapter-trimming fraction;
        - per-run candidate CPM.

    This module does NOT:
        - perform approximate sequence matching;
        - allow mismatches or indels;
        - infer parent tRNA identity;
        - infer genomic origin;
        - quantify unresolved candidate tRFs;
        - use processed miRNA abundance as a tRF surrogate;
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
from pathlib import Path
from typing import Any, Iterator

import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class ExactSequenceCountingResult:
    """Container for M7.4C.2 outputs."""

    run_counts: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# FASTQ reader
# ============================================================================


def iter_fastq_sequences(
    path: Path,
) -> Iterator[str]:
    """
    Yield sequence strings from an uncompressed FASTQ file.

    FASTQ records are validated structurally:
        line 1: identifier beginning with "@"
        line 2: sequence
        line 3: "+" line
        line 4: quality

    Sequence and quality lengths must match.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"FASTQ file not found: {path}"
        )

    if not path.is_file():

        raise RuntimeError(
            f"FASTQ path is not a file: {path}"
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
                    "Invalid FASTQ identifier at record "
                    f"{record_number} in {path}"
                )

            if not plus.startswith("+"):

                raise RuntimeError(
                    "Invalid FASTQ '+' line at record "
                    f"{record_number} in {path}"
                )

            sequence_value = sequence.strip()
            quality_value = quality.strip()

            if len(
                sequence_value
            ) != len(
                quality_value
            ):

                raise RuntimeError(
                    "FASTQ sequence/quality length mismatch at "
                    f"record {record_number} in {path}"
                )

            yield sequence_value


# ============================================================================
# Adapter trimming
# ============================================================================


def trim_suffix_adapter(
    sequence: str,
    *,
    adapter: str,
    minimum_overlap_nt: int,
) -> tuple[str, bool, int]:
    """
    Trim an exact 3' adapter or exact partial adapter suffix.

    No adapter mismatches are allowed.

    The function searches adapter prefixes from the maximum possible overlap
    down to minimum_overlap_nt.

    Returns
    -------
    tuple
        trimmed_sequence,
        adapter_trimmed,
        adapter_overlap_nt
    """

    if minimum_overlap_nt <= 0:

        raise ValueError(
            "minimum_overlap_nt must be > 0."
        )

    if not adapter:

        return (
            sequence,
            False,
            0,
        )

    maximum_overlap = min(
        len(
            adapter
        ),
        len(
            sequence
        ),
    )

    for overlap in range(
        maximum_overlap,
        minimum_overlap_nt - 1,
        -1,
    ):

        adapter_prefix = adapter[
            :overlap
        ]

        if sequence.endswith(
            adapter_prefix
        ):

            return (
                sequence[
                    :-overlap
                ],
                True,
                overlap,
            )

    return (
        sequence,
        False,
        0,
    )


# ============================================================================
# One-run counting
# ============================================================================


def count_candidate_in_fastq(
    *,
    fastq_path: Path,
    candidate_sequence: str,
    adapter_sequence: str,
    minimum_adapter_overlap_nt: int,
    minimum_read_length: int,
    maximum_read_length: int,
) -> dict[str, int | float | bool | str]:
    """
    Count exact candidate-sequence-positive reads in one FASTQ.

    A read can contribute at most one candidate-positive count.

    Candidate-positive reads are classified as:

        EXACT_FULL_LENGTH
            The complete usable read equals the candidate sequence.

        EXACT_SUBSTRING
            A longer usable read contains the candidate sequence exactly.

    Neither class permits mismatches or indels.
    """

    candidate = candidate_sequence.upper().strip()

    if not candidate:

        raise ValueError(
            "Candidate sequence must not be empty."
        )

    if minimum_read_length <= 0:

        raise ValueError(
            "minimum_read_length must be > 0."
        )

    if maximum_read_length < minimum_read_length:

        raise ValueError(
            "maximum_read_length must be >= minimum_read_length."
        )

    # ----------------------------------------------------------------------
    # Counters
    # ----------------------------------------------------------------------

    total_reads = 0

    reads_adapter_trimmed = 0
    total_adapter_bases_trimmed = 0

    empty_after_trim = 0

    reads_shorter_than_min = 0
    reads_longer_than_max = 0

    usable_reads = 0

    exact_positive_reads = 0
    exact_full_length_reads = 0
    exact_substring_reads = 0

    # ----------------------------------------------------------------------
    # Iterate FASTQ
    # ----------------------------------------------------------------------

    for raw_sequence in iter_fastq_sequences(
        fastq_path
    ):

        total_reads += 1

        sequence = raw_sequence.upper()

        (
            sequence,
            adapter_trimmed,
            adapter_overlap_nt,
        ) = trim_suffix_adapter(
            sequence,
            adapter=adapter_sequence,
            minimum_overlap_nt=minimum_adapter_overlap_nt,
        )

        if adapter_trimmed:

            reads_adapter_trimmed += 1

            total_adapter_bases_trimmed += int(
                adapter_overlap_nt
            )

        if not sequence:

            empty_after_trim += 1
            continue

        sequence_length = len(
            sequence
        )

        if sequence_length < minimum_read_length:

            reads_shorter_than_min += 1
            continue

        if sequence_length > maximum_read_length:

            reads_longer_than_max += 1
            continue

        usable_reads += 1

        # ------------------------------------------------------------------
        # Exact full-length candidate read
        # ------------------------------------------------------------------

        if sequence == candidate:

            exact_positive_reads += 1
            exact_full_length_reads += 1

            continue

        # ------------------------------------------------------------------
        # Exact candidate substring
        #
        # Still exact sequence evidence; no mismatches or indels.
        # ------------------------------------------------------------------

        if candidate in sequence:

            exact_positive_reads += 1
            exact_substring_reads += 1

    # ----------------------------------------------------------------------
    # Fractions / normalized metrics
    # ----------------------------------------------------------------------

    if total_reads > 0:

        adapter_trimmed_fraction = (
            reads_adapter_trimmed
            /
            total_reads
        )

        usable_read_fraction = (
            usable_reads
            /
            total_reads
        )

        filtered_read_fraction = (
            (
                empty_after_trim
                +
                reads_shorter_than_min
                +
                reads_longer_than_max
            )
            /
            total_reads
        )

    else:

        adapter_trimmed_fraction = 0.0
        usable_read_fraction = 0.0
        filtered_read_fraction = 0.0

    if usable_reads > 0:

        candidate_cpm = (
            exact_positive_reads
            /
            usable_reads
            *
            1_000_000
        )

        candidate_fraction = (
            exact_positive_reads
            /
            usable_reads
        )

        full_length_candidate_cpm = (
            exact_full_length_reads
            /
            usable_reads
            *
            1_000_000
        )

        substring_candidate_cpm = (
            exact_substring_reads
            /
            usable_reads
            *
            1_000_000
        )

    else:

        candidate_cpm = 0.0
        candidate_fraction = 0.0
        full_length_candidate_cpm = 0.0
        substring_candidate_cpm = 0.0

    # ----------------------------------------------------------------------
    # Internal consistency
    # ----------------------------------------------------------------------

    if (
        exact_full_length_reads
        +
        exact_substring_reads
        !=
        exact_positive_reads
    ):

        raise RuntimeError(
            "Exact candidate counting consistency failure."
        )

    expected_total = (
        usable_reads
        +
        empty_after_trim
        +
        reads_shorter_than_min
        +
        reads_longer_than_max
    )

    if expected_total != total_reads:

        raise RuntimeError(
            "FASTQ preprocessing accounting failure: "
            f"total={total_reads}, accounted={expected_total}."
        )

    return {
        "fastq_path":
            str(
                fastq_path
            ),

        # ------------------------------------------------------------------
        # Read QC
        # ------------------------------------------------------------------

        "total_reads":
            int(
                total_reads
            ),

        "reads_adapter_trimmed":
            int(
                reads_adapter_trimmed
            ),

        "adapter_trimmed_fraction":
            float(
                adapter_trimmed_fraction
            ),

        "total_adapter_bases_trimmed":
            int(
                total_adapter_bases_trimmed
            ),

        "empty_after_trim":
            int(
                empty_after_trim
            ),

        "reads_shorter_than_min":
            int(
                reads_shorter_than_min
            ),

        "reads_longer_than_max":
            int(
                reads_longer_than_max
            ),

        "usable_reads":
            int(
                usable_reads
            ),

        "usable_read_fraction":
            float(
                usable_read_fraction
            ),

        "filtered_read_fraction":
            float(
                filtered_read_fraction
            ),

        # ------------------------------------------------------------------
        # Candidate counts
        # ------------------------------------------------------------------

        "exact_candidate_positive_reads":
            int(
                exact_positive_reads
            ),

        "exact_full_length_candidate_reads":
            int(
                exact_full_length_reads
            ),

        "exact_candidate_substring_reads":
            int(
                exact_substring_reads
            ),

        # ------------------------------------------------------------------
        # Candidate normalization
        # ------------------------------------------------------------------

        "candidate_cpm":
            float(
                candidate_cpm
            ),

        "full_length_candidate_cpm":
            float(
                full_length_candidate_cpm
            ),

        "substring_candidate_cpm":
            float(
                substring_candidate_cpm
            ),

        "candidate_fraction":
            float(
                candidate_fraction
            ),

        # ------------------------------------------------------------------
        # Observation
        # ------------------------------------------------------------------

        "candidate_exact_match_observed":
            bool(
                exact_positive_reads
                >
                0
            ),
    }


# ============================================================================
# Multi-run counting
# ============================================================================


def build_run_counts(
    *,
    run_fastq_paths: dict[str, Path],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Count the candidate sequence across selected external runs."""

    if not run_fastq_paths:

        raise RuntimeError(
            "No FASTQ files supplied for M7.4C.2."
        )

    candidate_sequence = str(
        config[
            "candidate"
        ][
            "exact_sequence"
        ]
    )

    preprocessing = config[
        "preprocessing"
    ]

    counting = config[
        "counting"
    ]

    adapter_config = preprocessing[
        "adapter_trimming"
    ]

    adapter_enabled = bool(
        adapter_config[
            "enabled"
        ]
    )

    if adapter_enabled:

        adapter_sequence = str(
            adapter_config[
                "adapter_sequence"
            ]
        )

        minimum_adapter_overlap_nt = int(
            adapter_config[
                "minimum_adapter_overlap_nt"
            ]
        )

    else:

        adapter_sequence = ""

        minimum_adapter_overlap_nt = 1

    minimum_read_length = int(
        counting[
            "minimum_read_length_after_trim"
        ]
    )

    maximum_read_length = int(
        counting[
            "maximum_read_length_after_trim"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for (
        run_accession,
        fastq_path,
    ) in sorted(
        run_fastq_paths.items()
    ):

        result = count_candidate_in_fastq(
            fastq_path=fastq_path,
            candidate_sequence=candidate_sequence,
            adapter_sequence=adapter_sequence,
            minimum_adapter_overlap_nt=minimum_adapter_overlap_nt,
            minimum_read_length=minimum_read_length,
            maximum_read_length=maximum_read_length,
        )

        observed = bool(
            result[
                "candidate_exact_match_observed"
            ]
        )

        status = (
            config[
                "statuses"
            ][
                "run"
            ][
                "observed"
            ]
            if observed
            else
            config[
                "statuses"
            ][
                "run"
            ][
                "not_observed"
            ]
        )

        rows.append(
            {
                "run_accession":
                    run_accession,

                "candidate_trf":
                    config[
                        "candidate"
                    ][
                        "trf_id"
                    ],

                "candidate_sequence":
                    candidate_sequence,

                "candidate_length_nt":
                    int(
                        config[
                            "candidate"
                        ][
                            "sequence_length_nt"
                        ]
                    ),

                **result,

                "run_status":
                    status,

                # ----------------------------------------------------------
                # Safeguards
                # ----------------------------------------------------------

                "approximate_matching_performed":
                    False,

                "mismatches_allowed":
                    False,

                "indels_allowed":
                    False,

                "reverse_complement_search_performed":
                    False,

                "parent_trna_inference_performed":
                    False,

                "genomic_origin_inference_performed":
                    False,

                "clinical_association_performed":
                    False,

                "causal_inference_performed":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "M7.4C.2 produced an empty run-count table."
        )

    if result[
        "run_accession"
    ].duplicated().any():

        duplicates = (
            result.loc[
                result[
                    "run_accession"
                ].duplicated(
                    keep=False
                ),
                "run_accession",
            ]
            .astype(
                str
            )
            .unique()
            .tolist()
        )

        raise RuntimeError(
            "Duplicate run accessions in M7.4C.2 output: "
            f"{duplicates}"
        )

    return result.sort_values(
        "run_accession",
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    run_counts: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.4C.2 aggregate summary."""

    mode = str(
        config[
            "execution"
        ][
            "mode"
        ]
    )

    runs_processed = int(
        len(
            run_counts
        )
    )

    runs_with_exact_match = int(
        run_counts[
            "candidate_exact_match_observed"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    runs_without_exact_match = int(
        runs_processed
        -
        runs_with_exact_match
    )

    # ----------------------------------------------------------------------
    # Aggregate preprocessing metrics
    # ----------------------------------------------------------------------

    total_reads = int(
        run_counts[
            "total_reads"
        ].sum()
    )

    reads_adapter_trimmed = int(
        run_counts[
            "reads_adapter_trimmed"
        ].sum()
    )

    total_adapter_bases_trimmed = int(
        run_counts[
            "total_adapter_bases_trimmed"
        ].sum()
    )

    empty_after_trim = int(
        run_counts[
            "empty_after_trim"
        ].sum()
    )

    reads_shorter_than_min = int(
        run_counts[
            "reads_shorter_than_min"
        ].sum()
    )

    reads_longer_than_max = int(
        run_counts[
            "reads_longer_than_max"
        ].sum()
    )

    usable_reads = int(
        run_counts[
            "usable_reads"
        ].sum()
    )

    # ----------------------------------------------------------------------
    # Candidate metrics
    # ----------------------------------------------------------------------

    exact_positive_reads = int(
        run_counts[
            "exact_candidate_positive_reads"
        ].sum()
    )

    exact_full_length_reads = int(
        run_counts[
            "exact_full_length_candidate_reads"
        ].sum()
    )

    exact_substring_reads = int(
        run_counts[
            "exact_candidate_substring_reads"
        ].sum()
    )

    # ----------------------------------------------------------------------
    # Aggregate fractions
    # ----------------------------------------------------------------------

    if total_reads > 0:

        aggregate_adapter_trimmed_fraction = (
            reads_adapter_trimmed
            /
            total_reads
        )

        aggregate_usable_read_fraction = (
            usable_reads
            /
            total_reads
        )

    else:

        aggregate_adapter_trimmed_fraction = 0.0
        aggregate_usable_read_fraction = 0.0

    if usable_reads > 0:

        combined_candidate_cpm = (
            exact_positive_reads
            /
            usable_reads
            *
            1_000_000
        )

        combined_full_length_candidate_cpm = (
            exact_full_length_reads
            /
            usable_reads
            *
            1_000_000
        )

        combined_substring_candidate_cpm = (
            exact_substring_reads
            /
            usable_reads
            *
            1_000_000
        )

        combined_candidate_fraction = (
            exact_positive_reads
            /
            usable_reads
        )

    else:

        combined_candidate_cpm = 0.0
        combined_full_length_candidate_cpm = 0.0
        combined_substring_candidate_cpm = 0.0
        combined_candidate_fraction = 0.0

    # ----------------------------------------------------------------------
    # Run-level distribution
    # ----------------------------------------------------------------------

    run_cpm_min = float(
        run_counts[
            "candidate_cpm"
        ].min()
    )

    run_cpm_max = float(
        run_counts[
            "candidate_cpm"
        ].max()
    )

    run_cpm_median = float(
        run_counts[
            "candidate_cpm"
        ].median()
    )

    run_cpm_mean = float(
        run_counts[
            "candidate_cpm"
        ].mean()
    )

    # ----------------------------------------------------------------------
    # Observation status
    # ----------------------------------------------------------------------

    candidate_observed = bool(
        exact_positive_reads
        >
        0
    )

    # ----------------------------------------------------------------------
    # Determine whether selected run set is complete
    # ----------------------------------------------------------------------

    if mode == "full":

        expected_runs = {
            str(
                accession
            )
            for accession in config[
                "execution"
            ][
                "full_run_accessions"
            ]
        }

        observed_runs = set(
            run_counts[
                "run_accession"
            ].astype(
                str
            )
        )

        full_run_set_complete = bool(
            observed_runs
            ==
            expected_runs
        )

    else:

        expected_runs = {
            str(
                accession
            )
            for accession in config[
                "execution"
            ][
                "pilot_run_accessions"
            ]
        }

        observed_runs = set(
            run_counts[
                "run_accession"
            ].astype(
                str
            )
        )

        full_run_set_complete = False

    selected_run_set_complete = bool(
        observed_runs
        ==
        expected_runs
    )

    # ----------------------------------------------------------------------
    # Overall status
    # ----------------------------------------------------------------------

    if not selected_run_set_complete:

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "partial"
            ]
        )

        next_stage = (
            config[
                "next_stage"
            ][
                "if_failed"
            ]
        )

    elif candidate_observed:

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "observed"
            ]
        )

        if mode == "pilot":

            next_stage = (
                config[
                    "next_stage"
                ][
                    "if_pilot_observed"
                ]
            )

        else:

            next_stage = (
                config[
                    "next_stage"
                ][
                    "if_full_complete"
                ]
            )

    else:

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "not_observed"
            ]
        )

        if mode == "pilot":

            next_stage = (
                config[
                    "next_stage"
                ][
                    "if_pilot_not_observed"
                ]
            )

        else:

            next_stage = (
                config[
                    "next_stage"
                ][
                    "if_full_complete"
                ]
            )

    # ----------------------------------------------------------------------
    # Internal consistency
    # ----------------------------------------------------------------------

    if (
        exact_full_length_reads
        +
        exact_substring_reads
        !=
        exact_positive_reads
    ):

        raise RuntimeError(
            "Aggregate candidate count consistency failure."
        )

    accounted_reads = (
        usable_reads
        +
        empty_after_trim
        +
        reads_shorter_than_min
        +
        reads_longer_than_max
    )

    if accounted_reads != total_reads:

        raise RuntimeError(
            "Aggregate preprocessing accounting failure: "
            f"total={total_reads}, accounted={accounted_reads}."
        )

    return pd.DataFrame(
        [
            {
                # ----------------------------------------------------------
                # Identity
                # ----------------------------------------------------------

                "milestone":
                    "M7.4C.2",

                "execution_mode":
                    mode,

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

                "sra_study":
                    config[
                        "dataset"
                    ][
                        "sra_study"
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

                "candidate_length_nt":
                    int(
                        config[
                            "candidate"
                        ][
                            "sequence_length_nt"
                        ]
                    ),

                # ----------------------------------------------------------
                # Run completion
                # ----------------------------------------------------------

                "runs_processed":
                    runs_processed,

                "selected_run_set_complete":
                    selected_run_set_complete,

                "full_run_set_complete":
                    full_run_set_complete,

                "runs_with_exact_match":
                    runs_with_exact_match,

                "runs_without_exact_match":
                    runs_without_exact_match,

                # ----------------------------------------------------------
                # Read QC
                # ----------------------------------------------------------

                "total_reads":
                    total_reads,

                "reads_adapter_trimmed":
                    reads_adapter_trimmed,

                "aggregate_adapter_trimmed_fraction":
                    float(
                        aggregate_adapter_trimmed_fraction
                    ),

                "total_adapter_bases_trimmed":
                    total_adapter_bases_trimmed,

                "empty_after_trim":
                    empty_after_trim,

                "reads_shorter_than_min":
                    reads_shorter_than_min,

                "reads_longer_than_max":
                    reads_longer_than_max,

                "usable_reads_after_preprocessing":
                    usable_reads,

                "aggregate_usable_read_fraction":
                    float(
                        aggregate_usable_read_fraction
                    ),

                # ----------------------------------------------------------
                # Candidate counts
                # ----------------------------------------------------------

                "exact_candidate_positive_reads":
                    exact_positive_reads,

                "exact_full_length_candidate_reads":
                    exact_full_length_reads,

                "exact_candidate_substring_reads":
                    exact_substring_reads,

                # ----------------------------------------------------------
                # Candidate normalization
                # ----------------------------------------------------------

                "combined_candidate_cpm":
                    float(
                        combined_candidate_cpm
                    ),

                "combined_full_length_candidate_cpm":
                    float(
                        combined_full_length_candidate_cpm
                    ),

                "combined_substring_candidate_cpm":
                    float(
                        combined_substring_candidate_cpm
                    ),

                "combined_candidate_fraction":
                    float(
                        combined_candidate_fraction
                    ),

                # ----------------------------------------------------------
                # Per-run CPM distribution
                # ----------------------------------------------------------

                "run_candidate_cpm_min":
                    run_cpm_min,

                "run_candidate_cpm_max":
                    run_cpm_max,

                "run_candidate_cpm_median":
                    run_cpm_median,

                "run_candidate_cpm_mean":
                    run_cpm_mean,

                # ----------------------------------------------------------
                # Interpretation
                # ----------------------------------------------------------

                "candidate_exact_match_observed":
                    candidate_observed,

                "exact_sequence_evidence_only":
                    True,

                "zero_count_interpreted_as_biological_absence":
                    False,

                # ----------------------------------------------------------
                # Safeguards
                # ----------------------------------------------------------

                "approximate_matching_performed":
                    False,

                "mismatches_allowed":
                    False,

                "indels_allowed":
                    False,

                "reverse_complement_search_performed":
                    False,

                "parent_trna_inference_performed":
                    False,

                "genomic_origin_inference_performed":
                    False,

                "clinical_association_performed":
                    False,

                "causal_inference_performed":
                    False,

                # ----------------------------------------------------------
                # Resolution
                # ----------------------------------------------------------

                "overall_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Final consistency validation
# ============================================================================


def _validate_result(
    *,
    run_counts: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """Validate M7.4C.2 results before returning them."""

    if run_counts.empty:

        raise RuntimeError(
            "M7.4C.2 run-count output is empty."
        )

    if summary.empty or len(
        summary
    ) != 1:

        raise RuntimeError(
            "M7.4C.2 summary must contain exactly one row."
        )

    row = summary.iloc[
        0
    ]

    if int(
        row[
            "runs_processed"
        ]
    ) != len(
        run_counts
    ):

        raise RuntimeError(
            "M7.4C.2 run count does not match summary."
        )

    # ----------------------------------------------------------------------
    # Candidate decomposition
    # ----------------------------------------------------------------------

    if int(
        row[
            "exact_candidate_positive_reads"
        ]
    ) != (
        int(
            row[
                "exact_full_length_candidate_reads"
            ]
        )
        +
        int(
            row[
                "exact_candidate_substring_reads"
            ]
        )
    ):

        raise RuntimeError(
            "M7.4C.2 candidate count decomposition mismatch."
        )

    # ----------------------------------------------------------------------
    # Scientific safeguards
    # ----------------------------------------------------------------------

    required_false_fields = [
        "approximate_matching_performed",
        "mismatches_allowed",
        "indels_allowed",
        "reverse_complement_search_performed",
        "parent_trna_inference_performed",
        "genomic_origin_inference_performed",
        "clinical_association_performed",
        "causal_inference_performed",
    ]

    for field in required_false_fields:

        if bool(
            row[
                field
            ]
        ):

            raise RuntimeError(
                "M7.4C.2 safeguard violation: "
                f"{field}=True."
            )

    run_required_false_fields = [
        "approximate_matching_performed",
        "mismatches_allowed",
        "indels_allowed",
        "reverse_complement_search_performed",
        "parent_trna_inference_performed",
        "genomic_origin_inference_performed",
        "clinical_association_performed",
        "causal_inference_performed",
    ]

    for field in run_required_false_fields:

        if (
            run_counts[
                field
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
                "M7.4C.2 run-level safeguard violation: "
                f"{field}=True."
            )


# ============================================================================
# Public API
# ============================================================================


def count_external_candidate_sequence(
    *,
    run_fastq_paths: dict[str, Path],
    config: dict[str, Any],
) -> ExactSequenceCountingResult:
    """Execute M7.4C.2 exact candidate-sequence counting."""

    run_counts = build_run_counts(
        run_fastq_paths=run_fastq_paths,
        config=config,
    )

    summary = build_summary(
        run_counts=run_counts,
        config=config,
    )

    _validate_result(
        run_counts=run_counts,
        summary=summary,
    )

    return ExactSequenceCountingResult(
        run_counts=run_counts,
        summary=summary,
    )