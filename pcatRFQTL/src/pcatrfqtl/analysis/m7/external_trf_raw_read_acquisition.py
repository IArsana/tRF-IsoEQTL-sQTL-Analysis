"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/external_trf_raw_read_acquisition.py

Description:
    Core logic for M7.4C.1 External tRF Raw-Read Acquisition Audit.

    This stage consumes run-level SRA metadata for the pilot dataset
    GSE80400 and creates an acquisition manifest for downstream exact
    candidate-sequence counting.

    M7.4C.1 does NOT download raw reads and does NOT quantify tRFs.

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
class ExternalTrfRawReadAcquisitionResult:
    """Container for M7.4C.1 output tables."""

    run_inventory: pd.DataFrame
    acquisition_manifest: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional text."""

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


def _as_int(
    value: Any,
    *,
    default: int = 0,
) -> int:
    """Convert scalar to integer safely."""

    if value is None:
        return default

    try:
        if pd.isna(value):
            return default
    except (
        TypeError,
        ValueError,
    ):
        pass

    try:
        return int(float(value))
    except (
        TypeError,
        ValueError,
    ):
        return default


def _normalize_column_names(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Strip whitespace from SRA metadata column names."""

    result = dataframe.copy()

    result.columns = [
        str(column).strip()
        for column in result.columns
    ]

    return result


def _get_column(
    dataframe: pd.DataFrame,
    *candidates: str,
) -> str | None:
    """
    Find the first available column from a list of possible SRA names.
    """

    lookup = {
        str(column).lower():
            str(column)
        for column in dataframe.columns
    }

    for candidate in candidates:

        hit = lookup.get(
            candidate.lower()
        )

        if hit is not None:
            return hit

    return None


# ============================================================================
# Run inventory
# ============================================================================


def build_run_inventory(
    *,
    runinfo: pd.DataFrame,
) -> pd.DataFrame:
    """
    Normalize NCBI SRA RunInfo metadata into an M7.4C.1 run inventory.
    """

    runinfo = _normalize_column_names(
        runinfo
    )

    if runinfo.empty:
        raise RuntimeError(
            "SRA RunInfo metadata is empty."
        )

    run_col = _get_column(
        runinfo,
        "Run",
    )

    if run_col is None:
        raise RuntimeError(
            "SRA RunInfo does not contain a Run accession column."
        )

    experiment_col = _get_column(
        runinfo,
        "Experiment",
    )

    sample_col = _get_column(
        runinfo,
        "Sample",
    )

    biosample_col = _get_column(
        runinfo,
        "BioSample",
    )

    spots_col = _get_column(
        runinfo,
        "spots",
        "Spots",
    )

    bases_col = _get_column(
        runinfo,
        "bases",
        "Bases",
    )

    size_col = _get_column(
        runinfo,
        "size_MB",
        "Size_MB",
    )

    strategy_col = _get_column(
        runinfo,
        "LibraryStrategy",
    )

    source_col = _get_column(
        runinfo,
        "LibrarySource",
    )

    selection_col = _get_column(
        runinfo,
        "LibrarySelection",
    )

    layout_col = _get_column(
        runinfo,
        "LibraryLayout",
    )

    platform_col = _get_column(
        runinfo,
        "Platform",
    )

    model_col = _get_column(
        runinfo,
        "Model",
    )

    organism_col = _get_column(
        runinfo,
        "ScientificName",
    )

    project_col = _get_column(
        runinfo,
        "BioProject",
    )

    study_col = _get_column(
        runinfo,
        "SRAStudy",
    )

    release_col = _get_column(
        runinfo,
        "ReleaseDate",
    )

    load_col = _get_column(
        runinfo,
        "LoadDate",
    )

    rows: list[
        dict[str, Any]
    ] = []

    for _, row in runinfo.iterrows():

        run_accession = _normalize_text(
            row[
                run_col
            ]
        )

        if run_accession is None:
            continue

        rows.append(
            {
                "run_accession":
                    run_accession,

                "experiment_accession":
                    (
                        _normalize_text(
                            row[
                                experiment_col
                            ]
                        )
                        if experiment_col
                        else None
                    ),

                "sample_accession":
                    (
                        _normalize_text(
                            row[
                                sample_col
                            ]
                        )
                        if sample_col
                        else None
                    ),

                "biosample_accession":
                    (
                        _normalize_text(
                            row[
                                biosample_col
                            ]
                        )
                        if biosample_col
                        else None
                    ),

                "spots":
                    (
                        _as_int(
                            row[
                                spots_col
                            ]
                        )
                        if spots_col
                        else 0
                    ),

                "bases":
                    (
                        _as_int(
                            row[
                                bases_col
                            ]
                        )
                        if bases_col
                        else 0
                    ),

                "size_mb":
                    (
                        float(
                            row[
                                size_col
                            ]
                        )
                        if (
                            size_col
                            and
                            pd.notna(
                                row[
                                    size_col
                                ]
                            )
                        )
                        else None
                    ),

                "library_strategy":
                    (
                        _normalize_text(
                            row[
                                strategy_col
                            ]
                        )
                        if strategy_col
                        else None
                    ),

                "library_source":
                    (
                        _normalize_text(
                            row[
                                source_col
                            ]
                        )
                        if source_col
                        else None
                    ),

                "library_selection":
                    (
                        _normalize_text(
                            row[
                                selection_col
                            ]
                        )
                        if selection_col
                        else None
                    ),

                "library_layout":
                    (
                        _normalize_text(
                            row[
                                layout_col
                            ]
                        )
                        if layout_col
                        else None
                    ),

                "platform":
                    (
                        _normalize_text(
                            row[
                                platform_col
                            ]
                        )
                        if platform_col
                        else None
                    ),

                "instrument_model":
                    (
                        _normalize_text(
                            row[
                                model_col
                            ]
                        )
                        if model_col
                        else None
                    ),

                "organism":
                    (
                        _normalize_text(
                            row[
                                organism_col
                            ]
                        )
                        if organism_col
                        else None
                    ),

                "bioproject":
                    (
                        _normalize_text(
                            row[
                                project_col
                            ]
                        )
                        if project_col
                        else None
                    ),

                "sra_study":
                    (
                        _normalize_text(
                            row[
                                study_col
                            ]
                        )
                        if study_col
                        else None
                    ),

                "release_date":
                    (
                        _normalize_text(
                            row[
                                release_col
                            ]
                        )
                        if release_col
                        else None
                    ),

                "load_date":
                    (
                        _normalize_text(
                            row[
                                load_col
                            ]
                        )
                        if load_col
                        else None
                    ),
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        raise RuntimeError(
            "No SRA run accessions were recovered."
        )

    if result[
        "run_accession"
    ].duplicated().any():

        duplicated = (
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
            "Duplicate SRA Run accessions detected: "
            f"{duplicated}"
        )

    return result.sort_values(
        "run_accession",
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Metadata validation
# ============================================================================


def validate_run_inventory(
    *,
    run_inventory: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Validate pilot SRA metadata against expected technical properties."""

    expected = config[
        "expected_metadata"
    ]

    minimum_runs = int(
        expected[
            "minimum_runs"
        ]
    )

    run_count = int(
        len(
            run_inventory
        )
    )

    enough_runs = bool(
        run_count
        >=
        minimum_runs
    )

    expected_organism = str(
        expected[
            "organism"
        ]
    )

    organism_match = bool(
        (
            run_inventory[
                "organism"
            ]
            .fillna(
                ""
            )
            ==
            expected_organism
        ).all()
    )

    expected_strategy = str(
        expected[
            "library_strategy"
        ]
    ).lower()

    strategy_match = bool(
        (
            run_inventory[
                "library_strategy"
            ]
            .fillna(
                ""
            )
            .str.lower()
            ==
            expected_strategy
        ).all()
    )

    expected_layout = str(
        expected[
            "library_layout"
        ]
    ).lower()

    layout_match = bool(
        (
            run_inventory[
                "library_layout"
            ]
            .fillna(
                ""
            )
            .str.lower()
            ==
            expected_layout
        ).all()
    )

    expected_platform = str(
        expected[
            "platform"
        ]
    ).lower()

    platform_match = bool(
        (
            run_inventory[
                "platform"
            ]
            .fillna(
                ""
            )
            .str.lower()
            ==
            expected_platform
        ).all()
    )

    public_runs_present = bool(
        run_count
        >
        0
    )

    metadata_complete = bool(
        enough_runs
        and
        organism_match
        and
        strategy_match
        and
        layout_match
        and
        platform_match
    )

    return {
        "run_count":
            run_count,

        "minimum_run_requirement_met":
            enough_runs,

        "organism_match":
            organism_match,

        "library_strategy_match":
            strategy_match,

        "library_layout_match":
            layout_match,

        "platform_match":
            platform_match,

        "public_runs_present":
            public_runs_present,

        "metadata_complete":
            metadata_complete,
    }


# ============================================================================
# Acquisition manifest
# ============================================================================


def build_acquisition_manifest(
    *,
    run_inventory: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build raw-read acquisition manifest.

    This manifest is consumed by M7.4C.2 but does not itself download data.
    """

    manifest = run_inventory[
        [
            "run_accession",
            "experiment_accession",
            "sample_accession",
            "biosample_accession",
            "spots",
            "bases",
            "size_mb",
            "library_layout",
        ]
    ].copy()

    manifest[
        "prefetch_planned"
    ] = True

    manifest[
        "fasterq_dump_planned"
    ] = True

    manifest[
        "split_files_required"
    ] = (
        manifest[
            "library_layout"
        ]
        .fillna(
            ""
        )
        .str.upper()
        ==
        "PAIRED"
    )

    manifest[
        "download_performed"
    ] = False

    manifest[
        "fastq_generated"
    ] = False

    manifest[
        "fastq_checksum_verified"
    ] = False

    manifest[
        "read_qc_performed"
    ] = False

    manifest[
        "adapter_trimmed"
    ] = False

    manifest[
        "candidate_sequence_search_performed"
    ] = False

    manifest[
        "candidate_exact_count_performed"
    ] = False

    return manifest


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    run_inventory: pd.DataFrame,
    validation: dict[str, Any],
    toolchain: dict[str, bool],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.4C.1 summary."""

    public_runs_present = bool(
        validation[
            "public_runs_present"
        ]
    )

    metadata_complete = bool(
        validation[
            "metadata_complete"
        ]
    )

    toolchain_ready = bool(
        toolchain.get(
            "prefetch",
            False,
        )
        and
        toolchain.get(
            "fasterq-dump",
            False,
        )
    )

    if (
        public_runs_present
        and
        metadata_complete
    ):

        acquisition_ready = True

        if toolchain_ready:

            acquisition_status = (
                config[
                    "statuses"
                ][
                    "acquisition"
                ][
                    "ready"
                ]
            )

        else:

            acquisition_status = (
                config[
                    "statuses"
                ][
                    "acquisition"
                ][
                    "ready_toolchain_missing"
                ]
            )

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

    elif not public_runs_present:

        acquisition_ready = False

        acquisition_status = (
            config[
                "statuses"
            ][
                "acquisition"
            ][
                "no_public_runs"
            ]
        )

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "blocked"
            ]
        )

        next_stage = (
            config[
                "next_stage"
            ][
                "if_blocked"
            ]
        )

    else:

        acquisition_ready = False

        acquisition_status = (
            config[
                "statuses"
            ][
                "acquisition"
            ][
                "metadata_incomplete"
            ]
        )

        overall_status = (
            config[
                "statuses"
            ][
                "overall"
            ][
                "blocked"
            ]
        )

        next_stage = (
            config[
                "next_stage"
            ][
                "if_blocked"
            ]
        )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.4C.1",

                "dataset_key":
                    config[
                        "pilot_dataset"
                    ][
                        "dataset_key"
                    ],

                "geo_accession":
                    config[
                        "pilot_dataset"
                    ][
                        "geo_accession"
                    ],

                "sra_study":
                    config[
                        "pilot_dataset"
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

                "runs_identified":
                    int(
                        len(
                            run_inventory
                        )
                    ),

                "total_spots":
                    int(
                        run_inventory[
                            "spots"
                        ].sum()
                    ),

                "total_bases":
                    int(
                        run_inventory[
                            "bases"
                        ].sum()
                    ),

                "public_runs_present":
                    public_runs_present,

                "run_metadata_complete":
                    metadata_complete,

                "prefetch_available":
                    bool(
                        toolchain.get(
                            "prefetch",
                            False,
                        )
                    ),

                "fasterq_dump_available":
                    bool(
                        toolchain.get(
                            "fasterq-dump",
                            False,
                        )
                    ),

                "sra_toolkit_ready":
                    toolchain_ready,

                "raw_sequence_download_performed":
                    False,

                "fastq_generation_performed":
                    False,

                "candidate_sequence_search_performed":
                    False,

                "exact_match_counting_performed":
                    False,

                "trf_quantification_performed":
                    False,

                "candidate_detection_claimed":
                    False,

                "candidate_absence_claimed":
                    False,

                "acquisition_manifest_ready":
                    acquisition_ready,

                "acquisition_status":
                    acquisition_status,

                "overall_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def audit_external_trf_raw_read_acquisition(
    *,
    runinfo: pd.DataFrame,
    toolchain: dict[str, bool],
    config: dict[str, Any],
) -> ExternalTrfRawReadAcquisitionResult:
    """Execute M7.4C.1 core analysis."""

    run_inventory = build_run_inventory(
        runinfo=runinfo,
    )

    validation = validate_run_inventory(
        run_inventory=run_inventory,
        config=config,
    )

    acquisition_manifest = build_acquisition_manifest(
        run_inventory=run_inventory,
    )

    summary = build_summary(
        run_inventory=run_inventory,
        validation=validation,
        toolchain=toolchain,
        config=config,
    )

    return ExternalTrfRawReadAcquisitionResult(
        run_inventory=run_inventory,
        acquisition_manifest=acquisition_manifest,
        summary=summary,
    )