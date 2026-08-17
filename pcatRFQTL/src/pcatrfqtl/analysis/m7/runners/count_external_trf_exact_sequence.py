"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/count_external_trf_exact_sequence.py

Description:
    Runner for M7.4C.2 Exact Candidate tRF Sequence Counting.

    The runner validates locked M7.4C.1, validates the selected SRA runs
    against the upstream acquisition manifest, downloads missing public SRA
    runs, converts them to FASTQ, and performs conservative exact candidate
    sequence counting.

    Existing FASTQ files are reused.

    In full mode, all 11 GSE80400 runs defined in configuration must be
    successfully resolved before the stage is considered complete.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pcatrfqtl.analysis.m7.external_trf_exact_sequence_counting import (
    count_external_candidate_sequence,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


# ============================================================================
# JSON helpers
# ============================================================================


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert pandas row into JSON-safe scalar values."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:

            output[
                str(
                    key
                )
            ] = None

            continue

        try:

            if pd.isna(
                value
            ):

                output[
                    str(
                        key
                    )
                ] = None

                continue

        except (
            TypeError,
            ValueError,
        ):
            pass

        if hasattr(
            value,
            "item",
        ):

            value = value.item()

        output[
            str(
                key
            )
        ] = value

    return output


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert DataFrame into JSON-safe records."""

    if dataframe.empty:
        return []

    cleaned = dataframe.astype(
        object
    ).where(
        pd.notna(
            dataframe
        ),
        None,
    )

    records = cleaned.to_dict(
        orient="records"
    )

    for record in records:

        for key, value in list(
            record.items()
        ):

            if hasattr(
                value,
                "item",
            ):

                record[
                    key
                ] = value.item()

    return records


# ============================================================================
# Runner
# ============================================================================


class M74C2ExternalTrfExactSequenceRunner:
    """Execute M7.4C.2."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        )

        self.config_path = Path(
            config_path
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M7.4C.2 configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.4C.2 config not found: {self.config_path}"
            )

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            config = yaml.safe_load(
                handle
            )

        if not isinstance(
            config,
            dict,
        ):

            raise RuntimeError(
                "M7.4C.2 config root must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.4C.2":

            raise RuntimeError(
                "Unexpected M7.4C.2 milestone."
            )

        if config.get(
            "stage"
        ) != "external_trf_exact_candidate_sequence_counting":

            raise RuntimeError(
                "Unexpected M7.4C.2 stage."
            )

        mode = (
            config
            .get(
                "execution",
                {},
            )
            .get(
                "mode"
            )
        )

        if mode not in {
            "pilot",
            "full",
        }:

            raise RuntimeError(
                f"Unsupported M7.4C.2 execution mode: {mode!r}"
            )

        return config

    # ======================================================================
    # Locked upstream
    # ======================================================================

    def _load_and_validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Load and validate locked M7.4C.1."""

        upstream = (
            config[
                "upstream_qc"
            ][
                "raw_read_acquisition"
            ]
        )

        path = (
            self.project_root
            /
            str(
                upstream[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M7.4C.1 QC not found: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        if payload.get(
            "milestone"
        ) != upstream[
            "expected_milestone"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.1 milestone."
            )

        if payload.get(
            "stage"
        ) != upstream[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected M7.4C.1 stage."
            )

        observed_status = (
            payload
            .get(
                "summary",
                {},
            )
            .get(
                "overall_status"
            )
        )

        expected_status = upstream[
            "expected_overall_status"
        ]

        if observed_status != expected_status:

            raise RuntimeError(
                "Unexpected M7.4C.1 overall status: "
                f"{observed_status!r}; "
                f"expected {expected_status!r}."
            )

        if not bool(
            payload
            .get(
                "summary",
                {},
            )
            .get(
                "acquisition_manifest_ready",
                False,
            )
        ):

            raise RuntimeError(
                "M7.4C.1 acquisition manifest is not ready."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Toolchain
    # ======================================================================

    def _require_toolchain(
        self,
    ) -> dict[str, str | None]:
        """Require SRA Toolkit commands needed for M7.4C.2."""

        required_commands = [
            "prefetch",
            "fasterq-dump",
        ]

        result: dict[
            str,
            str | None
        ] = {}

        missing: list[
            str
        ] = []

        for command in required_commands:

            executable = shutil.which(
                command
            )

            result[
                command
            ] = executable

            if executable is None:

                missing.append(
                    command
                )

        result[
            "vdb-validate"
        ] = shutil.which(
            "vdb-validate"
        )

        if missing:

            raise RuntimeError(
                "M7.4C.2 requires SRA Toolkit. "
                f"Missing commands: {missing}"
            )

        return result

    # ======================================================================
    # Selected runs
    # ======================================================================

    def _selected_runs(
        self,
        *,
        config: dict[str, Any],
    ) -> list[str]:
        """Resolve pilot or full execution run set."""

        execution = config[
            "execution"
        ]

        mode = execution[
            "mode"
        ]

        if mode == "pilot":

            values = execution[
                "pilot_run_accessions"
            ]

        else:

            values = execution[
                "full_run_accessions"
            ]

        if not isinstance(
            values,
            list,
        ) or not values:

            raise RuntimeError(
                f"No SRA runs configured for mode={mode!r}."
            )

        runs = [
            str(
                value
            ).strip()
            for value in values
        ]

        if any(
            not run
            for run in runs
        ):

            raise RuntimeError(
                "Selected run accessions contain empty values."
            )

        if len(
            runs
        ) != len(
            set(
                runs
            )
        ):

            raise RuntimeError(
                "Selected run accessions contain duplicates."
            )

        return runs

    # ======================================================================
    # Upstream manifest consistency
    # ======================================================================

    def _validate_selected_runs_against_manifest(
        self,
        *,
        selected_runs: list[str],
        upstream_qc: dict[str, Any],
    ) -> None:
        """Require every selected run to exist in locked M7.4C.1 manifest."""

        manifest = upstream_qc.get(
            "acquisition_manifest",
            [],
        )

        if not isinstance(
            manifest,
            list,
        ):

            raise RuntimeError(
                "M7.4C.1 acquisition_manifest must be a list."
            )

        upstream_runs = {
            str(
                row.get(
                    "run_accession"
                )
            )
            for row in manifest
            if isinstance(
                row,
                dict,
            )
            and row.get(
                "run_accession"
            )
        }

        missing = [
            run
            for run in selected_runs
            if run not in upstream_runs
        ]

        if missing:

            raise RuntimeError(
                "M7.4C.2 selected runs missing from locked "
                f"M7.4C.1 manifest: {missing}"
            )

    # ======================================================================
    # Command execution
    # ======================================================================

    def _run_command(
        self,
        command: list[str],
        *,
        cwd: Path | None = None,
    ) -> None:
        """Execute one external command and fail on non-zero exit."""

        logger.info(
            "Executing: %s",
            " ".join(
                command
            ),
        )

        result = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False,
        )

        if result.returncode != 0:

            raise RuntimeError(
                "Command failed.\n"
                f"Command: {' '.join(command)}\n"
                f"Return code: {result.returncode}\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

    # ======================================================================
    # FASTQ preparation
    # ======================================================================

    def _prepare_fastq(
        self,
        *,
        run_accession: str,
        config: dict[str, Any],
        toolchain: dict[str, str | None],
    ) -> tuple[
        Path,
        bool,
    ]:
        """
        Resolve or generate one FASTQ.

        Returns
        -------
        tuple
            fastq_path,
            reused_existing_fastq
        """

        paths = config[
            "paths"
        ]

        cache_dir = (
            self.project_root
            /
            str(
                paths[
                    "sra_cache_directory"
                ]
            )
        )

        fastq_dir = (
            self.project_root
            /
            str(
                paths[
                    "fastq_directory"
                ]
            )
        )

        cache_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        fastq_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        expected_fastq = (
            fastq_dir
            /
            f"{run_accession}.fastq"
        )

        # ------------------------------------------------------------------
        # Reuse existing FASTQ
        # ------------------------------------------------------------------

        if expected_fastq.exists():

            if expected_fastq.stat().st_size <= 0:

                raise RuntimeError(
                    f"Existing FASTQ is empty: {expected_fastq}"
                )

            logger.info(
                "FASTQ already exists, reusing: %s",
                expected_fastq,
            )

            return (
                expected_fastq,
                True,
            )

        # ------------------------------------------------------------------
        # Prefetch
        # ------------------------------------------------------------------

        self._run_command(
            [
                str(
                    toolchain[
                        "prefetch"
                    ]
                ),
                run_accession,
                "--output-directory",
                str(
                    cache_dir
                ),
            ]
        )

        prefetched_run = (
            cache_dir
            /
            run_accession
        )

        if not prefetched_run.exists():

            raise RuntimeError(
                "Prefetch output directory not found: "
                f"{prefetched_run}"
            )

        # ------------------------------------------------------------------
        # Validate SRA archive when possible
        # ------------------------------------------------------------------

        validator = toolchain.get(
            "vdb-validate"
        )

        if (
            config[
                "acquisition"
            ].get(
                "validate_prefetched_run",
                False,
            )
            and
            validator is not None
        ):

            self._run_command(
                [
                    validator,
                    str(
                        prefetched_run
                    ),
                ]
            )

        # ------------------------------------------------------------------
        # fasterq-dump
        # ------------------------------------------------------------------

        threads = int(
            config[
                "acquisition"
            ][
                "fasterq_threads"
            ]
        )

        self._run_command(
            [
                str(
                    toolchain[
                        "fasterq-dump"
                    ]
                ),
                str(
                    prefetched_run
                ),
                "--outdir",
                str(
                    fastq_dir
                ),
                "--threads",
                str(
                    threads
                ),
            ]
        )

        # ------------------------------------------------------------------
        # Resolve output
        # ------------------------------------------------------------------

        if expected_fastq.exists():

            fastq_path = expected_fastq

        else:

            candidates = sorted(
                fastq_dir.glob(
                    f"{run_accession}*.fastq"
                )
            )

            if len(
                candidates
            ) != 1:

                raise RuntimeError(
                    "Expected exactly one SINGLE-end FASTQ for "
                    f"{run_accession}, found: {candidates}"
                )

            fastq_path = candidates[
                0
            ]

        if fastq_path.stat().st_size <= 0:

            raise RuntimeError(
                f"Generated FASTQ is empty: {fastq_path}"
            )

        return (
            fastq_path,
            False,
        )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.4C.2."""

        config = self._load_config()

        (
            upstream_qc,
            upstream_path,
        ) = self._load_and_validate_upstream(
            config=config,
        )

        toolchain = self._require_toolchain()

        selected_runs = self._selected_runs(
            config=config,
        )

        self._validate_selected_runs_against_manifest(
            selected_runs=selected_runs,
            upstream_qc=upstream_qc,
        )

        logger.info(
            "M7.4C.2 mode=%s | selected runs=%d.",
            config[
                "execution"
            ][
                "mode"
            ],
            len(
                selected_runs
            ),
        )

        # ------------------------------------------------------------------
        # FASTQ acquisition / reuse
        # ------------------------------------------------------------------

        run_fastq_paths: dict[
            str,
            Path
        ] = {}

        reused_runs: list[
            str
        ] = []

        generated_runs: list[
            str
        ] = []

        for index, run_accession in enumerate(
            selected_runs,
            start=1,
        ):

            logger.info(
                "Preparing run %d/%d: %s.",
                index,
                len(
                    selected_runs
                ),
                run_accession,
            )

            (
                fastq_path,
                reused,
            ) = self._prepare_fastq(
                run_accession=run_accession,
                config=config,
                toolchain=toolchain,
            )

            run_fastq_paths[
                run_accession
            ] = fastq_path

            if reused:

                reused_runs.append(
                    run_accession
                )

            else:

                generated_runs.append(
                    run_accession
                )

        # ------------------------------------------------------------------
        # Exact counting
        # ------------------------------------------------------------------

        result = count_external_candidate_sequence(
            run_fastq_paths=run_fastq_paths,
            config=config,
        )

        # ------------------------------------------------------------------
        # Output paths
        # ------------------------------------------------------------------

        processed_dir = (
            self.project_root
            /
            str(
                config[
                    "paths"
                ][
                    "processed_directory"
                ]
            )
        )

        processed_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        run_counts_path = (
            processed_dir
            /
            str(
                config[
                    "paths"
                ][
                    "run_count_filename"
                ]
            )
        )

        summary_path = (
            processed_dir
            /
            str(
                config[
                    "paths"
                ][
                    "summary_filename"
                ]
            )
        )

        qc_path = (
            self.project_root
            /
            str(
                config[
                    "paths"
                ][
                    "qc_relative_path"
                ]
            )
        )

        qc_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ------------------------------------------------------------------
        # Persist processed artifacts
        # ------------------------------------------------------------------

        write_parquet(
            result.run_counts,
            run_counts_path,
            index=False,
        )

        write_parquet(
            result.summary,
            summary_path,
            index=False,
        )

        summary = _json_safe_row(
            result.summary.iloc[
                0
            ]
        )

        # ------------------------------------------------------------------
        # QC report
        # ------------------------------------------------------------------

        report = {
            "milestone":
                "M7.4C.2",

            "stage":
                "external_trf_exact_candidate_sequence_counting",

            "execution":
                config[
                    "execution"
                ],

            "candidate":
                config[
                    "candidate"
                ],

            "dataset":
                config[
                    "dataset"
                ],

            "preprocessing":
                config[
                    "preprocessing"
                ],

            "counting":
                config[
                    "counting"
                ],

            "normalization":
                config[
                    "normalization"
                ],

            "interpretation":
                config[
                    "interpretation"
                ],

            "upstream_validation": {
                "m7_4c1_loaded":
                    True,

                "m7_4c1_path":
                    str(
                        upstream_path
                    ),

                "m7_4c1_status":
                    (
                        upstream_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "overall_status"
                        )
                    ),

                "acquisition_manifest_ready":
                    bool(
                        upstream_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "acquisition_manifest_ready",
                            False,
                        )
                    ),
            },

            "toolchain":
                toolchain,

            "fastq_resolution": {
                "selected_runs":
                    selected_runs,

                "selected_run_count":
                    len(
                        selected_runs
                    ),

                "reused_fastq_runs":
                    reused_runs,

                "reused_fastq_count":
                    len(
                        reused_runs
                    ),

                "generated_fastq_runs":
                    generated_runs,

                "generated_fastq_count":
                    len(
                        generated_runs
                    ),

                "all_selected_fastq_resolved":
                    bool(
                        len(
                            run_fastq_paths
                        )
                        ==
                        len(
                            selected_runs
                        )
                    ),
            },

            "summary":
                summary,

            "run_counts":
                _records_with_json_nulls(
                    result.run_counts
                ),

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs": {
                "run_counts":
                    str(
                        run_counts_path
                    ),

                "summary":
                    str(
                        summary_path
                    ),
            },
        }

        with qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            )

        # ------------------------------------------------------------------
        # Logging
        # ------------------------------------------------------------------

        logger.info(
            "M7.4C.2 complete."
        )

        logger.info(
            "Runs=%d | exact-match runs=%d | "
            "total reads=%d | usable=%d.",
            summary[
                "runs_processed"
            ],
            summary[
                "runs_with_exact_match"
            ],
            summary[
                "total_reads"
            ],
            summary[
                "usable_reads_after_preprocessing"
            ],
        )

        logger.info(
            "Adapter-trimmed=%d (%.6f) | usable fraction=%.6f.",
            summary[
                "reads_adapter_trimmed"
            ],
            summary[
                "aggregate_adapter_trimmed_fraction"
            ],
            summary[
                "aggregate_usable_read_fraction"
            ],
        )

        logger.info(
            "Exact-positive=%d | full-length=%d | substring=%d.",
            summary[
                "exact_candidate_positive_reads"
            ],
            summary[
                "exact_full_length_candidate_reads"
            ],
            summary[
                "exact_candidate_substring_reads"
            ],
        )

        logger.info(
            "Combined CPM=%.6f | full-length CPM=%.6f | "
            "substring CPM=%.6f.",
            summary[
                "combined_candidate_cpm"
            ],
            summary[
                "combined_full_length_candidate_cpm"
            ],
            summary[
                "combined_substring_candidate_cpm"
            ],
        )

        logger.info(
            "Run CPM min=%.6f | median=%.6f | mean=%.6f | max=%.6f.",
            summary[
                "run_candidate_cpm_min"
            ],
            summary[
                "run_candidate_cpm_median"
            ],
            summary[
                "run_candidate_cpm_mean"
            ],
            summary[
                "run_candidate_cpm_max"
            ],
        )

        logger.info(
            "Candidate exact match observed=%s.",
            summary[
                "candidate_exact_match_observed"
            ],
        )

        logger.info(
            "Selected run set complete=%s | full run set complete=%s.",
            summary[
                "selected_run_set_complete"
            ],
            summary[
                "full_run_set_complete"
            ],
        )

        logger.info(
            "Overall=%s | next=%s.",
            summary[
                "overall_status"
            ],
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC output=%s.",
            qc_path,
        )

        return report