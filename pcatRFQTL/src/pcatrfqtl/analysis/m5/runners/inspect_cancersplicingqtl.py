"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/inspect_cancersplicingqtl.py

Description:
    Runner for M5.6C.1 CancerSplicingQTL PRAD download inspection.

    The current official PRAD artifact is expected at:

        data/raw/m5/cancersplicingqtl/prad/PRAD_sQTLs.xlsx

    This runner:
        - discovers supported raw files;
        - reads Excel workbooks and text archives;
        - preserves file/sheet provenance;
        - validates table-schema compatibility;
        - executes the M5.6C.1 core inspection;
        - writes deterministic Parquet and JSON QC outputs.

    Important safeguards:
        - Raw files are never modified.
        - Incompatible sheets are never silently merged.
        - Significant-hit density is not interpreted as complete test density.
        - No formal colocalization is performed.
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

import gzip
import json
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pcatrfqtl.analysis.m5.cancersplicingqtl_inspection import (
    DATASET_NOT_PRESENT,
    inspect_cancersplicingqtl,
)
from pcatrfqtl.io.parquet import (
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


# ============================================================================
# Runner
# ============================================================================


class M56C1CancerSplicingQTLRunner:
    """Execute M5.6C.1 CancerSplicingQTL PRAD dataset inspection."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        )

        self.config_path = Path(
            config_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ======================================================================
    # Output paths
    # ======================================================================

    @property
    def schema_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prad_dataset_schema.parquet"
        )

    @property
    def feature_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prad_feature_coverage.parquet"
        )

    @property
    def readiness_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prad_coloc_readiness.parquet"
        )

    @property
    def candidate_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "candidate_cancersplicingqtl_resolution.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_6c_1_cancersplicingqtl_inspection.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M5.6C.1 config not found: "
                f"{self.config_path}"
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

            raise ValueError(
                "M5.6C.1 configuration must be a YAML mapping."
            )

        required_sections = {
            "resource",
            "source_scope",
            "candidate_leads",
            "expected_schema",
            "standard_error_audit",
            "inspection_policy",
        }

        missing = (
            required_sections
            -
            set(
                config
            )
        )

        if missing:

            raise ValueError(
                "M5.6C.1 config missing sections: "
                f"{sorted(missing)}"
            )

        return config

    # ======================================================================
    # Raw directory and discovery
    # ======================================================================

    def _raw_directory(
        self,
        config: dict[str, Any],
    ) -> Path:

        relative = (
            config[
                "resource"
            ][
                "raw_directory"
            ][
                "relative_path"
            ]
        )

        return (
            self.project_root
            /
            str(
                relative
            )
        )

    def _discover_files(
        self,
        config: dict[str, Any],
    ) -> list[Path]:

        directory = self._raw_directory(
            config
        )

        if not directory.exists():

            return []

        extensions = {
            str(
                extension
            ).lower()
            for extension
            in config[
                "resource"
            ][
                "accepted_extensions"
            ]
        }

        recursive = bool(
            config[
                "resource"
            ].get(
                "recursive_file_discovery",
                True,
            )
        )

        iterator = (
            directory.rglob(
                "*"
            )
            if recursive
            else directory.glob(
                "*"
            )
        )

        return sorted(
            path
            for path
            in iterator
            if (
                path.is_file()
                and
                any(
                    str(
                        path
                    ).lower().endswith(
                        extension
                    )
                    for extension
                    in extensions
                )
            )
        )

    # ======================================================================
    # Readers
    # ======================================================================

    @staticmethod
    def _detect_separator(
        filename: str,
    ) -> str:

        lower = str(
            filename
        ).lower()

        if (
            lower.endswith(
                ".csv"
            )
            or
            ".csv." in lower
        ):

            return ","

        return "\t"

    def _read_plain_file(
        self,
        path: Path,
    ) -> pd.DataFrame:

        return pd.read_csv(
            path,
            sep=self._detect_separator(
                path.name
            ),
            low_memory=False,
        )

    def _read_gzip_file(
        self,
        path: Path,
    ) -> pd.DataFrame:

        with gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace",
        ) as handle:

            return pd.read_csv(
                handle,
                sep=self._detect_separator(
                    path.name
                ),
                low_memory=False,
            )

    def _read_excel_file(
        self,
        path: Path,
        *,
        config: dict[str, Any],
    ) -> list[
        tuple[
            str,
            pd.DataFrame,
        ]
    ]:

        workbook_policy = config.get(
            "workbook_policy",
            {},
        )

        inspect_all = bool(
            workbook_policy.get(
                "inspect_all_sheets",
                True,
            )
        )

        skip_empty = bool(
            workbook_policy.get(
                "skip_empty_sheets",
                True,
            )
        )

        workbook = pd.ExcelFile(
            path
        )

        if inspect_all:

            sheets = list(
                workbook.sheet_names
            )

        else:

            sheets = (
                [
                    workbook.sheet_names[
                        0
                    ]
                ]
                if workbook.sheet_names
                else []
            )

        logger.info(
            "Excel workbook %s sheets: %s",
            path.name,
            sheets,
        )

        results: list[
            tuple[
                str,
                pd.DataFrame,
            ]
        ] = []

        for sheet_name in sheets:

            frame = pd.read_excel(
                workbook,
                sheet_name=sheet_name,
            )

            if (
                frame.empty
                and
                skip_empty
            ):

                logger.info(
                    "Skipping empty sheet %s [%s].",
                    path.name,
                    sheet_name,
                )

                continue

            results.append(
                (
                    str(
                        sheet_name
                    ),
                    frame,
                )
            )

        return results

    def _read_zip_file(
        self,
        path: Path,
    ) -> list[
        tuple[
            str,
            pd.DataFrame,
        ]
    ]:

        results: list[
            tuple[
                str,
                pd.DataFrame,
            ]
        ] = []

        with zipfile.ZipFile(
            path,
            "r",
        ) as archive:

            for member in sorted(
                archive.namelist()
            ):

                lower = member.lower()

                if not (
                    lower.endswith(
                        ".txt"
                    )
                    or
                    lower.endswith(
                        ".tsv"
                    )
                    or
                    lower.endswith(
                        ".csv"
                    )
                ):

                    continue

                separator = (
                    ","
                    if lower.endswith(
                        ".csv"
                    )
                    else "\t"
                )

                with archive.open(
                    member
                ) as handle:

                    frame = pd.read_csv(
                        handle,
                        sep=separator,
                        low_memory=False,
                    )

                results.append(
                    (
                        member,
                        frame,
                    )
                )

        return results

    def _read_source_file(
        self,
        path: Path,
        *,
        config: dict[str, Any],
    ) -> list[
        tuple[
            str,
            pd.DataFrame,
        ]
    ]:

        lower = (
            path.name.lower()
        )

        if lower.endswith(
            ".zip"
        ):

            return self._read_zip_file(
                path
            )

        if (
            lower.endswith(
                ".xlsx"
            )
            or
            lower.endswith(
                ".xls"
            )
        ):

            return self._read_excel_file(
                path,
                config=config,
            )

        if lower.endswith(
            ".gz"
        ):

            return [
                (
                    path.name,
                    self._read_gzip_file(
                        path
                    ),
                )
            ]

        return [
            (
                path.name,
                self._read_plain_file(
                    path
                ),
            )
        ]

    # ======================================================================
    # Dataset loading
    # ======================================================================

    def _load_dataset(
        self,
        paths: list[Path],
        *,
        config: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        list[
            dict[
                str,
                Any,
            ]
        ],
    ]:

        frames: list[
            pd.DataFrame
        ] = []

        source_metadata: list[
            dict[
                str,
                Any,
            ]
        ] = []

        for path in paths:

            members = self._read_source_file(
                path,
                config=config,
            )

            for member_name, frame in members:

                analytical_columns = [
                    str(
                        column
                    )
                    for column
                    in frame.columns
                ]

                enriched = frame.copy()

                enriched[
                    "_m56c1_source_file"
                ] = str(
                    path
                )

                enriched[
                    "_m56c1_source_member"
                ] = str(
                    member_name
                )

                frames.append(
                    enriched
                )

                source_metadata.append(
                    {
                        "path":
                            str(
                                path
                            ),

                        "member":
                            str(
                                member_name
                            ),

                        "rows":
                            int(
                                len(
                                    frame
                                )
                            ),

                        "column_count":
                            int(
                                len(
                                    analytical_columns
                                )
                            ),

                        "columns":
                            analytical_columns,
                    }
                )

                logger.info(
                    "Loaded %s [%s]: %d rows, %d columns.",
                    path.name,
                    member_name,
                    len(
                        frame
                    ),
                    len(
                        analytical_columns
                    ),
                )

        if not frames:

            return (
                pd.DataFrame(),
                source_metadata,
            )

        # ------------------------------------------------------------------
        # Schema compatibility
        # ------------------------------------------------------------------

        schema_groups: dict[
            tuple[
                str,
                ...,
            ],
            list[
                pd.DataFrame
            ],
        ] = {}

        for frame in frames:

            schema = tuple(
                sorted(
                    str(
                        column
                    )
                    for column
                    in frame.columns
                    if not str(
                        column
                    ).startswith(
                        "_m56c1_"
                    )
                )
            )

            schema_groups.setdefault(
                schema,
                [],
            ).append(
                frame
            )

        fail_on_incompatible = bool(
            config.get(
                "workbook_policy",
                {},
            ).get(
                "fail_on_incompatible_sheet_schema",
                True,
            )
        )

        if (
            len(
                schema_groups
            )
            > 1
            and
            fail_on_incompatible
        ):

            schema_summary = [
                {
                    "table_count":
                        len(
                            grouped_frames
                        ),

                    "column_count":
                        len(
                            schema
                        ),

                    "columns":
                        list(
                            schema
                        ),
                }
                for schema, grouped_frames
                in schema_groups.items()
            ]

            raise RuntimeError(
                "CancerSplicingQTL input contains incompatible "
                "analytical schemas. "
                f"Schema groups: {schema_summary}"
            )

        combined = pd.concat(
            frames,
            ignore_index=True,
            sort=False,
        )

        logger.info(
            "Combined CancerSplicingQTL PRAD dataset: "
            "%d rows, %d columns.",
            len(
                combined
            ),
            len(
                combined.columns
            ),
        )

        return (
            combined,
            source_metadata,
        )

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        raw_directory = self._raw_directory(
            config
        )

        source_files = self._discover_files(
            config
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        if not source_files:

            report = {
                "milestone":
                    "M5.6C.1",

                "stage":
                    "cancersplicingqtl_prad_download_inspection",

                "policy": {
                    "dataset_directly_inspected":
                        False,

                    "formal_colocalization_performed":
                        False,

                    "fine_mapping_performed":
                        False,

                    "causal_inference_performed":
                        False,
                },

                "summary": {
                    "raw_directory":
                        str(
                            raw_directory
                        ),

                    "source_files_discovered":
                        0,

                    "classification":
                        DATASET_NOT_PRESENT,

                    "formal_coloc_ready":
                        False,

                    "coloc_abf_ready":
                        False,

                    "coloc_susie_ready":
                        False,
                },
            }

            with self.qc_output.open(
                "w",
                encoding="utf-8",
            ) as handle:

                json.dump(
                    report,
                    handle,
                    indent=2,
                    ensure_ascii=False,
                )

            logger.warning(
                "No CancerSplicingQTL PRAD input found under %s.",
                raw_directory,
            )

            return report

        logger.info(
            "Discovered %d source file(s).",
            len(
                source_files
            ),
        )

        (
            dataframe,
            source_metadata,
        ) = self._load_dataset(
            source_files,
            config=config,
        )

        result = inspect_cancersplicingqtl(
            dataframe,
            config=config,
            source_file_count=len(
                source_files
            ),
        )

        write_parquet(
            result.schema,
            self.schema_output,
            index=False,
        )

        write_parquet(
            result.feature_coverage,
            self.feature_output,
            index=False,
        )

        write_parquet(
            result.readiness,
            self.readiness_output,
            index=False,
        )

        write_parquet(
            result.candidates,
            self.candidate_output,
            index=False,
        )

        readiness = result.readiness.iloc[
            0
        ]

        report = {
            "milestone":
                "M5.6C.1",

            "stage":
                "cancersplicingqtl_prad_download_inspection",

            "policy": {
                "dataset_directly_inspected":
                    True,

                "raw_dataset_modified":
                    False,

                "significant_only_qtl_allowed_for_coloc":
                    False,

                "missing_variants_interpreted_as_null":
                    False,

                "allele_orientation_inferred":
                    False,

                "diagnostic_se_derivation_allowed":
                    True,

                "derived_se_used_for_formal_coloc":
                    False,

                "cross_build_coordinate_join_performed":
                    False,

                "canonical_rsid_cross_build_matching_allowed":
                    True,

                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "inputs": {
                "raw_directory":
                    str(
                        raw_directory
                    ),

                "source_files":
                    [
                        str(
                            path
                        )
                        for path
                        in source_files
                    ],
            },

            "summary": {
                "source_files_discovered":
                    int(
                        len(
                            source_files
                        )
                    ),

                "source_tables_loaded":
                    int(
                        len(
                            source_metadata
                        )
                    ),

                "association_rows":
                    int(
                        readiness[
                            "association_rows"
                        ]
                    ),

                "genome_build":
                    str(
                        readiness[
                            "genome_build"
                        ]
                    ),

                "canonical_rsids":
                    int(
                        readiness[
                            "canonical_rsid_count"
                        ]
                    ),

                "canonical_rsid_fraction":
                    float(
                        readiness[
                            "canonical_rsid_fraction"
                        ]
                    ),

                "features":
                    int(
                        readiness[
                            "feature_count"
                        ]
                    ),

                "source_significance_filtered":
                    bool(
                        readiness[
                            "source_significance_filtered"
                        ]
                    ),

                "complete_tested_variant_feature_matrix_available":
                    bool(
                        readiness[
                            "complete_tested_variant_feature_matrix_available"
                        ]
                    ),

                "dense_unfiltered_summary_statistics_available":
                    bool(
                        readiness[
                            "dense_unfiltered_summary_statistics_available"
                        ]
                    ),

                "minimum_p_value":
                    readiness[
                        "minimum_p_value"
                    ],

                "maximum_p_value":
                    readiness[
                        "maximum_p_value"
                    ],

                "p_values_above_0_05":
                    int(
                        readiness[
                            "p_values_above_0_05"
                        ]
                    ),

                "beta_directly_reported":
                    bool(
                        readiness[
                            "beta_directly_reported"
                        ]
                    ),

                "t_statistic_directly_reported":
                    bool(
                        readiness[
                            "t_statistic_directly_reported"
                        ]
                    ),

                "standard_error_directly_reported":
                    bool(
                        readiness[
                            "standard_error_directly_reported"
                        ]
                    ),

                "derived_se_diagnostic_valid_fraction":
                    float(
                        readiness[
                            "derived_se_diagnostic_valid_fraction"
                        ]
                    ),

                "sample_size_column_available":
                    bool(
                        readiness[
                            "sample_size_column_available"
                        ]
                    ),

                "allele_frequency_available":
                    bool(
                        readiness[
                            "allele_frequency_available"
                        ]
                    ),

                "effect_allele_orientation_verified":
                    bool(
                        readiness[
                            "effect_allele_orientation_verified"
                        ]
                    ),

                "classification":
                    str(
                        readiness[
                            "classification"
                        ]
                    ),

                "formal_coloc_ready":
                    bool(
                        readiness[
                            "formal_coloc_ready"
                        ]
                    ),

                "coloc_abf_ready":
                    bool(
                        readiness[
                            "coloc_abf_ready"
                        ]
                    ),

                "coloc_susie_ready":
                    bool(
                        readiness[
                            "coloc_susie_ready"
                        ]
                    ),
            },

            "source_tables":
                source_metadata,

            "candidate_status":
                result.candidates.to_dict(
                    orient="records"
                ),

            "outputs": {
                "schema":
                    str(
                        self.schema_output
                    ),

                "feature_coverage":
                    str(
                        self.feature_output
                    ),

                "readiness":
                    str(
                        self.readiness_output
                    ),

                "candidate_resolution":
                    str(
                        self.candidate_output
                    ),
            },
        }

        with self.qc_output.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M5.6C.1 complete."
        )

        logger.info(
            "Rows=%d | canonical rsIDs=%d | features=%d.",
            report[
                "summary"
            ][
                "association_rows"
            ],
            report[
                "summary"
            ][
                "canonical_rsids"
            ],
            report[
                "summary"
            ][
                "features"
            ],
        )

        logger.info(
            "Source significance-filtered=%s | "
            "full tested matrix=%s | dense unfiltered=%s.",
            report[
                "summary"
            ][
                "source_significance_filtered"
            ],
            report[
                "summary"
            ][
                "complete_tested_variant_feature_matrix_available"
            ],
            report[
                "summary"
            ][
                "dense_unfiltered_summary_statistics_available"
            ],
        )

        logger.info(
            "Diagnostic derived-SE valid fraction=%.4f.",
            report[
                "summary"
            ][
                "derived_se_diagnostic_valid_fraction"
            ],
        )

        logger.info(
            "Classification=%s | formal coloc=%s | "
            "ABF=%s | SuSiE=%s.",
            report[
                "summary"
            ][
                "classification"
            ],
            report[
                "summary"
            ][
                "formal_coloc_ready"
            ],
            report[
                "summary"
            ][
                "coloc_abf_ready"
            ],
            report[
                "summary"
            ][
                "coloc_susie_ready"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report