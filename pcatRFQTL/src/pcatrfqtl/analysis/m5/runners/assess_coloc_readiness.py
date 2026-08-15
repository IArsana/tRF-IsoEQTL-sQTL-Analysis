"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/assess_coloc_readiness.py

Description:
    Runner for M5.5 colocalization readiness assessment.

    The runner evaluates whether disease-side GWAS summary statistics and
    Moradi regulatory-QTL evidence are suitable for formal colocalization.

    Moradi QTL input:
        M4.1 is stored as partitioned Parquet artifacts under:

            data/processed/m4/moradi_qtl_index/

        All partitions are loaded and concatenated in memory while retaining
        explicit source_partition provenance.

    Outputs:
        gwas_coloc_readiness.parquet
        qtl_coloc_readiness.parquet
        coloc_pair_readiness.parquet
        m5_5_coloc_readiness.json

    Scientific safeguards:
        - M5.5 evaluates readiness only.
        - No colocalization is executed.
        - No fine-mapping is executed.
        - No causal inference is performed.
        - Cross-build coordinate joins are not performed.
        - Significant-only QTL evidence is not treated as dense locus-level
          summary statistics.
        - Missing QTL variants are not interpreted as null associations.

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
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
import yaml

from pcatrfqtl.analysis.m5.coloc_readiness import (
    assess_gwas_readiness,
    assess_pair_readiness,
    assess_qtl_readiness,
)
from pcatrfqtl.io.parquet import (
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


# ============================================================================
# Logging
# ============================================================================


logger = get_logger(
    __name__
)


# ============================================================================
# Parquet readers
# ============================================================================


def _read_parquet_compat(
    path: Path,
) -> pd.DataFrame:
    """
    Read one Parquet artifact without pandas metadata reconstruction.

    M5 uses the PyArrow compatibility reader because nested/list columns can
    fail when reconstructed through pandas metadata.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"Input Parquet not found: {path}"
        )

    if not path.is_file():

        raise FileNotFoundError(
            f"Expected Parquet file, received: {path}"
        )

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    dataframe = table.to_pandas(
        ignore_metadata=True,
    )

    return dataframe


def _read_parquet_directory_compat(
    directory: Path,
) -> tuple[
    pd.DataFrame,
    list[dict[str, Any]],
]:
    """
    Load all M4.1 Moradi QTL Parquet partitions.

    Each input partition receives an explicit source_partition column using
    the filename stem.

    Returns:
        dataframe:
            Concatenated Moradi QTL dataframe.

        partition_metadata:
            Row counts and column information for every source partition.
    """

    if not directory.exists():

        raise FileNotFoundError(
            "Moradi QTL index directory not found: "
            f"{directory}"
        )

    if not directory.is_dir():

        raise NotADirectoryError(
            "Expected Moradi QTL index directory, received: "
            f"{directory}"
        )

    paths = sorted(
        directory.glob(
            "*.parquet"
        )
    )

    if not paths:

        raise FileNotFoundError(
            "No Parquet partitions found under Moradi QTL index directory: "
            f"{directory}"
        )

    frames: list[
        pd.DataFrame
    ] = []

    partition_metadata: list[
        dict[str, Any]
    ] = []

    for path in paths:

        frame = _read_parquet_compat(
            path
        )

        frame = frame.copy()

        partition_name = (
            path.stem
        )

        # --------------------------------------------------------------
        # Preserve partition provenance explicitly.
        #
        # If a previous source_partition column somehow exists, avoid
        # silently overwriting it.
        # --------------------------------------------------------------

        if "source_partition" in frame.columns:

            existing = (
                frame[
                    "source_partition"
                ]
                .astype("string")
                .dropna()
                .unique()
                .tolist()
            )

            if existing:

                logger.warning(
                    "Partition %s already contains source_partition=%s; "
                    "the runner will preserve the original values.",
                    partition_name,
                    existing,
                )

            else:

                frame[
                    "source_partition"
                ] = partition_name

        else:

            frame[
                "source_partition"
            ] = partition_name

        frames.append(
            frame
        )

        partition_metadata.append(
            {
                "source_partition":
                    partition_name,

                "path":
                    str(
                        path
                    ),

                "rows":
                    int(
                        len(
                            frame
                        )
                    ),

                "columns":
                    [
                        str(
                            column
                        )
                        for column
                        in frame.columns
                    ],
            }
        )

        logger.info(
            "Loaded Moradi QTL partition %s: %d rows, %d columns.",
            partition_name,
            len(
                frame
            ),
            len(
                frame.columns
            ),
        )

    dataframe = pd.concat(
        frames,
        ignore_index=True,
        sort=False,
    )

    logger.info(
        "Loaded %d Moradi QTL partitions with %d total rows.",
        len(
            paths
        ),
        len(
            dataframe
        ),
    )

    return (
        dataframe,
        partition_metadata,
    )


# ============================================================================
# Utility helpers
# ============================================================================


def _load_yaml(
    path: Path,
) -> dict[str, Any]:
    """Load YAML config."""

    if not path.exists():

        raise FileNotFoundError(
            f"M5.5 config not found: {path}"
        )

    with path.open(
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
            "M5.5 configuration must contain a YAML mapping."
        )

    return config


def _value_counts(
    dataframe: pd.DataFrame,
    column: str,
) -> dict[str, int]:
    """Return safe value counts for one dataframe column."""

    if dataframe.empty:

        return {}

    if column not in dataframe.columns:

        return {}

    counts = (
        dataframe[
            column
        ]
        .astype("string")
        .fillna("MISSING")
        .value_counts(
            dropna=False
        )
    )

    return {
        str(
            key
        ):
            int(
                value
            )
        for key, value
        in counts.items()
    }


def _boolean_count(
    dataframe: pd.DataFrame,
    column: str,
) -> int:
    """Count truthy values safely."""

    if dataframe.empty:

        return 0

    if column not in dataframe.columns:

        return 0

    return int(
        dataframe[
            column
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )


# ============================================================================
# Runner
# ============================================================================


class M55ColocReadinessRunner:
    """Execute M5.5 colocalization readiness assessment."""

    def __init__(
        self,
        *,
        gwas_path: str | Path,
        qtl_directory: str | Path,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.gwas_path = Path(
            gwas_path
        )

        self.qtl_directory = Path(
            qtl_directory
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
    def gwas_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_coloc_readiness.parquet"
        )

    @property
    def qtl_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "qtl_coloc_readiness.parquet"
        )

    @property
    def pair_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "coloc_pair_readiness.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_5_coloc_readiness.json"
        )

    # ======================================================================
    # Configuration validation
    # ======================================================================

    def _validate_config(
        self,
        config: dict[str, Any],
    ) -> None:
        """Validate minimum M5.5 configuration structure."""

        required_top_level = {
            "milestone",
            "gwas",
            "qtl",
            "readiness_policy",
        }

        missing = (
            required_top_level
            - set(
                config
            )
        )

        if missing:

            raise ValueError(
                "M5.5 config missing required sections: "
                f"{sorted(missing)}"
            )

        required_qtl = {
            "source_scope",
            "full_summary_statistics_available",
        }

        missing_qtl = (
            required_qtl
            - set(
                config[
                    "qtl"
                ]
            )
        )

        if missing_qtl:

            raise ValueError(
                "M5.5 QTL config missing required fields: "
                f"{sorted(missing_qtl)}"
            )

        required_policy = {
            "minimum_shared_variants_descriptive",
        }

        missing_policy = (
            required_policy
            - set(
                config[
                    "readiness_policy"
                ]
            )
        )

        if missing_policy:

            raise ValueError(
                "M5.5 readiness policy missing required fields: "
                f"{sorted(missing_policy)}"
            )

    # ======================================================================
    # Main execution
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute the full M5.5 readiness assessment."""

        # ------------------------------------------------------------------
        # Configuration.
        # ------------------------------------------------------------------

        config = _load_yaml(
            self.config_path
        )

        self._validate_config(
            config
        )

        # ------------------------------------------------------------------
        # Source loading.
        # ------------------------------------------------------------------

        gwas = _read_parquet_compat(
            self.gwas_path
        )

        (
            qtl,
            qtl_partition_metadata,
        ) = _read_parquet_directory_compat(
            self.qtl_directory
        )

        logger.info(
            "Loaded GWAS standardized rows: %d.",
            len(
                gwas
            ),
        )

        logger.info(
            "Loaded Moradi QTL index rows: %d.",
            len(
                qtl
            ),
        )

        if gwas.empty:

            raise RuntimeError(
                "M5.5 GWAS input is empty."
            )

        if qtl.empty:

            raise RuntimeError(
                "M5.5 Moradi QTL input is empty."
            )

        # ------------------------------------------------------------------
        # Output directories.
        # ------------------------------------------------------------------

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ==================================================================
        # Disease-side readiness
        # ==================================================================

        logger.info(
            "Assessing disease-side GWAS colocalization readiness."
        )

        gwas_result = assess_gwas_readiness(
            gwas
        )

        gwas_readiness = (
            gwas_result.dataframe
        )

        # ==================================================================
        # Regulatory-side readiness
        # ==================================================================

        logger.info(
            "Assessing Moradi QTL colocalization readiness."
        )

        qtl_result = assess_qtl_readiness(
            qtl,
            source_scope=str(
                config[
                    "qtl"
                ][
                    "source_scope"
                ]
            ),
            full_summary_statistics_available=bool(
                config[
                    "qtl"
                ][
                    "full_summary_statistics_available"
                ]
            ),
        )

        qtl_readiness = (
            qtl_result.dataframe
        )

        # ==================================================================
        # GWAS × QTL pair-level readiness
        # ==================================================================

        logger.info(
            "Assessing GWAS × QTL candidate-level readiness."
        )

        pair_result = assess_pair_readiness(
            gwas_locus=gwas,
            qtl_rows=qtl,
            minimum_shared_variants=int(
                config[
                    "readiness_policy"
                ][
                    "minimum_shared_variants_descriptive"
                ]
            ),
        )

        pair_readiness = (
            pair_result.dataframe
        )

        # ==================================================================
        # Persist outputs
        # ==================================================================

        write_parquet(
            gwas_readiness,
            self.gwas_output,
            index=False,
        )

        write_parquet(
            qtl_readiness,
            self.qtl_output,
            index=False,
        )

        write_parquet(
            pair_readiness,
            self.pair_output,
            index=False,
        )

        # ==================================================================
        # Summary calculations
        # ==================================================================

        formal_ready_pairs = (
            _boolean_count(
                pair_readiness,
                "formal_coloc_ready",
            )
        )

        coloc_abf_ready_pairs = (
            _boolean_count(
                pair_readiness,
                "coloc_abf_ready",
            )
        )

        coloc_susie_ready_pairs = (
            _boolean_count(
                pair_readiness,
                "coloc_susie_ready",
            )
        )

        gwas_basic_ready = (
            _boolean_count(
                gwas_readiness,
                "basic_coloc_vector_ready",
            )
        )

        qtl_basic_ready = (
            _boolean_count(
                qtl_readiness,
                "basic_coloc_vector_ready",
            )
        )

        qtl_full_locus_ready = (
            _boolean_count(
                qtl_readiness,
                "full_locus_summary_available",
            )
        )

        qtl_dense_locus_ready = (
            _boolean_count(
                qtl_readiness,
                "dense_locus_available",
            )
        )

        shared_variant_total = 0

        shared_variant_max = 0

        if (
            not pair_readiness.empty
            and
            "shared_variant_count"
            in pair_readiness.columns
        ):

            shared_series = pd.to_numeric(
                pair_readiness[
                    "shared_variant_count"
                ],
                errors="coerce",
            ).fillna(
                0
            )

            shared_variant_total = int(
                shared_series.sum()
            )

            shared_variant_max = int(
                shared_series.max()
            )

        # ==================================================================
        # QC report
        # ==================================================================

        report = {
            "milestone":
                "M5.5",

            "stage":
                "colocalization_readiness_assessment",

            "policy": {
                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,

                "cross_build_coordinate_join_performed":
                    False,

                "cross_build_matching_method":
                    (
                        config[
                            "readiness_policy"
                        ].get(
                            "cross_build_variant_matching",
                            "canonical_rsid_only",
                        )
                    ),

                "gwas_source_scope":
                    config[
                        "gwas"
                    ].get(
                        "source_scope"
                    ),

                "qtl_source":
                    config[
                        "qtl"
                    ].get(
                        "source"
                    ),

                "qtl_source_scope":
                    config[
                        "qtl"
                    ][
                        "source_scope"
                    ],

                "qtl_full_summary_statistics_available":
                    bool(
                        config[
                            "qtl"
                        ][
                            "full_summary_statistics_available"
                        ]
                    ),

                "dense_variant_coverage_required":
                    bool(
                        config[
                            "readiness_policy"
                        ].get(
                            "require_dense_variant_coverage",
                            True,
                        )
                    ),

                "variant_overlap_required":
                    bool(
                        config[
                            "readiness_policy"
                        ].get(
                            "require_variant_overlap",
                            True,
                        )
                    ),

                "minimum_shared_variants_descriptive":
                    int(
                        config[
                            "readiness_policy"
                        ][
                            "minimum_shared_variants_descriptive"
                        ]
                    ),

                "ld_matrix_required_for_susie":
                    bool(
                        config[
                            "readiness_policy"
                        ].get(
                            "require_ld_matrix_for_susie",
                            True,
                        )
                    ),

                "significant_only_qtl_treated_as_dense_locus":
                    False,
            },

            "inputs": {
                "gwas":
                    str(
                        self.gwas_path
                    ),

                "qtl_directory":
                    str(
                        self.qtl_directory
                    ),

                "qtl_partition_count":
                    int(
                        len(
                            qtl_partition_metadata
                        )
                    ),
            },

            "summary": {
                "gwas_rows":
                    int(
                        len(
                            gwas
                        )
                    ),

                "qtl_rows":
                    int(
                        len(
                            qtl
                        )
                    ),

                "gwas_loci_assessed":
                    int(
                        len(
                            gwas_readiness
                        )
                    ),

                "qtl_units_assessed":
                    int(
                        len(
                            qtl_readiness
                        )
                    ),

                "study_lead_pairs_assessed":
                    int(
                        len(
                            pair_readiness
                        )
                    ),

                "gwas_basic_coloc_vector_ready_units":
                    gwas_basic_ready,

                "qtl_basic_coloc_vector_ready_units":
                    qtl_basic_ready,

                "qtl_full_locus_summary_units":
                    qtl_full_locus_ready,

                "qtl_dense_locus_units":
                    qtl_dense_locus_ready,

                "formal_coloc_ready_pairs":
                    formal_ready_pairs,

                "formal_coloc_blocked_pairs":
                    int(
                        len(
                            pair_readiness
                        )
                        -
                        formal_ready_pairs
                    ),

                "coloc_abf_ready_pairs":
                    coloc_abf_ready_pairs,

                "coloc_susie_ready_pairs":
                    coloc_susie_ready_pairs,

                "shared_variant_observations_total":
                    shared_variant_total,

                "maximum_shared_variants_per_pair":
                    shared_variant_max,
            },

            "gwas_status_counts":
                _value_counts(
                    gwas_readiness,
                    "formal_coloc_data_status",
                ),

            "qtl_status_counts":
                _value_counts(
                    qtl_readiness,
                    "formal_coloc_data_status",
                ),

            "pair_status_counts":
                _value_counts(
                    pair_readiness,
                    "readiness_status",
                ),

            "qtl_partitions":
                qtl_partition_metadata,

            "outputs": {
                "gwas_readiness":
                    str(
                        self.gwas_output
                    ),

                "qtl_readiness":
                    str(
                        self.qtl_output
                    ),

                "pair_readiness":
                    str(
                        self.pair_output
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

        # ==================================================================
        # Logging
        # ==================================================================

        logger.info(
            "M5.5 complete."
        )

        logger.info(
            "GWAS readiness units: %d | basic vector-ready: %d.",
            len(
                gwas_readiness
            ),
            gwas_basic_ready,
        )

        logger.info(
            "QTL readiness units: %d | basic vector-ready: %d | "
            "full-locus: %d | dense-locus: %d.",
            len(
                qtl_readiness
            ),
            qtl_basic_ready,
            qtl_full_locus_ready,
            qtl_dense_locus_ready,
        )

        logger.info(
            "Formal coloc-ready pairs: %d/%d.",
            formal_ready_pairs,
            len(
                pair_readiness
            ),
        )

        logger.info(
            "coloc.abf-ready pairs: %d | coloc.susie-ready pairs: %d.",
            coloc_abf_ready_pairs,
            coloc_susie_ready_pairs,
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report