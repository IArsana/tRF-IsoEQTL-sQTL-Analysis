"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/prioritize_gwas_sumstats.py

Description:
    Runner for M5.3C.2 prostate cancer GWAS summary-statistics study
    prioritization.

    The M5.3C.1 discovery artifact contains physical Arrow list columns.
    Some pandas/PyArrow combinations cannot reconstruct pandas dtype
    metadata such as:

        list<item: string>[pyarrow]

    Therefore this runner uses a compatibility reader based directly on
    pyarrow.parquet.read_table() and converts the Arrow table to pandas while
    ignoring incompatible pandas reconstruction metadata.

    Scientific safeguards:
        - No study is downloaded automatically.
        - No sample size is inferred.
        - No ancestry is inferred.
        - Association-row counts are not interpreted as sample size.
        - Automatic keyword classification is triage only.
        - Generic prostate cancer studies require manual metadata review.
        - No locus association analysis is performed.
        - No colocalization is performed.
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

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

from pcatrfqtl.analysis.m5.gwas_sumstats_prioritization import (
    M53CProstateGWASStudyPrioritizer,
    REVIEW,
    TIER_1,
    TIER_2,
    TIER_3,
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
# Compatibility reader
# ============================================================================


def _read_discovery_parquet(
    path: Path,
) -> pd.DataFrame:
    """
    Read the M5.3C.1 discovery Parquet without incompatible pandas metadata.

    The Arrow physical schema and stored values are preserved. Only pandas
    reconstruction metadata is ignored.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"M5.3C.1 discovery Parquet not found: {path}"
        )

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    dataframe = table.to_pandas(
        ignore_metadata=True,
    )

    if dataframe.empty:

        raise RuntimeError(
            "M5.3C.1 discovery dataset contains zero rows."
        )

    return dataframe


# ============================================================================
# Runner
# ============================================================================


class M53CProstateGWASPrioritizationRunner:
    """Execute M5.3C.2 study prioritization."""

    def __init__(
        self,
        *,
        discovery_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.discovery_path = Path(
            discovery_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ----------------------------------------------------------------------
    # Output paths
    # ----------------------------------------------------------------------

    @property
    def ranked_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prostate_gwas_sumstats_prioritized.parquet"
        )

    @property
    def recommended_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prostate_gwas_sumstats_recommended.parquet"
        )

    @property
    def qc_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_2_prostate_gwas_sumstats_prioritization.json"
        )

    # ----------------------------------------------------------------------
    # Generic counters
    # ----------------------------------------------------------------------

    @staticmethod
    def _count(
        dataframe: pd.DataFrame,
        column: str,
        value: Any,
    ) -> int:
        """Count one categorical value safely."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[
                column
            ]
            .eq(
                value
            )
            .fillna(
                False
            )
            .sum()
        )

    @staticmethod
    def _true_count(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count True values safely."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[
                column
            ]
            .fillna(
                False
            )
            .eq(
                True
            )
            .sum()
        )

    # ----------------------------------------------------------------------
    # Main
    # ----------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M5.3C.2 prioritization."""

        if not self.discovery_path.exists():

            raise FileNotFoundError(
                "M5.3C.1 discovery table not found: "
                f"{self.discovery_path}"
            )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ------------------------------------------------------------------
        # Load discovery table
        # ------------------------------------------------------------------

        studies = _read_discovery_parquet(
            self.discovery_path
        )

        logger.info(
            "Loaded M5.3C.1 discovery dataset: %d studies, %d columns.",
            len(
                studies
            ),
            len(
                studies.columns
            ),
        )

        if (
            "study_accession"
            not in studies.columns
        ):

            raise ValueError(
                "M5.3C.1 discovery dataset does not contain "
                "'study_accession'."
            )

        duplicated_accessions = int(
            studies[
                "study_accession"
            ]
            .dropna()
            .duplicated()
            .sum()
        )

        if duplicated_accessions:

            logger.warning(
                "M5.3C.1 discovery table contains %d duplicated "
                "study accession rows.",
                duplicated_accessions,
            )

        # ------------------------------------------------------------------
        # Prioritize
        # ------------------------------------------------------------------

        logger.info(
            "Starting M5.3C.2 prioritization."
        )

        prioritized = (
            M53CProstateGWASStudyPrioritizer
            .prioritize(
                studies
            )
        )

        recommended = (
            prioritized.loc[
                prioritized[
                    "recommended_for_locus_retrieval"
                ]
                .fillna(
                    False
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        # ------------------------------------------------------------------
        # Persist outputs
        # ------------------------------------------------------------------

        write_parquet(
            prioritized,
            self.ranked_output_path,
            index=False,
        )

        write_parquet(
            recommended,
            self.recommended_output_path,
            index=False,
        )

        # ------------------------------------------------------------------
        # QC report
        # ------------------------------------------------------------------

        report = {
            "milestone":
                "M5.3C.2",

            "stage":
                "prostate_gwas_sumstats_prioritization",

            "input":
                str(
                    self.discovery_path
                ),

            "policy": {
                "automatic_download_performed":
                    False,

                "sample_size_inferred":
                    False,

                "ancestry_inferred":
                    False,

                "catalog_association_rows_used_as_sample_size":
                    False,

                "harmonised_sumstats_preferred":
                    True,

                "keyword_scoring_is_triage_only":
                    True,

                "generic_prostate_cancer_auto_tier1":
                    False,

                "secondary_phenotype_checked_before_tier1":
                    True,

                "manual_validation_before_download_required":
                    True,

                "nested_arrow_metadata_compatibility_reader_used":
                    True,

                "shared_parquet_reader_modified":
                    False,

                "colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "input_studies":
                    len(
                        studies
                    ),

                "unique_study_accessions":
                    int(
                        studies[
                            "study_accession"
                        ]
                        .dropna()
                        .nunique()
                    ),

                "duplicated_study_accession_rows":
                    duplicated_accessions,

                "tier_1_primary_susceptibility":
                    self._count(
                        prioritized,
                        "priority_tier",
                        TIER_1,
                    ),

                "tier_2_secondary_disease":
                    self._count(
                        prioritized,
                        "priority_tier",
                        TIER_2,
                    ),

                "tier_3_non_susceptibility":
                    self._count(
                        prioritized,
                        "priority_tier",
                        TIER_3,
                    ),

                "review_required":
                    self._count(
                        prioritized,
                        "priority_tier",
                        REVIEW,
                    ),

                "generic_prostate_cancer_review":
                    int(
                        (
                            prioritized[
                                "generic_prostate_cancer_signal"
                            ]
                            .fillna(
                                False
                            )
                            &
                            prioritized[
                                "priority_tier"
                            ]
                            .eq(
                                REVIEW
                            )
                        )
                        .sum()
                    ),

                "recommended_for_locus_retrieval":
                    self._true_count(
                        prioritized,
                        "recommended_for_locus_retrieval",
                    ),

                "harmonised_sumstats":
                    int(
                        prioritized[
                            "preferred_file_type"
                        ]
                        .astype("string")
                        .str.upper()
                        .eq(
                            "HARMONISED"
                        )
                        .fillna(
                            False
                        )
                        .sum()
                    ),
            },

            "tier_1_accessions":
                prioritized.loc[
                    prioritized[
                        "priority_tier"
                    ]
                    .eq(
                        TIER_1
                    ),
                    "study_accession",
                ]
                .dropna()
                .astype(str)
                .tolist(),

            "tier_2_accessions":
                prioritized.loc[
                    prioritized[
                        "priority_tier"
                    ]
                    .eq(
                        TIER_2
                    ),
                    "study_accession",
                ]
                .dropna()
                .astype(str)
                .tolist(),

            "review_required_accessions":
                prioritized.loc[
                    prioritized[
                        "priority_tier"
                    ]
                    .eq(
                        REVIEW
                    ),
                    "study_accession",
                ]
                .dropna()
                .astype(str)
                .tolist(),

            "generic_prostate_cancer_review_accessions":
                prioritized.loc[
                    (
                        prioritized[
                            "generic_prostate_cancer_signal"
                        ]
                        .fillna(
                            False
                        )
                        &
                        prioritized[
                            "priority_tier"
                        ]
                        .eq(
                            REVIEW
                        )
                    ),
                    "study_accession",
                ]
                .dropna()
                .astype(str)
                .tolist(),

            "recommended_accessions":
                recommended[
                    "study_accession"
                ]
                .dropna()
                .astype(str)
                .tolist(),

            "ranked_output":
                str(
                    self.ranked_output_path
                ),

            "recommended_output":
                str(
                    self.recommended_output_path
                ),
        }

        with self.qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        # ------------------------------------------------------------------
        # Logging
        # ------------------------------------------------------------------

        logger.info(
            "M5.3C.2 complete."
        )

        logger.info(
            "Tier 1: %d | Tier 2: %d | Tier 3: %d | "
            "Review: %d | Recommended: %d",
            report[
                "summary"
            ][
                "tier_1_primary_susceptibility"
            ],
            report[
                "summary"
            ][
                "tier_2_secondary_disease"
            ],
            report[
                "summary"
            ][
                "tier_3_non_susceptibility"
            ],
            report[
                "summary"
            ][
                "review_required"
            ],
            report[
                "summary"
            ][
                "recommended_for_locus_retrieval"
            ],
        )

        logger.info(
            "Generic prostate-cancer studies requiring review: %d",
            report[
                "summary"
            ][
                "generic_prostate_cancer_review"
            ],
        )

        logger.info(
            "Ranked output: %s",
            self.ranked_output_path,
        )

        logger.info(
            "Recommended output: %s",
            self.recommended_output_path,
        )

        logger.info(
            "QC output: %s",
            self.qc_path,
        )

        return report