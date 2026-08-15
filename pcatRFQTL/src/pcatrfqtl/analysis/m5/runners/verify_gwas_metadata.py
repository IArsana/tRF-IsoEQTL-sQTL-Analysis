"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/verify_gwas_metadata.py

Description:
    Runner for M5.3C.2B official GWAS Catalog metadata verification.

    This runner enriches M5.3C.2 REVIEW_REQUIRED prostate cancer studies
    using official GWAS Catalog All Studies and All Ancestry metadata.

    The verification stage evaluates:
        - official study metadata availability,
        - ancestry metadata availability,
        - case-control support,
        - prostate cancer phenotype compatibility,
        - primary susceptibility eligibility,
        - secondary disease phenotype status,
        - non-susceptibility phenotype status,
        - remaining manual-review requirements.

    Important safeguards:
        - No GWAS summary-statistics file is downloaded.
        - Sample size is never inferred from association counts.
        - Case/control evidence is derived from official ancestry metadata
          when available, with sample-description text used only as fallback.
        - Generic prostate cancer wording alone does not establish a primary
          susceptibility GWAS.
        - No candidate ranking based solely on sample size is performed.
        - No locus-level association analysis is performed.
        - No LD is recalculated.
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

from pcatrfqtl.analysis.m5.gwas_metadata_verification import (
    MANUAL_REVIEW,
    NON_SUSCEPTIBILITY,
    PRIMARY,
    SECONDARY,
    M53CMetadataVerifier,
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
# Input readers
# ============================================================================


def _read_nested_parquet(
    path: Path,
) -> pd.DataFrame:
    """
    Read M5.3C Parquet artifacts containing nested Arrow list columns.

    Some pandas/PyArrow combinations cannot reconstruct pandas extension
    metadata such as:

        list<item: string>[pyarrow]

    Therefore this reader ignores pandas reconstruction metadata while
    preserving the physical Arrow schema and values.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"Parquet file not found: {path}"
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
            f"Parquet input contains zero rows: {path}"
        )

    return dataframe


def _read_tsv(
    path: Path,
) -> pd.DataFrame:
    """
    Read one official GWAS Catalog TSV metadata artifact.

    All columns are initially retained as nullable strings so the verifier
    controls downstream numerical parsing explicitly.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"GWAS Catalog metadata file not found: {path}"
        )

    dataframe = pd.read_csv(
        path,
        sep="\t",
        dtype="string",
        low_memory=False,
    )

    if dataframe.empty:

        raise RuntimeError(
            f"GWAS Catalog metadata contains zero rows: {path}"
        )

    return dataframe


# ============================================================================
# Runner
# ============================================================================


class M53CMetadataVerificationRunner:
    """
    Execute M5.3C.2B official GWAS Catalog metadata verification.
    """

    def __init__(
        self,
        *,
        prioritized_path: str | Path,
        studies_metadata_path: str | Path,
        ancestry_metadata_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.prioritized_path = Path(
            prioritized_path
        )

        self.studies_metadata_path = Path(
            studies_metadata_path
        )

        self.ancestry_metadata_path = Path(
            ancestry_metadata_path
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
    def verified_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prostate_gwas_metadata_verified.parquet"
        )

    @property
    def primary_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prostate_gwas_primary_candidates.parquet"
        )

    @property
    def secondary_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prostate_gwas_secondary_candidates.parquet"
        )

    @property
    def manual_review_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "prostate_gwas_metadata_manual_review.parquet"
        )

    @property
    def qc_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_2b_gwas_metadata_verification.json"
        )

    # ----------------------------------------------------------------------
    # QC helpers
    # ----------------------------------------------------------------------

    @staticmethod
    def _status_count(
        dataframe: pd.DataFrame,
        status: str,
    ) -> int:
        """Count one metadata verification status."""

        if (
            "metadata_verification_status"
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                "metadata_verification_status"
            ]
            .eq(
                status
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
        """Count True values from one boolean-compatible column."""

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

    @staticmethod
    def _accessions_for_status(
        dataframe: pd.DataFrame,
        status: str,
    ) -> list[str]:
        """Return study accessions belonging to one verification status."""

        if (
            "metadata_verification_status"
            not in dataframe.columns
        ):
            return []

        return (
            dataframe.loc[
                dataframe[
                    "metadata_verification_status"
                ]
                .eq(
                    status
                ),
                "study_accession",
            ]
            .dropna()
            .astype(str)
            .tolist()
        )

    @staticmethod
    def _accessions_for_boolean(
        dataframe: pd.DataFrame,
        column: str,
    ) -> list[str]:
        """Return study accessions where one boolean flag is True."""

        if column not in dataframe.columns:
            return []

        return (
            dataframe.loc[
                dataframe[
                    column
                ]
                .fillna(
                    False
                )
                .eq(
                    True
                ),
                "study_accession",
            ]
            .dropna()
            .astype(str)
            .tolist()
        )

    # ----------------------------------------------------------------------
    # Main
    # ----------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M5.3C.2B metadata verification."""

        # ------------------------------------------------------------------
        # Validate inputs
        # ------------------------------------------------------------------

        for path in (
            self.prioritized_path,
            self.studies_metadata_path,
            self.ancestry_metadata_path,
        ):

            if not path.exists():

                raise FileNotFoundError(
                    f"Required M5.3C.2B input not found: {path}"
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
        # Load M5.3C.2 prioritization
        # ------------------------------------------------------------------

        prioritized = _read_nested_parquet(
            self.prioritized_path
        )

        required_prioritized_columns = {
            "study_accession",
            "priority_tier",
        }

        missing_prioritized_columns = (
            required_prioritized_columns
            - set(
                prioritized.columns
            )
        )

        if missing_prioritized_columns:

            raise ValueError(
                "M5.3C.2 prioritized table is missing required columns: "
                f"{sorted(missing_prioritized_columns)}"
            )

        candidates = (
            prioritized.loc[
                prioritized[
                    "priority_tier"
                ]
                .eq(
                    "REVIEW_REQUIRED"
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        logger.info(
            "M5.3C.2B review candidates: %d",
            len(
                candidates
            ),
        )

        if candidates.empty:

            raise RuntimeError(
                "M5.3C.2B found no REVIEW_REQUIRED studies to verify."
            )

        candidate_accessions = (
            candidates[
                "study_accession"
            ]
            .dropna()
            .astype("string")
            .str.upper()
        )

        duplicated_candidate_accessions = int(
            candidate_accessions
            .duplicated()
            .sum()
        )

        if duplicated_candidate_accessions:

            raise RuntimeError(
                "M5.3C.2B candidate table contains duplicated study "
                f"accessions: {duplicated_candidate_accessions}"
            )

        # ------------------------------------------------------------------
        # Load official metadata
        # ------------------------------------------------------------------

        studies_metadata = _read_tsv(
            self.studies_metadata_path
        )

        ancestry_metadata = _read_tsv(
            self.ancestry_metadata_path
        )

        logger.info(
            "Official All Studies metadata loaded: %d rows, %d columns.",
            len(
                studies_metadata
            ),
            len(
                studies_metadata.columns
            ),
        )

        logger.info(
            "Official All Ancestry metadata loaded: %d rows, %d columns.",
            len(
                ancestry_metadata
            ),
            len(
                ancestry_metadata.columns
            ),
        )

        # ------------------------------------------------------------------
        # Verify
        # ------------------------------------------------------------------

        result = M53CMetadataVerifier.verify(
            candidates=candidates,
            studies_metadata=studies_metadata,
            ancestry_metadata=ancestry_metadata,
        )

        if len(
            result
        ) != len(
            candidates
        ):

            raise RuntimeError(
                "M5.3C.2B cardinality changed during verification: "
                f"input={len(candidates)}, output={len(result)}"
            )

        logger.info(
            "Metadata verification completed for %d candidate studies.",
            len(
                result
            ),
        )

        # ------------------------------------------------------------------
        # Partition verified results
        # ------------------------------------------------------------------

        primary = (
            result.loc[
                result[
                    "metadata_verification_status"
                ]
                .eq(
                    PRIMARY
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        secondary = (
            result.loc[
                result[
                    "metadata_verification_status"
                ]
                .eq(
                    SECONDARY
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        non_susceptibility = (
            result.loc[
                result[
                    "metadata_verification_status"
                ]
                .eq(
                    NON_SUSCEPTIBILITY
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        manual_review = (
            result.loc[
                result[
                    "metadata_verification_status"
                ]
                .eq(
                    MANUAL_REVIEW
                )
            ]
            .copy()
            .reset_index(
                drop=True
            )
        )

        # ------------------------------------------------------------------
        # Persist tables
        # ------------------------------------------------------------------

        write_parquet(
            result,
            self.verified_output_path,
            index=False,
        )

        write_parquet(
            primary,
            self.primary_output_path,
            index=False,
        )

        write_parquet(
            secondary,
            self.secondary_output_path,
            index=False,
        )

        write_parquet(
            manual_review,
            self.manual_review_output_path,
            index=False,
        )

        # ------------------------------------------------------------------
        # Additional study-level metadata summaries
        # ------------------------------------------------------------------

        studies_with_ancestry = self._true_count(
            result,
            "official_ancestry_metadata_found",
        )

        studies_with_study_metadata = self._true_count(
            result,
            "official_study_metadata_found",
        )

        studies_with_any_metadata = self._true_count(
            result,
            "official_metadata_found",
        )

        studies_with_case_control_support = self._true_count(
            result,
            "case_control_supported",
        )

        # Numeric ancestry totals are descriptive provenance only.
        cases_total = None
        controls_total = None

        if (
            "ancestry_reported_cases"
            in result.columns
        ):

            cases_total_value = pd.to_numeric(
                result[
                    "ancestry_reported_cases"
                ],
                errors="coerce",
            ).sum(
                min_count=1
            )

            if pd.notna(
                cases_total_value
            ):

                cases_total = float(
                    cases_total_value
                )

        if (
            "ancestry_reported_controls"
            in result.columns
        ):

            controls_total_value = pd.to_numeric(
                result[
                    "ancestry_reported_controls"
                ],
                errors="coerce",
            ).sum(
                min_count=1
            )

            if pd.notna(
                controls_total_value
            ):

                controls_total = float(
                    controls_total_value
                )

        # ------------------------------------------------------------------
        # QC report
        # ------------------------------------------------------------------

        report = {
            "milestone":
                "M5.3C.2B",

            "stage":
                "official_gwas_metadata_verification",

            "inputs": {
                "prioritized":
                    str(
                        self.prioritized_path
                    ),

                "all_studies_metadata":
                    str(
                        self.studies_metadata_path
                    ),

                "all_ancestry_metadata":
                    str(
                        self.ancestry_metadata_path
                    ),
            },

            "policy": {
                "official_gwas_catalog_study_metadata_used":
                    True,

                "official_gwas_catalog_ancestry_metadata_used":
                    True,

                "sample_size_inferred":
                    False,

                "association_count_used_as_sample_size":
                    False,

                "case_control_numeric_metadata_preferred":
                    True,

                "sample_description_case_control_fallback_allowed":
                    True,

                "generic_prostate_cancer_auto_primary":
                    False,

                "secondary_phenotype_overrides_primary":
                    True,

                "non_susceptibility_outcome_overrides_primary":
                    True,

                "automatic_sumstats_download_performed":
                    False,

                "candidate_ranking_performed":
                    False,

                "locus_association_analysis_performed":
                    False,

                "ld_recalculated":
                    False,

                "colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "input_review_candidates":
                    len(
                        candidates
                    ),

                "output_verified_studies":
                    len(
                        result
                    ),

                "official_study_metadata_found":
                    studies_with_study_metadata,

                "official_ancestry_metadata_found":
                    studies_with_ancestry,

                "official_any_metadata_found":
                    studies_with_any_metadata,

                "studies_with_case_control_support":
                    studies_with_case_control_support,

                "primary_susceptibility_candidates":
                    len(
                        primary
                    ),

                "secondary_disease_phenotypes":
                    len(
                        secondary
                    ),

                "non_susceptibility":
                    len(
                        non_susceptibility
                    ),

                "manual_review":
                    len(
                        manual_review
                    ),

                "descriptive_sum_ancestry_reported_cases":
                    cases_total,

                "descriptive_sum_ancestry_reported_controls":
                    controls_total,
            },

            "primary_candidate_accessions":
                self._accessions_for_status(
                    result,
                    PRIMARY,
                ),

            "secondary_candidate_accessions":
                self._accessions_for_status(
                    result,
                    SECONDARY,
                ),

            "non_susceptibility_accessions":
                self._accessions_for_status(
                    result,
                    NON_SUSCEPTIBILITY,
                ),

            "manual_review_accessions":
                self._accessions_for_status(
                    result,
                    MANUAL_REVIEW,
                ),

            "case_control_supported_accessions":
                self._accessions_for_boolean(
                    result,
                    "case_control_supported",
                ),

            "missing_study_metadata_accessions":
                (
                    result.loc[
                        ~result[
                            "official_study_metadata_found"
                        ]
                        .fillna(
                            False
                        ),
                        "study_accession",
                    ]
                    .dropna()
                    .astype(str)
                    .tolist()
                ),

            "missing_ancestry_metadata_accessions":
                (
                    result.loc[
                        ~result[
                            "official_ancestry_metadata_found"
                        ]
                        .fillna(
                            False
                        ),
                        "study_accession",
                    ]
                    .dropna()
                    .astype(str)
                    .tolist()
                ),

            "outputs": {
                "verified":
                    str(
                        self.verified_output_path
                    ),

                "primary_candidates":
                    str(
                        self.primary_output_path
                    ),

                "secondary_candidates":
                    str(
                        self.secondary_output_path
                    ),

                "manual_review":
                    str(
                        self.manual_review_output_path
                    ),

                "qc":
                    str(
                        self.qc_path
                    ),
            },
        }

        # ------------------------------------------------------------------
        # Cardinality safeguard
        # ------------------------------------------------------------------

        classified_total = (
            len(
                primary
            )
            + len(
                secondary
            )
            + len(
                non_susceptibility
            )
            + len(
                manual_review
            )
        )

        if classified_total != len(
            result
        ):

            raise RuntimeError(
                "M5.3C.2B verification statuses do not partition all rows: "
                f"classified={classified_total}, result={len(result)}"
            )

        # ------------------------------------------------------------------
        # Write QC
        # ------------------------------------------------------------------

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
            "M5.3C.2B complete."
        )

        logger.info(
            "Primary: %d | Secondary: %d | "
            "Non-susceptibility: %d | Manual review: %d",
            len(
                primary
            ),
            len(
                secondary
            ),
            len(
                non_susceptibility
            ),
            len(
                manual_review
            ),
        )

        logger.info(
            "Case-control supported studies: %d",
            studies_with_case_control_support,
        )

        logger.info(
            "Official study metadata coverage: %d/%d",
            studies_with_study_metadata,
            len(
                result
            ),
        )

        logger.info(
            "Official ancestry metadata coverage: %d/%d",
            studies_with_ancestry,
            len(
                result
            ),
        )

        logger.info(
            "Verified output: %s",
            self.verified_output_path,
        )

        logger.info(
            "Primary candidates: %s",
            self.primary_output_path,
        )

        logger.info(
            "QC output: %s",
            self.qc_path,
        )

        return report