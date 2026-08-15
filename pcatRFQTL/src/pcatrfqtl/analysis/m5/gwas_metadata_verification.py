"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_metadata_verification.py

Description:
    M5.3C.2B study-level metadata verification for prostate cancer GWAS
    summary-statistics candidates.

    This module enriches M5.3C.2 review candidates using official
    GWAS Catalog study and ancestry metadata.

    The goal is to distinguish primary prostate cancer susceptibility GWAS
    candidates from secondary disease phenotypes, non-susceptibility
    phenotypes, and studies requiring manual review.

    Important safeguards:
        - Sample size is read from official GWAS Catalog metadata.
        - Association counts are never interpreted as sample size.
        - Case/control evidence preferentially uses official NUMBER OF CASES
          and NUMBER OF CONTROLS fields from ancestry metadata when available.
        - Sample descriptions are used only as supporting evidence.
        - Generic "prostate cancer" wording alone does not establish a
          susceptibility study.
        - No summary-statistics file is downloaded automatically.
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

import re
from typing import Any

import numpy as np
import pandas as pd


PRIMARY = "PRIMARY_SUSCEPTIBILITY_CANDIDATE"
SECONDARY = "SECONDARY_DISEASE_PHENOTYPE"
NON_SUSCEPTIBILITY = "NON_SUSCEPTIBILITY"
MANUAL_REVIEW = "MANUAL_REVIEW"


# ============================================================================
# Patterns
# ============================================================================

PROSTATE_CANCER_PATTERN = re.compile(
    r"\bprostate cancer\b|\bprostate carcinoma\b",
    flags=re.IGNORECASE,
)

SECONDARY_PATTERN = re.compile(
    (
        r"\baggressive\b"
        r"|\badvanced\b"
        r"|\bmetastatic\b"
        r"|\blethal\b"
        r"|\bfatal\b"
        r"|\bhigh[- ]risk\b"
        r"|\bgleason\b"
        r"|\bstage\b"
        r"|\bclinically significant\b"
    ),
    flags=re.IGNORECASE,
)

NON_SUSCEPTIBILITY_PATTERN = re.compile(
    (
        r"\bsurvival\b"
        r"|\bprognos"
        r"|\brecurrence\b"
        r"|\bmortality\b"
        r"|\btreatment response\b"
        r"|\btherapy response\b"
        r"|\bresistance\b"
        r"|\btoxicity\b"
    ),
    flags=re.IGNORECASE,
)

PSA_PATTERN = re.compile(
    r"\bpsa\b|\bprostate[- ]specific antigen\b",
    flags=re.IGNORECASE,
)

CASE_PATTERN = re.compile(
    r"\bcase[s]?\b",
    flags=re.IGNORECASE,
)

CONTROL_PATTERN = re.compile(
    r"\bcontrol[s]?\b",
    flags=re.IGNORECASE,
)


# ============================================================================
# Generic helpers
# ============================================================================


def _is_missing(
    value: Any,
) -> bool:
    """Safely identify missing scalar metadata."""

    if value is None:
        return True

    if isinstance(
        value,
        (
            list,
            tuple,
            set,
            np.ndarray,
        ),
    ):
        return False

    try:
        missing = pd.isna(
            value
        )

    except (
        TypeError,
        ValueError,
    ):
        return False

    if isinstance(
        missing,
        (
            bool,
            np.bool_,
        ),
    ):
        return bool(
            missing
        )

    return False


def _text(
    value: Any,
) -> str:
    """Convert scalar/list metadata into searchable text."""

    if _is_missing(
        value
    ):
        return ""

    if isinstance(
        value,
        np.ndarray,
    ):
        values = value.tolist()

    elif isinstance(
        value,
        (
            list,
            tuple,
            set,
        ),
    ):
        values = list(
            value
        )

    else:
        return str(
            value
        )

    return " | ".join(
        str(item)
        for item
        in values
        if not _is_missing(
            item
        )
    )


def _normalize_column_name(
    value: Any,
) -> str:
    """Normalize metadata column names for robust matching."""

    text = str(
        value
    )

    text = text.replace(
        "\ufeff",
        "",
    )

    text = text.replace(
        "_",
        " ",
    )

    text = re.sub(
        r"\s+",
        " ",
        text,
    )

    return (
        text
        .strip()
        .lower()
    )


def _find_column(
    dataframe: pd.DataFrame,
    candidates: tuple[str, ...],
) -> str | None:
    """Resolve a metadata column using normalized names."""

    lookup = {
        _normalize_column_name(
            column
        ):
            column
        for column
        in dataframe.columns
    }

    for candidate in candidates:

        resolved = lookup.get(
            _normalize_column_name(
                candidate
            )
        )

        if resolved is not None:
            return resolved

    return None


def _numeric(
    value: Any,
) -> float | None:
    """Convert one metadata value to numeric when possible."""

    if _is_missing(
        value
    ):
        return None

    numeric = pd.to_numeric(
        pd.Series(
            [
                value
            ]
        ),
        errors="coerce",
    ).iloc[0]

    if pd.isna(
        numeric
    ):
        return None

    return float(
        numeric
    )


# ============================================================================
# Verifier
# ============================================================================


class M53CMetadataVerifier:
    """Verify M5.3C study candidates using official GWAS metadata."""

    # ------------------------------------------------------------------
    # Official All Studies schema
    # ------------------------------------------------------------------

    STUDY_ACCESSION_COLUMNS = (
        "STUDY ACCESSION",
        "study_accession",
    )

    REPORTED_TRAIT_COLUMNS = (
        "DISEASE/TRAIT",
        "reported_trait",
    )

    TITLE_COLUMNS = (
        "STUDY",
        "study_title",
    )

    INITIAL_SAMPLE_SIZE_COLUMNS = (
        "INITIAL SAMPLE SIZE",
        "initial_sample_size",
    )

    REPLICATION_SAMPLE_SIZE_COLUMNS = (
        "REPLICATION SAMPLE SIZE",
        "replication_sample_size",
    )

    PMID_COLUMNS = (
        "PUBMED ID",
        "PUBMEDID",
        "pmid",
    )

    FULL_SUMSTATS_COLUMNS = (
        "FULL SUMMARY STATISTICS",
        "full_summary_statistics",
    )

    SUMSTATS_LOCATION_COLUMNS = (
        "SUMMARY STATS LOCATION",
        "summary_stats_location",
    )

    STATISTICAL_MODEL_COLUMNS = (
        "STATISTICAL MODEL",
        "statistical_model",
    )

    BACKGROUND_TRAIT_COLUMNS = (
        "BACKGROUND TRAIT",
        "background_trait",
    )

    COHORT_COLUMNS = (
        "COHORT",
        "cohort",
    )

    # ------------------------------------------------------------------
    # Normalize All Studies metadata
    # ------------------------------------------------------------------

    @classmethod
    def normalize_study_metadata(
        cls,
        metadata: pd.DataFrame,
    ) -> pd.DataFrame:
        """Normalize official GWAS Catalog All Studies metadata."""

        accession_col = _find_column(
            metadata,
            cls.STUDY_ACCESSION_COLUMNS,
        )

        if accession_col is None:

            raise ValueError(
                "Could not identify GWAS Catalog study accession column. "
                f"Observed columns: {list(metadata.columns)}"
            )

        mappings = {
            "official_reported_trait":
                cls.REPORTED_TRAIT_COLUMNS,

            "official_study_title":
                cls.TITLE_COLUMNS,

            "official_initial_sample_size":
                cls.INITIAL_SAMPLE_SIZE_COLUMNS,

            "official_replication_sample_size":
                cls.REPLICATION_SAMPLE_SIZE_COLUMNS,

            "official_pmid":
                cls.PMID_COLUMNS,

            "official_full_summary_statistics":
                cls.FULL_SUMSTATS_COLUMNS,

            "official_summary_stats_location":
                cls.SUMSTATS_LOCATION_COLUMNS,

            "official_statistical_model":
                cls.STATISTICAL_MODEL_COLUMNS,

            "official_background_trait":
                cls.BACKGROUND_TRAIT_COLUMNS,

            "official_cohort":
                cls.COHORT_COLUMNS,
        }

        result = pd.DataFrame(
            index=metadata.index
        )

        result[
            "study_accession"
        ] = (
            metadata[
                accession_col
            ]
            .astype("string")
            .str.strip()
            .str.upper()
        )

        for target, candidates in mappings.items():

            source = _find_column(
                metadata,
                candidates,
            )

            if source is None:

                result[
                    target
                ] = pd.Series(
                    pd.NA,
                    index=metadata.index,
                    dtype="string",
                )

            else:

                result[
                    target
                ] = metadata[
                    source
                ].astype(
                    "string"
                )

        return result

    # ------------------------------------------------------------------
    # Ancestry metadata
    # ------------------------------------------------------------------

    @staticmethod
    def build_ancestry_summary(
        ancestry: pd.DataFrame,
    ) -> pd.DataFrame:
        """Aggregate official ancestry records by study accession."""

        accession_col = _find_column(
            ancestry,
            (
                "STUDY ACCESSION",
                "study_accession",
            ),
        )

        if accession_col is None:

            raise ValueError(
                "Could not identify study accession in ancestry metadata. "
                f"Observed columns: {list(ancestry.columns)}"
            )

        initial_description_col = _find_column(
            ancestry,
            (
                "INITIAL SAMPLE DESCRIPTION",
                "initial_sample_description",
            ),
        )

        replication_description_col = _find_column(
            ancestry,
            (
                "REPLICATION SAMPLE DESCRIPTION",
                "replication_sample_description",
            ),
        )

        stage_col = _find_column(
            ancestry,
            (
                "STAGE",
                "stage",
            ),
        )

        individuals_col = _find_column(
            ancestry,
            (
                "NUMBER OF INDIVIDUALS",
                "number_of_individuals",
            ),
        )

        ancestry_col = _find_column(
            ancestry,
            (
                "BROAD ANCESTRAL CATEGORY",
                "broad_ancestral_category",
            ),
        )

        cases_col = _find_column(
            ancestry,
            (
                "NUMBER OF CASES",
                "number_of_cases",
            ),
        )

        controls_col = _find_column(
            ancestry,
            (
                "NUMBER OF CONTROLS",
                "number_of_controls",
            ),
        )

        sample_description_col = _find_column(
            ancestry,
            (
                "SAMPLE DESCRIPTION",
                "sample_description",
            ),
        )

        working = pd.DataFrame(
            index=ancestry.index
        )

        working[
            "study_accession"
        ] = (
            ancestry[
                accession_col
            ]
            .astype("string")
            .str.strip()
            .str.upper()
        )

        working[
            "stage"
        ] = (
            ancestry[
                stage_col
            ].astype("string")
            if stage_col
            else pd.Series(
                pd.NA,
                index=ancestry.index,
                dtype="string",
            )
        )

        working[
            "ancestry"
        ] = (
            ancestry[
                ancestry_col
            ].astype("string")
            if ancestry_col
            else pd.Series(
                pd.NA,
                index=ancestry.index,
                dtype="string",
            )
        )

        working[
            "number_of_individuals"
        ] = (
            pd.to_numeric(
                ancestry[
                    individuals_col
                ],
                errors="coerce",
            )
            if individuals_col
            else np.nan
        )

        working[
            "number_of_cases"
        ] = (
            pd.to_numeric(
                ancestry[
                    cases_col
                ],
                errors="coerce",
            )
            if cases_col
            else np.nan
        )

        working[
            "number_of_controls"
        ] = (
            pd.to_numeric(
                ancestry[
                    controls_col
                ],
                errors="coerce",
            )
            if controls_col
            else np.nan
        )

        working[
            "initial_sample_description"
        ] = (
            ancestry[
                initial_description_col
            ].astype("string")
            if initial_description_col
            else pd.Series(
                pd.NA,
                index=ancestry.index,
                dtype="string",
            )
        )

        working[
            "replication_sample_description"
        ] = (
            ancestry[
                replication_description_col
            ].astype("string")
            if replication_description_col
            else pd.Series(
                pd.NA,
                index=ancestry.index,
                dtype="string",
            )
        )

        working[
            "sample_description"
        ] = (
            ancestry[
                sample_description_col
            ].astype("string")
            if sample_description_col
            else pd.Series(
                pd.NA,
                index=ancestry.index,
                dtype="string",
            )
        )

        grouped = (
            working
            .groupby(
                "study_accession",
                dropna=False,
            )
            .agg(
                ancestry_categories=(
                    "ancestry",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                ancestry_stages=(
                    "stage",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                ancestry_metadata_rows=(
                    "study_accession",
                    "size",
                ),
                ancestry_reported_n=(
                    "number_of_individuals",
                    "sum",
                ),
                ancestry_reported_cases=(
                    "number_of_cases",
                    "sum",
                ),
                ancestry_reported_controls=(
                    "number_of_controls",
                    "sum",
                ),
                ancestry_initial_descriptions=(
                    "initial_sample_description",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                ancestry_replication_descriptions=(
                    "replication_sample_description",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                ancestry_sample_descriptions=(
                    "sample_description",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
            )
            .reset_index()
        )

        return grouped

    # ------------------------------------------------------------------
    # Classification
    # ------------------------------------------------------------------

    @classmethod
    def classify_row(
        cls,
        row: pd.Series,
    ) -> tuple[
        str,
        bool,
        bool,
        str,
    ]:
        """Classify one enriched GWAS study."""

        trait = _text(
            row.get(
                "official_reported_trait"
            )
        )

        title = _text(
            row.get(
                "official_study_title"
            )
        )

        initial_description = _text(
            row.get(
                "ancestry_initial_descriptions"
            )
        )

        replication_description = _text(
            row.get(
                "ancestry_replication_descriptions"
            )
        )

        sample_description = _text(
            row.get(
                "ancestry_sample_descriptions"
            )
        )

        combined = " | ".join(
            (
                trait,
                title,
                initial_description,
                replication_description,
                sample_description,
            )
        )

        prostate = bool(
            PROSTATE_CANCER_PATTERN.search(
                combined
            )
        )

        secondary = bool(
            SECONDARY_PATTERN.search(
                combined
            )
        )

        non_susceptibility = bool(
            NON_SUSCEPTIBILITY_PATTERN.search(
                combined
            )
        )

        psa = bool(
            PSA_PATTERN.search(
                combined
            )
        )

        cases = _numeric(
            row.get(
                "ancestry_reported_cases"
            )
        )

        controls = _numeric(
            row.get(
                "ancestry_reported_controls"
            )
        )

        numeric_case_control_supported = (
            cases is not None
            and controls is not None
            and cases > 0
            and controls > 0
        )

        text_case_control_supported = (
            bool(
                CASE_PATTERN.search(
                    initial_description
                )
            )
            and bool(
                CONTROL_PATTERN.search(
                    initial_description
                )
            )
        )

        case_control_supported = (
            numeric_case_control_supported
            or text_case_control_supported
        )

        rationale: list[str] = []

        if prostate:

            rationale.append(
                "prostate cancer phenotype supported"
            )

        if numeric_case_control_supported:

            rationale.append(
                "official ancestry metadata reports cases and controls"
            )

        elif text_case_control_supported:

            rationale.append(
                "sample description contains cases and controls"
            )

        if secondary:

            rationale.append(
                "secondary disease phenotype detected"
            )

        if non_susceptibility:

            rationale.append(
                "non-susceptibility outcome detected"
            )

        if psa:

            rationale.append(
                "PSA-related phenotype detected"
            )

        if non_susceptibility:

            status = NON_SUSCEPTIBILITY

        elif secondary:

            status = SECONDARY

        elif (
            prostate
            and case_control_supported
            and not psa
        ):

            status = PRIMARY

        else:

            status = MANUAL_REVIEW

        primary_candidate = (
            status
            == PRIMARY
        )

        return (
            status,
            case_control_supported,
            primary_candidate,
            "; ".join(
                rationale
            ),
        )

    # ------------------------------------------------------------------
    # Main verification
    # ------------------------------------------------------------------

    @classmethod
    def verify(
        cls,
        *,
        candidates: pd.DataFrame,
        studies_metadata: pd.DataFrame,
        ancestry_metadata: pd.DataFrame,
    ) -> pd.DataFrame:
        """Enrich and verify M5.3C.2 review candidates."""

        review_accessions = set(
            candidates[
                "study_accession"
            ]
            .dropna()
            .astype("string")
            .str.upper()
            .tolist()
        )

        study_metadata = cls.normalize_study_metadata(
            studies_metadata
        )

        study_metadata = (
            study_metadata.loc[
                study_metadata[
                    "study_accession"
                ]
                .isin(
                    review_accessions
                )
            ]
            .copy()
            .drop_duplicates(
                subset=[
                    "study_accession"
                ],
                keep="first",
            )
        )

        ancestry_summary = cls.build_ancestry_summary(
            ancestry_metadata
        )

        ancestry_summary = (
            ancestry_summary.loc[
                ancestry_summary[
                    "study_accession"
                ]
                .isin(
                    review_accessions
                )
            ]
            .copy()
        )

        result = (
            candidates
            .merge(
                study_metadata,
                on="study_accession",
                how="left",
                validate="one_to_one",
            )
            .merge(
                ancestry_summary,
                on="study_accession",
                how="left",
                validate="one_to_one",
            )
        )

        classifications = [
            cls.classify_row(
                row
            )
            for _, row
            in result.iterrows()
        ]

        result[
            "metadata_verification_status"
        ] = [
            value[0]
            for value
            in classifications
        ]

        result[
            "case_control_supported"
        ] = [
            value[1]
            for value
            in classifications
        ]

        result[
            "primary_susceptibility_candidate"
        ] = [
            value[2]
            for value
            in classifications
        ]

        result[
            "metadata_verification_rationale"
        ] = [
            value[3]
            for value
            in classifications
        ]

        result[
            "official_study_metadata_found"
        ] = (
            result[
                "official_reported_trait"
            ]
            .notna()
            |
            result[
                "official_study_title"
            ]
            .notna()
        )

        result[
            "official_ancestry_metadata_found"
        ] = (
            result[
                "ancestry_metadata_rows"
            ]
            .fillna(
                0
            )
            .gt(
                0
            )
        )

        result[
            "official_metadata_found"
        ] = (
            result[
                "official_study_metadata_found"
            ]
            |
            result[
                "official_ancestry_metadata_found"
            ]
        )

        return result