"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/clinical_cohort.py

Description:
    Core logic for M7.2 TCGA-PRAD Clinical Cohort Construction.

    This stage harmonizes public GDC TCGA-PRAD clinical metadata into a
    one-row-per-case cohort and audits candidate clinical endpoints.

    Important design principles:
        - one output row represents one GDC case;
        - the primary-disease diagnosis is preferred over prior or secondary
          malignancy diagnoses;
        - demographic.sex_at_birth is preserved as sex_at_birth;
        - missing recurrence information is never converted to "No";
        - recurrence/progression events used for time-to-event readiness must
          have a valid event/censor time;
        - raw event labels and analyzable events are reported separately;
        - endpoint selection is not based on statistical significance.

    M7.2 does NOT:
        - perform survival modeling;
        - perform Kaplan-Meier analysis;
        - perform tRF clinical association analysis;
        - join candidate tRF expression;
        - join germline genotype;
        - impute missing clinical variables;
        - select an endpoint based on significance;
        - perform causal inference.

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
class ClinicalCohortResult:
    """Container for M7.2 outputs."""

    cohort: pd.DataFrame
    endpoint_readiness: pd.DataFrame
    completeness: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Generic helpers
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


def _normalize_lower(
    value: Any,
) -> str | None:
    """Normalize optional text to lowercase."""

    normalized = _normalize_text(
        value
    )

    if normalized is None:
        return None

    return normalized.lower()


def _to_float(
    value: Any,
) -> float | None:
    """Convert optional scalar value to float."""

    if value is None:
        return None

    try:

        converted = float(
            value
        )

    except (
        TypeError,
        ValueError,
    ):

        return None

    if pd.isna(
        converted
    ):

        return None

    return converted


def _positive_days(
    value: Any,
) -> float | None:
    """
    Return strictly positive time in days.

    Zero and negative values are not considered usable follow-up/event times.
    """

    converted = _to_float(
        value
    )

    if converted is None:

        return None

    if converted <= 0:

        return None

    return converted


def _max_numeric(
    values: list[Any],
) -> float | None:
    """Return maximum valid numeric value."""

    valid: list[
        float
    ] = []

    for value in values:

        converted = _to_float(
            value
        )

        if converted is not None:

            valid.append(
                converted
            )

    if not valid:

        return None

    return max(
        valid
    )


def _minimum_positive(
    values: list[Any],
) -> float | None:
    """Return minimum strictly positive numeric value."""

    valid: list[
        float
    ] = []

    for value in values:

        converted = _positive_days(
            value
        )

        if converted is not None:

            valid.append(
                converted
            )

    if not valid:

        return None

    return min(
        valid
    )


def _as_dict_list(
    value: Any,
) -> list[dict[str, Any]]:
    """Normalize arbitrary nested value to list of dictionaries."""

    if not isinstance(
        value,
        list,
    ):

        return []

    return [
        item
        for item in value
        if isinstance(
            item,
            dict,
        )
    ]


# ============================================================================
# Diagnosis selection
# ============================================================================


def _diagnosis_is_primary_disease(
    diagnosis: dict[str, Any],
) -> bool:
    """
    Determine whether GDC explicitly marks a diagnosis as primary disease.
    """

    value = diagnosis.get(
        "diagnosis_is_primary_disease"
    )

    if isinstance(
        value,
        bool,
    ):

        return value

    normalized = _normalize_lower(
        value
    )

    return normalized in {
        "true",
        "yes",
        "1",
    }


def _classification_is_primary(
    diagnosis: dict[str, Any],
) -> bool:
    """
    Detect source-reported primary tumor classification.

    This is only a fallback when diagnosis_is_primary_disease is unavailable.
    """

    classification = _normalize_lower(
        diagnosis.get(
            "classification_of_tumor"
        )
    )

    if classification is None:

        return False

    return classification in {
        "primary",
        "primary tumor",
        "primary tumour",
    }


def _diagnosis_looks_prostate_origin(
    diagnosis: dict[str, Any],
) -> bool:
    """
    Conservative source-text prostate-origin fallback.

    This is not a gene/variant inference. It is used only after explicit
    primary-disease fields fail to identify a diagnosis.

    The TCGA project itself is TCGA-PRAD, but patients may contain prior
    primary diagnoses from other organs.
    """

    fields = [
        diagnosis.get(
            "tissue_or_organ_of_origin"
        ),
        diagnosis.get(
            "site_of_resection_or_biopsy"
        ),
        diagnosis.get(
            "primary_diagnosis"
        ),
    ]

    text = " ".join(
        value.lower()
        for value in (
            _normalize_text(
                field
            )
            for field in fields
        )
        if value is not None
    )

    return (
        "prostate"
        in text
        or
        "prostatic"
        in text
    )


def _diagnosis_age_sort_value(
    diagnosis: dict[str, Any],
) -> float:
    """Return deterministic age-based sort value."""

    age = _to_float(
        diagnosis.get(
            "age_at_diagnosis"
        )
    )

    if age is None:

        return float(
            "inf"
        )

    return age


def _select_diagnosis(
    diagnoses: list[dict[str, Any]],
) -> tuple[
    dict[str, Any] | None,
    str,
]:
    """
    Select the most defensible primary-disease diagnosis.

    Selection hierarchy:

        1. diagnosis_is_primary_disease == True
        2. classification_of_tumor explicitly indicates Primary
        3. source fields indicate prostate origin
        4. otherwise deterministic fallback

    Within each tier:
        - prefer diagnosis with primary_diagnosis populated;
        - prefer earliest age_at_diagnosis;
        - preserve stable source order.

    This prevents unrelated prior malignancies from being selected merely
    because they happen to appear first in GDC diagnoses.
    """

    if not diagnoses:

        return (
            None,
            "NO_DIAGNOSIS_AVAILABLE",
        )

    valid = [
        diagnosis
        for diagnosis in diagnoses
        if isinstance(
            diagnosis,
            dict,
        )
    ]

    if not valid:

        return (
            None,
            "NO_DIAGNOSIS_AVAILABLE",
        )

    primary_disease = [
        diagnosis
        for diagnosis in valid
        if _diagnosis_is_primary_disease(
            diagnosis
        )
    ]

    if primary_disease:

        pool = primary_disease

        selection_method = (
            "DIAGNOSIS_IS_PRIMARY_DISEASE"
        )

    else:

        explicitly_primary = [
            diagnosis
            for diagnosis in valid
            if _classification_is_primary(
                diagnosis
            )
        ]

        if explicitly_primary:

            pool = explicitly_primary

            selection_method = (
                "CLASSIFICATION_OF_TUMOR_PRIMARY"
            )

        else:

            prostate_origin = [
                diagnosis
                for diagnosis in valid
                if _diagnosis_looks_prostate_origin(
                    diagnosis
                )
            ]

            if prostate_origin:

                pool = prostate_origin

                selection_method = (
                    "PROSTATE_ORIGIN_SOURCE_TEXT_FALLBACK"
                )

            else:

                pool = valid

                selection_method = (
                    "DETERMINISTIC_DIAGNOSIS_FALLBACK"
                )

    indexed: list[
        tuple[
            int,
            float,
            int,
            dict[str, Any],
        ]
    ] = []

    for index, diagnosis in enumerate(
        pool
    ):

        has_primary_diagnosis = (
            _normalize_text(
                diagnosis.get(
                    "primary_diagnosis"
                )
            )
            is not None
        )

        indexed.append(
            (
                0
                if has_primary_diagnosis
                else 1,

                _diagnosis_age_sort_value(
                    diagnosis
                ),

                index,

                diagnosis,
            )
        )

    indexed.sort(
        key=lambda item: (
            item[
                0
            ],
            item[
                1
            ],
            item[
                2
            ],
        )
    )

    return (
        indexed[
            0
        ][
            3
        ],
        selection_method,
    )


# ============================================================================
# Follow-up extraction
# ============================================================================


def _follow_up_days(
    follow_ups: list[dict[str, Any]],
) -> list[float]:
    """Collect valid follow-up durations."""

    values: list[
        float
    ] = []

    for follow_up in follow_ups:

        value = _positive_days(
            follow_up.get(
                "days_to_follow_up"
            )
        )

        if value is not None:

            values.append(
                value
            )

    return values


def _collect_follow_up_statuses(
    follow_ups: list[dict[str, Any]],
) -> list[str]:
    """
    Collect explicit recurrence/progression status observations.

    Only explicit source values are used.
    """

    values: list[
        str
    ] = []

    for follow_up in follow_ups:

        status = _normalize_lower(
            follow_up.get(
                "progression_or_recurrence"
            )
        )

        if status is not None:

            values.append(
                status
            )

    return values


def _collect_follow_up_event_times(
    follow_ups: list[dict[str, Any]],
) -> list[float]:
    """
    Collect candidate recurrence/progression event times.

    These fields are used only when explicitly present in GDC.
    """

    fields = [
        "days_to_recurrence",
        "days_to_progression",
        "days_to_progression_free",
        "days_to_first_event",
    ]

    values: list[
        float
    ] = []

    for follow_up in follow_ups:

        for field in fields:

            value = _positive_days(
                follow_up.get(
                    field
                )
            )

            if value is not None:

                values.append(
                    value
                )

    return values


def _collect_disease_response_values(
    follow_ups: list[dict[str, Any]],
) -> list[str]:
    """
    Preserve disease-response values for descriptive QC only.

    disease_response is NOT converted automatically into a recurrence event.
    """

    values: list[
        str
    ] = []

    for follow_up in follow_ups:

        value = _normalize_text(
            follow_up.get(
                "disease_response"
            )
        )

        if value is not None:

            values.append(
                value
            )

    return values


# ============================================================================
# Overall survival derivation
# ============================================================================


def _derive_os(
    *,
    vital_status: str | None,
    days_to_death: Any,
    diagnosis_last_follow_up: Any,
    follow_up_days: list[float],
) -> tuple[
    int | None,
    float | None,
    bool,
]:
    """
    Derive conservative overall-survival endpoint.

    Dead:
        event = 1
        requires positive days_to_death.

    Alive:
        event = 0
        requires at least one positive follow-up time.

    Unknown/missing:
        unusable.
    """

    status = _normalize_lower(
        vital_status
    )

    death_time = _positive_days(
        days_to_death
    )

    diagnosis_follow_up = _positive_days(
        diagnosis_last_follow_up
    )

    censor_candidates = list(
        follow_up_days
    )

    if diagnosis_follow_up is not None:

        censor_candidates.append(
            diagnosis_follow_up
        )

    censor_time = (
        max(
            censor_candidates
        )
        if censor_candidates
        else None
    )

    if status in {
        "dead",
        "deceased",
    }:

        if death_time is None:

            return (
                1,
                None,
                False,
            )

        return (
            1,
            death_time,
            True,
        )

    if status == "alive":

        if censor_time is None:

            return (
                0,
                None,
                False,
            )

        return (
            0,
            censor_time,
            True,
        )

    return (
        None,
        None,
        False,
    )


# ============================================================================
# Recurrence/progression derivation
# ============================================================================


def _derive_recurrence(
    *,
    diagnosis: dict[str, Any] | None,
    follow_ups: list[dict[str, Any]],
) -> tuple[
    int | None,
    float | None,
    bool,
    str,
]:
    """
    Derive conservative recurrence/progression endpoint.

    Important rules:
        - missing status is not interpreted as No;
        - disease_response is not automatically converted into recurrence;
        - explicit positive status requires event time to be analytically
          usable;
        - explicit negative status requires censoring follow-up time;
        - conflicting positive/negative longitudinal statuses are resolved
          in favor of an observed positive event, but this is recorded in the
          derivation status.
    """

    diagnosis = (
        diagnosis
        if isinstance(
            diagnosis,
            dict,
        )
        else {}
    )

    statuses: list[
        str
    ] = []

    diagnosis_status = _normalize_lower(
        diagnosis.get(
            "progression_or_recurrence"
        )
    )

    if diagnosis_status is not None:

        statuses.append(
            diagnosis_status
        )

    statuses.extend(
        _collect_follow_up_statuses(
            follow_ups
        )
    )

    positive_values = {
        "yes",
        "true",
        "1",
    }

    negative_values = {
        "no",
        "false",
        "0",
    }

    has_positive = any(
        value in positive_values
        for value in statuses
    )

    has_negative = any(
        value in negative_values
        for value in statuses
    )

    diagnosis_event_time = _positive_days(
        diagnosis.get(
            "days_to_recurrence"
        )
    )

    event_times = _collect_follow_up_event_times(
        follow_ups
    )

    if diagnosis_event_time is not None:

        event_times.append(
            diagnosis_event_time
        )

    event_time = (
        min(
            event_times
        )
        if event_times
        else None
    )

    # ----------------------------------------------------------------------
    # Explicit positive event
    # ----------------------------------------------------------------------

    if has_positive:

        if event_time is None:

            return (
                1,
                None,
                False,
                (
                    "POSITIVE_EVENT_LABEL_WITHOUT_USABLE_TIME"
                    if not has_negative
                    else
                    "CONFLICTING_STATUS_POSITIVE_EVENT_WITHOUT_USABLE_TIME"
                ),
            )

        return (
            1,
            event_time,
            True,
            (
                "POSITIVE_EVENT_WITH_USABLE_TIME"
                if not has_negative
                else
                "CONFLICTING_STATUS_POSITIVE_EVENT_WITH_USABLE_TIME"
            ),
        )

    # ----------------------------------------------------------------------
    # Explicit negative/censored status
    # ----------------------------------------------------------------------

    if has_negative:

        censor_candidates = _follow_up_days(
            follow_ups
        )

        diagnosis_follow_up = _positive_days(
            diagnosis.get(
                "days_to_last_follow_up"
            )
        )

        if diagnosis_follow_up is not None:

            censor_candidates.append(
                diagnosis_follow_up
            )

        censor_time = (
            max(
                censor_candidates
            )
            if censor_candidates
            else None
        )

        if censor_time is None:

            return (
                0,
                None,
                False,
                "NEGATIVE_EVENT_LABEL_WITHOUT_USABLE_CENSOR_TIME",
            )

        return (
            0,
            censor_time,
            True,
            "NEGATIVE_EVENT_WITH_USABLE_CENSOR_TIME",
        )

    # ----------------------------------------------------------------------
    # No explicit recurrence/progression status
    # ----------------------------------------------------------------------

    return (
        None,
        None,
        False,
        "RECURRENCE_STATUS_NOT_REPORTED",
    )


# ============================================================================
# Cohort construction
# ============================================================================


def build_clinical_cohort(
    *,
    cases: list[dict[str, Any]],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Build one-row-per-case TCGA-PRAD clinical cohort.
    """

    days_per_year = float(
        config[
            "harmonization"
        ][
            "age"
        ][
            "days_per_year"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for case in cases:

        if not isinstance(
            case,
            dict,
        ):

            continue

        # ------------------------------------------------------------------
        # Nested entities
        # ------------------------------------------------------------------

        demographic = case.get(
            "demographic"
        )

        if not isinstance(
            demographic,
            dict,
        ):

            demographic = {}

        diagnoses = _as_dict_list(
            case.get(
                "diagnoses"
            )
        )

        follow_ups = _as_dict_list(
            case.get(
                "follow_ups"
            )
        )

        samples = _as_dict_list(
            case.get(
                "samples"
            )
        )

        # ------------------------------------------------------------------
        # Primary disease diagnosis
        # ------------------------------------------------------------------

        (
            selected_diagnosis,
            diagnosis_selection_method,
        ) = _select_diagnosis(
            diagnoses
        )

        diagnosis = (
            selected_diagnosis
            if selected_diagnosis is not None
            else {}
        )

        diagnosis_available = bool(
            diagnosis
        )

        diagnosis_primary_disease_flag = (
            _diagnosis_is_primary_disease(
                diagnosis
            )
            if diagnosis
            else False
        )

        # ------------------------------------------------------------------
        # Age
        # ------------------------------------------------------------------

        age_at_diagnosis_days = _to_float(
            diagnosis.get(
                "age_at_diagnosis"
            )
        )

        if age_at_diagnosis_days is not None:

            age_at_diagnosis_years = (
                age_at_diagnosis_days
                /
                days_per_year
            )

            age_source = (
                "diagnoses.age_at_diagnosis"
            )

        else:

            age_at_index = _to_float(
                demographic.get(
                    "age_at_index"
                )
            )

            age_at_diagnosis_years = (
                age_at_index
            )

            age_source = (
                "demographic.age_at_index"
                if age_at_index is not None
                else None
            )

        # ------------------------------------------------------------------
        # Follow-up
        # ------------------------------------------------------------------

        follow_up_days = _follow_up_days(
            follow_ups
        )

        diagnosis_last_follow_up = (
            diagnosis.get(
                "days_to_last_follow_up"
            )
        )

        max_follow_up_days = _max_numeric(
            follow_up_days
            +
            [
                diagnosis_last_follow_up
            ]
        )

        disease_response_values = (
            _collect_disease_response_values(
                follow_ups
            )
        )

        # ------------------------------------------------------------------
        # Overall survival
        # ------------------------------------------------------------------

        vital_status = _normalize_text(
            demographic.get(
                "vital_status"
            )
        )

        (
            os_event,
            os_time_days,
            os_usable,
        ) = _derive_os(
            vital_status=vital_status,
            days_to_death=demographic.get(
                "days_to_death"
            ),
            diagnosis_last_follow_up=diagnosis_last_follow_up,
            follow_up_days=follow_up_days,
        )

        # ------------------------------------------------------------------
        # Recurrence / progression
        # ------------------------------------------------------------------

        (
            recurrence_event,
            recurrence_time_days,
            recurrence_usable,
            recurrence_derivation_status,
        ) = _derive_recurrence(
            diagnosis=diagnosis,
            follow_ups=follow_ups,
        )

        # ------------------------------------------------------------------
        # Identifiers
        # ------------------------------------------------------------------

        case_id = _normalize_text(
            case.get(
                "case_id"
            )
        )

        patient_id = _normalize_text(
            case.get(
                "submitter_id"
            )
        )

        project = case.get(
            "project"
        )

        if not isinstance(
            project,
            dict,
        ):

            project = {}

        project_id = (
            _normalize_text(
                project.get(
                    "project_id"
                )
            )
            or
            str(
                config[
                    "target_cohort"
                ][
                    "project_id"
                ]
            )
        )

        # ------------------------------------------------------------------
        # Sample summary
        # ------------------------------------------------------------------

        sample_types = sorted(
            {
                str(value)
                for value in (
                    _normalize_text(
                        sample.get(
                            "sample_type"
                        )
                    )
                    for sample in samples
                )
                if value is not None
            }
        )

        # ------------------------------------------------------------------
        # Row
        # ------------------------------------------------------------------

        rows.append(
            {
                # ----------------------------------------------------------
                # Identifiers
                # ----------------------------------------------------------

                "case_id":
                    case_id,

                "patient_id":
                    patient_id,

                "project_id":
                    project_id,

                "primary_site":
                    _normalize_text(
                        case.get(
                            "primary_site"
                        )
                    ),

                "disease_type":
                    _normalize_text(
                        case.get(
                            "disease_type"
                        )
                    ),

                # ----------------------------------------------------------
                # Demographic
                # ----------------------------------------------------------

                "age_at_diagnosis_years":
                    age_at_diagnosis_years,

                "age_source":
                    age_source,

                "sex_at_birth":
                    _normalize_text(
                        demographic.get(
                            "sex_at_birth"
                        )
                    ),

                "race":
                    _normalize_text(
                        demographic.get(
                            "race"
                        )
                    ),

                "ethnicity":
                    _normalize_text(
                        demographic.get(
                            "ethnicity"
                        )
                    ),

                "vital_status":
                    vital_status,

                # ----------------------------------------------------------
                # Selected diagnosis
                # ----------------------------------------------------------

                "selected_diagnosis_id":
                    _normalize_text(
                        diagnosis.get(
                            "diagnosis_id"
                        )
                    ),

                "selected_diagnosis_submitter_id":
                    _normalize_text(
                        diagnosis.get(
                            "submitter_id"
                        )
                    ),

                "diagnosis_selection_method":
                    diagnosis_selection_method,

                "diagnosis_is_primary_disease":
                    diagnosis_primary_disease_flag,

                "classification_of_tumor":
                    _normalize_text(
                        diagnosis.get(
                            "classification_of_tumor"
                        )
                    ),

                "primary_diagnosis":
                    _normalize_text(
                        diagnosis.get(
                            "primary_diagnosis"
                        )
                    ),

                "tissue_or_organ_of_origin":
                    _normalize_text(
                        diagnosis.get(
                            "tissue_or_organ_of_origin"
                        )
                    ),

                "morphology":
                    _normalize_text(
                        diagnosis.get(
                            "morphology"
                        )
                    ),

                "tumor_grade":
                    _normalize_text(
                        diagnosis.get(
                            "tumor_grade"
                        )
                    ),

                "ajcc_pathologic_stage":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_pathologic_stage"
                        )
                    ),

                "ajcc_pathologic_t":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_pathologic_t"
                        )
                    ),

                "ajcc_pathologic_n":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_pathologic_n"
                        )
                    ),

                "ajcc_pathologic_m":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_pathologic_m"
                        )
                    ),

                "ajcc_clinical_stage":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_clinical_stage"
                        )
                    ),

                "ajcc_clinical_t":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_clinical_t"
                        )
                    ),

                "ajcc_clinical_n":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_clinical_n"
                        )
                    ),

                "ajcc_clinical_m":
                    _normalize_text(
                        diagnosis.get(
                            "ajcc_clinical_m"
                        )
                    ),

                # ----------------------------------------------------------
                # Source recurrence fields
                # ----------------------------------------------------------

                "progression_or_recurrence":
                    _normalize_text(
                        diagnosis.get(
                            "progression_or_recurrence"
                        )
                    ),

                "days_to_recurrence":
                    _to_float(
                        diagnosis.get(
                            "days_to_recurrence"
                        )
                    ),

                # ----------------------------------------------------------
                # Nested record summaries
                # ----------------------------------------------------------

                "diagnosis_record_count":
                    int(
                        len(
                            diagnoses
                        )
                    ),

                "follow_up_record_count":
                    int(
                        len(
                            follow_ups
                        )
                    ),

                "sample_record_count":
                    int(
                        len(
                            samples
                        )
                    ),

                "sample_types":
                    sample_types,

                "disease_response_values":
                    disease_response_values,

                "max_follow_up_days":
                    max_follow_up_days,

                # ----------------------------------------------------------
                # OS
                # ----------------------------------------------------------

                "os_event":
                    os_event,

                "os_time_days":
                    os_time_days,

                "os_usable":
                    bool(
                        os_usable
                    ),

                # ----------------------------------------------------------
                # Recurrence / progression
                # ----------------------------------------------------------

                "recurrence_event":
                    recurrence_event,

                "recurrence_time_days":
                    recurrence_time_days,

                "recurrence_usable":
                    bool(
                        recurrence_usable
                    ),

                "recurrence_derivation_status":
                    recurrence_derivation_status,

                # ----------------------------------------------------------
                # QC
                # ----------------------------------------------------------

                "diagnosis_available":
                    diagnosis_available,

                "clinical_usable":
                    bool(
                        case_id is not None
                        and
                        patient_id is not None
                        and
                        diagnosis_available
                    ),
            }
        )

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "No TCGA-PRAD cases available after clinical harmonization."
        )

    if result[
        "case_id"
    ].isna().any():

        raise RuntimeError(
            "M7.2 clinical cohort contains missing case_id."
        )

    if result[
        "case_id"
    ].duplicated().any():

        duplicates = (
            result.loc[
                result[
                    "case_id"
                ].duplicated(
                    keep=False
                ),
                "case_id",
            ]
            .astype(
                str
            )
            .unique()
            .tolist()
        )

        raise RuntimeError(
            "M7.2 produced duplicate case_id rows: "
            f"{duplicates[:10]}"
        )

    return result.sort_values(
        [
            "patient_id",
            "case_id",
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Completeness audit
# ============================================================================


def build_completeness(
    *,
    cohort: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build clinical variable completeness table.
    """

    fields = [
        "age_at_diagnosis_years",
        "sex_at_birth",
        "race",
        "ethnicity",
        "vital_status",
        "primary_diagnosis",
        "diagnosis_is_primary_disease",
        "classification_of_tumor",
        "tissue_or_organ_of_origin",
        "tumor_grade",
        "ajcc_pathologic_stage",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "ajcc_pathologic_m",
        "ajcc_clinical_stage",
        "ajcc_clinical_t",
        "ajcc_clinical_n",
        "ajcc_clinical_m",
        "progression_or_recurrence",
        "days_to_recurrence",
        "max_follow_up_days",
        "os_event",
        "os_time_days",
        "recurrence_event",
        "recurrence_time_days",
        "recurrence_derivation_status",
    ]

    rows: list[
        dict[str, Any]
    ] = []

    total_cases = int(
        len(
            cohort
        )
    )

    for field in fields:

        if field not in cohort.columns:

            continue

        non_missing_cases = int(
            cohort[
                field
            ]
            .notna()
            .sum()
        )

        missing_cases = (
            total_cases
            -
            non_missing_cases
        )

        non_missing_fraction = (
            non_missing_cases
            /
            total_cases
            if total_cases
            else 0.0
        )

        rows.append(
            {
                "field":
                    field,

                "total_cases":
                    total_cases,

                "non_missing_cases":
                    non_missing_cases,

                "missing_cases":
                    missing_cases,

                "non_missing_fraction":
                    non_missing_fraction,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Endpoint readiness
# ============================================================================


def build_endpoint_readiness(
    *,
    cohort: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Assess clinical endpoint readiness without fitting models.

    IMPORTANT:
        Event counts used for downstream analytical readiness are restricted
        to valid time-event pairs.

        Raw event labels are retained separately for QC. Therefore an event
        label without an event/follow-up time does not count as an analyzable
        survival event.
    """

    policy = config[
        "endpoint_readiness"
    ]

    minimum_fraction = float(
        policy[
            "minimum_non_missing_fraction"
        ]
    )

    minimum_events = int(
        policy[
            "minimum_events"
        ]
    )

    minimum_pairs = int(
        policy[
            "minimum_time_event_pairs"
        ]
    )

    minimum_unique_cases = int(
        policy.get(
            "minimum_unique_cases",
            0,
        )
    )

    endpoint_specs = [
        {
            "endpoint_id":
                "OVERALL_SURVIVAL",

            "event_column":
                "os_event",

            "time_column":
                "os_time_days",

            "usable_column":
                "os_usable",
        },
        {
            "endpoint_id":
                "RECURRENCE_PROGRESSION",

            "event_column":
                "recurrence_event",

            "time_column":
                "recurrence_time_days",

            "usable_column":
                "recurrence_usable",
        },
    ]

    total_cases = int(
        len(
            cohort
        )
    )

    unique_cases = int(
        cohort[
            "case_id"
        ].nunique()
    )

    rows: list[
        dict[str, Any]
    ] = []

    for spec in endpoint_specs:

        event_column = str(
            spec[
                "event_column"
            ]
        )

        time_column = str(
            spec[
                "time_column"
            ]
        )

        usable_column = str(
            spec[
                "usable_column"
            ]
        )

        # ------------------------------------------------------------------
        # Strict usable time-event mask
        # ------------------------------------------------------------------

        declared_usable = (
            cohort[
                usable_column
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
        )

        valid_time = (
            pd.to_numeric(
                cohort[
                    time_column
                ],
                errors="coerce",
            )
            >
            0
        ).fillna(
            False
        )

        valid_event_label = (
            cohort[
                event_column
            ]
            .isin(
                [
                    0,
                    1,
                ]
            )
            .fillna(
                False
            )
        )

        usable_mask = (
            declared_usable
            &
            valid_time
            &
            valid_event_label
        )

        usable_cases = int(
            usable_mask.sum()
        )

        usable_fraction = (
            usable_cases
            /
            total_cases
            if total_cases
            else 0.0
        )

        # ------------------------------------------------------------------
        # Raw positive event labels
        # ------------------------------------------------------------------

        raw_positive_mask = (
            cohort[
                event_column
            ]
            .eq(
                1
            )
            .fillna(
                False
            )
        )

        raw_positive_event_labels = int(
            raw_positive_mask.sum()
        )

        # ------------------------------------------------------------------
        # Analyzable events
        # ------------------------------------------------------------------

        analyzable_event_mask = (
            usable_mask
            &
            raw_positive_mask
        )

        analyzable_events = int(
            analyzable_event_mask.sum()
        )

        # ------------------------------------------------------------------
        # Analyzable censored observations
        # ------------------------------------------------------------------

        analyzable_censored_mask = (
            usable_mask
            &
            cohort[
                event_column
            ]
            .eq(
                0
            )
            .fillna(
                False
            )
        )

        analyzable_censored_cases = int(
            analyzable_censored_mask.sum()
        )

        # ------------------------------------------------------------------
        # Positive labels lacking usable time
        # ------------------------------------------------------------------

        events_without_usable_time = int(
            (
                raw_positive_mask
                &
                ~usable_mask
            ).sum()
        )

        # ------------------------------------------------------------------
        # Defensive invariants
        # ------------------------------------------------------------------

        if analyzable_events > usable_cases:

            raise RuntimeError(
                f"{spec['endpoint_id']}: analyzable events "
                f"({analyzable_events}) exceed usable cases "
                f"({usable_cases})."
            )

        if (
            analyzable_events
            +
            analyzable_censored_cases
            !=
            usable_cases
        ):

            raise RuntimeError(
                f"{spec['endpoint_id']}: usable cases are not fully "
                "partitioned into event and censored observations."
            )

        # ------------------------------------------------------------------
        # Readiness
        # ------------------------------------------------------------------

        endpoint_ready = bool(
            unique_cases
            >=
            minimum_unique_cases
            and
            usable_fraction
            >=
            minimum_fraction
            and
            analyzable_events
            >=
            minimum_events
            and
            usable_cases
            >=
            minimum_pairs
        )

        if unique_cases < minimum_unique_cases:

            endpoint_status = (
                "ENDPOINT_NOT_READY_INSUFFICIENT_COHORT_SIZE"
            )

        elif analyzable_events < minimum_events:

            endpoint_status = (
                "ENDPOINT_NOT_READY_INSUFFICIENT_EVENTS"
            )

        elif usable_fraction < minimum_fraction:

            endpoint_status = (
                "ENDPOINT_NOT_READY_INSUFFICIENT_COMPLETENESS"
            )

        elif usable_cases < minimum_pairs:

            endpoint_status = (
                "ENDPOINT_NOT_READY_INSUFFICIENT_TIME_EVENT_PAIRS"
            )

        else:

            endpoint_status = (
                "ENDPOINT_READY_FOR_DOWNSTREAM_ASSOCIATION"
            )

        rows.append(
            {
                "endpoint_id":
                    str(
                        spec[
                            "endpoint_id"
                        ]
                    ),

                "total_cases":
                    total_cases,

                "unique_cases":
                    unique_cases,

                "usable_cases":
                    usable_cases,

                "usable_fraction":
                    usable_fraction,

                # Main downstream event count
                "events":
                    analyzable_events,

                # QC event counts
                "raw_positive_event_labels":
                    raw_positive_event_labels,

                "events_without_usable_time":
                    events_without_usable_time,

                "analyzable_censored_cases":
                    analyzable_censored_cases,

                "minimum_fraction_required":
                    minimum_fraction,

                "minimum_events_required":
                    minimum_events,

                "minimum_time_event_pairs_required":
                    minimum_pairs,

                "minimum_unique_cases_required":
                    minimum_unique_cases,

                "endpoint_ready":
                    endpoint_ready,

                "endpoint_status":
                    endpoint_status,

                "primary_endpoint_selected":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    cohort: pd.DataFrame,
    endpoint_readiness: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build M7.2 summary from analyzable endpoint counts.
    """

    os_rows = endpoint_readiness.loc[
        endpoint_readiness[
            "endpoint_id"
        ]
        ==
        "OVERALL_SURVIVAL"
    ]

    recurrence_rows = endpoint_readiness.loc[
        endpoint_readiness[
            "endpoint_id"
        ]
        ==
        "RECURRENCE_PROGRESSION"
    ]

    if len(
        os_rows
    ) != 1:

        raise RuntimeError(
            "Expected exactly one OVERALL_SURVIVAL readiness row."
        )

    if len(
        recurrence_rows
    ) != 1:

        raise RuntimeError(
            "Expected exactly one RECURRENCE_PROGRESSION readiness row."
        )

    os_row = os_rows.iloc[
        0
    ]

    recurrence_row = recurrence_rows.iloc[
        0
    ]

    primary_disease_selected = int(
        cohort[
            "diagnosis_is_primary_disease"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.2",

                "project_id":
                    "TCGA-PRAD",

                "cases_acquired":
                    int(
                        len(
                            cohort
                        )
                    ),

                "unique_cases":
                    int(
                        cohort[
                            "case_id"
                        ].nunique()
                    ),

                "clinically_usable_cases":
                    int(
                        cohort[
                            "clinical_usable"
                        ]
                        .astype(
                            bool
                        )
                        .sum()
                    ),

                "cases_with_explicit_primary_disease_diagnosis":
                    primary_disease_selected,

                "cases_with_sex_at_birth":
                    int(
                        cohort[
                            "sex_at_birth"
                        ]
                        .notna()
                        .sum()
                    ),

                # ----------------------------------------------------------
                # Overall survival
                # ----------------------------------------------------------

                "os_usable_cases":
                    int(
                        os_row[
                            "usable_cases"
                        ]
                    ),

                "os_events":
                    int(
                        os_row[
                            "events"
                        ]
                    ),

                "os_raw_positive_event_labels":
                    int(
                        os_row[
                            "raw_positive_event_labels"
                        ]
                    ),

                "os_events_without_usable_time":
                    int(
                        os_row[
                            "events_without_usable_time"
                        ]
                    ),

                "os_endpoint_ready":
                    bool(
                        os_row[
                            "endpoint_ready"
                        ]
                    ),

                # ----------------------------------------------------------
                # Recurrence / progression
                # ----------------------------------------------------------

                "recurrence_usable_cases":
                    int(
                        recurrence_row[
                            "usable_cases"
                        ]
                    ),

                "recurrence_events":
                    int(
                        recurrence_row[
                            "events"
                        ]
                    ),

                "recurrence_raw_positive_event_labels":
                    int(
                        recurrence_row[
                            "raw_positive_event_labels"
                        ]
                    ),

                "recurrence_events_without_usable_time":
                    int(
                        recurrence_row[
                            "events_without_usable_time"
                        ]
                    ),

                "recurrence_endpoint_ready":
                    bool(
                        recurrence_row[
                            "endpoint_ready"
                        ]
                    ),

                # ----------------------------------------------------------
                # Policy
                # ----------------------------------------------------------

                "primary_endpoint_selected":
                    False,

                "clinical_association_performed":
                    False,

                "survival_model_fitted":
                    False,

                "next_stage":
                    "M7.2B_ENDPOINT_READINESS_RESOLUTION",
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def construct_clinical_cohort(
    *,
    cases: list[dict[str, Any]],
    config: dict[str, Any],
) -> ClinicalCohortResult:
    """
    Execute M7.2 TCGA-PRAD clinical cohort construction.
    """

    cohort = build_clinical_cohort(
        cases=cases,
        config=config,
    )

    completeness = build_completeness(
        cohort=cohort,
    )

    endpoint_readiness = build_endpoint_readiness(
        cohort=cohort,
        config=config,
    )

    summary = build_summary(
        cohort=cohort,
        endpoint_readiness=endpoint_readiness,
    )

    return ClinicalCohortResult(
        cohort=cohort,
        endpoint_readiness=endpoint_readiness,
        completeness=completeness,
        summary=summary,
    )