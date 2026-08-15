"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_sumstats_prioritization.py

Description:
    M5.3C.2 deterministic prioritization of prostate cancer GWAS studies
    with available full summary statistics.

    The objective is to identify studies most appropriate for locus-level
    prostate cancer susceptibility analysis before downloading large
    genome-wide summary-statistics datasets.

    Priority classes:

        TIER_1_PRIMARY_SUSCEPTIBILITY
            Studies with explicit evidence that the phenotype represents
            prostate cancer susceptibility, risk, incidence, or a
            case-control disease analysis.

        TIER_2_SECONDARY_DISEASE_PHENOTYPE
            Studies of aggressive, advanced, metastatic, lethal, high-risk,
            subtype, stage, or other secondary prostate cancer phenotypes.

        TIER_3_NON_SUSCEPTIBILITY
            Studies dominated by prognosis, survival, recurrence, treatment
            response, PSA-only, toxicity, or other phenotypes not suitable
            as the primary prostate cancer susceptibility GWAS.

        REVIEW_REQUIRED
            Studies whose available text identifies prostate cancer but does
            not provide sufficient evidence to classify the study as a
            susceptibility, secondary disease, or non-susceptibility study.

    Important safeguards:

        - Generic "prostate cancer" wording alone does not establish that a
          study is a susceptibility GWAS.
        - This module does not infer unavailable ancestry metadata.
        - This module does not infer sample size from Catalog association
          counts.
        - Catalog association row count is retained only as provenance and
          is not treated as sample size.
        - Harmonised summary statistics availability is favored.
        - Keyword scoring is used only for study triage.
        - Automatic prioritization does not constitute phenotype validation.
        - Manual metadata validation is required before downloading studies.
        - No summary-statistics files are downloaded.
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
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


# ============================================================================
# Priority labels
# ============================================================================

TIER_1 = "TIER_1_PRIMARY_SUSCEPTIBILITY"

TIER_2 = "TIER_2_SECONDARY_DISEASE_PHENOTYPE"

TIER_3 = "TIER_3_NON_SUSCEPTIBILITY"

REVIEW = "REVIEW_REQUIRED"


# ============================================================================
# Phenotype keyword groups
# ============================================================================

PROSTATE_CANCER_PATTERNS = (
    r"\bprostate cancer\b",
    r"\bprostatic cancer\b",
    r"\bprostate carcinoma\b",
    r"\bprostatic carcinoma\b",
    r"\bprostate neoplasm",
)


# Generic "prostate cancer" is deliberately NOT included here.
#
# A study must contain additional susceptibility-compatible wording before it
# can automatically enter Tier 1.
SUSCEPTIBILITY_PATTERNS = (
    r"\bsusceptib",
    r"\bprostate cancer risk\b",
    r"\brisk of prostate cancer\b",
    r"\bprostate cancer incidence\b",
    r"\bincidence of prostate cancer\b",
    r"\bcase[- ]?control\b",
    r"\bdisease status\b",
    r"\boverall prostate cancer\b",
)


SECONDARY_DISEASE_PATTERNS = (
    r"\baggressive\b",
    r"\badvanced\b",
    r"\bmetastatic\b",
    r"\blethal\b",
    r"\bfatal\b",
    r"\bhigh[- ]risk\b",
    r"\blocali[sz]ed\b",
    r"\bsubtype\b",
    r"\bgleason\b",
    r"\bstage\b",
    r"\bclinically significant\b",
)


NON_SUSCEPTIBILITY_PATTERNS = (
    r"\bsurvival\b",
    r"\boverall survival\b",
    r"\bprogression[- ]free\b",
    r"\bprogression free\b",
    r"\bmortality\b",
    r"\bprognos",
    r"\brecurrence\b",
    r"\bbiochemical recurrence\b",
    r"\btreatment response\b",
    r"\btherapy response\b",
    r"\bdrug response\b",
    r"\bresistance\b",
    r"\btoxicity\b",
    r"\badverse event",
)


PSA_PATTERNS = (
    r"\bprostate[- ]specific antigen\b",
    r"\bpsa\b",
)


GENERAL_CANCER_PATTERNS = (
    r"\bcancer\b",
    r"\bcarcinoma\b",
    r"\bneoplasm\b",
    r"\bmalignan",
)


# ============================================================================
# Data model
# ============================================================================


@dataclass(frozen=True)
class StudyPriorityResult:
    """One deterministic study-prioritization result."""

    priority_score: int

    priority_tier: str

    susceptibility_signal: bool

    secondary_disease_signal: bool

    non_susceptibility_signal: bool

    psa_signal: bool

    prostate_cancer_signal: bool

    generic_prostate_cancer_signal: bool

    phenotype_ambiguous: bool

    harmonised_bonus: bool

    recommended_for_locus_retrieval: bool

    rationale: str


# ============================================================================
# Generic helpers
# ============================================================================


def _is_missing_scalar(
    value: Any,
) -> bool:
    """
    Detect missing scalar values safely.

    Nested arrays are checked before pandas.isna because pd.isna(array)
    returns an array of booleans rather than a scalar.
    """

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


def _flatten_nested_value(
    value: Any,
) -> list[str]:
    """
    Flatten scalar or nested study metadata into strings.

    PyArrow may reconstruct physical list columns as numpy.ndarray, therefore
    lists, tuples, sets, and numpy arrays are handled equivalently.
    """

    if _is_missing_scalar(
        value
    ):
        return []

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
        return [
            str(
                value
            )
        ]

    flattened: list[str] = []

    for item in values:

        if _is_missing_scalar(
            item
        ):
            continue

        if isinstance(
            item,
            (
                list,
                tuple,
                set,
                np.ndarray,
            ),
        ):
            flattened.extend(
                _flatten_nested_value(
                    item
                )
            )

        else:
            flattened.append(
                str(
                    item
                )
            )

    return flattened


def _as_text(
    value: Any,
) -> str:
    """Convert scalar or nested metadata into searchable text."""

    values = _flatten_nested_value(
        value
    )

    return " | ".join(
        item
        for item
        in values
        if item.strip()
    )


def _match_any(
    text: str,
    patterns: tuple[str, ...],
) -> bool:
    """Return True when at least one regex pattern matches."""

    return any(
        re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )
        is not None
        for pattern
        in patterns
    )


def _build_search_text(
    row: pd.Series,
) -> str:
    """
    Build one searchable phenotype string.

    Only phenotype-relevant discovery fields are used.
    """

    fields = (
        "reported_traits",
        "mapped_traits",
        "study_titles",
    )

    parts = [
        _as_text(
            row.get(
                field
            )
        )
        for field
        in fields
    ]

    return (
        " | ".join(
            part
            for part
            in parts
            if part
        )
        .strip()
        .lower()
    )


# ============================================================================
# Prioritizer
# ============================================================================


class M53CProstateGWASStudyPrioritizer:
    """
    Prioritize prostate GWAS studies for locus-level summary-statistics use.
    """

    @classmethod
    def prioritize_row(
        cls,
        row: pd.Series,
    ) -> StudyPriorityResult:
        """Prioritize one discovery-table row."""

        text = _build_search_text(
            row
        )

        # ------------------------------------------------------------------
        # Summary-statistics availability
        # ------------------------------------------------------------------

        sumstats_value = row.get(
            "sumstats_available",
            False,
        )

        sumstats_available = (
            False
            if _is_missing_scalar(
                sumstats_value
            )
            else bool(
                sumstats_value
            )
        )

        preferred_file_type = row.get(
            "preferred_file_type"
        )

        preferred_type = (
            ""
            if _is_missing_scalar(
                preferred_file_type
            )
            else str(
                preferred_file_type
            )
            .strip()
            .upper()
        )

        harmonised = (
            preferred_type
            == "HARMONISED"
        )

        # ------------------------------------------------------------------
        # Phenotype signals
        # ------------------------------------------------------------------

        prostate_cancer_signal = _match_any(
            text,
            PROSTATE_CANCER_PATTERNS,
        )

        susceptibility_signal = _match_any(
            text,
            SUSCEPTIBILITY_PATTERNS,
        )

        secondary_signal = _match_any(
            text,
            SECONDARY_DISEASE_PATTERNS,
        )

        non_susceptibility_signal = _match_any(
            text,
            NON_SUSCEPTIBILITY_PATTERNS,
        )

        psa_signal = _match_any(
            text,
            PSA_PATTERNS,
        )

        general_cancer_signal = _match_any(
            text,
            GENERAL_CANCER_PATTERNS,
        )

        psa_only = (
            psa_signal
            and not prostate_cancer_signal
        )

        generic_prostate_cancer_signal = (
            prostate_cancer_signal
            and not susceptibility_signal
            and not secondary_signal
            and not non_susceptibility_signal
        )

        phenotype_ambiguous = (
            not prostate_cancer_signal
            and not general_cancer_signal
        )

        # ------------------------------------------------------------------
        # Score
        # ------------------------------------------------------------------

        score = 0

        rationale: list[str] = []

        if sumstats_available:
            score += 2

            rationale.append(
                "full summary statistics available"
            )

        else:
            score -= 20

            rationale.append(
                "full summary statistics unavailable"
            )

        if harmonised:
            score += 4

            rationale.append(
                "harmonised summary statistics available"
            )

        if prostate_cancer_signal:
            score += 4

            rationale.append(
                "explicit prostate cancer phenotype"
            )

        if susceptibility_signal:
            score += 5

            rationale.append(
                "explicit susceptibility/risk-compatible phenotype"
            )

        if secondary_signal:
            score += 2

            rationale.append(
                "secondary prostate cancer phenotype detected"
            )

        if non_susceptibility_signal:
            score -= 8

            rationale.append(
                "prognosis/survival/recurrence/treatment phenotype detected"
            )

        if psa_signal:
            score -= 2

            rationale.append(
                "PSA-related phenotype detected"
            )

        if psa_only:
            score -= 6

            rationale.append(
                "PSA signal without explicit prostate cancer phenotype"
            )

        if generic_prostate_cancer_signal:
            rationale.append(
                "generic prostate cancer phenotype requires metadata review"
            )

        if phenotype_ambiguous:
            score -= 3

            rationale.append(
                "available text does not clearly identify prostate cancer"
            )

        # ------------------------------------------------------------------
        # Tier assignment
        #
        # Ordering is deliberate:
        #
        #   1. unusable/non-susceptibility
        #   2. secondary disease phenotypes
        #   3. explicit susceptibility
        #   4. generic/ambiguous prostate cancer
        #
        # This prevents aggressive/metastatic studies from being absorbed
        # into Tier 1 merely because their title also contains
        # "prostate cancer".
        # ------------------------------------------------------------------

        if not sumstats_available:

            tier = TIER_3

        elif non_susceptibility_signal:

            tier = TIER_3

        elif psa_only:

            tier = TIER_3

        elif (
            prostate_cancer_signal
            and secondary_signal
        ):

            tier = TIER_2

        elif (
            prostate_cancer_signal
            and susceptibility_signal
        ):

            tier = TIER_1

        elif generic_prostate_cancer_signal:

            tier = REVIEW

        else:

            tier = REVIEW

        # ------------------------------------------------------------------
        # Retrieval recommendation
        # ------------------------------------------------------------------

        recommended = (
            tier
            == TIER_1
            and harmonised
            and sumstats_available
        )

        return StudyPriorityResult(
            priority_score=score,
            priority_tier=tier,
            susceptibility_signal=(
                susceptibility_signal
            ),
            secondary_disease_signal=(
                secondary_signal
            ),
            non_susceptibility_signal=(
                non_susceptibility_signal
            ),
            psa_signal=(
                psa_signal
            ),
            prostate_cancer_signal=(
                prostate_cancer_signal
            ),
            generic_prostate_cancer_signal=(
                generic_prostate_cancer_signal
            ),
            phenotype_ambiguous=(
                phenotype_ambiguous
            ),
            harmonised_bonus=(
                harmonised
            ),
            recommended_for_locus_retrieval=(
                recommended
            ),
            rationale="; ".join(
                rationale
            ),
        )

    # ----------------------------------------------------------------------
    # Full table
    # ----------------------------------------------------------------------

    @classmethod
    def prioritize(
        cls,
        studies: pd.DataFrame,
    ) -> pd.DataFrame:
        """Prioritize all discovered prostate-related studies."""

        required = {
            "study_accession",
            "sumstats_available",
            "preferred_file_type",
        }

        missing = (
            required
            - set(
                studies.columns
            )
        )

        if missing:

            raise ValueError(
                "M5.3C.2 input is missing required columns: "
                f"{sorted(missing)}"
            )

        result = studies.copy()

        decisions = [
            cls.prioritize_row(
                row
            )
            for _, row
            in result.iterrows()
        ]

        result[
            "priority_score"
        ] = [
            decision.priority_score
            for decision
            in decisions
        ]

        result[
            "priority_tier"
        ] = [
            decision.priority_tier
            for decision
            in decisions
        ]

        result[
            "susceptibility_signal"
        ] = [
            decision.susceptibility_signal
            for decision
            in decisions
        ]

        result[
            "secondary_disease_signal"
        ] = [
            decision.secondary_disease_signal
            for decision
            in decisions
        ]

        result[
            "non_susceptibility_signal"
        ] = [
            decision.non_susceptibility_signal
            for decision
            in decisions
        ]

        result[
            "psa_signal"
        ] = [
            decision.psa_signal
            for decision
            in decisions
        ]

        result[
            "prostate_cancer_signal"
        ] = [
            decision.prostate_cancer_signal
            for decision
            in decisions
        ]

        result[
            "generic_prostate_cancer_signal"
        ] = [
            decision.generic_prostate_cancer_signal
            for decision
            in decisions
        ]

        result[
            "phenotype_ambiguous"
        ] = [
            decision.phenotype_ambiguous
            for decision
            in decisions
        ]

        result[
            "harmonised_bonus"
        ] = [
            decision.harmonised_bonus
            for decision
            in decisions
        ]

        result[
            "recommended_for_locus_retrieval"
        ] = [
            decision.recommended_for_locus_retrieval
            for decision
            in decisions
        ]

        result[
            "priority_rationale"
        ] = [
            decision.rationale
            for decision
            in decisions
        ]

        tier_order = {
            TIER_1:
                1,

            TIER_2:
                2,

            REVIEW:
                3,

            TIER_3:
                4,
        }

        result[
            "_priority_tier_order"
        ] = (
            result[
                "priority_tier"
            ]
            .map(
                tier_order
            )
            .astype(
                "Int64"
            )
        )

        result = (
            result
            .sort_values(
                by=[
                    "_priority_tier_order",
                    "priority_score",
                    "study_accession",
                ],
                ascending=[
                    True,
                    False,
                    True,
                ],
                kind="stable",
            )
            .drop(
                columns=[
                    "_priority_tier_order",
                ]
            )
            .reset_index(
                drop=True
            )
        )

        result[
            "priority_rank"
        ] = pd.Series(
            range(
                1,
                len(
                    result
                )
                + 1,
            ),
            dtype="Int64",
        )

        return result