"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/study_independence_harmonization.py

Description:
    Core logic for M7.5C Study Independence & Harmonization Audit.

    This module:
        - groups GWAS analyses into publication families;
        - identifies explicit publication/cohort relatedness;
        - attaches phenotype and ancestry metadata to exact variant hits;
        - harmonizes effect direction using exact allele orientation;
        - aggregates directions within publication families;
        - distinguishes within-publication from cross-publication
          directional heterogeneity;
        - retains study independence as unresolved.

    No network access, external allele lookup, meta-analysis, or causal
    inference occurs in this module.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
from typing import Any

import pandas as pd


# ============================================================================
# Result container
# ============================================================================


@dataclass(frozen=True)
class M75CResult:
    """Container for M7.5C outputs."""

    publication_family_inventory: pd.DataFrame
    study_pair_overlap: pd.DataFrame
    candidate_study_harmonization: pd.DataFrame
    candidate_effect_harmonization: pd.DataFrame
    publication_family_direction: pd.DataFrame
    candidate_direction_resolution: pd.DataFrame
    candidate_resolution: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str:
    """Normalize scalar text."""

    if value is None:

        return ""

    return str(
        value
    ).strip()


def _normalize_string_list(
    value: Any,
    *,
    ignored_values: set[str] | None = None,
) -> list[str]:
    """Normalize list-like metadata without inference."""

    if not isinstance(
        value,
        list,
    ):

        return []

    ignored = {
        value.upper()
        for value in (
            ignored_values
            or set()
        )
    }

    output: list[str] = []

    for item in value:

        text = _normalize_text(
            item
        )

        if not text:

            continue

        if text.upper() in ignored:

            continue

        output.append(
            text
        )

    return list(
        dict.fromkeys(
            output
        )
    )


def _safe_float(
    value: Any,
) -> float | None:
    """Convert to finite float when possible."""

    if value is None:

        return None

    text = str(
        value
    ).strip()

    if not text:

        return None

    try:

        result = float(
            text
        )

    except (
        TypeError,
        ValueError,
    ):

        return None

    if not math.isfinite(
        result
    ):

        return None

    return result


def _normalized_allele(
    value: Any,
) -> str | None:
    """Normalize allele without complement inference."""

    text = _normalize_text(
        value
    ).upper()

    if not text:

        return None

    return text


# ============================================================================
# Publication-family inventory
# ============================================================================


def build_publication_family_inventory(
    *,
    study_inventory: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Group eligible studies by PubMed ID."""

    ignored = {
        str(value).upper()
        for value in config[
            "cohort_policy"
        ][
            "ignored_unknown_values"
        ]
    }

    eligible = study_inventory.loc[
        study_inventory[
            "eligible_for_exact_lookup"
        ]
        .fillna(False)
        .astype(bool)
    ].copy()

    rows: list[dict[str, Any]] = []

    for pubmed_id, group in eligible.groupby(
        "pubmed_id",
        dropna=False,
        sort=True,
    ):

        accessions = sorted(
            group[
                "accession_id"
            ]
            .astype(str)
            .unique()
            .tolist()
        )

        traits = sorted(
            set(
                group[
                    "disease_trait"
                ]
                .dropna()
                .astype(str)
            )
        )

        cohorts: set[str] = set()
        ancestries: set[str] = set()

        for value in group[
            "cohort"
        ]:

            cohorts.update(
                _normalize_string_list(
                    value,
                    ignored_values=ignored,
                )
            )

        for value in group[
            "discovery_ancestry"
        ]:

            ancestries.update(
                _normalize_string_list(
                    value
                )
            )

        family_id = (
            f"PMID_{pubmed_id}"
            if pd.notna(
                pubmed_id
            )
            else
            "PMID_UNRESOLVED"
        )

        rows.append(
            {
                "publication_family_id":
                    family_id,

                "pubmed_id":
                    pubmed_id,

                "study_count":
                    len(
                        accessions
                    ),

                "study_accessions":
                    accessions,

                "disease_traits":
                    traits,

                "cohort_labels":
                    sorted(
                        cohorts
                    ),

                "discovery_ancestry_labels":
                    sorted(
                        ancestries
                    ),

                "multiple_accessions_same_publication":
                    bool(
                        len(
                            accessions
                        )
                        >
                        1
                    ),

                "independent_replication_unit":
                    False,

                "sample_independence_verified":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Pairwise publication / cohort relatedness
# ============================================================================


def build_study_pair_overlap(
    *,
    study_inventory: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Audit explicit publication and cohort overlap."""

    ignored = {
        str(value).upper()
        for value in config[
            "cohort_policy"
        ][
            "ignored_unknown_values"
        ]
    }

    eligible = study_inventory.loc[
        study_inventory[
            "eligible_for_exact_lookup"
        ]
        .fillna(False)
        .astype(bool)
    ].copy()

    records = eligible.to_dict(
        orient="records"
    )

    rows: list[dict[str, Any]] = []

    for left, right in combinations(
        records,
        2,
    ):

        left_pubmed = left.get(
            "pubmed_id"
        )

        right_pubmed = right.get(
            "pubmed_id"
        )

        same_pubmed = bool(
            left_pubmed is not None
            and
            right_pubmed is not None
            and
            str(
                left_pubmed
            )
            ==
            str(
                right_pubmed
            )
        )

        left_cohorts = {
            value.lower()
            for value in _normalize_string_list(
                left.get(
                    "cohort"
                ),
                ignored_values=ignored,
            )
        }

        right_cohorts = {
            value.lower()
            for value in _normalize_string_list(
                right.get(
                    "cohort"
                ),
                ignored_values=ignored,
            )
        }

        overlap = sorted(
            left_cohorts
            &
            right_cohorts
        )

        if same_pubmed:

            status = (
                "RELATED_ANALYSIS_SAME_PUBLICATION"
            )

        elif overlap:

            status = (
                "POTENTIAL_SAMPLE_OVERLAP_SHARED_COHORT"
            )

        else:

            status = (
                "NO_EXPLICIT_METADATA_OVERLAP_"
                "INDEPENDENCE_UNRESOLVED"
            )

        rows.append(
            {
                "study_a":
                    str(
                        left[
                            "accession_id"
                        ]
                    ),

                "study_b":
                    str(
                        right[
                            "accession_id"
                        ]
                    ),

                "pubmed_a":
                    left_pubmed,

                "pubmed_b":
                    right_pubmed,

                "same_pubmed":
                    same_pubmed,

                "cohort_overlap":
                    overlap,

                "explicit_cohort_overlap":
                    bool(
                        overlap
                    ),

                "potential_sample_overlap":
                    bool(
                        same_pubmed
                        or
                        overlap
                    ),

                "sample_independence_verified":
                    False,

                "relationship_status":
                    status,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Phenotype harmonization
# ============================================================================


def _phenotype_status(
    trait: Any,
    config: dict[str, Any],
) -> str:
    """Classify prostate-cancer phenotype compatibility."""

    policy = config[
        "phenotype_policy"
    ]

    statuses = policy[
        "statuses"
    ]

    text = _normalize_text(
        trait
    )

    lower = text.lower()

    exact = {
        str(value).lower()
        for value in policy[
            "exact_general_traits"
        ]
    }

    subgroup = {
        str(value).lower()
        for value in policy[
            "compatible_subgroup_traits"
        ]
    }

    if lower in exact:

        return statuses[
            "general"
        ]

    if lower in subgroup:

        return statuses[
            "subgroup"
        ]

    for prefix in policy[
        "accepted_prefixes"
    ]:

        if lower.startswith(
            str(
                prefix
            ).lower()
        ):

            return statuses[
                "general"
            ]

    return statuses[
        "unresolved"
    ]


# ============================================================================
# Candidate × study metadata
# ============================================================================


def build_candidate_study_harmonization(
    *,
    exact_hits: pd.DataFrame,
    study_inventory: pd.DataFrame,
    file_manifest: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Attach study/source metadata to each exact variant hit."""

    study_columns = [
        "accession_id",
        "pubmed_id",
        "disease_trait",
        "cohort",
        "discovery_ancestry",
        "replication_ancestry",
    ]

    file_columns = [
        "accession_id",
        "file_representation",
        "harmonised",
        "compression",
        "filename",
    ]

    merged = exact_hits.merge(
        study_inventory[
            study_columns
        ],
        on="accession_id",
        how="left",
        validate="many_to_one",
    )

    merged = merged.merge(
        file_manifest[
            file_columns
        ],
        on="accession_id",
        how="left",
        validate="many_to_one",
    )

    merged[
        "publication_family_id"
    ] = merged[
        "pubmed_id"
    ].apply(
        lambda value:
            (
                f"PMID_{value}"
                if pd.notna(
                    value
                )
                else
                "PMID_UNRESOLVED"
            )
    )

    merged[
        "phenotype_harmonization_status"
    ] = merged[
        "disease_trait"
    ].apply(
        lambda value:
            _phenotype_status(
                value,
                config,
            )
    )

    compatible_statuses = {
        config[
            "phenotype_policy"
        ][
            "statuses"
        ][
            "general"
        ],
        config[
            "phenotype_policy"
        ][
            "statuses"
        ][
            "subgroup"
        ],
    }

    merged[
        "phenotype_compatible"
    ] = merged[
        "phenotype_harmonization_status"
    ].isin(
        compatible_statuses
    )

    merged[
        "study_independence_verified"
    ] = False

    merged[
        "ancestry_independence_verified"
    ] = False

    return merged


# ============================================================================
# Canonical allele orientation
# ============================================================================


def _choose_canonical_hit(
    group: pd.DataFrame,
) -> pd.Series | None:
    """
    Deterministically choose one harmonised record with complete alleles.

    This record defines allele orientation only.
    It is NOT treated as biological reference evidence.
    """

    usable = group.loc[
        group[
            "harmonised"
        ]
        .fillna(False)
        .astype(bool)
    ].copy()

    usable[
        "_effect"
    ] = usable[
        "effect_allele"
    ].apply(
        _normalized_allele
    )

    usable[
        "_other"
    ] = usable[
        "other_allele"
    ].apply(
        _normalized_allele
    )

    usable = usable.loc[
        usable[
            "_effect"
        ].notna()
        &
        usable[
            "_other"
        ].notna()
    ]

    if usable.empty:

        return None

    usable = usable.sort_values(
        [
            "accession_id",
            "source_file_url",
        ],
        kind="stable",
    )

    return usable.iloc[
        0
    ]


# ============================================================================
# Effect harmonization
# ============================================================================


def _effect_row(
    *,
    record: pd.Series,
    canonical_effect: str | None,
    canonical_other: str | None,
    allele_orientation: str | None,
    effect_type: str | None,
    aligned_effect: float | None,
    status: str,
) -> dict[str, Any]:

    if aligned_effect is None:

        direction = None

    elif aligned_effect > 0:

        direction = "POSITIVE"

    elif aligned_effect < 0:

        direction = "NEGATIVE"

    else:

        direction = "ZERO"

    return {
        "rsid":
            str(
                record[
                    "rs_id"
                ]
            ),

        "accession_id":
            str(
                record[
                    "accession_id"
                ]
            ),

        "pubmed_id":
            record.get(
                "pubmed_id"
            ),

        "publication_family_id":
            record.get(
                "publication_family_id"
            ),

        "harmonised_source":
            bool(
                record.get(
                    "harmonised",
                    False,
                )
            ),

        "effect_allele":
            _normalized_allele(
                record.get(
                    "effect_allele"
                )
            ),

        "other_allele":
            _normalized_allele(
                record.get(
                    "other_allele"
                )
            ),

        "canonical_effect_allele":
            canonical_effect,

        "canonical_other_allele":
            canonical_other,

        "allele_orientation":
            allele_orientation,

        "effect_type":
            effect_type,

        "aligned_effect":
            aligned_effect,

        "aligned_direction":
            direction,

        "p_value":
            _safe_float(
                record.get(
                    "p_value"
                )
            ),

        "effect_harmonization_status":
            status,

        "effect_direction_usable":
            bool(
                aligned_effect
                is not None
            ),

        "study_independence_verified":
            False,
    }


def build_candidate_effect_harmonization(
    *,
    candidate_study_harmonization: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Harmonize alleles and directional effect signs."""

    policy = config[
        "effect_harmonization"
    ]

    statuses = policy[
        "statuses"
    ]

    rows: list[dict[str, Any]] = []

    for rsid, group in candidate_study_harmonization.groupby(
        "rs_id",
        sort=False,
    ):

        canonical = _choose_canonical_hit(
            group
        )

        if canonical is None:

            canonical_effect = None
            canonical_other = None

        else:

            canonical_effect = _normalized_allele(
                canonical[
                    "effect_allele"
                ]
            )

            canonical_other = _normalized_allele(
                canonical[
                    "other_allele"
                ]
            )

        for _, record in group.iterrows():

            effect = _normalized_allele(
                record.get(
                    "effect_allele"
                )
            )

            other = _normalized_allele(
                record.get(
                    "other_allele"
                )
            )

            harmonised_source = bool(
                record.get(
                    "harmonised",
                    False,
                )
            )

            beta = _safe_float(
                record.get(
                    "beta"
                )
            )

            odds_ratio = _safe_float(
                record.get(
                    "odds_ratio"
                )
            )

            orientation: str | None = None
            effect_type: str | None = None
            aligned_effect: float | None = None

            if (
                canonical_effect is None
                or
                canonical_other is None
            ):

                status = statuses[
                    "unresolved"
                ]

            elif (
                effect is None
                or
                other is None
            ):

                status = statuses[
                    "missing_alleles"
                ]

            elif (
                policy[
                    "harmonised_sources_only_for_directional_comparison"
                ]
                and
                not harmonised_source
            ):

                status = statuses[
                    "nonharmonised_source"
                ]

            elif (
                effect
                ==
                canonical_effect
                and
                other
                ==
                canonical_other
            ):

                orientation = "SAME"

                if beta is not None:

                    aligned_effect = beta
                    effect_type = "BETA"

                    status = statuses[
                        "usable"
                    ]

                elif (
                    odds_ratio is not None
                    and
                    odds_ratio > 0
                ):

                    aligned_effect = math.log(
                        odds_ratio
                    )

                    effect_type = "LOG_OR"

                    status = statuses[
                        "usable"
                    ]

                else:

                    status = statuses[
                        "missing_effect"
                    ]

            elif (
                effect
                ==
                canonical_other
                and
                other
                ==
                canonical_effect
            ):

                orientation = "REVERSED"

                if beta is not None:

                    aligned_effect = -beta
                    effect_type = "BETA"

                    status = statuses[
                        "usable"
                    ]

                elif (
                    odds_ratio is not None
                    and
                    odds_ratio > 0
                ):

                    aligned_effect = -math.log(
                        odds_ratio
                    )

                    effect_type = "LOG_OR"

                    status = statuses[
                        "usable"
                    ]

                else:

                    status = statuses[
                        "missing_effect"
                    ]

            else:

                status = statuses[
                    "allele_pair_mismatch"
                ]

            rows.append(
                _effect_row(
                    record=record,
                    canonical_effect=canonical_effect,
                    canonical_other=canonical_other,
                    allele_orientation=orientation,
                    effect_type=effect_type,
                    aligned_effect=aligned_effect,
                    status=status,
                )
            )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Publication-family directional aggregation
# ============================================================================


def build_publication_family_direction(
    *,
    candidate_effect_harmonization: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Aggregate directions within candidate × publication family.

    Multiple accessions from one publication are not treated as multiple
    directional replication units.
    """

    statuses = config[
        "publication_family_direction"
    ][
        "statuses"
    ]

    usable = candidate_effect_harmonization.loc[
        candidate_effect_harmonization[
            "effect_direction_usable"
        ]
        .fillna(False)
        .astype(bool)
    ].copy()

    rows: list[dict[str, Any]] = []

    grouped = usable.groupby(
        [
            "rsid",
            "publication_family_id",
            "pubmed_id",
        ],
        dropna=False,
        sort=True,
    )

    for (
        rsid,
        family_id,
        pubmed_id,
    ), group in grouped:

        directions = sorted(
            set(
                group[
                    "aligned_direction"
                ]
                .dropna()
                .astype(str)
            )
        )

        accessions = sorted(
            set(
                group[
                    "accession_id"
                ].astype(str)
            )
        )

        if directions == [
            "POSITIVE"
        ]:

            family_direction = "POSITIVE"

            status = statuses[
                "positive"
            ]

        elif directions == [
            "NEGATIVE"
        ]:

            family_direction = "NEGATIVE"

            status = statuses[
                "negative"
            ]

        elif len(
            directions
        ) > 1:

            family_direction = "MIXED"

            status = statuses[
                "mixed"
            ]

        else:

            family_direction = None

            status = statuses[
                "unavailable"
            ]

        rows.append(
            {
                "rsid":
                    str(
                        rsid
                    ),

                "publication_family_id":
                    str(
                        family_id
                    ),

                "pubmed_id":
                    pubmed_id,

                "study_count":
                    len(
                        accessions
                    ),

                "study_accessions":
                    accessions,

                "usable_effect_rows":
                    int(
                        len(
                            group
                        )
                    ),

                "observed_directions":
                    directions,

                "publication_family_direction":
                    family_direction,

                "within_publication_direction_mixed":
                    bool(
                        family_direction
                        ==
                        "MIXED"
                    ),

                "publication_family_direction_status":
                    status,

                "independent_replication_unit":
                    False,

                "sample_independence_verified":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Candidate publication-direction resolution
# ============================================================================


def build_candidate_direction_resolution(
    *,
    publication_family_direction: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve within- and cross-publication directional heterogeneity."""

    statuses = config[
        "candidate_direction"
    ][
        "statuses"
    ]

    rows: list[dict[str, Any]] = []

    for candidate in config[
        "candidates"
    ]:

        rsid = str(
            candidate[
                "rsid"
            ]
        )

        subset = publication_family_direction.loc[
            publication_family_direction[
                "rsid"
            ]
            ==
            rsid
        ].copy()

        mixed_count = int(
            (
                subset[
                    "publication_family_direction"
                ]
                ==
                "MIXED"
            ).sum()
        )

        positive_count = int(
            (
                subset[
                    "publication_family_direction"
                ]
                ==
                "POSITIVE"
            ).sum()
        )

        negative_count = int(
            (
                subset[
                    "publication_family_direction"
                ]
                ==
                "NEGATIVE"
            ).sum()
        )

        clean_directions = set()

        if positive_count > 0:

            clean_directions.add(
                "POSITIVE"
            )

        if negative_count > 0:

            clean_directions.add(
                "NEGATIVE"
            )

        within_publication_mixed = bool(
            mixed_count
            >
            0
        )

        cross_publication_mixed = bool(
            len(
                clean_directions
            )
            >
            1
        )

        family_count = int(
            len(
                subset
            )
        )

        if family_count <= 1:

            status = statuses[
                "single_publication_family"
            ]

        elif (
            within_publication_mixed
            and
            cross_publication_mixed
        ):

            status = statuses[
                "within_and_cross_publication_heterogeneity"
            ]

        elif within_publication_mixed:

            status = statuses[
                "within_publication_heterogeneity"
            ]

        elif cross_publication_mixed:

            status = statuses[
                "cross_publication_heterogeneity"
            ]

        elif (
            len(
                clean_directions
            )
            ==
            1
        ):

            status = statuses[
                "multi_publication_concordant"
            ]

        else:

            status = statuses[
                "direction_unresolved"
            ]

        rows.append(
            {
                "rsid":
                    rsid,

                "priority_rank":
                    int(
                        candidate[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        candidate[
                            "priority_class"
                        ]
                    ),

                "publication_families_assessed":
                    family_count,

                "publication_families_direction_positive":
                    positive_count,

                "publication_families_direction_negative":
                    negative_count,

                "publication_families_direction_mixed":
                    mixed_count,

                "within_publication_directional_heterogeneity":
                    within_publication_mixed,

                "cross_publication_directional_heterogeneity":
                    cross_publication_mixed,

                "independent_replication_verified":
                    False,

                "study_independence_verified":
                    False,

                "direction_resolution_status":
                    status,
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            "priority_rank",
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


# ============================================================================
# Final candidate resolution
# ============================================================================


def build_candidate_resolution(
    *,
    candidate_study_harmonization: pd.DataFrame,
    candidate_effect_harmonization: pd.DataFrame,
    candidate_direction_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build final candidate-level M7.5C evidence table."""

    rows: list[dict[str, Any]] = []

    for candidate in config[
        "candidates"
    ]:

        rsid = str(
            candidate[
                "rsid"
            ]
        )

        study_rows = candidate_study_harmonization.loc[
            candidate_study_harmonization[
                "rs_id"
            ]
            ==
            rsid
        ]

        effect_rows = candidate_effect_harmonization.loc[
            candidate_effect_harmonization[
                "rsid"
            ]
            ==
            rsid
        ]

        direction_row = candidate_direction_resolution.loc[
            candidate_direction_resolution[
                "rsid"
            ]
            ==
            rsid
        ].iloc[
            0
        ]

        accessions = sorted(
            set(
                study_rows[
                    "accession_id"
                ].astype(str)
            )
        )

        families = sorted(
            set(
                study_rows[
                    "publication_family_id"
                ].astype(str)
            )
        )

        usable = effect_rows.loc[
            effect_rows[
                "effect_direction_usable"
            ]
            .fillna(False)
            .astype(bool)
        ]

        rows.append(
            {
                "rsid":
                    rsid,

                "priority_rank":
                    int(
                        candidate[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        candidate[
                            "priority_class"
                        ]
                    ),

                "studies_with_exact_variant":
                    len(
                        accessions
                    ),

                "publication_families_with_exact_variant":
                    len(
                        families
                    ),

                "study_accessions":
                    accessions,

                "publication_family_ids":
                    families,

                "directionally_usable_rows":
                    int(
                        len(
                            usable
                        )
                    ),

                "publication_families_direction_positive":
                    int(
                        direction_row[
                            "publication_families_direction_positive"
                        ]
                    ),

                "publication_families_direction_negative":
                    int(
                        direction_row[
                            "publication_families_direction_negative"
                        ]
                    ),

                "publication_families_direction_mixed":
                    int(
                        direction_row[
                            "publication_families_direction_mixed"
                        ]
                    ),

                "within_publication_directional_heterogeneity":
                    bool(
                        direction_row[
                            "within_publication_directional_heterogeneity"
                        ]
                    ),

                "cross_publication_directional_heterogeneity":
                    bool(
                        direction_row[
                            "cross_publication_directional_heterogeneity"
                        ]
                    ),

                "independent_replication_verified":
                    False,

                "study_independence_verified":
                    False,

                "candidate_status":
                    str(
                        direction_row[
                            "direction_resolution_status"
                        ]
                    ),
            }
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            "priority_rank",
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    publication_family_inventory: pd.DataFrame,
    study_pair_overlap: pd.DataFrame,
    candidate_effect_harmonization: pd.DataFrame,
    candidate_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.5C summary."""

    related_pairs = (
        int(
            study_pair_overlap[
                "potential_sample_overlap"
            ]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        if not study_pair_overlap.empty
        else 0
    )

    usable_effect_rows = int(
        candidate_effect_harmonization[
            "effect_direction_usable"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    within_heterogeneous = int(
        candidate_resolution[
            "within_publication_directional_heterogeneity"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    cross_heterogeneous = int(
        candidate_resolution[
            "cross_publication_directional_heterogeneity"
        ]
        .fillna(False)
        .astype(bool)
        .sum()
    )

    any_heterogeneous = bool(
        within_heterogeneous > 0
        or
        cross_heterogeneous > 0
    )

    statuses = config[
        "overall_statuses"
    ]

    if any_heterogeneous:

        overall_status = statuses[
            "heterogeneous"
        ]

    else:

        concordant_status = config[
            "candidate_direction"
        ][
            "statuses"
        ][
            "multi_publication_concordant"
        ]

        concordant_count = int(
            (
                candidate_resolution[
                    "candidate_status"
                ]
                ==
                concordant_status
            ).sum()
        )

        if concordant_count > 0:

            overall_status = statuses[
                "concordant"
            ]

        else:

            overall_status = statuses[
                "direction_unresolved"
            ]

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.5C",

                "candidate_variants_assessed":
                    int(
                        len(
                            candidate_resolution
                        )
                    ),

                "publication_families_identified":
                    int(
                        len(
                            publication_family_inventory
                        )
                    ),

                "study_pairs_assessed":
                    int(
                        len(
                            study_pair_overlap
                        )
                    ),

                "study_pairs_with_explicit_relatedness":
                    related_pairs,

                "effect_rows_assessed":
                    int(
                        len(
                            candidate_effect_harmonization
                        )
                    ),

                "directionally_usable_effect_rows":
                    usable_effect_rows,

                "directionally_unusable_effect_rows":
                    int(
                        len(
                            candidate_effect_harmonization
                        )
                        -
                        usable_effect_rows
                    ),

                "candidates_with_within_publication_heterogeneity":
                    within_heterogeneous,

                "candidates_with_cross_publication_heterogeneity":
                    cross_heterogeneous,

                "independent_replications_verified":
                    0,

                "independent_replication_claimed":
                    False,

                "study_independence_verified":
                    False,

                "different_accession_used_as_independence":
                    False,

                "different_pubmed_used_as_proof_of_independence":
                    False,

                "cohort_nonoverlap_used_as_proof_of_independence":
                    False,

                "nonharmonised_sources_used_for_directional_confirmation":
                    False,

                "p_value_used_to_select_direction":
                    False,

                "meta_analysis_performed":
                    False,

                "pooled_effect_estimated":
                    False,

                "overall_status":
                    overall_status,

                "next_stage":
                    config[
                        "next_stage"
                    ],
            }
        ]
    )


# ============================================================================
# Validation
# ============================================================================


def _validate_result(
    *,
    candidate_effect_harmonization: pd.DataFrame,
    candidate_resolution: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict[str, Any],
) -> None:
    """Enforce scientific safeguards."""

    expected_rsids = {
        str(value)
        for value in config[
            "validation"
        ][
            "expected_rsids"
        ]
    }

    observed_rsids = set(
        candidate_resolution[
            "rsid"
        ].astype(str)
    )

    if observed_rsids != expected_rsids:

        raise RuntimeError(
            "M7.5C candidate rsID mismatch."
        )

    if len(
        candidate_resolution
    ) != int(
        config[
            "validation"
        ][
            "expected_candidate_count"
        ]
    ):

        raise RuntimeError(
            "M7.5C candidate count mismatch."
        )

    if len(
        candidate_effect_harmonization
    ) != int(
        config[
            "validation"
        ][
            "expected_exact_hit_rows"
        ]
    ):

        raise RuntimeError(
            "M7.5C effect-row count mismatch."
        )

    if len(
        summary
    ) != 1:

        raise RuntimeError(
            "M7.5C summary must contain exactly one row."
        )

    row = summary.iloc[
        0
    ]

    forbidden_true = [
        "independent_replication_claimed",
        "study_independence_verified",
        "different_accession_used_as_independence",
        "different_pubmed_used_as_proof_of_independence",
        "cohort_nonoverlap_used_as_proof_of_independence",
        "nonharmonised_sources_used_for_directional_confirmation",
        "p_value_used_to_select_direction",
        "meta_analysis_performed",
        "pooled_effect_estimated",
    ]

    for field in forbidden_true:

        if bool(
            row[
                field
            ]
        ):

            raise RuntimeError(
                f"M7.5C safeguard violation: {field}=True."
            )

    if int(
        row[
            "independent_replications_verified"
        ]
    ) != 0:

        raise RuntimeError(
            "M7.5C cannot verify independent replication."
        )


# ============================================================================
# Public API
# ============================================================================


def assess_study_independence_harmonization(
    *,
    study_inventory: pd.DataFrame,
    file_manifest: pd.DataFrame,
    exact_hits: pd.DataFrame,
    config: dict[str, Any],
) -> M75CResult:
    """Execute complete M7.5C."""

    publication_family_inventory = (
        build_publication_family_inventory(
            study_inventory=study_inventory,
            config=config,
        )
    )

    study_pair_overlap = (
        build_study_pair_overlap(
            study_inventory=study_inventory,
            config=config,
        )
    )

    candidate_study_harmonization = (
        build_candidate_study_harmonization(
            exact_hits=exact_hits,
            study_inventory=study_inventory,
            file_manifest=file_manifest,
            config=config,
        )
    )

    candidate_effect_harmonization = (
        build_candidate_effect_harmonization(
            candidate_study_harmonization=(
                candidate_study_harmonization
            ),
            config=config,
        )
    )

    publication_family_direction = (
        build_publication_family_direction(
            candidate_effect_harmonization=(
                candidate_effect_harmonization
            ),
            config=config,
        )
    )

    candidate_direction_resolution = (
        build_candidate_direction_resolution(
            publication_family_direction=(
                publication_family_direction
            ),
            config=config,
        )
    )

    candidate_resolution = (
        build_candidate_resolution(
            candidate_study_harmonization=(
                candidate_study_harmonization
            ),
            candidate_effect_harmonization=(
                candidate_effect_harmonization
            ),
            candidate_direction_resolution=(
                candidate_direction_resolution
            ),
            config=config,
        )
    )

    summary = build_summary(
        publication_family_inventory=(
            publication_family_inventory
        ),
        study_pair_overlap=(
            study_pair_overlap
        ),
        candidate_effect_harmonization=(
            candidate_effect_harmonization
        ),
        candidate_resolution=(
            candidate_resolution
        ),
        config=config,
    )

    _validate_result(
        candidate_effect_harmonization=(
            candidate_effect_harmonization
        ),
        candidate_resolution=(
            candidate_resolution
        ),
        summary=summary,
        config=config,
    )

    return M75CResult(
        publication_family_inventory=(
            publication_family_inventory
        ),
        study_pair_overlap=(
            study_pair_overlap
        ),
        candidate_study_harmonization=(
            candidate_study_harmonization
        ),
        candidate_effect_harmonization=(
            candidate_effect_harmonization
        ),
        publication_family_direction=(
            publication_family_direction
        ),
        candidate_direction_resolution=(
            candidate_direction_resolution
        ),
        candidate_resolution=(
            candidate_resolution
        ),
        summary=summary,
    )