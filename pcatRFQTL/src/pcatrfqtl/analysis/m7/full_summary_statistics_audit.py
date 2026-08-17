"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/full_summary_statistics_audit.py

Description:
    Core evidence-resolution logic for M7.5B Full Summary Statistics Audit.

    This module performs no network operations.

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


@dataclass(frozen=True)
class FullSummaryStatisticsAuditResult:
    """Container for M7.5B processed results."""

    study_candidate_audit: pd.DataFrame
    candidate_resolution: pd.DataFrame
    summary: pd.DataFrame


STUDY_CANDIDATE_COLUMNS = [
    "accession_id",
    "rsid",
    "study_reused_from_m5",
    "study_eligible",
    "summary_statistics_file_auditable",
    "exact_variant_found",
    "exact_hit_rows",
    "independent_replication_established",
    "study_independence_verified",
    "phenotype_compatibility_verified",
    "ancestry_compatibility_verified",
    "effect_allele_harmonization_performed",
    "association_status",
]


def _candidate_rsids(
    config: dict[str, Any],
) -> list[str]:
    """Return candidate rsIDs ordered by priority."""

    return [
        str(
            row["rsid"]
        )
        for row in sorted(
            config["candidates"],
            key=lambda row: int(
                row["priority_rank"]
            ),
        )
    ]


def build_study_candidate_audit(
    *,
    study_inventory: pd.DataFrame,
    file_manifest: pd.DataFrame,
    exact_variant_hits: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one row per eligible study × candidate."""

    statuses = config[
        "statuses"
    ][
        "association"
    ]

    previously_used = {
        str(value)
        for value in config[
            "previously_used_gwas"
        ][
            "study_accessions"
        ]
    }

    if study_inventory.empty:

        return pd.DataFrame(
            columns=STUDY_CANDIDATE_COLUMNS
        )

    eligible = study_inventory.loc[
        study_inventory[
            "eligible_for_exact_lookup"
        ]
        .fillna(False)
        .astype(bool)
    ].copy()

    rows: list[dict[str, Any]] = []

    for _, study in eligible.iterrows():

        accession = str(
            study["accession_id"]
        )

        if file_manifest.empty:

            auditable = False

        else:

            selected = file_manifest.loc[
                (
                    file_manifest[
                        "accession_id"
                    ].astype(str)
                    ==
                    accession
                )
                &
                file_manifest[
                    "selected_for_scan"
                ]
                .fillna(False)
                .astype(bool)
            ]

            auditable = bool(
                not selected.empty
                and
                selected[
                    "scan_success"
                ]
                .fillna(False)
                .astype(bool)
                .any()
            )

        for rsid in _candidate_rsids(
            config
        ):

            if exact_variant_hits.empty:

                hits = exact_variant_hits

            else:

                hits = exact_variant_hits.loc[
                    (
                        exact_variant_hits[
                            "accession_id"
                        ].astype(str)
                        ==
                        accession
                    )
                    &
                    (
                        exact_variant_hits[
                            "rs_id"
                        ].astype(str)
                        ==
                        rsid
                    )
                ]

            found = bool(
                not hits.empty
            )

            reused = bool(
                accession in previously_used
            )

            if not auditable:

                status = statuses[
                    "not_auditable"
                ]

            elif found and reused:

                status = statuses[
                    "exact_reused"
                ]

            elif found:

                status = statuses[
                    "exact_nonreused"
                ]

            else:

                status = statuses[
                    "not_found"
                ]

            rows.append(
                {
                    "accession_id":
                        accession,

                    "rsid":
                        rsid,

                    "study_reused_from_m5":
                        reused,

                    "study_eligible":
                        True,

                    "summary_statistics_file_auditable":
                        auditable,

                    "exact_variant_found":
                        found,

                    "exact_hit_rows":
                        int(
                            len(hits)
                        ),

                    "independent_replication_established":
                        False,

                    "study_independence_verified":
                        False,

                    "phenotype_compatibility_verified":
                        False,

                    "ancestry_compatibility_verified":
                        False,

                    "effect_allele_harmonization_performed":
                        False,

                    "association_status":
                        status,
                }
            )

    return pd.DataFrame(
        rows,
        columns=STUDY_CANDIDATE_COLUMNS,
    )


def build_candidate_resolution(
    *,
    study_candidate_audit: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve exact-variant evidence for each candidate."""

    statuses = config[
        "statuses"
    ][
        "candidate"
    ]

    rows: list[dict[str, Any]] = []

    for candidate in config[
        "candidates"
    ]:

        rsid = str(
            candidate["rsid"]
        )

        subset = (
            study_candidate_audit.loc[
                study_candidate_audit[
                    "rsid"
                ]
                ==
                rsid
            ].copy()
            if not study_candidate_audit.empty
            else study_candidate_audit.copy()
        )

        if subset.empty:

            auditable = subset
            found = subset

        else:

            auditable = subset.loc[
                subset[
                    "summary_statistics_file_auditable"
                ]
                .fillna(False)
                .astype(bool)
            ]

            found = subset.loc[
                subset[
                    "exact_variant_found"
                ]
                .fillna(False)
                .astype(bool)
            ]

        eligible_count = (
            int(
                subset[
                    "accession_id"
                ].nunique()
            )
            if not subset.empty
            else 0
        )

        auditable_count = (
            int(
                auditable[
                    "accession_id"
                ].nunique()
            )
            if not auditable.empty
            else 0
        )

        found_accessions = (
            sorted(
                set(
                    found[
                        "accession_id"
                    ].astype(str)
                )
            )
            if not found.empty
            else []
        )

        reused_accessions = (
            sorted(
                set(
                    found.loc[
                        found[
                            "study_reused_from_m5"
                        ]
                        .fillna(False)
                        .astype(bool),
                        "accession_id",
                    ].astype(str)
                )
            )
            if not found.empty
            else []
        )

        nonreused_accessions = (
            sorted(
                set(
                    found.loc[
                        ~found[
                            "study_reused_from_m5"
                        ]
                        .fillna(False)
                        .astype(bool),
                        "accession_id",
                    ].astype(str)
                )
            )
            if not found.empty
            else []
        )

        if nonreused_accessions:

            candidate_status = statuses[
                "nonreused_found"
            ]

        elif reused_accessions:

            candidate_status = statuses[
                "reused_only"
            ]

        elif (
            eligible_count > 0
            and
            auditable_count == eligible_count
        ):

            candidate_status = statuses[
                "not_found"
            ]

        else:

            candidate_status = statuses[
                "incomplete"
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

                "eligible_studies_assessed":
                    eligible_count,

                "auditable_studies":
                    auditable_count,

                "studies_with_exact_variant":
                    len(
                        found_accessions
                    ),

                "reused_studies_with_exact_variant":
                    len(
                        reused_accessions
                    ),

                "nonreused_studies_with_exact_variant":
                    len(
                        nonreused_accessions
                    ),

                "reused_study_accessions":
                    reused_accessions,

                "nonreused_study_accessions":
                    nonreused_accessions,

                "exact_variant_evidence_found":
                    bool(
                        found_accessions
                    ),

                "independent_replication_established":
                    False,

                "study_independence_verified":
                    False,

                "phenotype_compatibility_verified":
                    False,

                "ancestry_compatibility_verified":
                    False,

                "effect_allele_harmonization_performed":
                    False,

                "candidate_status":
                    candidate_status,
            }
        )

    return (
        pd.DataFrame(rows)
        .sort_values(
            "priority_rank",
            kind="stable",
        )
        .reset_index(
            drop=True
        )
    )


def build_summary(
    *,
    study_inventory: pd.DataFrame,
    file_manifest: pd.DataFrame,
    exact_variant_hits: pd.DataFrame,
    candidate_resolution: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.5B summary."""

    studies_discovered = int(
        len(study_inventory)
    )

    full_available = (
        int(
            study_inventory[
                "full_summary_stats_available"
            ]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        if not study_inventory.empty
        else 0
    )

    direct_studies = (
        int(
            study_inventory[
                "direct_prostate_phenotype"
            ]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        if not study_inventory.empty
        else 0
    )

    eligible_studies = (
        int(
            study_inventory[
                "eligible_for_exact_lookup"
            ]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        if not study_inventory.empty
        else 0
    )

    excluded_or_unresolved = (
        full_available
        -
        eligible_studies
    )

    selected_files = (
        int(
            file_manifest[
                "selected_for_scan"
            ]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        if not file_manifest.empty
        else 0
    )

    successful_files = (
        int(
            file_manifest[
                "scan_success"
            ]
            .fillna(False)
            .astype(bool)
            .sum()
        )
        if not file_manifest.empty
        else 0
    )

    nonreused_candidates = int(
        (
            candidate_resolution[
                "nonreused_studies_with_exact_variant"
            ]
            >
            0
        ).sum()
    )

    reused_only_candidates = int(
        (
            (
                candidate_resolution[
                    "reused_studies_with_exact_variant"
                ]
                >
                0
            )
            &
            (
                candidate_resolution[
                    "nonreused_studies_with_exact_variant"
                ]
                ==
                0
            )
        ).sum()
    )

    candidate_incomplete = bool(
        (
            candidate_resolution[
                "candidate_status"
            ]
            ==
            config[
                "statuses"
            ][
                "candidate"
            ][
                "incomplete"
            ]
        ).any()
    )

    audit_complete = bool(
        eligible_studies > 0
        and
        selected_files == eligible_studies
        and
        successful_files == selected_files
        and
        not candidate_incomplete
    )

    statuses = config[
        "statuses"
    ][
        "overall"
    ]

    if not audit_complete:

        overall_status = statuses[
            "incomplete"
        ]

    elif nonreused_candidates > 0:

        overall_status = statuses[
            "nonreused_evidence"
        ]

    elif reused_only_candidates > 0:

        overall_status = statuses[
            "reused_only"
        ]

    else:

        overall_status = statuses[
            "no_exact_variants"
        ]

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.5B",

                "candidate_variants_assessed":
                    int(
                        len(
                            candidate_resolution
                        )
                    ),

                "prostate_related_studies_discovered":
                    studies_discovered,

                "full_summary_statistics_available_studies":
                    full_available,

                "direct_prostate_phenotype_studies":
                    direct_studies,

                "phenotype_excluded_or_unresolved_studies":
                    excluded_or_unresolved,

                "eligible_full_sumstats_studies":
                    eligible_studies,

                "summary_statistics_files_selected":
                    selected_files,

                "summary_statistics_files_scanned_successfully":
                    successful_files,

                "exact_variant_rows_found":
                    int(
                        len(
                            exact_variant_hits
                        )
                    ),

                "candidates_with_nonreused_exact_variant_evidence":
                    nonreused_candidates,

                "candidates_with_reused_only_exact_variant_evidence":
                    reused_only_candidates,

                "full_summary_statistics_audit_complete":
                    audit_complete,

                "independent_replications_verified":
                    0,

                "independent_replication_claimed":
                    False,

                "study_independence_verified":
                    False,

                "phenotype_compatibility_verified":
                    False,

                "ancestry_compatibility_verified":
                    False,

                "effect_allele_harmonization_performed":
                    False,

                "effect_direction_concordance_claimed":
                    False,

                "proxy_variant_lookup_performed":
                    False,

                "coordinate_substitution_performed":
                    False,

                "candidate_substitution_performed":
                    False,

                "excluded_phenotype_means_negative_evidence":
                    False,

                "candidate_not_found_means_no_association":
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


def _validate_result(
    *,
    study_inventory: pd.DataFrame,
    candidate_resolution: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict[str, Any],
) -> None:
    """Enforce M7.5B structural and scientific safeguards."""

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
            "M7.5B candidate rsID set mismatch."
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
            "M7.5B candidate count mismatch."
        )

    if len(summary) != 1:

        raise RuntimeError(
            "M7.5B summary must contain exactly one row."
        )

    summary_row = summary.iloc[0]

    forbidden_true = [
        "independent_replication_claimed",
        "study_independence_verified",
        "phenotype_compatibility_verified",
        "ancestry_compatibility_verified",
        "effect_allele_harmonization_performed",
        "effect_direction_concordance_claimed",
        "proxy_variant_lookup_performed",
        "coordinate_substitution_performed",
        "candidate_substitution_performed",
        "excluded_phenotype_means_negative_evidence",
        "candidate_not_found_means_no_association",
    ]

    for field in forbidden_true:

        if bool(
            summary_row[field]
        ):

            raise RuntimeError(
                f"M7.5B safeguard violation: {field}=True."
            )

    if int(
        summary_row[
            "independent_replications_verified"
        ]
    ) != 0:

        raise RuntimeError(
            "M7.5B cannot verify independent replication."
        )

    if not study_inventory.empty:

        eligible = study_inventory[
            "eligible_for_exact_lookup"
        ].fillna(False).astype(bool)

        full_available = study_inventory[
            "full_summary_stats_available"
        ].fillna(False).astype(bool)

        direct = study_inventory[
            "direct_prostate_phenotype"
        ].fillna(False).astype(bool)

        if (
            eligible
            &
            ~full_available
        ).any():

            raise RuntimeError(
                "Eligible study lacks full summary statistics."
            )

        if (
            eligible
            &
            ~direct
        ).any():

            raise RuntimeError(
                "Non-direct phenotype entered M7.5B exact lookup."
            )


def assess_full_summary_statistics_audit(
    *,
    study_inventory: pd.DataFrame,
    file_manifest: pd.DataFrame,
    exact_variant_hits: pd.DataFrame,
    config: dict[str, Any],
) -> FullSummaryStatisticsAuditResult:
    """Execute M7.5B core evidence resolution."""

    study_candidate_audit = build_study_candidate_audit(
        study_inventory=study_inventory,
        file_manifest=file_manifest,
        exact_variant_hits=exact_variant_hits,
        config=config,
    )

    candidate_resolution = build_candidate_resolution(
        study_candidate_audit=study_candidate_audit,
        config=config,
    )

    summary = build_summary(
        study_inventory=study_inventory,
        file_manifest=file_manifest,
        exact_variant_hits=exact_variant_hits,
        candidate_resolution=candidate_resolution,
        config=config,
    )

    _validate_result(
        study_inventory=study_inventory,
        candidate_resolution=candidate_resolution,
        summary=summary,
        config=config,
    )

    return FullSummaryStatisticsAuditResult(
        study_candidate_audit=study_candidate_audit,
        candidate_resolution=candidate_resolution,
        summary=summary,
    )