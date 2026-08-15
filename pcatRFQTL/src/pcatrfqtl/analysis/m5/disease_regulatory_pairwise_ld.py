"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/disease_regulatory_pairwise_ld.py

Description:
    M5.3C.3G disease-proxy ↔ regulatory-proxy pairwise LD utilities.

    This stage evaluates direct pairwise LD between:
        - suggestive prostate cancer GWAS variants identified in M5.3C.3E/F;
        - Moradi regulatory LD proxy variants identified previously in M5.4.

    The main use case is rs10216902, where disease-associated regional
    variants and regulatory proxies coexist but do not form a lead-centered
    LD bridge.

    Integration is canonical-rsID based.

    LD interpretation:
        HIGH:
            r² >= 0.80

        MODERATE:
            0.50 <= r² < 0.80

        LOW:
            r² < 0.50

        UNAVAILABLE:
            pairwise LD could not be retrieved or the queried variant was
            unavailable in the reference panel.

    Scientific safeguards:
        - Physical proximity is not interpreted as LD.
        - Shared broad locus is not interpreted as a shared causal signal.
        - Pairwise LD does not prove colocalization.
        - Cross-population consistency is sensitivity evidence, not
          independent replication.
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

from dataclasses import dataclass
from typing import Any

import pandas as pd


HIGH_LD_R2 = 0.80
MODERATE_LD_R2 = 0.50


@dataclass(frozen=True)
class PairwiseCandidateBuildResult:
    """Candidate disease × regulatory variant-pair scaffold."""

    dataframe: pd.DataFrame

    disease_variants: int

    regulatory_variants: int

    candidate_pairs: int


@dataclass(frozen=True)
class PairwiseLDSynthesisResult:
    """Cross-population pairwise LD synthesis."""

    dataframe: pd.DataFrame

    tested_pairs: int

    high_ld_pairs: int

    moderate_or_high_ld_pairs: int


def _normalize_rsid(
    series: pd.Series,
) -> pd.Series:
    """Normalize canonical rsIDs."""

    result = (
        series
        .astype("string")
        .str.strip()
        .str.lower()
    )

    valid = result.str.fullmatch(
        r"rs\d+",
        na=False,
    )

    return result.where(
        valid,
        pd.NA,
    )


def _resolve_column(
    dataframe: pd.DataFrame,
    aliases: tuple[str, ...],
    *,
    required: bool = True,
) -> str | None:
    """Resolve the first matching column alias."""

    normalized = {
        str(column).strip().lower():
            str(column)
        for column in dataframe.columns
    }

    for alias in aliases:

        match = normalized.get(
            alias.lower()
        )

        if match is not None:
            return match

    if required:

        raise ValueError(
            "Unable to resolve required column. "
            f"Expected one of {aliases}; "
            f"available={list(dataframe.columns)}"
        )

    return None


def prepare_disease_variants(
    recurrence: pd.DataFrame,
    *,
    lead_rsid: str = "rs10216902",
) -> pd.DataFrame:
    """
    Prepare suggestive disease variants for one candidate locus.

    This function keeps one row per disease variant identity.
    """

    required = {
        "lead_rsid",
        "variant_identity_key",
        "rsid",
        "minimum_p_value",
        "study_count",
        "genome_wide_in_any_study",
    }

    missing = (
        required
        - set(
            recurrence.columns
        )
    )

    if missing:

        raise ValueError(
            "Disease recurrence table missing columns: "
            f"{sorted(missing)}"
        )

    normalized_lead = (
        str(
            lead_rsid
        )
        .strip()
        .lower()
    )

    out = recurrence.copy()

    out[
        "lead_rsid"
    ] = _normalize_rsid(
        out[
            "lead_rsid"
        ]
    )

    out[
        "disease_proxy_rsid"
    ] = _normalize_rsid(
        out[
            "rsid"
        ]
    )

    out = out.loc[
        out[
            "lead_rsid"
        ].eq(
            normalized_lead
        )
        &
        out[
            "disease_proxy_rsid"
        ].notna()
    ].copy()

    return out


def prepare_regulatory_variants(
    regulatory_bridges: dict[
        str,
        pd.DataFrame,
    ],
    *,
    lead_rsid: str = "rs10216902",
) -> pd.DataFrame:
    """
    Collapse regulatory proxies observed across population-specific
    regulatory bridge tables.
    """

    normalized_lead = (
        str(
            lead_rsid
        )
        .strip()
        .lower()
    )

    records: list[
        dict[str, Any]
    ] = []

    for population, dataframe in regulatory_bridges.items():

        if dataframe.empty:
            continue

        lead_column = _resolve_column(
            dataframe,
            (
                "lead_rsid",
                "candidate_rsid",
                "candidate_lead_rsid",
            ),
        )

        proxy_column = _resolve_column(
            dataframe,
            (
                "proxy_rsid",
                "ld_proxy_rsid",
                "neighbor_rsid",
                "bridge_rsid",
                "variant_rsid",
                "regulatory_proxy_rsid",
            ),
        )

        feature_column = _resolve_column(
            dataframe,
            (
                "feature_id",
                "feature",
                "normalized_feature_id",
                "m4_feature_id",
            ),
            required=False,
        )

        source_column = _resolve_column(
            dataframe,
            (
                "source_partition",
                "partition",
                "moradi_partition",
                "source_table",
            ),
            required=False,
        )

        current = dataframe.copy()

        current[
            "_lead"
        ] = _normalize_rsid(
            current[
                lead_column
            ]
        )

        current[
            "_proxy"
        ] = _normalize_rsid(
            current[
                proxy_column
            ]
        )

        current = current.loc[
            current[
                "_lead"
            ].eq(
                normalized_lead
            )
            &
            current[
                "_proxy"
            ].notna()
        ]

        for _, row in current.iterrows():

            records.append(
                {
                    "lead_rsid":
                        normalized_lead,

                    "regulatory_proxy_rsid":
                        str(
                            row[
                                "_proxy"
                            ]
                        ),

                    "regulatory_source_population":
                        str(
                            population
                        ).upper(),

                    "regulatory_feature_id":
                        (
                            None
                            if (
                                feature_column is None
                                or pd.isna(
                                    row[
                                        feature_column
                                    ]
                                )
                            )
                            else str(
                                row[
                                    feature_column
                                ]
                            )
                        ),

                    "regulatory_source_partition":
                        (
                            None
                            if (
                                source_column is None
                                or pd.isna(
                                    row[
                                        source_column
                                    ]
                                )
                            )
                            else str(
                                row[
                                    source_column
                                ]
                            )
                        ),
                }
            )

    if not records:

        return pd.DataFrame(
            columns=[
                "lead_rsid",
                "regulatory_proxy_rsid",
                "regulatory_populations",
                "regulatory_population_count",
                "regulatory_feature_ids",
                "regulatory_source_partitions",
            ]
        )

    raw = pd.DataFrame(
        records
    )

    collapsed_records: list[
        dict[str, Any]
    ] = []

    for (
        lead,
        proxy,
    ), group in raw.groupby(
        [
            "lead_rsid",
            "regulatory_proxy_rsid",
        ],
        sort=True,
    ):

        populations = sorted(
            group[
                "regulatory_source_population"
            ]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        feature_ids = sorted(
            group[
                "regulatory_feature_id"
            ]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        partitions = sorted(
            group[
                "regulatory_source_partition"
            ]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        collapsed_records.append(
            {
                "lead_rsid":
                    str(
                        lead
                    ),

                "regulatory_proxy_rsid":
                    str(
                        proxy
                    ),

                "regulatory_populations":
                    "|".join(
                        populations
                    ),

                "regulatory_population_count":
                    int(
                        len(
                            populations
                        )
                    ),

                "regulatory_feature_ids":
                    "|".join(
                        feature_ids
                    ),

                "regulatory_source_partitions":
                    "|".join(
                        partitions
                    ),
            }
        )

    return pd.DataFrame(
        collapsed_records
    )


def build_pairwise_candidate_scaffold(
    *,
    recurrence: pd.DataFrame,
    regulatory_bridges: dict[
        str,
        pd.DataFrame,
    ],
    lead_rsid: str = "rs10216902",
) -> PairwiseCandidateBuildResult:
    """
    Generate every disease-proxy × regulatory-proxy pair for one lead.

    Self-pairs are excluded because same-proxy identity has already been
    evaluated in M5.3C.3F.
    """

    disease = prepare_disease_variants(
        recurrence,
        lead_rsid=lead_rsid,
    )

    regulatory = prepare_regulatory_variants(
        regulatory_bridges,
        lead_rsid=lead_rsid,
    )

    if disease.empty or regulatory.empty:

        scaffold = pd.DataFrame()

    else:

        disease = disease.copy()

        regulatory = regulatory.copy()

        disease[
            "_join_key"
        ] = 1

        regulatory[
            "_join_key"
        ] = 1

        scaffold = disease.merge(
            regulatory,
            on=[
                "_join_key",
                "lead_rsid",
            ],
            how="inner",
        )

        scaffold = scaffold.drop(
            columns=[
                "_join_key",
            ]
        )

        scaffold[
            "is_same_proxy"
        ] = (
            scaffold[
                "disease_proxy_rsid"
            ]
            .eq(
                scaffold[
                    "regulatory_proxy_rsid"
                ]
            )
            .fillna(
                False
            )
        )

        scaffold = scaffold.loc[
            ~scaffold[
                "is_same_proxy"
            ]
        ].copy()

        scaffold[
            "pair_id"
        ] = (
            scaffold[
                "disease_proxy_rsid"
            ].astype("string")
            + "__"
            + scaffold[
                "regulatory_proxy_rsid"
            ].astype("string")
        )

        scaffold = (
            scaffold
            .drop_duplicates(
                subset=[
                    "lead_rsid",
                    "disease_proxy_rsid",
                    "regulatory_proxy_rsid",
                ]
            )
            .reset_index(
                drop=True
            )
        )

    return PairwiseCandidateBuildResult(
        dataframe=scaffold,
        disease_variants=int(
            disease[
                "disease_proxy_rsid"
            ].nunique()
            if not disease.empty
            else 0
        ),
        regulatory_variants=int(
            regulatory[
                "regulatory_proxy_rsid"
            ].nunique()
            if not regulatory.empty
            else 0
        ),
        candidate_pairs=int(
            len(
                scaffold
            )
        ),
    )


def classify_pairwise_ld(
    r2: Any,
) -> str:
    """Classify one pairwise r² value."""

    value = pd.to_numeric(
        pd.Series(
            [
                r2,
            ]
        ),
        errors="coerce",
    ).iloc[
        0
    ]

    if pd.isna(
        value
    ):
        return "UNAVAILABLE"

    value = float(
        value
    )

    if value >= HIGH_LD_R2:
        return "HIGH"

    if value >= MODERATE_LD_R2:
        return "MODERATE"

    return "LOW"


def synthesize_pairwise_ld(
    population_results: dict[
        str,
        pd.DataFrame,
    ],
) -> PairwiseLDSynthesisResult:
    """
    Synthesize pairwise disease-regulatory LD across populations.

    Expected per-population columns:
        lead_rsid
        disease_proxy_rsid
        regulatory_proxy_rsid
        r2
        query_status
    """

    frames: list[
        pd.DataFrame
    ] = []

    for population, dataframe in population_results.items():

        if dataframe.empty:
            continue

        current = dataframe.copy()

        current[
            "population"
        ] = (
            str(
                population
            )
            .strip()
            .upper()
        )

        current[
            "r2"
        ] = pd.to_numeric(
            current[
                "r2"
            ],
            errors="coerce",
        )

        current[
            "ld_class"
        ] = current[
            "r2"
        ].map(
            classify_pairwise_ld
        )

        frames.append(
            current
        )

    if not frames:

        return PairwiseLDSynthesisResult(
            dataframe=pd.DataFrame(),
            tested_pairs=0,
            high_ld_pairs=0,
            moderate_or_high_ld_pairs=0,
        )

    combined = pd.concat(
        frames,
        ignore_index=True,
    )

    records: list[
        dict[str, Any]
    ] = []

    for (
        lead,
        disease_proxy,
        regulatory_proxy,
    ), group in combined.groupby(
        [
            "lead_rsid",
            "disease_proxy_rsid",
            "regulatory_proxy_rsid",
        ],
        sort=True,
    ):

        valid_r2 = pd.to_numeric(
            group[
                "r2"
            ],
            errors="coerce",
        ).dropna()

        populations_tested = sorted(
            group[
                "population"
            ]
            .astype(str)
            .unique()
            .tolist()
        )

        high_populations = sorted(
            group.loc[
                group[
                    "ld_class"
                ].eq(
                    "HIGH"
                ),
                "population",
            ]
            .astype(str)
            .unique()
            .tolist()
        )

        moderate_populations = sorted(
            group.loc[
                group[
                    "ld_class"
                ].isin(
                    [
                        "MODERATE",
                        "HIGH",
                    ]
                ),
                "population",
            ]
            .astype(str)
            .unique()
            .tolist()
        )

        unavailable_populations = sorted(
            group.loc[
                group[
                    "ld_class"
                ].eq(
                    "UNAVAILABLE"
                ),
                "population",
            ]
            .astype(str)
            .unique()
            .tolist()
        )

        records.append(
            {
                "lead_rsid":
                    str(
                        lead
                    ),

                "disease_proxy_rsid":
                    str(
                        disease_proxy
                    ),

                "regulatory_proxy_rsid":
                    str(
                        regulatory_proxy
                    ),

                "populations_tested":
                    "|".join(
                        populations_tested
                    ),

                "population_count":
                    int(
                        len(
                            populations_tested
                        )
                    ),

                "minimum_r2":
                    (
                        float(
                            valid_r2.min()
                        )
                        if not valid_r2.empty
                        else None
                    ),

                "maximum_r2":
                    (
                        float(
                            valid_r2.max()
                        )
                        if not valid_r2.empty
                        else None
                    ),

                "moderate_or_high_populations":
                    "|".join(
                        moderate_populations
                    ),

                "moderate_or_high_population_count":
                    int(
                        len(
                            moderate_populations
                        )
                    ),

                "high_ld_populations":
                    "|".join(
                        high_populations
                    ),

                "high_ld_population_count":
                    int(
                        len(
                            high_populations
                        )
                    ),

                "unavailable_populations":
                    "|".join(
                        unavailable_populations
                    ),

                "unavailable_population_count":
                    int(
                        len(
                            unavailable_populations
                        )
                    ),

                "moderate_or_high_ld_in_any_population":
                    bool(
                        len(
                            moderate_populations
                        )
                        > 0
                    ),

                "high_ld_in_any_population":
                    bool(
                        len(
                            high_populations
                        )
                        > 0
                    ),

                "moderate_or_high_ld_in_multiple_populations":
                    bool(
                        len(
                            moderate_populations
                        )
                        >= 2
                    ),

                "high_ld_in_multiple_populations":
                    bool(
                        len(
                            high_populations
                        )
                        >= 2
                    ),
            }
        )

    synthesis = pd.DataFrame(
        records
    )

    moderate_count = int(
        synthesis[
            "moderate_or_high_ld_in_any_population"
        ]
        .fillna(
            False
        )
        .sum()
    )

    high_count = int(
        synthesis[
            "high_ld_in_any_population"
        ]
        .fillna(
            False
        )
        .sum()
    )

    return PairwiseLDSynthesisResult(
        dataframe=synthesis,
        tested_pairs=int(
            len(
                synthesis
            )
        ),
        high_ld_pairs=high_count,
        moderate_or_high_ld_pairs=moderate_count,
    )