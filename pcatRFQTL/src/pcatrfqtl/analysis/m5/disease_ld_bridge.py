"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/disease_ld_bridge.py

Description:
    M5.3C.3F disease-signal ↔ LD bridge integration.

    This module connects suggestive prostate cancer GWAS variants identified
    in M5.3C.3E with candidate tRF-QTL leads using previously normalized
    population-specific LD evidence.

    It also evaluates whether disease-associated LD proxies correspond to
    variants participating in the existing Moradi regulatory bridge.

    Integration is rsID-based only.

    Important build policy:
        - Harmonised prostate cancer GWAS loci are GRCh38.
        - Existing LD evidence and Moradi QTL evidence are GRCh37/hg19.
        - No coordinate join is performed across assemblies.
        - No silent liftover is performed.
        - Cross-build integration is allowed only through canonical rsIDs.

    LD interpretation:
        HIGH:
            r² >= 0.80

        MODERATE:
            0.50 <= r² < 0.80

        BELOW_THRESHOLD:
            r² < 0.50

    Scientific safeguards:
        - Physical proximity is not interpreted as LD.
        - Suggestive disease association is not interpreted as causality.
        - LD bridge does not imply colocalization.
        - Regulatory bridge does not imply shared causal variant.
        - Cross-population recurrence is not independent replication.
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


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class DiseaseLDBridgeResult:
    """Population-specific disease ↔ LD bridge result."""

    population: str

    bridge: pd.DataFrame

    input_disease_variants: int

    disease_variants_with_rsid: int

    disease_variants_without_rsid: int

    matched_ld_rows: int

    matched_disease_variant_identities: int


@dataclass(frozen=True)
class DiseaseRegulatoryBridgeResult:
    """Disease ↔ LD ↔ regulatory bridge result."""

    population: str

    bridge: pd.DataFrame

    disease_ld_rows: int

    regulatory_rows: int

    integrated_rows: int

    integrated_disease_variants: int


# ============================================================================
# Normalization helpers
# ============================================================================


def _normalize_rsid(
    series: pd.Series,
) -> pd.Series:
    """Normalize canonical rsID values."""

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


def _numeric(
    series: pd.Series,
) -> pd.Series:
    """Convert source series to numeric."""

    return pd.to_numeric(
        series,
        errors="coerce",
    )


def _resolve_column(
    dataframe: pd.DataFrame,
    aliases: tuple[str, ...],
    *,
    required: bool = True,
) -> str | None:
    """Resolve first available source-column alias."""

    lower_map = {
        str(column).strip().lower():
            str(column)
        for column in dataframe.columns
    }

    for alias in aliases:

        match = lower_map.get(
            alias.lower()
        )

        if match is not None:
            return match

    if required:

        raise ValueError(
            "Could not resolve required column. "
            f"Expected one of: {aliases}; "
            f"available={list(dataframe.columns)}"
        )

    return None


# ============================================================================
# LD normalization
# ============================================================================


def normalize_ld_evidence(
    dataframe: pd.DataFrame,
    *,
    population: str,
) -> pd.DataFrame:
    """
    Normalize pre-existing LD evidence into the schema required by M5.3C.3F.

    Expected core source fields from the existing M5 LD preparation include:
        lead_rsid
        neighbor_rsid
        r2
        rsid_valid
        r2_valid
        is_self_pair
    """

    lead_column = _resolve_column(
        dataframe,
        (
            "lead_rsid",
            "lead",
            "query_rsid",
        ),
    )

    neighbor_column = _resolve_column(
        dataframe,
        (
            "neighbor_rsid",
            "proxy_rsid",
            "ld_proxy_rsid",
            "rsid",
        ),
    )

    r2_column = _resolve_column(
        dataframe,
        (
            "r2",
            "r_squared",
            "rsquared",
        ),
    )

    rsid_valid_column = _resolve_column(
        dataframe,
        (
            "rsid_valid",
        ),
        required=False,
    )

    r2_valid_column = _resolve_column(
        dataframe,
        (
            "r2_valid",
        ),
        required=False,
    )

    self_pair_column = _resolve_column(
        dataframe,
        (
            "is_self_pair",
        ),
        required=False,
    )

    out = pd.DataFrame(
        index=dataframe.index,
    )

    out["population"] = (
        str(population)
        .strip()
        .upper()
    )

    out["lead_rsid"] = _normalize_rsid(
        dataframe[
            lead_column
        ]
    )

    out["ld_variant_rsid"] = _normalize_rsid(
        dataframe[
            neighbor_column
        ]
    )

    out["r2"] = _numeric(
        dataframe[
            r2_column
        ]
    )

    if rsid_valid_column is not None:

        rsid_valid = (
            dataframe[
                rsid_valid_column
            ]
            .fillna(False)
            .astype(bool)
        )

    else:

        rsid_valid = (
            out["lead_rsid"].notna()
            &
            out["ld_variant_rsid"].notna()
        )

    if r2_valid_column is not None:

        r2_valid = (
            dataframe[
                r2_valid_column
            ]
            .fillna(False)
            .astype(bool)
        )

    else:

        r2_valid = (
            out["r2"].notna()
            &
            out["r2"].between(
                0,
                1,
                inclusive="both",
            )
        )

    if self_pair_column is not None:

        is_self_pair = (
            dataframe[
                self_pair_column
            ]
            .fillna(False)
            .astype(bool)
        )

    else:

        is_self_pair = (
            out["lead_rsid"]
            .eq(
                out["ld_variant_rsid"]
            )
            .fillna(False)
        )

    out["rsid_valid"] = rsid_valid

    out["r2_valid"] = r2_valid

    out["is_self_pair"] = is_self_pair

    out["ld_class"] = "BELOW_THRESHOLD"

    out.loc[
        out["r2"].ge(
            MODERATE_LD_R2
        ),
        "ld_class",
    ] = "MODERATE"

    out.loc[
        out["r2"].ge(
            HIGH_LD_R2
        ),
        "ld_class",
    ] = "HIGH"

    out["passes_moderate_ld"] = (
        out["r2"].ge(
            MODERATE_LD_R2
        )
        .fillna(False)
    )

    out["passes_high_ld"] = (
        out["r2"].ge(
            HIGH_LD_R2
        )
        .fillna(False)
    )

    return out


# ============================================================================
# Disease variant preparation
# ============================================================================


def prepare_disease_variants(
    recurrence: pd.DataFrame,
) -> pd.DataFrame:
    """
    Prepare M5.3C.3E suggestive disease variants for rsID-based LD matching.

    No association rows are removed from the source artifact. This function
    creates an analysis-specific bridge table.
    """

    required = {
        "lead_rsid",
        "variant_identity_key",
        "rsid",
        "variant_id",
        "chromosome",
        "base_pair_location",
        "distance_to_lead",
        "study_count",
        "minimum_p_value",
        "genome_wide_in_any_study",
        "effect_direction_consistent",
    }

    missing = (
        required
        - set(
            recurrence.columns
        )
    )

    if missing:

        raise ValueError(
            "M5.3C.3E recurrence table missing columns: "
            f"{sorted(missing)}"
        )

    out = recurrence.copy()

    out["lead_rsid"] = _normalize_rsid(
        out["lead_rsid"]
    )

    out["disease_variant_rsid"] = _normalize_rsid(
        out["rsid"]
    )

    out["minimum_p_value"] = _numeric(
        out["minimum_p_value"]
    )

    out["study_count"] = (
        _numeric(
            out["study_count"]
        )
        .astype("Int64")
    )

    out["disease_rsid_usable"] = (
        out["disease_variant_rsid"].notna()
    )

    out["disease_variant_is_candidate_lead"] = (
        out["disease_variant_rsid"]
        .eq(
            out["lead_rsid"]
        )
        .fillna(False)
    )

    return out


# ============================================================================
# Disease ↔ LD bridge
# ============================================================================


def build_disease_ld_bridge(
    *,
    recurrence: pd.DataFrame,
    ld_evidence: pd.DataFrame,
    population: str,
) -> DiseaseLDBridgeResult:
    """
    Match suggestive disease variants to candidate-lead LD proxies.

    Matching key:
        candidate lead rsID
        +
        disease variant rsID == LD neighbor rsID

    The candidate lead itself is handled separately from proxy matching.
    """

    disease = prepare_disease_variants(
        recurrence
    )

    ld = normalize_ld_evidence(
        ld_evidence,
        population=population,
    )

    usable_disease = disease.loc[
        disease[
            "disease_rsid_usable"
        ]
    ].copy()

    usable_ld = ld.loc[
        ld["rsid_valid"]
        &
        ld["r2_valid"]
        &
        ~ld["is_self_pair"]
    ].copy()

    matched = usable_disease.merge(
        usable_ld,
        left_on=[
            "lead_rsid",
            "disease_variant_rsid",
        ],
        right_on=[
            "lead_rsid",
            "ld_variant_rsid",
        ],
        how="inner",
        validate="many_to_many",
    )

    if not matched.empty:

        matched[
            "bridge_type"
        ] = "LD_MEDIATED_DISEASE_BRIDGE"

        matched[
            "bridge_eligible_primary"
        ] = matched[
            "passes_high_ld"
        ].fillna(
            False
        )

        matched[
            "bridge_eligible_sensitivity"
        ] = matched[
            "passes_moderate_ld"
        ].fillna(
            False
        )

    return DiseaseLDBridgeResult(
        population=(
            str(population)
            .strip()
            .upper()
        ),
        bridge=matched,
        input_disease_variants=int(
            len(
                disease
            )
        ),
        disease_variants_with_rsid=int(
            disease[
                "disease_rsid_usable"
            ].sum()
        ),
        disease_variants_without_rsid=int(
            (
                ~disease[
                    "disease_rsid_usable"
                ]
            ).sum()
        ),
        matched_ld_rows=int(
            len(
                matched
            )
        ),
        matched_disease_variant_identities=int(
            matched[
                "variant_identity_key"
            ].nunique()
            if not matched.empty
            else 0
        ),
    )


# ============================================================================
# Regulatory bridge normalization
# ============================================================================


def normalize_regulatory_bridge(
    dataframe: pd.DataFrame,
    *,
    population: str,
) -> pd.DataFrame:
    """
    Normalize the previously generated Moradi LD bridge.

    The function deliberately accepts a small alias set because earlier M5
    bridge artifacts may use proxy_rsid, neighbor_rsid, or ld_proxy_rsid.
    """

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
        ),
    )

    r2_column = _resolve_column(
        dataframe,
        (
            "r2",
            "ld_r2",
        ),
        required=False,
    )

    out = dataframe.copy()

    out["population"] = (
        str(population)
        .strip()
        .upper()
    )

    out["regulatory_lead_rsid"] = _normalize_rsid(
        dataframe[
            lead_column
        ]
    )

    out["regulatory_proxy_rsid"] = _normalize_rsid(
        dataframe[
            proxy_column
        ]
    )

    if r2_column is not None:

        out["regulatory_r2"] = _numeric(
            dataframe[
                r2_column
            ]
        )

    else:

        out["regulatory_r2"] = pd.NA

    return out


# ============================================================================
# Disease ↔ regulatory bridge
# ============================================================================


def build_disease_regulatory_bridge(
    *,
    disease_ld_bridge: pd.DataFrame,
    regulatory_bridge: pd.DataFrame,
    population: str,
) -> DiseaseRegulatoryBridgeResult:
    """
    Identify disease-associated LD proxies that are also regulatory proxies.

    This is a same-proxy rsID bridge:

        candidate lead
             ↓ LD
        disease-associated variant
             =
        regulatory bridge proxy

    It does NOT prove colocalization.
    """

    regulatory = normalize_regulatory_bridge(
        regulatory_bridge,
        population=population,
    )

    if disease_ld_bridge.empty:

        integrated = pd.DataFrame()

    else:

        integrated = (
            disease_ld_bridge.merge(
                regulatory,
                left_on=[
                    "lead_rsid",
                    "disease_variant_rsid",
                ],
                right_on=[
                    "regulatory_lead_rsid",
                    "regulatory_proxy_rsid",
                ],
                how="inner",
                suffixes=(
                    "_disease",
                    "_regulatory",
                ),
            )
        )

        if not integrated.empty:

            integrated[
                "integration_class"
            ] = (
                "DISEASE_REGULATORY_SHARED_LD_PROXY"
            )

            integrated[
                "colocalization_supported"
            ] = False

            integrated[
                "colocalization_not_assessed"
            ] = True

    return DiseaseRegulatoryBridgeResult(
        population=(
            str(population)
            .strip()
            .upper()
        ),
        bridge=integrated,
        disease_ld_rows=int(
            len(
                disease_ld_bridge
            )
        ),
        regulatory_rows=int(
            len(
                regulatory
            )
        ),
        integrated_rows=int(
            len(
                integrated
            )
        ),
        integrated_disease_variants=int(
            integrated[
                "disease_variant_rsid"
            ].nunique()
            if not integrated.empty
            else 0
        ),
    )


# ============================================================================
# Cross-population synthesis
# ============================================================================


def synthesize_disease_ld_bridges(
    bridges: dict[
        str,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    """
    Synthesize disease-LD bridges across reference populations.

    Population recurrence is interpreted as sensitivity consistency, not
    biological cohort replication.
    """

    frames = []

    for population, dataframe in bridges.items():

        if dataframe.empty:
            continue

        current = dataframe.copy()

        current[
            "population"
        ] = (
            str(population)
            .strip()
            .upper()
        )

        frames.append(
            current
        )

    if not frames:

        return pd.DataFrame()

    combined = pd.concat(
        frames,
        ignore_index=True,
    )

    records: list[
        dict[str, Any]
    ] = []

    for (
        lead_rsid,
        disease_variant_rsid,
    ), group in combined.groupby(
        [
            "lead_rsid",
            "disease_variant_rsid",
        ],
        sort=True,
    ):

        populations = sorted(
            group[
                "population"
            ]
            .astype(str)
            .unique()
            .tolist()
        )

        r2_values = pd.to_numeric(
            group["r2"],
            errors="coerce",
        ).dropna()

        first = group.sort_values(
            "minimum_p_value",
            ascending=True,
        ).iloc[0]

        high_population_count = int(
            group[
                "passes_high_ld"
            ]
            .fillna(False)
            .groupby(
                group["population"]
            )
            .any()
            .sum()
        )

        moderate_population_count = int(
            group[
                "passes_moderate_ld"
            ]
            .fillna(False)
            .groupby(
                group["population"]
            )
            .any()
            .sum()
        )

        records.append(
            {
                "lead_rsid":
                    str(
                        lead_rsid
                    ),

                "disease_variant_rsid":
                    str(
                        disease_variant_rsid
                    ),

                "variant_identity_key":
                    (
                        None
                        if pd.isna(
                            first[
                                "variant_identity_key"
                            ]
                        )
                        else str(
                            first[
                                "variant_identity_key"
                            ]
                        )
                    ),

                "variant_id":
                    (
                        None
                        if pd.isna(
                            first[
                                "variant_id"
                            ]
                        )
                        else str(
                            first[
                                "variant_id"
                            ]
                        )
                    ),

                "minimum_disease_p_value":
                    (
                        None
                        if pd.isna(
                            first[
                                "minimum_p_value"
                            ]
                        )
                        else float(
                            first[
                                "minimum_p_value"
                            ]
                        )
                    ),

                "disease_study_count":
                    (
                        None
                        if pd.isna(
                            first[
                                "study_count"
                            ]
                        )
                        else int(
                            first[
                                "study_count"
                            ]
                        )
                    ),

                "genome_wide_in_any_disease_study":
                    bool(
                        group[
                            "genome_wide_in_any_study"
                        ]
                        .fillna(False)
                        .any()
                    ),

                "ld_populations":
                    "|".join(
                        populations
                    ),

                "ld_population_count":
                    int(
                        len(
                            populations
                        )
                    ),

                "high_ld_population_count":
                    high_population_count,

                "moderate_or_high_ld_population_count":
                    moderate_population_count,

                "minimum_r2":
                    (
                        float(
                            r2_values.min()
                        )
                        if not r2_values.empty
                        else None
                    ),

                "maximum_r2":
                    (
                        float(
                            r2_values.max()
                        )
                        if not r2_values.empty
                        else None
                    ),

                "high_ld_in_any_population":
                    bool(
                        high_population_count
                        > 0
                    ),

                "high_ld_in_multiple_populations":
                    bool(
                        high_population_count
                        >= 2
                    ),

                "ld_sensitivity_consistent":
                    bool(
                        moderate_population_count
                        >= 2
                    ),
            }
        )

    result = pd.DataFrame(
        records
    )

    if result.empty:
        return result

    return result.sort_values(
        [
            "lead_rsid",
            "high_ld_population_count",
            "minimum_disease_p_value",
        ],
        ascending=[
            True,
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )