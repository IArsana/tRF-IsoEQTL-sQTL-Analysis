"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/coloc_readiness.py

Description:
    Core logic for M5.5 colocalization readiness assessment.

    This module deliberately separates:

        1. dataset-level GWAS readiness;
        2. dataset-level Moradi QTL readiness;
        3. candidate-level readiness.

    Important data-model rule:
        Raw M4.1 Moradi QTL index partitions do not contain candidate
        lead_rsid values. Therefore candidate association must not be inferred
        directly from the raw QTL index.

    Scientific safeguards:
        - Significant-only QTL resources are not treated as dense regional
          summary statistics.
        - Missing QTL variants are not interpreted as null associations.
        - Cross-build coordinate joins are not performed.
        - Canonical rsID matching is the only allowed cross-build identity
          bridge in this stage.
        - No formal colocalization is performed.
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

import numpy as np
import pandas as pd


# ============================================================================
# Result containers
# ============================================================================


@dataclass(frozen=True)
class DatasetReadiness:
    """Container for dataset-level readiness output."""

    dataframe: pd.DataFrame


@dataclass(frozen=True)
class PairReadiness:
    """Container for candidate-level readiness output."""

    dataframe: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_rsid(
    series: pd.Series,
) -> pd.Series:
    """Normalize canonical rsID values."""

    values = (
        series
        .astype("string")
        .str.strip()
        .str.lower()
    )

    return values.where(
        values.str.fullmatch(
            r"rs\d+",
            na=False,
        ),
        pd.NA,
    )


def _numeric(
    series: pd.Series,
) -> pd.Series:
    """Convert a Series to numeric with invalid values set missing."""

    return pd.to_numeric(
        series,
        errors="coerce",
    )


def _resolve_column(
    dataframe: pd.DataFrame,
    aliases: tuple[str, ...],
) -> str | None:
    """Resolve one logical field from a set of candidate column aliases."""

    lookup = {
        str(column).lower():
            str(column)
        for column
        in dataframe.columns
    }

    for alias in aliases:

        resolved = lookup.get(
            alias.lower()
        )

        if resolved is not None:

            return resolved

    return None


def _fraction_present(
    series: pd.Series,
) -> float:
    """Return non-missing fraction."""

    if len(series) == 0:

        return 0.0

    return float(
        series.notna().mean()
    )


def _derive_maf(
    allele_frequency: pd.Series,
) -> pd.Series:
    """Derive minor-allele frequency from allele frequency."""

    values = _numeric(
        allele_frequency
    )

    result = pd.Series(
        np.nan,
        index=values.index,
        dtype="float64",
    )

    valid = values.between(
        0.0,
        1.0,
        inclusive="both",
    )

    result.loc[
        valid
    ] = np.minimum(
        values.loc[
            valid
        ],
        1.0
        - values.loc[
            valid
        ],
    )

    return result


# ============================================================================
# GWAS readiness
# ============================================================================


def assess_gwas_readiness(
    dataframe: pd.DataFrame,
) -> DatasetReadiness:
    """
    Assess disease-side GWAS readiness by study × lead locus.
    """

    required = {
        "study_accession",
        "lead_rsid",
        "rsid",
        "beta",
        "standard_error",
        "effect_allele_frequency",
        "p_value",
    }

    missing = (
        required
        - set(
            dataframe.columns
        )
    )

    if missing:

        raise ValueError(
            "GWAS standardized locus dataset missing required columns: "
            f"{sorted(missing)}"
        )

    source = dataframe.copy()

    source["_rsid"] = _normalize_rsid(
        source["rsid"]
    )

    source["_beta"] = _numeric(
        source["beta"]
    )

    source["_se"] = _numeric(
        source["standard_error"]
    )

    source["_varbeta"] = (
        source["_se"]
        ** 2
    )

    source["_p"] = _numeric(
        source["p_value"]
    )

    source["_maf"] = _derive_maf(
        source["effect_allele_frequency"]
    )

    records: list[
        dict[str, Any]
    ] = []

    for (
        study,
        lead,
    ), group in source.groupby(
        [
            "study_accession",
            "lead_rsid",
        ],
        sort=True,
        dropna=False,
    ):

        rows = int(
            len(group)
        )

        unique_rsids = int(
            group["_rsid"]
            .dropna()
            .nunique()
        )

        rsid_fraction = _fraction_present(
            group["_rsid"]
        )

        beta_fraction = _fraction_present(
            group["_beta"]
        )

        varbeta_fraction = _fraction_present(
            group["_varbeta"]
        )

        p_fraction = _fraction_present(
            group["_p"]
        )

        maf_fraction = _fraction_present(
            group["_maf"]
        )

        vector_ready = bool(
            rsid_fraction == 1.0
            and
            beta_fraction == 1.0
            and
            varbeta_fraction == 1.0
        )

        dense_locus = bool(
            rows >= 100
        )

        if vector_ready and dense_locus:

            status = (
                "GWAS_VECTOR_READY_METADATA_REVIEW"
            )

        elif not dense_locus:

            status = (
                "GWAS_LOCUS_COVERAGE_INSUFFICIENT"
            )

        else:

            status = (
                "GWAS_VECTOR_INCOMPLETE"
            )

        records.append(
            {
                "study_accession":
                    str(study),

                "lead_rsid":
                    str(lead),

                "variant_rows":
                    rows,

                "unique_rsids":
                    unique_rsids,

                "rsid_complete_fraction":
                    rsid_fraction,

                "beta_complete_fraction":
                    beta_fraction,

                "varbeta_complete_fraction":
                    varbeta_fraction,

                "p_value_complete_fraction":
                    p_fraction,

                "maf_complete_fraction":
                    maf_fraction,

                "full_locus_summary_available":
                    True,

                "dense_locus_available":
                    dense_locus,

                "basic_coloc_vector_ready":
                    vector_ready,

                "trait_type":
                    "cc",

                "susie_ld_matrix_available":
                    False,

                "formal_coloc_data_status":
                    status,
            }
        )

    return DatasetReadiness(
        dataframe=pd.DataFrame(
            records
        )
    )


# ============================================================================
# Moradi QTL readiness
# ============================================================================


def assess_qtl_readiness(
    dataframe: pd.DataFrame,
    *,
    source_scope: str,
    full_summary_statistics_available: bool,
) -> DatasetReadiness:
    """
    Assess raw M4.1 Moradi QTL dataset readiness.

    This function intentionally does NOT require lead_rsid.

    Assessment unit:
        source_partition

    Optional feature-level grouping is used only when a stable feature column
    is available.
    """

    if "source_partition" not in dataframe.columns:

        raise ValueError(
            "Moradi QTL dataframe requires source_partition provenance."
        )

    variant_col = _resolve_column(
        dataframe,
        (
            "qtl_rsid",
            "variant_rsid",
            "rsid",
            "snp",
            "snp_id",
            "variant_id",
        ),
    )

    beta_col = _resolve_column(
        dataframe,
        (
            "beta",
            "effect_size",
            "effect",
            "qtl_beta",
        ),
    )

    se_col = _resolve_column(
        dataframe,
        (
            "standard_error",
            "se",
            "stderr",
            "qtl_se",
        ),
    )

    p_col = _resolve_column(
        dataframe,
        (
            "p_value",
            "pvalue",
            "pval",
            "p",
        ),
    )

    af_col = _resolve_column(
        dataframe,
        (
            "effect_allele_frequency",
            "eaf",
            "maf",
            "allele_frequency",
        ),
    )

    feature_col = _resolve_column(
        dataframe,
        (
            "feature_id",
            "normalized_feature_id",
            "gene_id",
            "exon_id",
            "intron_id",
            "isoform_id",
        ),
    )

    source = dataframe.copy()

    if variant_col is not None:

        source["_rsid"] = _normalize_rsid(
            source[
                variant_col
            ]
        )

    else:

        source["_rsid"] = pd.NA

    if beta_col is not None:

        source["_beta"] = _numeric(
            source[
                beta_col
            ]
        )

    else:

        source["_beta"] = np.nan

    if se_col is not None:

        source["_se"] = _numeric(
            source[
                se_col
            ]
        )

        source["_varbeta"] = (
            source["_se"]
            ** 2
        )

    else:

        source["_se"] = np.nan
        source["_varbeta"] = np.nan

    if p_col is not None:

        source["_p"] = _numeric(
            source[
                p_col
            ]
        )

    else:

        source["_p"] = np.nan

    if af_col is not None:

        source["_maf"] = _derive_maf(
            source[
                af_col
            ]
        )

    else:

        source["_maf"] = np.nan

    grouping = [
        "source_partition",
    ]

    if feature_col is not None:

        grouping.append(
            feature_col
        )

    records: list[
        dict[str, Any]
    ] = []

    for keys, group in source.groupby(
        grouping,
        sort=True,
        dropna=False,
    ):

        if not isinstance(
            keys,
            tuple,
        ):

            keys = (
                keys,
            )

        partition = keys[
            0
        ]

        feature = (
            keys[
                1
            ]
            if len(keys) > 1
            else None
        )

        rows = int(
            len(group)
        )

        unique_rsids = int(
            group["_rsid"]
            .dropna()
            .nunique()
        )

        rsid_fraction = _fraction_present(
            group["_rsid"]
        )

        beta_fraction = _fraction_present(
            group["_beta"]
        )

        varbeta_fraction = _fraction_present(
            group["_varbeta"]
        )

        p_fraction = _fraction_present(
            group["_p"]
        )

        maf_fraction = _fraction_present(
            group["_maf"]
        )

        abf_vector_ready = bool(
            rsid_fraction == 1.0
            and
            (
                (
                    beta_fraction == 1.0
                    and
                    varbeta_fraction == 1.0
                )
                or
                (
                    p_fraction == 1.0
                    and
                    maf_fraction == 1.0
                )
            )
        )

        dense_locus = bool(
            full_summary_statistics_available
            and
            rows >= 100
        )

        if not full_summary_statistics_available:

            status = (
                "QTL_SIGNIFICANT_ONLY_NOT_COLOC_READY"
            )

        elif not dense_locus:

            status = (
                "QTL_LOCUS_COVERAGE_INSUFFICIENT"
            )

        elif not abf_vector_ready:

            status = (
                "QTL_VECTOR_INCOMPLETE"
            )

        else:

            status = (
                "QTL_VECTOR_READY"
            )

        records.append(
            {
                "source_partition":
                    str(partition),

                "feature_id":
                    (
                        None
                        if feature is None
                        or pd.isna(
                            feature
                        )
                        else str(
                            feature
                        )
                    ),

                "source_scope":
                    str(
                        source_scope
                    ),

                "variant_rows":
                    rows,

                "unique_rsids":
                    unique_rsids,

                "rsid_complete_fraction":
                    rsid_fraction,

                "beta_complete_fraction":
                    beta_fraction,

                "varbeta_complete_fraction":
                    varbeta_fraction,

                "p_value_complete_fraction":
                    p_fraction,

                "maf_complete_fraction":
                    maf_fraction,

                "full_locus_summary_available":
                    bool(
                        full_summary_statistics_available
                    ),

                "dense_locus_available":
                    dense_locus,

                "basic_coloc_vector_ready":
                    abf_vector_ready,

                "trait_type":
                    "quant",

                "susie_ld_matrix_available":
                    False,

                "formal_coloc_data_status":
                    status,
            }
        )

    return DatasetReadiness(
        dataframe=pd.DataFrame(
            records
        )
    )


# ============================================================================
# Candidate-level readiness
# ============================================================================


def assess_pair_readiness(
    *,
    gwas_locus: pd.DataFrame,
    qtl_rows: pd.DataFrame,
    minimum_shared_variants: int = 10,
) -> PairReadiness:
    """
    Assess descriptive GWAS ↔ raw-QTL variant overlap.

    Important:
        Raw Moradi M4.1 does not contain candidate lead_rsid mappings.

        Therefore this function does NOT claim candidate-specific QTL
        colocalization readiness from raw QTL data alone.

        It only reports descriptive rsID overlap between each GWAS
        study-lead locus and the complete raw Moradi QTL index.

    Formal coloc remains blocked when the QTL resource is significant-only.
    """

    required_gwas = {
        "study_accession",
        "lead_rsid",
        "rsid",
    }

    missing = (
        required_gwas
        - set(
            gwas_locus.columns
        )
    )

    if missing:

        raise ValueError(
            "GWAS input missing required columns: "
            f"{sorted(missing)}"
        )

    qtl_variant_col = _resolve_column(
        qtl_rows,
        (
            "qtl_rsid",
            "variant_rsid",
            "rsid",
            "snp",
            "snp_id",
            "variant_id",
        ),
    )

    if qtl_variant_col is None:

        raise ValueError(
            "Unable to identify a canonical variant-rsID field "
            "in raw Moradi QTL input."
        )

    gwas = gwas_locus.copy()

    qtl = qtl_rows.copy()

    gwas["_rsid"] = _normalize_rsid(
        gwas["rsid"]
    )

    qtl["_rsid"] = _normalize_rsid(
        qtl[
            qtl_variant_col
        ]
    )

    qtl_rsids = set(
        qtl["_rsid"]
        .dropna()
        .astype(str)
    )

    records: list[
        dict[str, Any]
    ] = []

    for (
        study,
        lead,
    ), group in gwas.groupby(
        [
            "study_accession",
            "lead_rsid",
        ],
        sort=True,
        dropna=False,
    ):

        gwas_rsids = set(
            group["_rsid"]
            .dropna()
            .astype(str)
        )

        shared = (
            gwas_rsids
            &
            qtl_rsids
        )

        shared_count = int(
            len(shared)
        )

        if shared_count == 0:

            readiness_status = (
                "NO_SHARED_VARIANTS_WITH_RAW_QTL_INDEX"
            )

        elif shared_count < minimum_shared_variants:

            readiness_status = (
                "LIMITED_RAW_QTL_VARIANT_OVERLAP"
            )

        else:

            readiness_status = (
                "RAW_QTL_OVERLAP_PRESENT_FORMAL_COLOC_BLOCKED"
            )

        records.append(
            {
                "study_accession":
                    str(study),

                "lead_rsid":
                    str(lead),

                "gwas_unique_rsids":
                    int(
                        len(
                            gwas_rsids
                        )
                    ),

                "raw_qtl_unique_rsids":
                    int(
                        len(
                            qtl_rsids
                        )
                    ),

                "shared_variant_count":
                    shared_count,

                "minimum_shared_variants":
                    int(
                        minimum_shared_variants
                    ),

                "cross_build_match_method":
                    "canonical_rsid_only",

                "candidate_specific_qtl_mapping_available":
                    False,

                "formal_coloc_ready":
                    False,

                "coloc_abf_ready":
                    False,

                "coloc_susie_ready":
                    False,

                "readiness_status":
                    readiness_status,

                "blocking_reason":
                    (
                        "Raw Moradi M4.1 QTL index does not provide "
                        "candidate-specific lead mapping and is configured "
                        "as a significant-only QTL resource rather than "
                        "dense full-locus summary statistics."
                    ),
            }
        )

    return PairReadiness(
        dataframe=pd.DataFrame(
            records
        )
    )