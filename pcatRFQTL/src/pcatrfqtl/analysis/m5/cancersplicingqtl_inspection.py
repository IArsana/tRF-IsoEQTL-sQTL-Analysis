"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/cancersplicingqtl_inspection.py

Description:
    Core inspection logic for M5.6C.1 CancerSplicingQTL PRAD.

    This module evaluates the actual downloaded CancerSplicingQTL PRAD
    workbook and determines whether it can resolve the regulatory-QTL
    data gap identified during M5.5 and M5.6.

    The actual public PRAD table contains:

        Cancer type
        SNP ID
        SNP position
        Alleles
        Gene
        Splicing type
        Splicing exon
        AS ID
        Splice position
        Beta
        T-stat
        R
        P-value

    Critical source-scope rule:
        CancerSplicingQTL reports identified/significant sQTL associations.
        Therefore high row counts or hundreds of variants per feature do NOT
        establish availability of complete unfiltered locus-wide summary
        statistics.

    T-statistic policy:
        A diagnostic standard error may mathematically be calculated as:

            abs(beta / t)

        for finite nonzero t statistics.

        This diagnostic value is NOT promoted to formal colocalization input
        in M5.6C.1.

    Scientific safeguards:
        - Significant-only associations are not treated as a complete test set.
        - Absent variants are not interpreted as null associations.
        - Allele orientation is not inferred from the raw "Alleles" field.
        - Cross-build coordinate joins are prohibited.
        - Canonical rsID is the only allowed cross-build identity bridge.
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
# Classification constants
# ============================================================================


DATASET_NOT_PRESENT = (
    "DATASET_NOT_PRESENT"
)

EMPTY_DATASET = (
    "EMPTY_DATASET"
)

SCHEMA_MISMATCH = (
    "SCHEMA_MISMATCH"
)

NON_PRAD_CONTENT_DETECTED = (
    "NON_PRAD_CONTENT_DETECTED"
)

VARIANT_IDENTITY_INCOMPLETE = (
    "VARIANT_IDENTITY_INCOMPLETE"
)

FEATURE_IDENTITY_INCOMPLETE = (
    "FEATURE_IDENTITY_INCOMPLETE"
)

SIGNIFICANT_ONLY_NOT_COLOC_READY = (
    "SIGNIFICANT_ONLY_NOT_COLOC_READY"
)

FORMAL_ABF_READY_SUSIE_BLOCKED = (
    "FORMAL_ABF_READY_SUSIE_BLOCKED"
)

FORMAL_COLOC_READY = (
    "FORMAL_COLOC_READY"
)


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class CancerSplicingQTLInspectionResult:
    """Container for M5.6C.1 outputs."""

    schema: pd.DataFrame

    feature_coverage: pd.DataFrame

    readiness: pd.DataFrame

    candidates: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_rsid_series(
    series: pd.Series,
) -> pd.Series:
    """Normalize canonical dbSNP rs identifiers."""

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
    """Convert values to numeric while preserving invalid values as missing."""

    return pd.to_numeric(
        series,
        errors="coerce",
    )


def _fraction_present(
    series: pd.Series,
) -> float:
    """Return fraction of nonmissing values."""

    if len(series) == 0:

        return 0.0

    return float(
        series.notna().mean()
    )


def _safe_int(
    value: Any,
) -> int:
    """Safely cast count-like values to int."""

    if pd.isna(
        value
    ):

        return 0

    return int(
        value
    )


# ============================================================================
# Schema resolution
# ============================================================================


def resolve_expected_schema(
    dataframe: pd.DataFrame,
    *,
    config: dict[str, Any],
) -> dict[str, str | None]:
    """
    Resolve logical fields using the known actual CancerSplicingQTL schema.

    Unlike earlier exploratory versions, M5.6C.1 no longer performs broad
    alias guessing for the official PRAD workbook.
    """

    expected_schema = config[
        "expected_schema"
    ]

    actual_columns = {
        str(
            column
        )
        for column
        in dataframe.columns
    }

    resolved: dict[
        str,
        str | None,
    ] = {}

    for logical_field, properties in expected_schema.items():

        actual_column = str(
            properties[
                "actual_column"
            ]
        )

        if actual_column in actual_columns:

            resolved[
                str(
                    logical_field
                )
            ] = actual_column

        else:

            resolved[
                str(
                    logical_field
                )
            ] = None

    return resolved


def build_schema_table(
    dataframe: pd.DataFrame,
    *,
    resolved_schema: dict[str, str | None],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build schema-level audit table."""

    expected_schema = config[
        "expected_schema"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    for logical_field, properties in expected_schema.items():

        actual_column = resolved_schema.get(
            logical_field
        )

        required = bool(
            properties.get(
                "required",
                False,
            )
        )

        if actual_column is None:

            rows.append(
                {
                    "logical_field":
                        str(
                            logical_field
                        ),

                    "expected_column":
                        str(
                            properties[
                                "actual_column"
                            ]
                        ),

                    "actual_column":
                        None,

                    "required":
                        required,

                    "available":
                        False,

                    "nonmissing_fraction":
                        0.0,

                    "dtype":
                        None,

                    "unique_values":
                        0,
                }
            )

            continue

        series = dataframe[
            actual_column
        ]

        rows.append(
            {
                "logical_field":
                    str(
                        logical_field
                    ),

                "expected_column":
                    str(
                        properties[
                            "actual_column"
                        ]
                    ),

                "actual_column":
                    str(
                        actual_column
                    ),

                "required":
                    required,

                "available":
                    True,

                "nonmissing_fraction":
                    _fraction_present(
                        series
                    ),

                "dtype":
                    str(
                        series.dtype
                    ),

                "unique_values":
                    int(
                        series.nunique(
                            dropna=True
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Standard-error diagnostic
# ============================================================================


def derive_standard_error_diagnostic(
    dataframe: pd.DataFrame,
    *,
    beta_column: str,
    t_column: str,
) -> pd.DataFrame:
    """
    Derive diagnostic SE = abs(beta / t).

    This is diagnostic only.

    The derived value must not be treated as formal coloc input during
    M5.6C.1.
    """

    beta = _numeric(
        dataframe[
            beta_column
        ]
    )

    t_stat = _numeric(
        dataframe[
            t_column
        ]
    )

    valid = (
        beta.notna()
        &
        t_stat.notna()
        &
        np.isfinite(
            beta
        )
        &
        np.isfinite(
            t_stat
        )
        &
        t_stat.ne(
            0
        )
    )

    derived = pd.Series(
        np.nan,
        index=dataframe.index,
        dtype="float64",
    )

    derived.loc[
        valid
    ] = (
        beta.loc[
            valid
        ]
        /
        t_stat.loc[
            valid
        ]
    ).abs()

    return pd.DataFrame(
        {
            "beta_numeric":
                beta,

            "t_stat_numeric":
                t_stat,

            "derived_se_diagnostic":
                derived,

            "derived_se_valid":
                valid,
        },
        index=dataframe.index,
    )


# ============================================================================
# Feature coverage
# ============================================================================


def assess_feature_coverage(
    dataframe: pd.DataFrame,
    *,
    resolved_schema: dict[str, str | None],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Describe per-AS-feature coverage.

    IMPORTANT:
        Coverage is descriptive only.

        It does not establish dense/unfiltered summary-statistic availability
        because the underlying source is significance-filtered.
    """

    feature_column = resolved_schema.get(
        "feature_id"
    )

    variant_column = resolved_schema.get(
        "variant_rsid"
    )

    p_column = resolved_schema.get(
        "p_value"
    )

    if (
        feature_column is None
        or
        variant_column is None
    ):

        return pd.DataFrame(
            columns=[
                "feature_id",
                "association_rows",
                "unique_canonical_rsids",
                "minimum_dense_threshold_pass",
                "strict_dense_threshold_pass",
                "minimum_p_value",
                "maximum_p_value",
                "p_values_above_0_05",
                "coverage_interpretation",
            ]
        )

    policy = config[
        "inspection_policy"
    ]

    minimum_dense = int(
        policy[
            "minimum_dense_variants_per_feature"
        ]
    )

    strict_dense = int(
        policy[
            "minimum_dense_variants_per_feature_strict"
        ]
    )

    source = dataframe.copy()

    source[
        "_feature"
    ] = (
        source[
            feature_column
        ]
        .astype("string")
        .str.strip()
    )

    source[
        "_rsid"
    ] = _normalize_rsid_series(
        source[
            variant_column
        ]
    )

    if p_column is not None:

        source[
            "_p"
        ] = _numeric(
            source[
                p_column
            ]
        )

    else:

        source[
            "_p"
        ] = np.nan

    rows: list[
        dict[str, Any]
    ] = []

    for feature_id, group in source.groupby(
        "_feature",
        sort=False,
        dropna=True,
    ):

        p_values = (
            group[
                "_p"
            ]
            .dropna()
        )

        unique_rsids = int(
            group[
                "_rsid"
            ]
            .dropna()
            .nunique()
        )

        rows.append(
            {
                "feature_id":
                    str(
                        feature_id
                    ),

                "association_rows":
                    int(
                        len(
                            group
                        )
                    ),

                "unique_canonical_rsids":
                    unique_rsids,

                "minimum_dense_threshold_pass":
                    bool(
                        unique_rsids
                        >=
                        minimum_dense
                    ),

                "strict_dense_threshold_pass":
                    bool(
                        unique_rsids
                        >=
                        strict_dense
                    ),

                "minimum_p_value":
                    (
                        None
                        if p_values.empty
                        else float(
                            p_values.min()
                        )
                    ),

                "maximum_p_value":
                    (
                        None
                        if p_values.empty
                        else float(
                            p_values.max()
                        )
                    ),

                "p_values_above_0_05":
                    int(
                        p_values.gt(
                            0.05
                        ).sum()
                    ),

                "coverage_interpretation":
                    (
                        "DESCRIPTIVE_SIGNIFICANT_HIT_COVERAGE_ONLY"
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Dataset readiness
# ============================================================================


def assess_dataset_readiness(
    dataframe: pd.DataFrame,
    *,
    resolved_schema: dict[str, str | None],
    schema_table: pd.DataFrame,
    feature_coverage: pd.DataFrame,
    config: dict[str, Any],
    source_file_count: int,
) -> pd.DataFrame:
    """Assess actual CancerSplicingQTL PRAD dataset readiness."""

    if dataframe.empty:

        return pd.DataFrame(
            [
                {
                    "resource":
                        "CancerSplicingQTL",

                    "classification":
                        EMPTY_DATASET,

                    "formal_coloc_ready":
                        False,

                    "coloc_abf_ready":
                        False,

                    "coloc_susie_ready":
                        False,
                }
            ]
        )

    policy = config[
        "inspection_policy"
    ]

    source_scope = config[
        "source_scope"
    ]

    # ----------------------------------------------------------------------
    # Required schema
    # ----------------------------------------------------------------------

    required_missing = schema_table.loc[
        schema_table[
            "required"
        ].astype(
            bool
        )
        &
        ~schema_table[
            "available"
        ].astype(
            bool
        ),
        "logical_field",
    ].astype(
        str
    ).tolist()

    if required_missing:

        classification = (
            SCHEMA_MISMATCH
        )

    else:

        classification = None

    # ----------------------------------------------------------------------
    # Cancer-type validation
    # ----------------------------------------------------------------------

    cancer_column = resolved_schema.get(
        "cancer_type"
    )

    expected_cancer = str(
        policy[
            "expected_cancer_value"
        ]
    ).strip().upper()

    if cancer_column is not None:

        cancer_values = (
            dataframe[
                cancer_column
            ]
            .astype("string")
            .str.strip()
            .str.upper()
        )

        non_prad_rows = int(
            cancer_values.ne(
                expected_cancer
            )
            .fillna(
                True
            )
            .sum()
        )

    else:

        non_prad_rows = int(
            len(
                dataframe
            )
        )

    if (
        classification is None
        and
        bool(
            policy[
                "require_prad_rows_only"
            ]
        )
        and
        non_prad_rows > 0
    ):

        classification = (
            NON_PRAD_CONTENT_DETECTED
        )

    # ----------------------------------------------------------------------
    # Canonical rsIDs
    # ----------------------------------------------------------------------

    variant_column = resolved_schema.get(
        "variant_rsid"
    )

    if variant_column is not None:

        rsids = _normalize_rsid_series(
            dataframe[
                variant_column
            ]
        )

        canonical_rsid_count = int(
            rsids.dropna().nunique()
        )

        canonical_rsid_fraction = float(
            rsids.notna().mean()
        )

    else:

        canonical_rsid_count = 0
        canonical_rsid_fraction = 0.0

    if (
        classification is None
        and
        canonical_rsid_fraction
        <
        float(
            policy[
                "minimum_canonical_rsid_fraction"
            ]
        )
    ):

        classification = (
            VARIANT_IDENTITY_INCOMPLETE
        )

    # ----------------------------------------------------------------------
    # Feature identity
    # ----------------------------------------------------------------------

    feature_column = resolved_schema.get(
        "feature_id"
    )

    if feature_column is not None:

        feature_values = (
            dataframe[
                feature_column
            ]
            .astype("string")
            .str.strip()
        )

        feature_fraction = float(
            feature_values.notna().mean()
        )

        feature_count = int(
            feature_values.dropna().nunique()
        )

    else:

        feature_fraction = 0.0
        feature_count = 0

    if (
        classification is None
        and
        feature_fraction < 0.95
    ):

        classification = (
            FEATURE_IDENTITY_INCOMPLETE
        )

    # ----------------------------------------------------------------------
    # Numeric fields
    # ----------------------------------------------------------------------

    beta_column = resolved_schema.get(
        "beta"
    )

    t_column = resolved_schema.get(
        "t_statistic"
    )

    p_column = resolved_schema.get(
        "p_value"
    )

    if beta_column is not None:

        beta = _numeric(
            dataframe[
                beta_column
            ]
        )

        beta_fraction = float(
            beta.notna().mean()
        )

    else:

        beta_fraction = 0.0

    if t_column is not None:

        t_stat = _numeric(
            dataframe[
                t_column
            ]
        )

        t_fraction = float(
            t_stat.notna().mean()
        )

    else:

        t_fraction = 0.0

    if p_column is not None:

        p_values = (
            _numeric(
                dataframe[
                    p_column
                ]
            )
            .dropna()
        )

        p_fraction = float(
            _numeric(
                dataframe[
                    p_column
                ]
            )
            .notna()
            .mean()
        )

    else:

        p_values = pd.Series(
            dtype="float64"
        )

        p_fraction = 0.0

    # ----------------------------------------------------------------------
    # P-value diagnostics
    #
    # These do NOT determine source completeness.
    # ----------------------------------------------------------------------

    minimum_p = (
        None
        if p_values.empty
        else float(
            p_values.min()
        )
    )

    maximum_p = (
        None
        if p_values.empty
        else float(
            p_values.max()
        )
    )

    p_above_005 = int(
        p_values.gt(
            0.05
        ).sum()
    )

    p_at_or_below_005 = int(
        p_values.le(
            0.05
        ).sum()
    )

    # ----------------------------------------------------------------------
    # Feature-coverage diagnostics
    # ----------------------------------------------------------------------

    if feature_coverage.empty:

        dense_feature_count = 0
        strict_dense_feature_count = 0
        maximum_variants_per_feature = 0
        median_variants_per_feature = 0.0

    else:

        dense_feature_count = int(
            feature_coverage[
                "minimum_dense_threshold_pass"
            ]
            .fillna(
                False
            )
            .sum()
        )

        strict_dense_feature_count = int(
            feature_coverage[
                "strict_dense_threshold_pass"
            ]
            .fillna(
                False
            )
            .sum()
        )

        maximum_variants_per_feature = _safe_int(
            feature_coverage[
                "unique_canonical_rsids"
            ].max()
        )

        median_variants_per_feature = float(
            feature_coverage[
                "unique_canonical_rsids"
            ].median()
        )

    # ----------------------------------------------------------------------
    # Diagnostic standard-error derivation
    # ----------------------------------------------------------------------

    se_diagnostic_performed = False

    se_diagnostic_valid_rows = 0

    se_diagnostic_valid_fraction = 0.0

    derived_se_minimum = None

    derived_se_median = None

    derived_se_maximum = None

    if (
        beta_column is not None
        and
        t_column is not None
        and
        bool(
            config[
                "standard_error_audit"
            ][
                "perform_diagnostic_derivation"
            ]
        )
    ):

        se_diagnostic = (
            derive_standard_error_diagnostic(
                dataframe,
                beta_column=beta_column,
                t_column=t_column,
            )
        )

        valid_se = (
            se_diagnostic[
                "derived_se_diagnostic"
            ]
            .dropna()
        )

        se_diagnostic_performed = True

        se_diagnostic_valid_rows = int(
            len(
                valid_se
            )
        )

        se_diagnostic_valid_fraction = float(
            se_diagnostic[
                "derived_se_valid"
            ]
            .astype(
                bool
            )
            .mean()
        )

        if not valid_se.empty:

            derived_se_minimum = float(
                valid_se.min()
            )

            derived_se_median = float(
                valid_se.median()
            )

            derived_se_maximum = float(
                valid_se.max()
            )

    # ----------------------------------------------------------------------
    # Source completeness
    #
    # This is the decisive criterion.
    # ----------------------------------------------------------------------

    significant_only_source = bool(
        source_scope[
            "source_selection"
        ][
            "significance_filtered"
        ]
    )

    full_test_matrix_available = bool(
        source_scope[
            "complete_tested_variant_feature_matrix"
        ][
            "available"
        ]
    )

    nonsignificant_universe_available = bool(
        source_scope[
            "nonsignificant_test_universe"
        ][
            "available"
        ]
    )

    dense_unfiltered_available = bool(
        source_scope[
            "dense_unfiltered_summary_statistics"
        ][
            "available"
        ]
    )

    # ----------------------------------------------------------------------
    # Formal availability metadata
    # ----------------------------------------------------------------------

    direct_se_available = False

    varbeta_available = False

    sample_size_available = False

    allele_frequency_available = False

    effect_allele_orientation_verified = bool(
        config[
            "allele_policy"
        ][
            "effect_allele_orientation_verified"
        ]
    )

    # ----------------------------------------------------------------------
    # Final classification
    #
    # The significant-only provenance takes precedence over raw row count,
    # p-value distribution, or apparent per-feature SNP density.
    # ----------------------------------------------------------------------

    if classification is None:

        if (
            significant_only_source
            or
            not full_test_matrix_available
            or
            not dense_unfiltered_available
        ):

            classification = (
                SIGNIFICANT_ONLY_NOT_COLOC_READY
            )

        else:

            abf_requirements_met = bool(
                canonical_rsid_fraction
                >=
                float(
                    policy[
                        "minimum_canonical_rsid_fraction"
                    ]
                )
                and
                feature_fraction >= 0.95
                and
                beta_fraction
                >=
                float(
                    policy[
                        "minimum_numeric_beta_fraction"
                    ]
                )
                and
                p_fraction
                >=
                float(
                    policy[
                        "minimum_numeric_p_fraction"
                    ]
                )
                and
                (
                    direct_se_available
                    or
                    varbeta_available
                )
                and
                sample_size_available
            )

            if abf_requirements_met:

                classification = (
                    FORMAL_ABF_READY_SUSIE_BLOCKED
                )

            else:

                classification = (
                    SCHEMA_MISMATCH
                )

    formal_coloc_ready = bool(
        classification
        in {
            FORMAL_ABF_READY_SUSIE_BLOCKED,
            FORMAL_COLOC_READY,
        }
    )

    coloc_abf_ready = formal_coloc_ready

    coloc_susie_ready = bool(
        classification
        ==
        FORMAL_COLOC_READY
    )

    return pd.DataFrame(
        [
            {
                "resource":
                    str(
                        config[
                            "resource"
                        ][
                            "display_name"
                        ]
                    ),

                "cancer":
                    expected_cancer,

                "genome_build":
                    str(
                        config[
                            "resource"
                        ][
                            "genome_build"
                        ][
                            "reported_build"
                        ]
                    ),

                "genome_build_alias":
                    str(
                        config[
                            "resource"
                        ][
                            "genome_build"
                        ][
                            "alias"
                        ]
                    ),

                "source_file_count":
                    int(
                        source_file_count
                    ),

                "association_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "non_prad_rows":
                    non_prad_rows,

                "required_schema_missing":
                    required_missing,

                "canonical_rsid_count":
                    canonical_rsid_count,

                "canonical_rsid_fraction":
                    canonical_rsid_fraction,

                "feature_count":
                    feature_count,

                "feature_identity_fraction":
                    feature_fraction,

                "beta_numeric_fraction":
                    beta_fraction,

                "t_stat_numeric_fraction":
                    t_fraction,

                "p_value_numeric_fraction":
                    p_fraction,

                "minimum_p_value":
                    minimum_p,

                "maximum_p_value":
                    maximum_p,

                "p_values_at_or_below_0_05":
                    p_at_or_below_005,

                "p_values_above_0_05":
                    p_above_005,

                "dense_feature_count_descriptive":
                    dense_feature_count,

                "strict_dense_feature_count_descriptive":
                    strict_dense_feature_count,

                "maximum_variants_per_feature_descriptive":
                    maximum_variants_per_feature,

                "median_variants_per_feature_descriptive":
                    median_variants_per_feature,

                "source_significance_filtered":
                    significant_only_source,

                "complete_tested_variant_feature_matrix_available":
                    full_test_matrix_available,

                "nonsignificant_test_universe_available":
                    nonsignificant_universe_available,

                "dense_unfiltered_summary_statistics_available":
                    dense_unfiltered_available,

                "beta_directly_reported":
                    bool(
                        beta_column is not None
                    ),

                "t_statistic_directly_reported":
                    bool(
                        t_column is not None
                    ),

                "standard_error_directly_reported":
                    direct_se_available,

                "varbeta_directly_reported":
                    varbeta_available,

                "derived_se_diagnostic_performed":
                    se_diagnostic_performed,

                "derived_se_diagnostic_valid_rows":
                    se_diagnostic_valid_rows,

                "derived_se_diagnostic_valid_fraction":
                    se_diagnostic_valid_fraction,

                "derived_se_diagnostic_minimum":
                    derived_se_minimum,

                "derived_se_diagnostic_median":
                    derived_se_median,

                "derived_se_diagnostic_maximum":
                    derived_se_maximum,

                "derived_se_allowed_for_formal_coloc":
                    bool(
                        config[
                            "standard_error_audit"
                        ][
                            "use_derived_se_for_formal_coloc"
                        ]
                    ),

                "sample_size_column_available":
                    sample_size_available,

                "allele_frequency_available":
                    allele_frequency_available,

                "effect_allele_orientation_verified":
                    effect_allele_orientation_verified,

                "cross_build_coordinate_join_allowed":
                    False,

                "canonical_rsid_cross_build_matching_allowed":
                    True,

                "classification":
                    classification,

                "formal_coloc_ready":
                    formal_coloc_ready,

                "coloc_abf_ready":
                    coloc_abf_ready,

                "coloc_susie_ready":
                    coloc_susie_ready,

                "formal_colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            }
        ]
    )


# ============================================================================
# Candidate inspection
# ============================================================================


def assess_candidate_presence(
    dataframe: pd.DataFrame,
    *,
    resolved_schema: dict[str, str | None],
    readiness: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Assess direct candidate rsID presence.

    Candidate presence in a significance-filtered sQTL table is descriptive
    evidence only and cannot override dataset-level coloc blocking.
    """

    variant_column = resolved_schema.get(
        "variant_rsid"
    )

    if variant_column is None:

        available_rsids: set[
            str
        ] = set()

    else:

        available_rsids = set(
            _normalize_rsid_series(
                dataframe[
                    variant_column
                ]
            )
            .dropna()
            .astype(
                str
            )
        )

    readiness_row = readiness.iloc[
        0
    ]

    resource_classification = str(
        readiness_row[
            "classification"
        ]
    )

    resource_ready = bool(
        readiness_row[
            "formal_coloc_ready"
        ]
    )

    rows: list[
        dict[str, Any]
    ] = []

    for lead_rsid in config[
        "candidate_leads"
    ]:

        canonical = (
            str(
                lead_rsid
            )
            .strip()
            .lower()
        )

        present = bool(
            canonical
            in
            available_rsids
        )

        if present:

            evidence_status = (
                "DIRECT_SIGNIFICANT_SQTL_RECORD_PRESENT"
            )

        else:

            evidence_status = (
                "NO_DIRECT_RECORD_IN_SIGNIFICANCE_FILTERED_RESOURCE"
            )

        rows.append(
            {
                "lead_rsid":
                    canonical,

                "direct_rsid_present":
                    present,

                "candidate_evidence_status":
                    evidence_status,

                "absence_interpreted_as_null":
                    False,

                "resource_classification":
                    resource_classification,

                "resource_formal_coloc_ready":
                    resource_ready,

                "candidate_formal_coloc_ready":
                    bool(
                        resource_ready
                        and
                        present
                    ),

                "recommended_action":
                    (
                        "RETAIN_AS_DESCRIPTIVE_REGULATORY_EVIDENCE"
                        if present
                        else
                        "DO_NOT_INTERPRET_ABSENCE_AS_NO_QTL_EFFECT"
                    ),

                "cross_build_matching_method":
                    "canonical_rsid_only",

                "formal_colocalization_performed":
                    False,

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Public API
# ============================================================================


def inspect_cancersplicingqtl(
    dataframe: pd.DataFrame,
    *,
    config: dict[str, Any],
    source_file_count: int,
) -> CancerSplicingQTLInspectionResult:
    """Execute complete M5.6C.1 inspection."""

    resolved_schema = (
        resolve_expected_schema(
            dataframe,
            config=config,
        )
    )

    schema = (
        build_schema_table(
            dataframe,
            resolved_schema=resolved_schema,
            config=config,
        )
    )

    feature_coverage = (
        assess_feature_coverage(
            dataframe,
            resolved_schema=resolved_schema,
            config=config,
        )
    )

    readiness = (
        assess_dataset_readiness(
            dataframe,
            resolved_schema=resolved_schema,
            schema_table=schema,
            feature_coverage=feature_coverage,
            config=config,
            source_file_count=source_file_count,
        )
    )

    candidates = (
        assess_candidate_presence(
            dataframe,
            resolved_schema=resolved_schema,
            readiness=readiness,
            config=config,
        )
    )

    return CancerSplicingQTLInspectionResult(
        schema=schema,
        feature_coverage=feature_coverage,
        readiness=readiness,
        candidates=candidates,
    )