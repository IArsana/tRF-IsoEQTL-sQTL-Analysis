"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_locus_standardization.py

Description:
    M5.3C.3D standardization and descriptive disease-signal QC utilities
    for retrieved prostate cancer GWAS loci.

    The module preserves all retrieved harmonised GWAS rows while producing
    a common analysis schema and descriptive regional disease-signal metrics.

    Variant-representation QC distinguishes:

        repeated_variant_identity
            Two or more source rows share the same normalized variant
            identity within the same study × candidate locus.

        exact_duplicate
            Repeated source rows are identical across the relevant
            association and harmonisation fields.

        conflicting_variant_representation
            Repeated source rows share the same normalized variant identity
            but differ in one or more relevant fields such as reference
            allele, effect estimate, standard error, allele frequency,
            p-value, or harmonisation metadata.

    No rows are removed at this stage.

    Significance definitions:
        Genome-wide significant:
            p < 5e-8

        Suggestive:
            p < 1e-5

        Suggestive-only:
            5e-8 <= p < 1e-5

    Important scientific safeguards:
        - All source rows are preserved.
        - No significance filtering is performed.
        - No source representation is arbitrarily selected.
        - No deduplication is performed.
        - No LD pruning is performed.
        - Physical distance is not interpreted as LD.
        - Regional significance does not imply colocalization.
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


# ============================================================================
# Constants
# ============================================================================


GENOME_WIDE_SIGNIFICANCE = 5e-8

SUGGESTIVE_SIGNIFICANCE = 1e-5


TEXTUAL_MISSING_VALUES = {
    "",
    "NA",
    "N/A",
    "NAN",
    "NULL",
    "NONE",
    ".",
}


REPRESENTATION_COMPARISON_COLUMNS = [
    "chromosome",
    "base_pair_location",
    "rsid",
    "variant_id",
    "effect_allele",
    "other_allele",
    "reference_allele",
    "beta",
    "standard_error",
    "effect_allele_frequency",
    "p_value",
    "hm_coordinate_conversion",
    "hm_code",
]


STANDARD_COLUMNS = [
    "study_accession",
    "lead_rsid",
    "chromosome",
    "base_pair_location",
    "rsid",
    "variant_id",
    "effect_allele",
    "other_allele",
    "reference_allele",
    "beta",
    "standard_error",
    "effect_allele_frequency",
    "p_value",
    "hm_coordinate_conversion",
    "hm_code",
    "lead_position",
    "distance_to_lead",
    "is_lead_variant",
    "is_genome_wide_significant",
    "is_suggestive_significant",
    "is_suggestive_only",
    "variant_identity_key",
    "is_repeated_variant_identity",
    "is_exact_duplicate_representation",
    "is_conflicting_variant_representation",
]


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class LocusStandardizationResult:
    """Standardized GWAS locus plus QC counters."""

    dataframe: pd.DataFrame

    input_rows: int

    output_rows: int

    malformed_coordinate_rows: int

    missing_p_value_rows: int

    missing_beta_rows: int

    missing_standard_error_rows: int

    missing_effect_allele_frequency_rows: int

    repeated_variant_identity_rows: int

    repeated_variant_identity_keys: int

    exact_duplicate_rows: int

    exact_duplicate_keys: int

    conflicting_variant_representation_rows: int

    conflicting_variant_representation_keys: int


# ============================================================================
# Source-column helpers
# ============================================================================


def _series_or_na(
    dataframe: pd.DataFrame,
    column: str,
) -> pd.Series:
    """Return one source column or an NA-filled replacement."""

    if column in dataframe.columns:
        return dataframe[column]

    return pd.Series(
        pd.NA,
        index=dataframe.index,
        dtype="object",
    )


def _nullable_string(
    series: pd.Series,
) -> pd.Series:
    """Normalize source strings and textual missing-value representations."""

    result = (
        series
        .astype("string")
        .str.strip()
    )

    upper = result.str.upper()

    return result.mask(
        upper.isin(
            TEXTUAL_MISSING_VALUES
        )
    )


def _nullable_numeric(
    series: pd.Series,
) -> pd.Series:
    """Convert source values to numeric with invalid values mapped to NA."""

    return pd.to_numeric(
        _nullable_string(
            series
        ),
        errors="coerce",
    )


# ============================================================================
# Variant normalization
# ============================================================================


def _normalize_rsid_series(
    series: pd.Series,
) -> pd.Series:
    """Normalize canonical rsIDs."""

    result = (
        _nullable_string(
            series
        )
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


def _normalize_chromosome_series(
    series: pd.Series,
) -> pd.Series:
    """Normalize chromosome labels without chr prefix."""

    result = (
        _nullable_string(
            series
        )
        .str.replace(
            r"^chr",
            "",
            regex=True,
            case=False,
        )
        .str.upper()
    )

    return result.replace(
        {
            "M": "MT",
        }
    )


# ============================================================================
# Variant identity
# ============================================================================


def _build_variant_identity_key(
    dataframe: pd.DataFrame,
) -> pd.Series:
    """
    Build normalized within-locus variant identity.

    Priority:
        1. variant_id
        2. chromosome:position:effect_allele:other_allele
    """

    identity = pd.Series(
        pd.NA,
        index=dataframe.index,
        dtype="string",
    )

    variant_id = (
        dataframe[
            "variant_id"
        ]
        .astype("string")
    )

    has_variant_id = (
        variant_id.notna()
    )

    identity.loc[
        has_variant_id
    ] = (
        "VID:"
        + variant_id.loc[
            has_variant_id
        ]
    )

    fallback_usable = (
        dataframe[
            "chromosome"
        ].notna()
        &
        dataframe[
            "base_pair_location"
        ].notna()
        &
        dataframe[
            "effect_allele"
        ].notna()
        &
        dataframe[
            "other_allele"
        ].notna()
    )

    fallback_key = (
        dataframe[
            "chromosome"
        ].astype("string")
        + ":"
        + dataframe[
            "base_pair_location"
        ].astype("string")
        + ":"
        + dataframe[
            "effect_allele"
        ].astype("string")
        + ":"
        + dataframe[
            "other_allele"
        ].astype("string")
    )

    use_fallback = (
        ~has_variant_id
        &
        fallback_usable
    )

    identity.loc[
        use_fallback
    ] = (
        "COORD:"
        + fallback_key.loc[
            use_fallback
        ]
    )

    return identity


# ============================================================================
# Repeated representation QC
# ============================================================================


def _normalized_comparison_value(
    value: Any,
) -> Any:
    """Normalize scalar values for exact-representation comparison."""

    if pd.isna(
        value
    ):
        return None

    if isinstance(
        value,
        str,
    ):
        return value.strip()

    return value


def _classify_variant_representations(
    dataframe: pd.DataFrame,
) -> tuple[
    pd.Series,
    pd.Series,
    pd.Series,
    dict[str, int],
]:
    """
    Classify repeated source representations within one locus.

    Returns:
        repeated_mask
        exact_duplicate_mask
        conflicting_mask
        metrics
    """

    identity = dataframe[
        "variant_identity_key"
    ]

    repeated_mask = (
        identity.notna()
        &
        identity.duplicated(
            keep=False
        )
    )

    exact_mask = pd.Series(
        False,
        index=dataframe.index,
        dtype=bool,
    )

    conflicting_mask = pd.Series(
        False,
        index=dataframe.index,
        dtype=bool,
    )

    repeated = dataframe.loc[
        repeated_mask
    ]

    repeated_keys = int(
        repeated[
            "variant_identity_key"
        ].nunique()
    )

    exact_keys = 0
    conflicting_keys = 0

    for _, group in repeated.groupby(
        "variant_identity_key",
        sort=False,
    ):

        normalized_rows = []

        for _, row in group.iterrows():

            normalized_rows.append(
                tuple(
                    _normalized_comparison_value(
                        row[
                            column
                        ]
                    )
                    for column
                    in REPRESENTATION_COMPARISON_COLUMNS
                )
            )

        unique_representations = set(
            normalized_rows
        )

        if len(
            unique_representations
        ) == 1:

            exact_keys += 1

            exact_mask.loc[
                group.index
            ] = True

        else:

            conflicting_keys += 1

            conflicting_mask.loc[
                group.index
            ] = True

    metrics = {
        "repeated_variant_identity_rows":
            int(
                repeated_mask.sum()
            ),

        "repeated_variant_identity_keys":
            repeated_keys,

        "exact_duplicate_rows":
            int(
                exact_mask.sum()
            ),

        "exact_duplicate_keys":
            exact_keys,

        "conflicting_variant_representation_rows":
            int(
                conflicting_mask.sum()
            ),

        "conflicting_variant_representation_keys":
            conflicting_keys,
    }

    return (
        repeated_mask,
        exact_mask,
        conflicting_mask,
        metrics,
    )


# ============================================================================
# Standardization
# ============================================================================


def standardize_gwas_locus(
    dataframe: pd.DataFrame,
    *,
    study_accession: str,
    lead_rsid: str,
    lead_position: int,
    expected_chromosome: str,
) -> LocusStandardizationResult:
    """Standardize one retrieved harmonised GWAS locus."""

    input_rows = int(
        len(
            dataframe
        )
    )

    out = pd.DataFrame(
        index=dataframe.index,
    )

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    out[
        "study_accession"
    ] = (
        str(
            study_accession
        )
        .strip()
        .upper()
    )

    out[
        "lead_rsid"
    ] = (
        str(
            lead_rsid
        )
        .strip()
        .lower()
    )

    # ------------------------------------------------------------------
    # Variant identity
    # ------------------------------------------------------------------

    out[
        "chromosome"
    ] = _normalize_chromosome_series(
        _series_or_na(
            dataframe,
            "chromosome",
        )
    )

    out[
        "base_pair_location"
    ] = (
        _nullable_numeric(
            _series_or_na(
                dataframe,
                "base_pair_location",
            )
        )
        .astype(
            "Int64"
        )
    )

    out[
        "rsid"
    ] = _normalize_rsid_series(
        _series_or_na(
            dataframe,
            "rsid",
        )
    )

    out[
        "variant_id"
    ] = _nullable_string(
        _series_or_na(
            dataframe,
            "variant_id",
        )
    )

    # ------------------------------------------------------------------
    # Alleles
    # ------------------------------------------------------------------

    for column in (
        "effect_allele",
        "other_allele",
        "reference_allele",
    ):

        out[
            column
        ] = (
            _nullable_string(
                _series_or_na(
                    dataframe,
                    column,
                )
            )
            .str.upper()
        )

    # ------------------------------------------------------------------
    # Association statistics
    # ------------------------------------------------------------------

    for column in (
        "beta",
        "standard_error",
        "effect_allele_frequency",
        "p_value",
    ):

        out[
            column
        ] = _nullable_numeric(
            _series_or_na(
                dataframe,
                column,
            )
        )

    # ------------------------------------------------------------------
    # Harmonisation provenance
    # ------------------------------------------------------------------

    out[
        "hm_coordinate_conversion"
    ] = _nullable_string(
        _series_or_na(
            dataframe,
            "hm_coordinate_conversion",
        )
    )

    out[
        "hm_code"
    ] = _nullable_string(
        _series_or_na(
            dataframe,
            "hm_code",
        )
    )

    # ------------------------------------------------------------------
    # Candidate lead relationship
    # ------------------------------------------------------------------

    out[
        "lead_position"
    ] = int(
        lead_position
    )

    out[
        "distance_to_lead"
    ] = (
        out[
            "base_pair_location"
        ]
        - int(
            lead_position
        )
    ).astype(
        "Int64"
    )

    normalized_lead = (
        str(
            lead_rsid
        )
        .strip()
        .lower()
    )

    out[
        "is_lead_variant"
    ] = (
        out[
            "rsid"
        ]
        .eq(
            normalized_lead
        )
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    # ------------------------------------------------------------------
    # Significance flags
    # ------------------------------------------------------------------

    out[
        "is_genome_wide_significant"
    ] = (
        out[
            "p_value"
        ]
        .lt(
            GENOME_WIDE_SIGNIFICANCE
        )
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    out[
        "is_suggestive_significant"
    ] = (
        out[
            "p_value"
        ]
        .lt(
            SUGGESTIVE_SIGNIFICANCE
        )
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    out[
        "is_suggestive_only"
    ] = (
        out[
            "p_value"
        ]
        .ge(
            GENOME_WIDE_SIGNIFICANCE
        )
        &
        out[
            "p_value"
        ]
        .lt(
            SUGGESTIVE_SIGNIFICANCE
        )
    ).fillna(
        False
    ).astype(
        bool
    )

    # ------------------------------------------------------------------
    # Variant identity + source representation QC
    # ------------------------------------------------------------------

    out[
        "variant_identity_key"
    ] = _build_variant_identity_key(
        out
    )

    (
        repeated_mask,
        exact_mask,
        conflicting_mask,
        representation_metrics,
    ) = _classify_variant_representations(
        out
    )

    out[
        "is_repeated_variant_identity"
    ] = repeated_mask.astype(
        bool
    )

    out[
        "is_exact_duplicate_representation"
    ] = exact_mask.astype(
        bool
    )

    out[
        "is_conflicting_variant_representation"
    ] = conflicting_mask.astype(
        bool
    )

    # ------------------------------------------------------------------
    # Coordinate QC
    # ------------------------------------------------------------------

    expected_chr = (
        str(
            expected_chromosome
        )
        .strip()
        .upper()
    )

    if expected_chr.startswith(
        "CHR"
    ):

        expected_chr = expected_chr[
            3:
        ]

    if expected_chr == "M":

        expected_chr = "MT"

    malformed_coordinate_mask = (
        out[
            "chromosome"
        ].isna()
        |
        out[
            "base_pair_location"
        ].isna()
        |
        ~out[
            "chromosome"
        ].eq(
            expected_chr
        )
    )

    # ------------------------------------------------------------------
    # Final schema + cardinality
    # ------------------------------------------------------------------

    out = out[
        STANDARD_COLUMNS
    ].copy()

    output_rows = int(
        len(
            out
        )
    )

    if output_rows != input_rows:

        raise RuntimeError(
            "GWAS locus standardization changed row cardinality: "
            f"input={input_rows}, output={output_rows}"
        )

    return LocusStandardizationResult(
        dataframe=out,
        input_rows=input_rows,
        output_rows=output_rows,
        malformed_coordinate_rows=int(
            malformed_coordinate_mask.sum()
        ),
        missing_p_value_rows=int(
            out[
                "p_value"
            ].isna().sum()
        ),
        missing_beta_rows=int(
            out[
                "beta"
            ].isna().sum()
        ),
        missing_standard_error_rows=int(
            out[
                "standard_error"
            ].isna().sum()
        ),
        missing_effect_allele_frequency_rows=int(
            out[
                "effect_allele_frequency"
            ].isna().sum()
        ),
        repeated_variant_identity_rows=(
            representation_metrics[
                "repeated_variant_identity_rows"
            ]
        ),
        repeated_variant_identity_keys=(
            representation_metrics[
                "repeated_variant_identity_keys"
            ]
        ),
        exact_duplicate_rows=(
            representation_metrics[
                "exact_duplicate_rows"
            ]
        ),
        exact_duplicate_keys=(
            representation_metrics[
                "exact_duplicate_keys"
            ]
        ),
        conflicting_variant_representation_rows=(
            representation_metrics[
                "conflicting_variant_representation_rows"
            ]
        ),
        conflicting_variant_representation_keys=(
            representation_metrics[
                "conflicting_variant_representation_keys"
            ]
        ),
    )


# ============================================================================
# Scalar converters
# ============================================================================


def _float_or_none(
    value: Any,
) -> float | None:

    if pd.isna(
        value
    ):
        return None

    return float(
        value
    )


def _int_or_none(
    value: Any,
) -> int | None:

    if pd.isna(
        value
    ):
        return None

    return int(
        value
    )


def _str_or_none(
    value: Any,
) -> str | None:

    if pd.isna(
        value
    ):
        return None

    return str(
        value
    )


# ============================================================================
# Disease signal summary
# ============================================================================


def summarize_disease_signal(
    dataframe: pd.DataFrame,
    *,
    study_accession: str,
    lead_rsid: str,
) -> dict[str, Any]:
    """Summarize descriptive disease-association evidence for one locus."""

    if dataframe.empty:

        return {
            "study_accession":
                study_accession,

            "lead_rsid":
                lead_rsid,

            "variant_rows":
                0,

            "variants_with_p_value":
                0,

            "lead_present":
                False,

            "lead_p_value":
                None,

            "lead_beta":
                None,

            "lead_standard_error":
                None,

            "lead_effect_allele_frequency":
                None,

            "minimum_p_value":
                None,

            "top_variant_rsid":
                None,

            "top_variant_id":
                None,

            "top_variant_position":
                None,

            "top_variant_distance_to_lead":
                None,

            "top_variant_effect_allele":
                None,

            "top_variant_other_allele":
                None,

            "top_variant_beta":
                None,

            "top_variant_standard_error":
                None,

            "top_variant_eaf":
                None,

            "genome_wide_significant_variants":
                0,

            "suggestive_significant_variants":
                0,

            "suggestive_only_variants":
                0,

            "has_genome_wide_signal":
                False,

            "has_suggestive_signal":
                False,

            "has_suggestive_only_signal":
                False,
        }

    # ------------------------------------------------------------------
    # Lead
    # ------------------------------------------------------------------

    lead_rows = dataframe.loc[
        dataframe[
            "is_lead_variant"
        ]
    ]

    lead_present = (
        not lead_rows.empty
    )

    lead_p = None
    lead_beta = None
    lead_se = None
    lead_eaf = None

    if lead_present:

        usable_lead = lead_rows.loc[
            lead_rows[
                "p_value"
            ].notna()
        ]

        if not usable_lead.empty:

            lead_row = dataframe.loc[
                usable_lead[
                    "p_value"
                ].idxmin()
            ]

        else:

            lead_row = lead_rows.iloc[
                0
            ]

        lead_p = _float_or_none(
            lead_row[
                "p_value"
            ]
        )

        lead_beta = _float_or_none(
            lead_row[
                "beta"
            ]
        )

        lead_se = _float_or_none(
            lead_row[
                "standard_error"
            ]
        )

        lead_eaf = _float_or_none(
            lead_row[
                "effect_allele_frequency"
            ]
        )

    # ------------------------------------------------------------------
    # Regional top variant
    # ------------------------------------------------------------------

    top_values = {
        "minimum_p_value":
            None,

        "top_variant_rsid":
            None,

        "top_variant_id":
            None,

        "top_variant_position":
            None,

        "top_variant_distance_to_lead":
            None,

        "top_variant_effect_allele":
            None,

        "top_variant_other_allele":
            None,

        "top_variant_beta":
            None,

        "top_variant_standard_error":
            None,

        "top_variant_eaf":
            None,
    }

    p_valid = dataframe.loc[
        dataframe[
            "p_value"
        ].notna()
    ]

    if not p_valid.empty:

        top = dataframe.loc[
            p_valid[
                "p_value"
            ].idxmin()
        ]

        top_values = {
            "minimum_p_value":
                _float_or_none(
                    top[
                        "p_value"
                    ]
                ),

            "top_variant_rsid":
                _str_or_none(
                    top[
                        "rsid"
                    ]
                ),

            "top_variant_id":
                _str_or_none(
                    top[
                        "variant_id"
                    ]
                ),

            "top_variant_position":
                _int_or_none(
                    top[
                        "base_pair_location"
                    ]
                ),

            "top_variant_distance_to_lead":
                _int_or_none(
                    top[
                        "distance_to_lead"
                    ]
                ),

            "top_variant_effect_allele":
                _str_or_none(
                    top[
                        "effect_allele"
                    ]
                ),

            "top_variant_other_allele":
                _str_or_none(
                    top[
                        "other_allele"
                    ]
                ),

            "top_variant_beta":
                _float_or_none(
                    top[
                        "beta"
                    ]
                ),

            "top_variant_standard_error":
                _float_or_none(
                    top[
                        "standard_error"
                    ]
                ),

            "top_variant_eaf":
                _float_or_none(
                    top[
                        "effect_allele_frequency"
                    ]
                ),
        }

    genome_wide_count = int(
        dataframe[
            "is_genome_wide_significant"
        ].sum()
    )

    suggestive_count = int(
        dataframe[
            "is_suggestive_significant"
        ].sum()
    )

    suggestive_only_count = int(
        dataframe[
            "is_suggestive_only"
        ].sum()
    )

    return {
        "study_accession":
            str(
                study_accession
            ),

        "lead_rsid":
            str(
                lead_rsid
            ),

        "variant_rows":
            int(
                len(
                    dataframe
                )
            ),

        "variants_with_p_value":
            int(
                dataframe[
                    "p_value"
                ].notna().sum()
            ),

        "lead_present":
            bool(
                lead_present
            ),

        "lead_p_value":
            lead_p,

        "lead_beta":
            lead_beta,

        "lead_standard_error":
            lead_se,

        "lead_effect_allele_frequency":
            lead_eaf,

        **top_values,

        "genome_wide_significant_variants":
            genome_wide_count,

        "suggestive_significant_variants":
            suggestive_count,

        "suggestive_only_variants":
            suggestive_only_count,

        "has_genome_wide_signal":
            bool(
                genome_wide_count
                > 0
            ),

        "has_suggestive_signal":
            bool(
                suggestive_count
                > 0
            ),

        "has_suggestive_only_signal":
            bool(
                suggestive_only_count
                > 0
            ),
    }