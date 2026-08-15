"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/whole_catalog_proxy_audit.py

Description:
    M5.3B whole-GWAS-Catalog proxy audit.

    This module audits LD proxy variants derived from M5 against the
    complete standardized GWAS Catalog association dataset, without
    restricting associations to the previously harmonized prostate
    cancer subset.

    Scientific objective:
        Determine whether M5 LD proxies are present among any GWAS Catalog
        association records and, when present, inspect the reported and
        mapped phenotypes associated with those variants.

    Important safeguards:
        - No GWAS significance filtering is performed.
        - No disease filtering is performed before rsID matching.
        - No GWAS association rows are deduplicated.
        - No LD is recalculated.
        - No coordinate matching is performed.
        - No colocalization is performed.
        - No causal inference is performed.
        - Prostate-related text annotations are exploratory audit flags
          and do not replace the harmonized disease classification used
          in M3.

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
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


RSID_PATTERN = re.compile(
    r"\brs\d+\b",
    flags=re.IGNORECASE,
)


# ============================================================================
# Helpers
# ============================================================================


def _first_existing_column(
    dataframe: pd.DataFrame,
    candidates: tuple[str, ...],
) -> str | None:
    """Resolve the first matching column using case-insensitive names."""

    lookup = {
        str(column).strip().lower():
            column
        for column
        in dataframe.columns
    }

    for candidate in candidates:

        resolved = lookup.get(
            candidate.strip().lower()
        )

        if resolved is not None:
            return resolved

    return None


def _normalize_rsid(
    value: Any,
) -> str | None:
    """Normalize one canonical rsID."""

    if value is None:
        return None

    text = str(
        value
    ).strip().lower()

    if re.fullmatch(
        r"rs\d+",
        text,
    ) is None:
        return None

    return text


def _extract_rsids(
    value: Any,
) -> list[str]:
    """
    Extract every canonical rsID from one GWAS variant field.

    This intentionally allows association rows containing multiple
    reported variants to contribute multiple rsID audit identities.
    """

    if value is None:
        return []

    try:

        if pd.isna(
            value
        ):
            return []

    except (
        TypeError,
        ValueError,
    ):
        pass

    matches = RSID_PATTERN.findall(
        str(
            value
        )
    )

    normalized = {
        match.lower()
        for match
        in matches
    }

    return sorted(
        normalized
    )


def _combine_text_columns(
    dataframe: pd.DataFrame,
    columns: Iterable[str],
) -> pd.Series:
    """Combine selected columns into one lower-case audit text."""

    available = [
        column
        for column
        in columns
        if column
        in dataframe.columns
    ]

    if not available:

        return pd.Series(
            "",
            index=dataframe.index,
            dtype="string",
        )

    text = (
        dataframe[
            available
        ]
        .fillna("")
        .astype("string")
        .agg(
            " | ".join,
            axis=1,
        )
        .str.lower()
    )

    return text


# ============================================================================
# Main audit
# ============================================================================


class M53BWholeCatalogProxyAudit:
    """
    Audit M5 LD proxies against all standardized GWAS Catalog associations.
    """

    GWAS_VARIANT_COLUMNS = (
        "canonical_rsid",
        "gwas_rsid",
        "rsid",
        "rs_id",
        "SNPS",
        "snps",
        "SNP",
        "snp",
        "variant",
        "variant_id",
    )

    REPORTED_TRAIT_COLUMNS = (
        "DISEASE/TRAIT",
        "disease_trait",
        "reported_trait",
        "trait",
    )

    MAPPED_TRAIT_COLUMNS = (
        "MAPPED_TRAIT",
        "mapped_trait",
    )

    MAPPED_TRAIT_URI_COLUMNS = (
        "MAPPED_TRAIT_URI",
        "mapped_trait_uri",
    )

    STUDY_COLUMNS = (
        "STUDY ACCESSION",
        "study_accession",
        "study_id",
    )

    PMID_COLUMNS = (
        "PUBMEDID",
        "pubmedid",
        "pmid",
    )

    PVALUE_COLUMNS = (
        "P-VALUE",
        "p_value",
        "pvalue",
    )

    PVALUE_MLOG_COLUMNS = (
        "PVALUE_MLOG",
        "pvalue_mlog",
        "p_value_mlog",
    )

    # ------------------------------------------------------------------
    # LD preparation
    # ------------------------------------------------------------------

    @classmethod
    def prepare_ld_proxy_set(
        cls,
        ld_evidence: pd.DataFrame,
        *,
        secondary_r2_threshold: float = 0.5,
    ) -> pd.DataFrame:
        """
        Prepare eligible LD proxy identities for whole-Catalog audit.
        """

        required = {
            "lead_rsid",
            "neighbor_rsid",
            "r2",
        }

        missing = (
            required
            - set(
                ld_evidence.columns
            )
        )

        if missing:

            raise ValueError(
                "LD evidence is missing required columns: "
                f"{sorted(missing)}"
            )

        working = (
            ld_evidence
            .copy()
        )

        working[
            "lead_rsid"
        ] = (
            working[
                "lead_rsid"
            ]
            .astype("string")
            .str.strip()
            .str.lower()
        )

        working[
            "neighbor_rsid"
        ] = (
            working[
                "neighbor_rsid"
            ]
            .astype("string")
            .str.strip()
            .str.lower()
        )

        working[
            "r2"
        ] = pd.to_numeric(
            working[
                "r2"
            ],
            errors="coerce",
        )

        usable = (
            working[
                "neighbor_rsid"
            ]
            .str.match(
                r"^rs\d+$",
                na=False,
            )
            &
            working[
                "r2"
            ]
            .ge(
                secondary_r2_threshold
            )
        )

        if (
            "is_self_pair"
            in working.columns
        ):

            usable &= (
                ~working[
                    "is_self_pair"
                ]
                .fillna(
                    False
                )
            )

        elif (
            "m5_bridge_eligible"
            in working.columns
        ):

            usable &= (
                working[
                    "m5_bridge_eligible"
                ]
                .fillna(
                    False
                )
            )

        proxies = (
            working.loc[
                usable
            ]
            .copy()
        )

        if proxies.empty:

            return proxies

        # Preserve each lead-proxy identity once for the audit.
        proxies = (
            proxies
            .sort_values(
                by=[
                    "lead_rsid",
                    "neighbor_rsid",
                    "r2",
                ],
                ascending=[
                    True,
                    True,
                    False,
                ],
                kind="stable",
            )
            .drop_duplicates(
                subset=[
                    "lead_rsid",
                    "neighbor_rsid",
                ],
                keep="first",
            )
            .reset_index(
                drop=True
            )
        )

        return proxies

    # ------------------------------------------------------------------
    # GWAS variant extraction
    # ------------------------------------------------------------------

    @classmethod
    def _resolve_variant_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str:
        """Resolve the source GWAS variant field."""

        column = _first_existing_column(
            dataframe,
            cls.GWAS_VARIANT_COLUMNS,
        )

        if column is None:

            raise ValueError(
                "Could not resolve a GWAS rsID/variant column. "
                f"Observed columns: {list(dataframe.columns)}"
            )

        return column

    @classmethod
    def _build_gwas_rsid_index(
        cls,
        dataframe: pd.DataFrame,
        *,
        source_file: str,
        proxy_rsids: set[str],
    ) -> pd.DataFrame:
        """
        Extract proxy-matching rsIDs from one standardized GWAS part.

        Only rows containing at least one requested proxy rsID are retained,
        which avoids exploding the complete Catalog unnecessarily.
        """

        variant_column = (
            cls._resolve_variant_column(
                dataframe
            )
        )

        working = (
            dataframe
            .copy()
        )

        working[
            "_m53b_source_file"
        ] = source_file

        working[
            "_m53b_source_row"
        ] = pd.Series(
            range(
                1,
                len(
                    working
                )
                + 1,
            ),
            index=working.index,
            dtype="Int64",
        )

        extracted = (
            working[
                variant_column
            ]
            .map(
                _extract_rsids
            )
        )

        contains_proxy = (
            extracted.map(
                lambda values:
                    bool(
                        proxy_rsids.intersection(
                            values
                        )
                    )
            )
        )

        if not contains_proxy.any():

            return pd.DataFrame()

        matched = (
            working.loc[
                contains_proxy
            ]
            .copy()
        )

        matched[
            "_m53b_all_extracted_rsids"
        ] = (
            extracted.loc[
                contains_proxy
            ]
        )

        matched[
            "gwas_audit_rsid"
        ] = (
            matched[
                "_m53b_all_extracted_rsids"
            ]
        )

        matched = (
            matched
            .explode(
                "gwas_audit_rsid",
                ignore_index=False,
            )
        )

        matched[
            "gwas_audit_rsid"
        ] = (
            matched[
                "gwas_audit_rsid"
            ]
            .astype("string")
            .str.lower()
        )

        matched = (
            matched.loc[
                matched[
                    "gwas_audit_rsid"
                ]
                .isin(
                    proxy_rsids
                )
            ]
            .copy()
        )

        matched.drop(
            columns=[
                "_m53b_all_extracted_rsids",
            ],
            inplace=True,
        )

        return matched.reset_index(
            drop=True
        )

    # ------------------------------------------------------------------
    # Phenotype annotation
    # ------------------------------------------------------------------

    @classmethod
    def _annotate_phenotypes(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Add exploratory phenotype audit flags.

        These flags are deliberately broad and must not replace the disease
        harmonization already used by M3.
        """

        result = (
            dataframe
            .copy()
        )

        reported_column = (
            _first_existing_column(
                result,
                cls.REPORTED_TRAIT_COLUMNS,
            )
        )

        mapped_column = (
            _first_existing_column(
                result,
                cls.MAPPED_TRAIT_COLUMNS,
            )
        )

        mapped_uri_column = (
            _first_existing_column(
                result,
                cls.MAPPED_TRAIT_URI_COLUMNS,
            )
        )

        study_column = (
            _first_existing_column(
                result,
                cls.STUDY_COLUMNS,
            )
        )

        pmid_column = (
            _first_existing_column(
                result,
                cls.PMID_COLUMNS,
            )
        )

        pvalue_column = (
            _first_existing_column(
                result,
                cls.PVALUE_COLUMNS,
            )
        )

        pvalue_mlog_column = (
            _first_existing_column(
                result,
                cls.PVALUE_MLOG_COLUMNS,
            )
        )

        result[
            "gwas_reported_trait"
        ] = (
            result[
                reported_column
            ].astype(
                "string"
            )
            if reported_column
            else pd.Series(
                pd.NA,
                index=result.index,
                dtype="string",
            )
        )

        result[
            "gwas_mapped_trait"
        ] = (
            result[
                mapped_column
            ].astype(
                "string"
            )
            if mapped_column
            else pd.Series(
                pd.NA,
                index=result.index,
                dtype="string",
            )
        )

        result[
            "gwas_mapped_trait_uri"
        ] = (
            result[
                mapped_uri_column
            ].astype(
                "string"
            )
            if mapped_uri_column
            else pd.Series(
                pd.NA,
                index=result.index,
                dtype="string",
            )
        )

        result[
            "gwas_study_accession"
        ] = (
            result[
                study_column
            ].astype(
                "string"
            )
            if study_column
            else pd.Series(
                pd.NA,
                index=result.index,
                dtype="string",
            )
        )

        result[
            "gwas_pmid"
        ] = (
            result[
                pmid_column
            ].astype(
                "string"
            )
            if pmid_column
            else pd.Series(
                pd.NA,
                index=result.index,
                dtype="string",
            )
        )

        result[
            "gwas_p_value"
        ] = (
            pd.to_numeric(
                result[
                    pvalue_column
                ],
                errors="coerce",
            )
            if pvalue_column
            else pd.Series(
                float("nan"),
                index=result.index,
                dtype="float64",
            )
        )

        result[
            "gwas_pvalue_mlog"
        ] = (
            pd.to_numeric(
                result[
                    pvalue_mlog_column
                ],
                errors="coerce",
            )
            if pvalue_mlog_column
            else pd.Series(
                float("nan"),
                index=result.index,
                dtype="float64",
            )
        )

        phenotype_text = (
            _combine_text_columns(
                result,
                (
                    "gwas_reported_trait",
                    "gwas_mapped_trait",
                ),
            )
        )

        result[
            "prostate_text_hit"
        ] = (
            phenotype_text
            .str.contains(
                r"\bprostate\b",
                regex=True,
                na=False,
            )
            .astype(
                "boolean"
            )
        )

        result[
            "prostate_cancer_text_hit"
        ] = (
            (
                phenotype_text
                .str.contains(
                    r"\bprostate\b",
                    regex=True,
                    na=False,
                )
            )
            &
            (
                phenotype_text
                .str.contains(
                    (
                        r"\bcancer\b"
                        r"|\bcarcinoma\b"
                        r"|\bneoplasm\b"
                        r"|\bmalignan"
                    ),
                    regex=True,
                    na=False,
                )
            )
        ).astype(
            "boolean"
        )

        return result

    # ------------------------------------------------------------------
    # Main scan
    # ------------------------------------------------------------------

    @classmethod
    def build(
        cls,
        *,
        ld_evidence: pd.DataFrame,
        gwas_paths: Iterable[Path],
        population: str,
        secondary_r2_threshold: float = 0.5,
    ) -> tuple[
        pd.DataFrame,
        dict[str, Any],
    ]:
        """
        Run whole-Catalog audit for one LD population.
        """

        proxies = (
            cls.prepare_ld_proxy_set(
                ld_evidence,
                secondary_r2_threshold=(
                    secondary_r2_threshold
                ),
            )
        )

        proxy_rsids = set(
            proxies[
                "neighbor_rsid"
            ]
            .dropna()
            .astype(str)
        )

        matched_parts: list[
            pd.DataFrame
        ] = []

        files_scanned = 0
        rows_scanned = 0

        for path in gwas_paths:

            path = Path(
                path
            )

            if not path.exists():

                raise FileNotFoundError(
                    f"GWAS Parquet file not found: {path}"
                )

            source = pd.read_parquet(
                path
            )

            files_scanned += 1
            rows_scanned += len(
                source
            )

            matched = (
                cls._build_gwas_rsid_index(
                    source,
                    source_file=path.name,
                    proxy_rsids=proxy_rsids,
                )
            )

            if not matched.empty:

                matched_parts.append(
                    matched
                )

            del source

        if matched_parts:

            catalog_matches = pd.concat(
                matched_parts,
                ignore_index=True,
                sort=False,
            )

            catalog_matches = (
                cls._annotate_phenotypes(
                    catalog_matches
                )
            )

            result = proxies.merge(
                catalog_matches,
                how="inner",
                left_on="neighbor_rsid",
                right_on="gwas_audit_rsid",
                validate="many_to_many",
            )

            result[
                "m5_ld_analysis_population"
            ] = (
                population
                .strip()
                .upper()
            )

            result[
                "m53b_match_type"
            ] = (
                "LD_PROXY_TO_WHOLE_GWAS_CATALOG"
            )

            result = (
                result
                .sort_values(
                    by=[
                        "lead_rsid",
                        "neighbor_rsid",
                        "r2",
                        "prostate_cancer_text_hit",
                    ],
                    ascending=[
                        True,
                        True,
                        False,
                        False,
                    ],
                    kind="stable",
                )
                .reset_index(
                    drop=True
                )
            )

            result[
                "m53b_audit_id"
            ] = [
                (
                    f"M53B-"
                    f"{population.upper()}-"
                    f"{index:08d}"
                )
                for index
                in range(
                    1,
                    len(
                        result
                    )
                    + 1,
                )
            ]

        else:

            result = pd.DataFrame()

        matched_proxy_rsids = (
            set(
                result[
                    "neighbor_rsid"
                ]
                .dropna()
                .astype(str)
            )
            if not result.empty
            else set()
        )

        prostate_matches = (
            result.loc[
                result[
                    "prostate_text_hit"
                ].fillna(
                    False
                )
            ]
            if (
                not result.empty
                and "prostate_text_hit"
                in result.columns
            )
            else pd.DataFrame()
        )

        prostate_cancer_matches = (
            result.loc[
                result[
                    "prostate_cancer_text_hit"
                ].fillna(
                    False
                )
            ]
            if (
                not result.empty
                and "prostate_cancer_text_hit"
                in result.columns
            )
            else pd.DataFrame()
        )

        report = {
            "milestone":
                "M5.3B",

            "stage":
                "whole_catalog_proxy_audit",

            "population":
                population.upper(),

            "policy": {
                "whole_catalog_scanned":
                    True,

                "disease_filter_before_matching":
                    False,

                "matching_key":
                    "canonical_rsid_extracted_from_gwas_variant_field",

                "secondary_r2_threshold":
                    secondary_r2_threshold,

                "physical_coordinate_matching":
                    False,

                "gwas_significance_filtering":
                    False,

                "association_deduplication":
                    False,

                "prostate_text_annotation_is_exploratory":
                    True,

                "m3_disease_harmonization_replaced":
                    False,

                "colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "gwas_files_scanned":
                    files_scanned,

                "gwas_rows_scanned":
                    rows_scanned,

                "eligible_proxy_rows":
                    len(
                        proxies
                    ),

                "unique_eligible_proxy_rsids":
                    len(
                        proxy_rsids
                    ),

                "whole_catalog_match_rows":
                    len(
                        result
                    ),

                "unique_proxy_rsids_found_in_catalog":
                    len(
                        matched_proxy_rsids
                    ),

                "unique_proxy_rsids_not_found_in_catalog":
                    len(
                        proxy_rsids
                        - matched_proxy_rsids
                    ),

                "prostate_text_match_rows":
                    len(
                        prostate_matches
                    ),

                "prostate_cancer_text_match_rows":
                    len(
                        prostate_cancer_matches
                    ),

                "unique_prostate_related_proxy_rsids":
                    (
                        int(
                            prostate_matches[
                                "neighbor_rsid"
                            ]
                            .nunique()
                        )
                        if not prostate_matches.empty
                        else 0
                    ),

                "unique_prostate_cancer_proxy_rsids":
                    (
                        int(
                            prostate_cancer_matches[
                                "neighbor_rsid"
                            ]
                            .nunique()
                        )
                        if not prostate_cancer_matches.empty
                        else 0
                    ),
            },

            "proxy_rsids_not_found_in_catalog":
                sorted(
                    proxy_rsids
                    - matched_proxy_rsids
                ),

            "matched_lead_rsids":
                (
                    sorted(
                        result[
                            "lead_rsid"
                        ]
                        .dropna()
                        .astype(str)
                        .unique()
                        .tolist()
                    )
                    if not result.empty
                    else []
                ),
        }

        return (
            result,
            report,
        )