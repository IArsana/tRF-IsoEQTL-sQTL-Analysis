"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/trfqtl_supplementary.py

Description:
    Standardization utilities for supplementary Cancer-tRFQTL workbook
    tables that are not part of the existing core S2/S10 standardizer.

    Supported supplementary sheets:

        S1  - survival-associated tRFQTLs
        S3  - cross-ancestry GWAS study metadata
        S4  - tRF immune-infiltration associations
        S5  - tRF drug-response associations
        S6  - Chinese cohort characteristics
        S7  - probe/primer sequence metadata
        S8  - UK Biobank cohort characteristics
        S9  - PLCO cohort characteristics
        S11 - POU2F1-regulated differentially expressed genes

    This module intentionally preserves source semantics rather than
    coercing all supplementary tables into one molecular-QTL schema.

    Standardization performs:

        - source-column preservation
        - source-row provenance
        - sheet and record-type annotation
        - conservative column-name normalization
        - known numeric-field conversion
        - genome-build annotation where applicable
        - common semantic aliases for downstream harmonization

    Standardization does not perform:

        - disease filtering
        - PRAD-only extraction
        - genome liftover
        - statistical significance filtering
        - identifier harmonization
        - candidate ranking
        - cross-source joins

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
import re
from typing import Any

import pandas as pd


DATASET_ID = "cancer_trfqtl"
SOURCE_STUDY = "Cancer-tRFQTL"
SOURCE_GENOME_BUILD = "hg19"


@dataclass(frozen=True)
class SupplementarySheetSpec:
    """Standardization specification for one Cancer-tRFQTL sheet."""

    sheet: str
    record_type: str
    category: str

    disease_column: str | None = None
    variant_column: str | None = None
    coordinate_column: str | None = None
    allele_column: str | None = None

    trf_column: str | None = None
    gene_column: str | None = None

    numeric_columns: tuple[str, ...] = ()
    p_value_columns: tuple[str, ...] = ()

    genome_build: str | None = None


SHEET_REGISTRY: dict[str, SupplementarySheetSpec] = {
    "S1": SupplementarySheetSpec(
        sheet="S1",
        record_type="survival_trfqtl",
        category="variant_feature",
        disease_column="Cancer type",
        variant_column="SNP ID",
        coordinate_column="SNP position (hg19)",
        allele_column="Alleles",
        trf_column="tRF",
        numeric_columns=(
            "Median survival time AA",
            "Median survival time Aa",
            "Median survival time aa",
            "P-value (survival)",
        ),
        p_value_columns=(
            "P-value (survival)",
        ),
        genome_build=SOURCE_GENOME_BUILD,
    ),
    "S3": SupplementarySheetSpec(
        sheet="S3",
        record_type="cross_ancestry_gwas_metadata",
        category="metadata",
        disease_column="Cancer type",
        numeric_columns=(
            "No. of cases in study",
            "No. of controls in study",
        ),
    ),
    "S4": SupplementarySheetSpec(
        sheet="S4",
        record_type="immune_infiltration",
        category="feature_annotation",
        disease_column="Cancer type",
        trf_column="tRF",
        numeric_columns=(
            "Correlation coefficient",
            "Correlation P-value",
        ),
        p_value_columns=(
            "Correlation P-value",
        ),
    ),
    "S5": SupplementarySheetSpec(
        sheet="S5",
        record_type="drug_response",
        category="feature_annotation",
        disease_column="Cancer type",
        trf_column="tRF",
        numeric_columns=(
            "Correlation coefficient",
            "Correlation P-value",
        ),
        p_value_columns=(
            "Correlation P-value",
        ),
    ),
    "S6": SupplementarySheetSpec(
        sheet="S6",
        record_type="chinese_cohort_metadata",
        category="metadata",
    ),
    "S7": SupplementarySheetSpec(
        sheet="S7",
        record_type="experimental_sequence",
        category="metadata",
    ),
    "S8": SupplementarySheetSpec(
        sheet="S8",
        record_type="uk_biobank_metadata",
        category="metadata",
    ),
    "S9": SupplementarySheetSpec(
        sheet="S9",
        record_type="plco_cohort_metadata",
        category="metadata",
    ),
    "S11": SupplementarySheetSpec(
        sheet="S11",
        record_type="pou2f1_de_gene",
        category="gene_annotation",
        gene_column="Gene ID",
        numeric_columns=(
            "log2 (Fold Change)",
            "P-value",
        ),
        p_value_columns=(
            "P-value",
        ),
    ),
}


class CancerTRFQTLSupplementaryStandardizer:
    """
    Standardize non-core Cancer-tRFQTL supplementary tables.

    S2 and S10 are intentionally excluded because they already have an
    existing dedicated standardizer and downstream harmonized outputs.
    """

    SUPPORTED_SHEETS = frozenset(
        SHEET_REGISTRY
    )

    CORE_LOCKED_SHEETS = frozenset(
        {
            "S2",
            "S10",
        }
    )

    # ------------------------------------------------------------------
    # Cleaning helpers
    # ------------------------------------------------------------------

    @staticmethod
    def clean_source_dataframe(
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Remove completely empty source rows/columns.

        Non-empty source values are never corrected or rewritten.
        """

        cleaned = dataframe.dropna(
            axis=0,
            how="all",
        )

        cleaned = cleaned.dropna(
            axis=1,
            how="all",
        )

        cleaned = cleaned.copy()

        cleaned.columns = [
            str(column).strip()
            for column in cleaned.columns
        ]

        return cleaned.reset_index(
            drop=True
        )

    @staticmethod
    def normalize_column_name(
        column: Any,
    ) -> str:
        """
        Convert a source column name to conservative snake_case.

        The original source column remains present separately using the
        ``source_*`` representation policy.
        """

        text = str(
            column
        ).strip()

        text = text.replace(
            "²",
            "2",
        )

        text = re.sub(
            r"[^\w]+",
            "_",
            text,
            flags=re.UNICODE,
        )

        text = re.sub(
            r"_+",
            "_",
            text,
        )

        return text.strip(
            "_"
        ).lower()

    @staticmethod
    def _numeric(
        series: pd.Series,
    ) -> pd.Series:
        """Conservatively convert a semantically numeric source column."""

        return pd.to_numeric(
            series,
            errors="coerce",
        )

    @staticmethod
    def _text(
        series: pd.Series,
    ) -> pd.Series:
        """Return nullable textual representation."""

        return series.astype(
            "string"
        )

    @staticmethod
    def _existing_column(
        dataframe: pd.DataFrame,
        column: str | None,
    ) -> pd.Series | None:
        """Return a source column only when explicitly defined and present."""

        if (
            column is None
            or column not in dataframe.columns
        ):
            return None

        return dataframe[
            column
        ]

    # ------------------------------------------------------------------
    # Schema verification
    # ------------------------------------------------------------------

    @classmethod
    def require_supported_sheet(
        cls,
        sheet: str,
    ) -> None:
        """Require a supplementary sheet handled by this standardizer."""

        if sheet in cls.CORE_LOCKED_SHEETS:
            raise ValueError(
                f"{sheet} is a locked core Cancer-tRFQTL table and "
                "must be handled by the existing S2/S10 standardizer."
            )

        if sheet not in cls.SUPPORTED_SHEETS:
            raise ValueError(
                "Unsupported Cancer-tRFQTL supplementary sheet: "
                f"{sheet}. Supported sheets: "
                f"{sorted(cls.SUPPORTED_SHEETS)}"
            )

    @staticmethod
    def _required_semantic_columns(
        spec: SupplementarySheetSpec,
    ) -> set[str]:
        """Return semantic source fields required for a sheet."""

        required: set[str] = set()

        for column in (
            spec.disease_column,
            spec.variant_column,
            spec.coordinate_column,
            spec.allele_column,
            spec.trf_column,
            spec.gene_column,
        ):
            if column is not None:
                required.add(
                    column
                )

        required.update(
            spec.numeric_columns
        )

        return required

    @classmethod
    def validate_schema(
        cls,
        dataframe: pd.DataFrame,
        sheet: str,
    ) -> None:
        """
        Validate known semantic fields.

        Metadata-only sheets without a documented fixed schema are
        intentionally accepted after basic non-empty validation.
        """

        cls.require_supported_sheet(
            sheet
        )

        spec = SHEET_REGISTRY[
            sheet
        ]

        if dataframe.empty:
            raise ValueError(
                f"Cancer-tRFQTL {sheet} is empty after source cleaning."
            )

        required = (
            cls._required_semantic_columns(
                spec
            )
        )

        missing = (
            required
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                f"Cancer-tRFQTL {sheet} is missing expected columns: "
                f"{sorted(missing)}"
            )

    # ------------------------------------------------------------------
    # Generic provenance
    # ------------------------------------------------------------------

    @classmethod
    def _base_standardized_dataframe(
        cls,
        dataframe: pd.DataFrame,
        *,
        sheet: str,
    ) -> pd.DataFrame:
        """Create provenance-preserving base standardized table."""

        spec = SHEET_REGISTRY[
            sheet
        ]

        result = pd.DataFrame(
            index=dataframe.index
        )

        result[
            "dataset_id"
        ] = DATASET_ID

        result[
            "source_study"
        ] = SOURCE_STUDY

        result[
            "source_sheet"
        ] = sheet

        # Excel actual data begin after:
        # row 1 = descriptive title
        # row 2 = column header
        #
        # Therefore first dataframe record corresponds to Excel row 3.
        result[
            "source_row"
        ] = (
            dataframe.index
            + 3
        ).astype(
            "int64"
        )

        result[
            "record_type"
        ] = spec.record_type

        result[
            "record_category"
        ] = spec.category

        result[
            "source_genome_build"
        ] = (
            spec.genome_build
            if spec.genome_build
            is not None
            else pd.NA
        )

        result[
            "standardization_status"
        ] = "STANDARDIZED"

        return result

    @classmethod
    def _add_source_columns(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """
        Preserve every source field.

        Source fields are prefixed to guarantee they cannot silently
        collide with standardized semantic fields.
        """

        used_names: set[str] = set(
            result.columns
        )

        for original_column in source.columns:

            normalized = (
                cls.normalize_column_name(
                    original_column
                )
            )

            base_name = (
                f"source_{normalized}"
                if normalized
                else "source_column"
            )

            target = base_name
            suffix = 2

            while target in used_names:
                target = (
                    f"{base_name}_{suffix}"
                )

                suffix += 1

            used_names.add(
                target
            )

            result[
                target
            ] = source[
                original_column
            ].copy()

    # ------------------------------------------------------------------
    # Common semantic aliases
    # ------------------------------------------------------------------

    @classmethod
    def _add_common_semantics(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
        *,
        spec: SupplementarySheetSpec,
    ) -> None:
        """Add semantic fields shared by downstream harmonization."""

        disease = cls._existing_column(
            source,
            spec.disease_column,
        )

        if disease is not None:
            result[
                "disease_raw"
            ] = cls._text(
                disease
            )

        variant = cls._existing_column(
            source,
            spec.variant_column,
        )

        if variant is not None:
            result[
                "variant_raw"
            ] = cls._text(
                variant
            )

        coordinate = cls._existing_column(
            source,
            spec.coordinate_column,
        )

        if coordinate is not None:
            result[
                "coordinate_raw"
            ] = cls._text(
                coordinate
            )

        allele = cls._existing_column(
            source,
            spec.allele_column,
        )

        if allele is not None:
            result[
                "alleles_raw"
            ] = cls._text(
                allele
            )

        trf = cls._existing_column(
            source,
            spec.trf_column,
        )

        if trf is not None:
            result[
                "feature_raw"
            ] = cls._text(
                trf
            )

            result[
                "feature_type"
            ] = "TRF"

            result[
                "trf_raw"
            ] = cls._text(
                trf
            )

        gene = cls._existing_column(
            source,
            spec.gene_column,
        )

        if gene is not None:
            result[
                "feature_raw"
            ] = cls._text(
                gene
            )

            result[
                "feature_type"
            ] = "GENE"

            result[
                "gene_raw"
            ] = cls._text(
                gene
            )

    # ------------------------------------------------------------------
    # Sheet-specific semantics
    # ------------------------------------------------------------------

    @classmethod
    def _standardize_s1(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Standardize survival-associated tRFQTL fields."""

        mappings = {
            "Median survival time AA":
                "median_survival_aa_homozygous_major",
            "Median survival time Aa":
                "median_survival_heterozygous",
            "Median survival time aa":
                "median_survival_aa_homozygous_minor",
            "P-value (survival)":
                "survival_p_value",
        }

        for (
            source_column,
            target_column,
        ) in mappings.items():

            result[
                target_column
            ] = cls._numeric(
                source[
                    source_column
                ]
            )

    @classmethod
    def _standardize_s3(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Standardize cross-ancestry GWAS metadata."""

        text_mappings = {
            "Sub-study name":
                "substudy_name",
            "Individual or summary level data":
                "data_level",
            "Design, location":
                "design_location",
            "Study reference":
                "study_reference",
        }

        numeric_mappings = {
            "No. of cases in study":
                "number_of_cases",
            "No. of controls in study":
                "number_of_controls",
        }

        for (
            source_column,
            target_column,
        ) in text_mappings.items():

            if (
                source_column
                in source.columns
            ):
                result[
                    target_column
                ] = cls._text(
                    source[
                        source_column
                    ]
                )

        for (
            source_column,
            target_column,
        ) in numeric_mappings.items():

            if (
                source_column
                in source.columns
            ):
                result[
                    target_column
                ] = cls._numeric(
                    source[
                        source_column
                    ]
                )

    @classmethod
    def _standardize_s4(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Standardize immune-infiltration association fields."""

        result[
            "cell_type_raw"
        ] = cls._text(
            source[
                "Cell type"
            ]
        )

        result[
            "correlation_coefficient"
        ] = cls._numeric(
            source[
                "Correlation coefficient"
            ]
        )

        result[
            "correlation_p_value"
        ] = cls._numeric(
            source[
                "Correlation P-value"
            ]
        )

        result[
            "method_raw"
        ] = cls._text(
            source[
                "Method"
            ]
        )

        result[
            "effect_type"
        ] = "correlation"

        result[
            "effect_value"
        ] = result[
            "correlation_coefficient"
        ]

        result[
            "p_value"
        ] = result[
            "correlation_p_value"
        ]

    @classmethod
    def _standardize_s5(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Standardize drug-response association fields."""

        result[
            "drug_raw"
        ] = cls._text(
            source[
                "Drug"
            ]
        )

        result[
            "correlation_coefficient"
        ] = cls._numeric(
            source[
                "Correlation coefficient"
            ]
        )

        result[
            "correlation_p_value"
        ] = cls._numeric(
            source[
                "Correlation P-value"
            ]
        )

        result[
            "effect_type"
        ] = "correlation"

        result[
            "effect_value"
        ] = result[
            "correlation_coefficient"
        ]

        result[
            "p_value"
        ] = result[
            "correlation_p_value"
        ]

    @classmethod
    def _standardize_metadata_sheet(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """
        Preserve metadata-only sheets without imposing invented semantics.

        S6, S8, and S9 may contain complex cohort-level layouts, so their
        source columns are preserved while only generic provenance fields
        are added.
        """

        del result
        del source

    @classmethod
    def _standardize_s7(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Standardize experimental probe/primer metadata."""

        preferred_column = (
            "qRT-PCR (5'-3')"
        )

        if (
            preferred_column
            in source.columns
        ):
            result[
                "experimental_sequence_raw"
            ] = cls._text(
                source[
                    preferred_column
                ]
            )

    @classmethod
    def _standardize_s11(
        cls,
        result: pd.DataFrame,
        source: pd.DataFrame,
    ) -> None:
        """Standardize POU2F1-regulated differential-expression fields."""

        result[
            "log2_fold_change"
        ] = cls._numeric(
            source[
                "log2 (Fold Change)"
            ]
        )

        result[
            "p_value"
        ] = cls._numeric(
            source[
                "P-value"
            ]
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def standardize(
        cls,
        dataframe: pd.DataFrame,
        *,
        sheet: str,
    ) -> pd.DataFrame:
        """
        Standardize one supported supplementary sheet.

        Row cardinality is preserved after the initial removal of wholly
        empty spreadsheet rows.
        """

        cls.require_supported_sheet(
            sheet
        )

        source = (
            cls.clean_source_dataframe(
                dataframe
            )
        )

        cls.validate_schema(
            source,
            sheet,
        )

        spec = SHEET_REGISTRY[
            sheet
        ]

        result = (
            cls._base_standardized_dataframe(
                source,
                sheet=sheet,
            )
        )

        cls._add_common_semantics(
            result,
            source,
            spec=spec,
        )

        if sheet == "S1":
            cls._standardize_s1(
                result,
                source,
            )

        elif sheet == "S3":
            cls._standardize_s3(
                result,
                source,
            )

        elif sheet == "S4":
            cls._standardize_s4(
                result,
                source,
            )

        elif sheet == "S5":
            cls._standardize_s5(
                result,
                source,
            )

        elif sheet in {
            "S6",
            "S8",
            "S9",
        }:
            cls._standardize_metadata_sheet(
                result,
                source,
            )

        elif sheet == "S7":
            cls._standardize_s7(
                result,
                source,
            )

        elif sheet == "S11":
            cls._standardize_s11(
                result,
                source,
            )

        cls._add_source_columns(
            result,
            source,
        )

        if (
            len(
                result
            )
            != len(
                source
            )
        ):
            raise RuntimeError(
                f"Cancer-tRFQTL {sheet} standardization changed "
                f"cardinality: {len(source)} -> {len(result)}"
            )

        return result