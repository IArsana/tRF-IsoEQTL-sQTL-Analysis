"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/trfqtl.py

Description:
    Source-specific standardization utilities for Cancer-tRFQTL
    supplementary datasets.

    Supported core source structures:

        S2:
            GWAS-associated tRFQTL associations.

        S10:
            GWAS effect records.

    Standardization responsibilities:

        - canonical dbSNP representation
        - hg19 genomic-coordinate parsing
        - chromosome and position standardization
        - allele parsing
        - tRF identifier preservation
        - cancer-type preservation
        - statistical-field conversion
        - LD conversion
        - usability flags
        - source provenance

    Important:
        Genome build is preserved as hg19.

        No liftover is performed in this module.

        Raw source values are never modified.

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
from typing import Any

import pandas as pd

from pcatrfqtl.standardization.common import (
    StandardizationUtils,
)


class TRFQTLStandardizer:
    """Standardize validated Cancer-tRFQTL records."""

    GENOME_BUILD = "hg19"

    COORDINATE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?P<chromosome>"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
        r")"
        r":"
        r"(?P<position>\d+)$",
        re.IGNORECASE,
    )

    ALLELE_PAIR_PATTERN = re.compile(
        r"^(?P<reference>[ACGT]+)"
        r"/"
        r"(?P<alternate>[ACGT]+)$",
        re.IGNORECASE,
    )

    SINGLE_ALLELE_PATTERN = re.compile(
        r"^[ACGT]+$",
        re.IGNORECASE,
    )

    # ==================================================================
    # rsID
    # ==================================================================

    @classmethod
    def standardize_rsid(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Standardize a Cancer-tRFQTL SNP identifier.

        Non-rsID source values are preserved but are not considered
        usable for rsID-dependent downstream analysis.
        """

        result = (
            StandardizationUtils
            .standardize_rsid(
                value
            )
        )

        return {
            "canonical_rsid":
                result.standardized_value,

            "rsid_usable":
                result.usable,

            "rsid_source":
                result.source,
        }

    # ==================================================================
    # Coordinate
    # ==================================================================

    @classmethod
    def standardize_coordinate(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Standardize a Cancer-tRFQTL hg19 coordinate.

        Expected source format:
            chromosome:position

        Examples:
            6:28958399
            chr6:28958399
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "coordinate_raw":
                    value,
                "chromosome":
                    None,
                "position":
                    None,
                "coordinate_usable":
                    False,
                "coordinate_source":
                    "missing",
                "genome_build":
                    cls.GENOME_BUILD,
            }

        text = str(
            value
        ).strip()

        match = (
            cls.COORDINATE_PATTERN
            .fullmatch(
                text
            )
        )

        if match is None:
            return {
                "coordinate_raw":
                    value,
                "chromosome":
                    None,
                "position":
                    None,
                "coordinate_usable":
                    False,
                "coordinate_source":
                    "unsupported_coordinate",
                "genome_build":
                    cls.GENOME_BUILD,
            }

        chromosome_result = (
            StandardizationUtils
            .standardize_chromosome(
                match.group(
                    "chromosome"
                )
            )
        )

        position_result = (
            StandardizationUtils
            .standardize_position(
                match.group(
                    "position"
                )
            )
        )

        usable = (
            chromosome_result.usable
            and position_result.usable
        )

        return {
            "coordinate_raw":
                value,

            "chromosome":
                chromosome_result
                .standardized_value,

            "position":
                position_result
                .standardized_value,

            "coordinate_usable":
                bool(
                    usable
                ),

            "coordinate_source":
                (
                    "parsed_hg19_coordinate"
                    if usable
                    else "unsupported_coordinate"
                ),

            "genome_build":
                cls.GENOME_BUILD,
        }

    # ==================================================================
    # Alleles
    # ==================================================================

    @classmethod
    def standardize_alleles(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Standardize the Cancer-tRFQTL allele field.

        Supported representations:
            A
            A/T
            AT/GC

        A single allele is preserved but cannot be interpreted as a
        reference/alternate pair.
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "alleles_raw":
                    value,
                "reference_allele":
                    None,
                "alternate_allele":
                    None,
                "allele_value":
                    None,
                "allele_pair_usable":
                    False,
                "allele_source":
                    "missing",
            }

        text = (
            str(
                value
            )
            .strip()
            .upper()
        )

        pair_match = (
            cls.ALLELE_PAIR_PATTERN
            .fullmatch(
                text
            )
        )

        if pair_match is not None:
            return {
                "alleles_raw":
                    value,

                "reference_allele":
                    pair_match.group(
                        "reference"
                    ),

                "alternate_allele":
                    pair_match.group(
                        "alternate"
                    ),

                "allele_value":
                    None,

                "allele_pair_usable":
                    True,

                "allele_source":
                    "reference_alternate_pair",
            }

        if (
            cls.SINGLE_ALLELE_PATTERN
            .fullmatch(
                text
            )
            is not None
        ):
            return {
                "alleles_raw":
                    value,

                "reference_allele":
                    None,

                "alternate_allele":
                    None,

                "allele_value":
                    text,

                "allele_pair_usable":
                    False,

                "allele_source":
                    "single_allele",
            }

        return {
            "alleles_raw":
                value,

            "reference_allele":
                None,

            "alternate_allele":
                None,

            "allele_value":
                None,

            "allele_pair_usable":
                False,

            "allele_source":
                "unsupported",
        }

    # ==================================================================
    # Cancer type
    # ==================================================================

    @classmethod
    def standardize_cancer_type(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """Standardize cancer-type source text without reclassification."""

        result = (
            StandardizationUtils
            .standardize_string(
                value
            )
        )

        return {
            "cancer_type":
                result.standardized_value,

            "cancer_type_usable":
                result.usable,

            "cancer_type_source":
                result.source,
        }

    # ==================================================================
    # tRF identifier
    # ==================================================================

    @classmethod
    def standardize_trf(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Preserve the source tRF identifier.

        No attempt is made here to translate between different tRF
        nomenclature systems.
        """

        result = (
            StandardizationUtils
            .standardize_string(
                value
            )
        )

        return {
            "trf_id":
                result.standardized_value,

            "trf_usable":
                result.usable,

            "trf_source":
                (
                    "source_identifier"
                    if result.usable
                    else result.source
                ),
        }

    # ==================================================================
    # Numeric fields
    # ==================================================================

    @classmethod
    def standardize_numeric_field(
        cls,
        value: Any,
    ) -> Any:
        """Return a finite numeric value or None."""

        result = (
            StandardizationUtils
            .standardize_numeric(
                value,
                allow_missing=True,
            )
        )

        return (
            result.standardized_value
            if result.usable
            else None
        )

    @classmethod
    def standardize_probability_field(
        cls,
        value: Any,
    ) -> Any:
        """Return a standardized probability value or None."""

        result = (
            StandardizationUtils
            .standardize_probability(
                value,
                allow_zero=True,
                allow_missing=True,
            )
        )

        return (
            result.standardized_value
            if result.usable
            else None
        )

    # ==================================================================
    # S2
    # ==================================================================

    @classmethod
    def standardize_s2_record(
        cls,
        record: dict[str, Any],
    ) -> dict[str, Any]:
        """
        Standardize one Cancer-tRFQTL S2 association record.
        """

        cancer = (
            cls
            .standardize_cancer_type(
                record.get(
                    "Cancer type"
                )
            )
        )

        rsid = (
            cls
            .standardize_rsid(
                record.get(
                    "SNP ID"
                )
            )
        )

        coordinate = (
            cls
            .standardize_coordinate(
                record.get(
                    "SNP position (hg19)"
                )
            )
        )

        alleles = (
            cls
            .standardize_alleles(
                record.get(
                    "Alleles"
                )
            )
        )

        trf = (
            cls
            .standardize_trf(
                record.get(
                    "tRF"
                )
            )
        )

        gwas_tag = (
            cls
            .standardize_rsid(
                record.get(
                    "GWAS tagSnp"
                )
            )
        )

        statistic = (
            cls
            .standardize_numeric_field(
                record.get(
                    "statistic"
                )
            )
        )

        p_value = (
            cls
            .standardize_probability_field(
                record.get(
                    "P-value"
                )
            )
        )

        fdr = (
            cls
            .standardize_probability_field(
                record.get(
                    "FDR"
                )
            )
        )

        ld_r2 = (
            cls
            .standardize_probability_field(
                record.get(
                    "LD (r2)"
                )
            )
        )

        variant_usable = (
            rsid[
                "rsid_usable"
            ]
            or coordinate[
                "coordinate_usable"
            ]
        )

        qtl_usable = (
            bool(
                variant_usable
            )
            and trf[
                "trf_usable"
            ]
            and p_value is not None
        )

        if qtl_usable:
            standardization_status = (
                "STANDARDIZED"
            )

        elif (
            variant_usable
            or trf[
                "trf_usable"
            ]
        ):
            standardization_status = (
                "PARTIAL"
            )

        else:
            standardization_status = (
                "PRESERVED"
            )

        return {
            "source_table":
                "S2",

            **cancer,

            "snp_id_raw":
                record.get(
                    "SNP ID"
                ),

            **rsid,

            **coordinate,

            **alleles,

            "trf_raw":
                record.get(
                    "tRF"
                ),

            **trf,

            "statistic":
                statistic,

            "p_value":
                p_value,

            "fdr":
                fdr,

            "gwas_tag_snp_raw":
                record.get(
                    "GWAS tagSnp"
                ),

            "gwas_tag_rsid":
                gwas_tag[
                    "canonical_rsid"
                ],

            "gwas_tag_rsid_usable":
                gwas_tag[
                    "rsid_usable"
                ],

            "ld_r2":
                ld_r2,

            "variant_usable":
                bool(
                    variant_usable
                ),

            "qtl_usable":
                bool(
                    qtl_usable
                ),

            "standardization_status":
                standardization_status,
        }

    # ==================================================================
    # S10
    # ==================================================================

    @classmethod
    def standardize_s10_record(
        cls,
        record: dict[str, Any],
    ) -> dict[str, Any]:
        """
        Standardize one Cancer-tRFQTL S10 GWAS-effect record.

        Known non-rsID source values, such as "6", remain preserved and
        are not converted into artificial rsIDs.
        """

        cancer = (
            cls
            .standardize_cancer_type(
                record.get(
                    "Cancer type"
                )
            )
        )

        rsid = (
            cls
            .standardize_rsid(
                record.get(
                    "SNP ID"
                )
            )
        )

        coordinate = (
            cls
            .standardize_coordinate(
                record.get(
                    "Position (hg19)"
                )
            )
        )

        effect_allele_result = (
            StandardizationUtils
            .standardize_string(
                record.get(
                    "A1 (effect allele)"
                )
            )
        )

        effect_allele = (
            effect_allele_result
            .standardized_value
        )

        if (
            effect_allele
            is not None
        ):
            effect_allele = (
                effect_allele.upper()
            )

        effect = (
            cls
            .standardize_numeric_field(
                record.get(
                    "Effect"
                )
            )
        )

        gwas_p_value = (
            cls
            .standardize_probability_field(
                record.get(
                    "GWAS P-value"
                )
            )
        )

        variant_usable = (
            rsid[
                "rsid_usable"
            ]
            or coordinate[
                "coordinate_usable"
            ]
        )

        effect_usable = (
            bool(
                variant_usable
            )
            and effect is not None
            and gwas_p_value is not None
        )

        if effect_usable:
            standardization_status = (
                "STANDARDIZED"
            )

        elif variant_usable:
            standardization_status = (
                "PARTIAL"
            )

        else:
            standardization_status = (
                "PRESERVED"
            )

        return {
            "source_table":
                "S10",

            **cancer,

            "snp_id_raw":
                record.get(
                    "SNP ID"
                ),

            **rsid,

            **coordinate,

            "effect_allele_raw":
                record.get(
                    "A1 (effect allele)"
                ),

            "effect_allele":
                effect_allele,

            "effect":
                effect,

            "gwas_p_value":
                gwas_p_value,

            "variant_usable":
                bool(
                    variant_usable
                ),

            "effect_usable":
                bool(
                    effect_usable
                ),

            "standardization_status":
                standardization_status,
        }

    # ==================================================================
    # DataFrame
    # ==================================================================

    @classmethod
    def standardize_s2_dataframe(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Standardize a validated S2 DataFrame."""

        records = (
            dataframe
            .to_dict(
                orient="records"
            )
        )

        return pd.DataFrame(
            [
                cls
                .standardize_s2_record(
                    record
                )
                for record
                in records
            ]
        )

    @classmethod
    def standardize_s10_dataframe(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Standardize a validated S10 DataFrame."""

        records = (
            dataframe
            .to_dict(
                orient="records"
            )
        )

        return pd.DataFrame(
            [
                cls
                .standardize_s10_record(
                    record
                )
                for record
                in records
            ]
        )