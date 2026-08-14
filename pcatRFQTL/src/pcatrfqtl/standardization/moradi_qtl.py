"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/moradi_qtl.py

Description:
    Source-specific standardization utilities for Moradi QTL
    supplementary datasets S1-S4.

    Supported source structures:

        S1:
            Significant cis-sQTL and cis-iso-eQTL records.

        S2:
            Cis-QTL records linked to GWAS variants through LD.

        S3:
            Trans-sQTL and trans-iso-eQTL records.

        S4:
            Trans-QTL records linked to GWAS variants through LD.

    Responsibilities:

        - canonical dbSNP rsID standardization
        - composite rsID decomposition
        - genomic-coordinate parsing
        - allele standardization
        - splicing-event identifier standardization
        - Ensembl transcript identifier standardization
        - numeric-field standardization
        - probability-field standardization
        - QTL usability flags
        - explicit provenance

    Important source conventions:

        S2 and S4 contain semantically inverted column names:

            sQTL_SNP
                genomic coordinate

            sQTL_SNP-pos
                rsID

            tag_SNP
                genomic coordinate

            tag_SNP_pos
                single or composite rsID

        S3 trans-exon uses source column "T" rather than "SNP".

        Known malformed source identifiers such as "INT1e+05" are
        preserved but are not converted to canonical INT identifiers.

        The Moradi genome build is not assigned in this module because
        it has not yet been independently verified.

    This module performs no file I/O and never modifies raw data.

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
from typing import Any, Literal

import pandas as pd

from pcatrfqtl.standardization.common import (
    StandardizationUtils,
)


FeatureType = Literal[
    "exon",
    "intron",
    "transcript",
]


class MoradiQTLStandardizer:
    """Standardize validated Moradi QTL records from S1-S4."""

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

    COMPOSITE_RSID_PATTERN = re.compile(
        r"^rs\d+(?::rs\d+)*$",
        re.IGNORECASE,
    )

    EXON_PATTERN = re.compile(
        r"^EX\d+$",
        re.IGNORECASE,
    )

    INTRON_PATTERN = re.compile(
        r"^INT\d+$",
        re.IGNORECASE,
    )

    TRANSCRIPT_PATTERN = re.compile(
        r"^ENST\d+(?:\.\d+)?$",
        re.IGNORECASE,
    )

    ALLELE_PATTERN = re.compile(
        r"^[ACGT]+(?:,[ACGT]+)*$",
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
        """Standardize one canonical dbSNP rsID."""

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
                bool(
                    result.usable
                ),

            "rsid_source":
                result.source,
        }

    @classmethod
    def standardize_composite_rsid(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Standardize a single or colon-separated composite rsID field.

        Examples:
            rs123
            rs123:rs456
            rs123:rs456:rs789
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "rsids":
                    None,

                "primary_rsid":
                    None,

                "composite_rsid":
                    False,

                "rsid_usable":
                    False,

                "rsid_source":
                    "missing",
            }

        text = (
            str(
                value
            )
            .strip()
        )

        if (
            cls.COMPOSITE_RSID_PATTERN
            .fullmatch(
                text
            )
            is None
        ):
            return {
                "rsids":
                    None,

                "primary_rsid":
                    None,

                "composite_rsid":
                    False,

                "rsid_usable":
                    False,

                "rsid_source":
                    "non_canonical_composite_rsid",
            }

        parts = (
            text.split(
                ":"
            )
        )

        standardized: list[str] = []

        for part in parts:

            result = (
                StandardizationUtils
                .standardize_rsid(
                    part
                )
            )

            if not result.usable:
                return {
                    "rsids":
                        None,

                    "primary_rsid":
                        None,

                    "composite_rsid":
                        False,

                    "rsid_usable":
                        False,

                    "rsid_source":
                        "non_canonical_composite_rsid",
                }

            standardized.append(
                result.standardized_value
            )

        return {
            "rsids":
                standardized,

            "primary_rsid":
                (
                    standardized[0]
                    if len(
                        standardized
                    )
                    == 1
                    else None
                ),

            "composite_rsid":
                (
                    len(
                        standardized
                    )
                    > 1
                ),

            "rsid_usable":
                True,

            "rsid_source":
                (
                    "composite_rsid"
                    if len(
                        standardized
                    )
                    > 1
                    else "single_rsid"
                ),
        }

    # ==================================================================
    # Coordinates
    # ==================================================================

    @classmethod
    def standardize_coordinate(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """Standardize chromosome:position representations."""

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

        text = (
            str(
                value
            )
            .strip()
        )

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

        chromosome = (
            StandardizationUtils
            .standardize_chromosome(
                match.group(
                    "chromosome"
                )
            )
        )

        position = (
            StandardizationUtils
            .standardize_position(
                match.group(
                    "position"
                )
            )
        )

        usable = (
            chromosome.usable
            and position.usable
        )

        return {
            "coordinate_raw":
                value,

            "chromosome":
                chromosome.standardized_value,

            "position":
                position.standardized_value,

            "coordinate_usable":
                bool(
                    usable
                ),

            "coordinate_source":
                (
                    "parsed_coordinate"
                    if usable
                    else "unsupported_coordinate"
                ),

            "genome_build":
                cls.GENOME_BUILD,
        }

    # ==================================================================
    # Feature identifiers
    # ==================================================================

    @classmethod
    def standardize_feature(
        cls,
        value: Any,
        feature_type: FeatureType,
    ) -> dict[str, Any]:
        """
        Standardize exon, intron, or transcript identifiers.

        Canonical examples:
            EX123
            INT123
            ENST00000318325

        Malformed identifiers such as INT1e+05 remain unconverted.
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "feature_id":
                    None,

                "feature_type":
                    feature_type,

                "feature_usable":
                    False,

                "feature_source":
                    "missing",
            }

        text = (
            str(
                value
            )
            .strip()
        )

        if feature_type == "exon":
            pattern = (
                cls.EXON_PATTERN
            )

        elif feature_type == "intron":
            pattern = (
                cls.INTRON_PATTERN
            )

        elif feature_type == "transcript":
            pattern = (
                cls.TRANSCRIPT_PATTERN
            )

        else:
            raise ValueError(
                "Unsupported Moradi feature type: "
                f"{feature_type!r}"
            )

        if (
            pattern.fullmatch(
                text
            )
            is None
        ):
            return {
                "feature_id":
                    None,

                "feature_type":
                    feature_type,

                "feature_usable":
                    False,

                "feature_source":
                    "non_canonical_feature_id",
            }

        if (
            feature_type
            == "transcript"
        ):
            standardized = (
                text.upper()
            )

        else:
            prefix = (
                "EX"
                if feature_type
                == "exon"
                else "INT"
            )

            digit_part = (
                text[
                    len(
                        prefix
                    ):
                ]
            )

            standardized = (
                prefix
                + digit_part
            )

        return {
            "feature_id":
                standardized,

            "feature_type":
                feature_type,

            "feature_usable":
                True,

            "feature_source":
                "canonical_source_identifier",
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
        Standardize Moradi allele representations.

        Supports:
            A
            G
            A,G
            A,C,G
            AT,G
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

                "alleles":
                    None,

                "allele_count":
                    0,

                "alleles_usable":
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

        if (
            cls.ALLELE_PATTERN
            .fullmatch(
                text
            )
            is None
        ):
            return {
                "alleles_raw":
                    value,

                "alleles":
                    None,

                "allele_count":
                    0,

                "alleles_usable":
                    False,

                "allele_source":
                    "unsupported",
            }

        alleles = (
            text.split(
                ","
            )
        )

        return {
            "alleles_raw":
                value,

            "alleles":
                alleles,

            "allele_count":
                len(
                    alleles
                ),

            "alleles_usable":
                True,

            "allele_source":
                (
                    "multi_allelic"
                    if len(
                        alleles
                    )
                    > 2
                    else (
                        "biallelic"
                        if len(
                            alleles
                        )
                        == 2
                        else "single_allele"
                    )
                ),
        }

    # ==================================================================
    # Numeric helpers
    # ==================================================================

    @classmethod
    def standardize_numeric(
        cls,
        value: Any,
    ) -> float | None:
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
    def standardize_probability(
        cls,
        value: Any,
    ) -> float | None:
        """Return a standardized probability or None."""

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
    # S1
    # ==================================================================

    @classmethod
    def standardize_s1_record(
        cls,
        record: dict[str, Any],
        *,
        feature_type: FeatureType,
    ) -> dict[str, Any]:
        """
        Standardize one Moradi S1 cis-QTL record.

        Source schema:
            SNP
            R
            A
            SNP_pos
            splicing_event
            Stat
            p_value
            FDR
            beta
            MAF
            AvgCall
            Rsq

        For iso-eQTL records, the feature column may still be represented
        through the source splicing_event field.
        """

        rsid = (
            cls.standardize_rsid(
                record.get(
                    "SNP"
                )
            )
        )

        coordinate = (
            cls.standardize_coordinate(
                record.get(
                    "SNP_pos"
                )
            )
        )

        feature = (
            cls.standardize_feature(
                record.get(
                    "splicing_event"
                ),
                feature_type,
            )
        )

        reference_allele = (
            cls.standardize_alleles(
                record.get(
                    "R"
                )
            )
        )

        alternate_allele = (
            cls.standardize_alleles(
                record.get(
                    "A"
                )
            )
        )

        p_value = (
            cls.standardize_probability(
                record.get(
                    "p_value"
                )
            )
        )

        fdr = (
            cls.standardize_probability(
                record.get(
                    "FDR"
                )
            )
        )

        maf = (
            cls.standardize_probability(
                record.get(
                    "MAF"
                )
            )
        )

        avg_call = (
            cls.standardize_probability(
                record.get(
                    "AvgCall"
                )
            )
        )

        rsq = (
            cls.standardize_probability(
                record.get(
                    "Rsq"
                )
            )
        )

        qtl_usable = (
            (
                rsid[
                    "rsid_usable"
                ]
                or coordinate[
                    "coordinate_usable"
                ]
            )
            and feature[
                "feature_usable"
            ]
            and p_value is not None
        )

        standardization_status = (
            "STANDARDIZED"
            if qtl_usable
            else (
                "PARTIAL"
                if (
                    rsid[
                        "rsid_usable"
                    ]
                    or coordinate[
                        "coordinate_usable"
                    ]
                    or feature[
                        "feature_usable"
                    ]
                )
                else "PRESERVED"
            )
        )

        return {
            "source_table":
                "S1",

            "qtl_scope":
                "cis",

            "snp_id_raw":
                record.get(
                    "SNP"
                ),

            **rsid,

            **coordinate,

            "reference_allele_raw":
                record.get(
                    "R"
                ),

            "reference_allele":
                (
                    reference_allele[
                        "alleles"
                    ][0]
                    if reference_allele[
                        "alleles_usable"
                    ]
                    else None
                ),

            "alternate_allele_raw":
                record.get(
                    "A"
                ),

            "alternate_alleles":
                (
                    alternate_allele[
                        "alleles"
                    ]
                    if alternate_allele[
                        "alleles_usable"
                    ]
                    else None
                ),

            "feature_raw":
                record.get(
                    "splicing_event"
                ),

            **feature,

            "statistic":
                cls.standardize_numeric(
                    record.get(
                        "Stat"
                    )
                ),

            "p_value":
                p_value,

            "fdr":
                fdr,

            "beta":
                cls.standardize_numeric(
                    record.get(
                        "beta"
                    )
                ),

            "maf":
                maf,

            "average_call_rate":
                avg_call,

            "rsq":
                rsq,

            "variant_usable":
                bool(
                    rsid[
                        "rsid_usable"
                    ]
                    or coordinate[
                        "coordinate_usable"
                    ]
                ),

            "qtl_usable":
                bool(
                    qtl_usable
                ),

            "standardization_status":
                standardization_status,
        }

    # ==================================================================
    # S2 / S4 LD-linked records
    # ==================================================================

    @classmethod
    def standardize_ld_record(
        cls,
        record: dict[str, Any],
        *,
        source_table: Literal[
            "S2",
            "S4",
        ],
        feature_type: FeatureType,
        qtl_scope: Literal[
            "cis",
            "trans",
        ],
        feature_column: str,
    ) -> dict[str, Any]:
        """
        Standardize one S2 or S4 LD-linked QTL record.

        Source semantic inversion is handled explicitly:

            sQTL_SNP
                QTL coordinate

            sQTL_SNP-pos
                QTL rsID

            tag_SNP
                GWAS tag coordinate

            tag_SNP_pos
                GWAS tag rsID/composite rsID
        """

        qtl_coordinate = (
            cls.standardize_coordinate(
                record.get(
                    "sQTL_SNP"
                )
            )
        )

        qtl_rsid = (
            cls.standardize_rsid(
                record.get(
                    "sQTL_SNP-pos"
                )
            )
        )

        tag_coordinate = (
            cls.standardize_coordinate(
                record.get(
                    "tag_SNP"
                )
            )
        )

        tag_rsid = (
            cls.standardize_composite_rsid(
                record.get(
                    "tag_SNP_pos"
                )
            )
        )

        feature = (
            cls.standardize_feature(
                record.get(
                    feature_column
                ),
                feature_type,
            )
        )

        ld_r2 = (
            cls.standardize_probability(
                record.get(
                    "LD"
                )
            )
        )

        qtl_variant_usable = (
            qtl_rsid[
                "rsid_usable"
            ]
            or qtl_coordinate[
                "coordinate_usable"
            ]
        )

        tag_variant_usable = (
            tag_rsid[
                "rsid_usable"
            ]
            or tag_coordinate[
                "coordinate_usable"
            ]
        )

        qtl_usable = (
            qtl_variant_usable
            and tag_variant_usable
            and feature[
                "feature_usable"
            ]
            and ld_r2 is not None
        )

        standardization_status = (
            "STANDARDIZED"
            if qtl_usable
            else (
                "PARTIAL"
                if (
                    qtl_variant_usable
                    or tag_variant_usable
                    or feature[
                        "feature_usable"
                    ]
                )
                else "PRESERVED"
            )
        )

        return {
            "source_table":
                source_table,

            "qtl_scope":
                qtl_scope,

            "qtl_coordinate_raw":
                record.get(
                    "sQTL_SNP"
                ),

            "qtl_chromosome":
                qtl_coordinate[
                    "chromosome"
                ],

            "qtl_position":
                qtl_coordinate[
                    "position"
                ],

            "qtl_coordinate_usable":
                qtl_coordinate[
                    "coordinate_usable"
                ],

            "qtl_rsid_raw":
                record.get(
                    "sQTL_SNP-pos"
                ),

            "qtl_rsid":
                qtl_rsid[
                    "canonical_rsid"
                ],

            "qtl_rsid_usable":
                qtl_rsid[
                    "rsid_usable"
                ],

            "feature_raw":
                record.get(
                    feature_column
                ),

            **feature,

            "tag_coordinate_raw":
                record.get(
                    "tag_SNP"
                ),

            "tag_chromosome":
                tag_coordinate[
                    "chromosome"
                ],

            "tag_position":
                tag_coordinate[
                    "position"
                ],

            "tag_coordinate_usable":
                tag_coordinate[
                    "coordinate_usable"
                ],

            "tag_rsid_raw":
                record.get(
                    "tag_SNP_pos"
                ),

            "tag_rsids":
                tag_rsid[
                    "rsids"
                ],

            "tag_primary_rsid":
                tag_rsid[
                    "primary_rsid"
                ],

            "tag_is_composite_rsid":
                tag_rsid[
                    "composite_rsid"
                ],

            "tag_rsid_usable":
                tag_rsid[
                    "rsid_usable"
                ],

            "ld_r2":
                ld_r2,

            "gwas_cancer":
                (
                    StandardizationUtils
                    .standardize_string(
                        record.get(
                            "gwas_cancer"
                        )
                    )
                    .standardized_value
                ),

            "genome_build":
                cls.GENOME_BUILD,

            "qtl_variant_usable":
                bool(
                    qtl_variant_usable
                ),

            "tag_variant_usable":
                bool(
                    tag_variant_usable
                ),

            "qtl_usable":
                bool(
                    qtl_usable
                ),

            "standardization_status":
                standardization_status,
        }

    # ==================================================================
    # S3
    # ==================================================================

    @classmethod
    def standardize_s3_record(
        cls,
        record: dict[str, Any],
        *,
        feature_type: FeatureType,
        feature_column: str,
        variant_column: str = "SNP",
    ) -> dict[str, Any]:
        """
        Standardize one Moradi S3 trans-QTL record.

        The variant column is configurable because the trans-exon source
        uses "T" instead of "SNP".
        """

        rsid = (
            cls.standardize_rsid(
                record.get(
                    variant_column
                )
            )
        )

        coordinate_value = (
            record.get(
                "SNP_pos"
            )
        )

        coordinate = (
            cls.standardize_coordinate(
                coordinate_value
            )
        )

        feature = (
            cls.standardize_feature(
                record.get(
                    feature_column
                ),
                feature_type,
            )
        )

        p_value = (
            cls.standardize_probability(
                record.get(
                    "p_value"
                )
            )
        )

        fdr = (
            cls.standardize_probability(
                record.get(
                    "FDR"
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
            variant_usable
            and feature[
                "feature_usable"
            ]
            and p_value is not None
        )

        standardization_status = (
            "STANDARDIZED"
            if qtl_usable
            else (
                "PARTIAL"
                if (
                    variant_usable
                    or feature[
                        "feature_usable"
                    ]
                )
                else "PRESERVED"
            )
        )

        result = {
            "source_table":
                "S3",

            "qtl_scope":
                "trans",

            "snp_id_raw":
                record.get(
                    variant_column
                ),

            "source_variant_column":
                variant_column,

            **rsid,

            **coordinate,

            "feature_raw":
                record.get(
                    feature_column
                ),

            **feature,

            "statistic":
                cls.standardize_numeric(
                    record.get(
                        "Stat"
                    )
                ),

            "p_value":
                p_value,

            "fdr":
                fdr,

            "beta":
                cls.standardize_numeric(
                    record.get(
                        "beta"
                    )
                ),

            "maf":
                cls.standardize_probability(
                    record.get(
                        "MAF"
                    )
                ),

            "average_call_rate":
                cls.standardize_probability(
                    record.get(
                        "AvgCall"
                    )
                ),

            "rsq":
                cls.standardize_probability(
                    record.get(
                        "Rsq"
                    )
                ),

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

        return result

    # ==================================================================
    # DataFrames
    # ==================================================================

    @classmethod
    def standardize_s1_dataframe(
        cls,
        dataframe: pd.DataFrame,
        *,
        feature_type: FeatureType,
    ) -> pd.DataFrame:
        """Standardize a Moradi S1 DataFrame."""

        return pd.DataFrame(
            [
                cls.standardize_s1_record(
                    record,
                    feature_type=(
                        feature_type
                    ),
                )
                for record
                in dataframe.to_dict(
                    orient="records"
                )
            ]
        )

    @classmethod
    def standardize_ld_dataframe(
        cls,
        dataframe: pd.DataFrame,
        *,
        source_table: Literal[
            "S2",
            "S4",
        ],
        feature_type: FeatureType,
        qtl_scope: Literal[
            "cis",
            "trans",
        ],
        feature_column: str,
    ) -> pd.DataFrame:
        """Standardize a Moradi S2 or S4 DataFrame."""

        return pd.DataFrame(
            [
                cls.standardize_ld_record(
                    record,
                    source_table=(
                        source_table
                    ),
                    feature_type=(
                        feature_type
                    ),
                    qtl_scope=(
                        qtl_scope
                    ),
                    feature_column=(
                        feature_column
                    ),
                )
                for record
                in dataframe.to_dict(
                    orient="records"
                )
            ]
        )

    @classmethod
    def standardize_s3_dataframe(
        cls,
        dataframe: pd.DataFrame,
        *,
        feature_type: FeatureType,
        feature_column: str,
        variant_column: str = "SNP",
    ) -> pd.DataFrame:
        """Standardize a Moradi S3 DataFrame."""

        return pd.DataFrame(
            [
                cls.standardize_s3_record(
                    record,
                    feature_type=(
                        feature_type
                    ),
                    feature_column=(
                        feature_column
                    ),
                    variant_column=(
                        variant_column
                    ),
                )
                for record
                in dataframe.to_dict(
                    orient="records"
                )
            ]
        )