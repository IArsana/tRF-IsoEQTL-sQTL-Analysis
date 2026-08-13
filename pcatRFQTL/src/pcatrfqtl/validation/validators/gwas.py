"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/validators/gwas.py

Description:
    Field-level validation and classification utilities for GWAS Catalog
    association data.

    GWAS Catalog "SNPS" values are treated as heterogeneous variant
    descriptors rather than strict dbSNP identifiers because source
    association records may contain:

        - single dbSNP rsIDs
        - multiple rsIDs
        - SNP-by-SNP interaction descriptors
        - coordinate-based variant descriptors
        - allele-resolved coordinate descriptors
        - insertion/deletion coordinate descriptors
        - HLA allele descriptors
        - array/probe identifiers
        - structural-variant identifiers
        - gene or source-specific variant labels
        - other non-empty source descriptors

    The validator also supports:

        - strict rsID validation when required
        - chromosome validation
        - interaction-aware chromosome fields
        - genomic-position validation
        - interaction-aware position fields
        - P-value validation
        - configurable zero-P-value policy
        - configurable missing-P-value policy
        - generic finite numeric validation

    Important:
        Validation never modifies raw source values.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import math
import re
from typing import Any

import pandas as pd


class GWASValidator:
    """
    Field-level validator for GWAS Catalog association records.

    Notes
    -----
    GWAS Catalog association exports use heterogeneous representations
    in the "SNPS" field. Therefore:

        SNPS != strict dbSNP identifier field

    Use:

        validate_snp_ids()

    only when a field is semantically required to contain rsIDs.

    Use:

        classify_variant_descriptor()
        validate_variant_descriptors()

    for the heterogeneous GWAS Catalog association "SNPS" field.
    """

    # ==================================================================
    # Stable rule identifiers
    # ==================================================================

    ANOMALY_IDS = {
        "snp": "GWAS-SNP-001",
        "variant": "GWAS-VARIANT-001",
        "chromosome": "GWAS-CHR-001",
        "position": "GWAS-POS-001",
        "pvalue": "GWAS-PVALUE-001",
        "numeric": "GWAS-NUMERIC-001",
        "missing": "GWAS-MISSING-001",
    }

    # ==================================================================
    # Core patterns
    # ==================================================================

    SNP_PATTERN = re.compile(
        r"^rs\d+$",
        re.IGNORECASE,
    )

    MULTI_RSID_PATTERN = re.compile(
        r"^rs\d+(?:[;,\s]+rs\d+)+$",
        re.IGNORECASE,
    )

    INTERACTION_PATTERN = re.compile(
        r"^rs\d+\s+[xX]\s+rs\d+$",
        re.IGNORECASE,
    )

    CHROMOSOME_PATTERN = re.compile(
        r"^(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)$",
        re.IGNORECASE,
    )

    # Examples:
    #   chr12:9098995
    #   12:9098995
    VARIANT_COORDINATE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
        r":\d+$",
        re.IGNORECASE,
    )

    # Examples:
    #   9:139737088:G:A
    #   chr1:12345:A:G
    VARIANT_ALLELIC_COLON_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
        r":\d+"
        r":[ACGTN]+"
        r":[ACGTN]+$",
        re.IGNORECASE,
    )

    # Examples:
    #   15:79372468:A_ATG
    #   12:46182832:CA_C
    #   7:85738582:TTA_T
    VARIANT_ALLELIC_UNDERSCORE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
        r":\d+"
        r":[ACGTN]+_[ACGTN]+$",
        re.IGNORECASE,
    )

    # Example:
    #   7:120812727_G_C
    COORDINATE_UNDERSCORE_ALLELE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
        r":\d+_[ACGTN]+_[ACGTN]+$",
        re.IGNORECASE,
    )

    # Examples:
    #   chr2:134471945:D
    #   chr7:19622838:I
    COORDINATE_INDEL_CODE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
        r":\d+"
        r":[ID]$",
        re.IGNORECASE,
    )

    # Example:
    #   rs11639303_G_A
    RSID_ALLELE_PATTERN = re.compile(
        r"^rs\d+_[ACGTN]+_[ACGTN]+$",
        re.IGNORECASE,
    )

    # Examples:
    #   HLA-DRB1*14:04
    #   HLA-B*07:02
    #   B*08:01
    #   DQA1*03:01
    HLA_PATTERN = re.compile(
        r"^(?:HLA-)?"
        r"[A-Z0-9]+(?:-[A-Z0-9]+)?"
        r"\*\d+(?::\d+)?$",
        re.IGNORECASE,
    )

    # Historical compact representation:
    #   HLA-A*3101
    HLA_COMPACT_PATTERN = re.compile(
        r"^(?:HLA-)?"
        r"[A-Z0-9]+(?:-[A-Z0-9]+)?"
        r"\*\d{4,}$",
        re.IGNORECASE,
    )

    # Examples:
    #   exm130158
    #   kgp4136779
    #   i4000416
    ARRAY_PROBE_PATTERN = re.compile(
        r"^(?:exm|kgp|i)\d+$",
        re.IGNORECASE,
    )

    # Examples:
    #   nsv831124
    #   esv3596105
    STRUCTURAL_VARIANT_PATTERN = re.compile(
        r"^(?:nsv|esv)\d+$",
        re.IGNORECASE,
    )

    # Examples:
    #   SLC22A8
    #   LYPD6B
    #   PPP5C
    #   downstreamRASGEF1B
    GENE_LABEL_PATTERN = re.compile(
        r"^[A-Za-z][A-Za-z0-9_-]*$"
    )

    # ==================================================================
    # Missing values
    # ==================================================================

    @staticmethod
    def _is_missing(
        value: Any,
    ) -> bool:
        """Return True when a value should be treated as missing."""

        if value is None:
            return True

        try:
            if pd.isna(value):
                return True
        except (
            TypeError,
            ValueError,
        ):
            pass

        if isinstance(
            value,
            str,
        ):
            return (
                value.strip().lower()
                in {
                    "",
                    "na",
                    "nan",
                    "none",
                    "null",
                }
            )

        return False

    # ==================================================================
    # Error formatting
    # ==================================================================

    @classmethod
    def _format_error(
        cls,
        rule_id: str,
        field_name: str,
        index: int,
        value: Any,
        message: str,
    ) -> str:
        """Create a standardized field-validation message."""

        return (
            f"[{rule_id}] "
            f"{field_name}[{index}] "
            f"{message}: {value!r}"
        )

    # ==================================================================
    # Strict rsID validation
    # ==================================================================

    @classmethod
    def _split_snp_field(
        cls,
        value: Any,
    ) -> list[str]:
        """
        Split a conventional multi-rsID field.

        Supported separators:
            semicolon
            comma
            whitespace
        """

        text = str(
            value
        ).strip()

        return [
            token
            for token
            in re.split(
                r"[;,\s]+",
                text,
            )
            if token
        ]

    @classmethod
    def validate_snp_ids(
        cls,
        values: list[Any],
        field_name: str = "SNPS",
        *,
        allow_missing: bool = True,
    ) -> list[str]:
        """
        Strictly validate dbSNP identifiers.

        This method intentionally rejects non-rsID GWAS representations.
        """

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):
            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            "missing SNP identifier",
                        )
                    )

                continue

            tokens = (
                cls._split_snp_field(
                    value
                )
            )

            if not tokens:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "snp"
                        ],
                        field_name,
                        index,
                        value,
                        "invalid SNP identifier",
                    )
                )

                continue

            invalid = any(
                cls.SNP_PATTERN
                .fullmatch(
                    token
                )
                is None
                for token
                in tokens
            )

            if invalid:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "snp"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "invalid dbSNP identifier "
                            "or multi-rsID representation"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # Variant descriptor classification
    # ==================================================================

    @classmethod
    def classify_variant_descriptor(
        cls,
        value: Any,
    ) -> str:
        """
        Classify a GWAS Catalog SNPS source descriptor.

        Returns
        -------
        str
            One of:

                missing
                single_rsid
                multi_rsid
                snp_interaction
                coordinate_variant
                allelic_variant
                coordinate_indel
                rsid_allelic_variant
                hla_allele
                array_probe
                structural_variant
                gene_or_source_label
                other
        """

        if cls._is_missing(
            value
        ):
            return "missing"

        text = str(
            value
        ).strip()

        if (
            cls.SNP_PATTERN
            .fullmatch(text)
        ):
            return "single_rsid"

        if (
            cls.INTERACTION_PATTERN
            .fullmatch(text)
        ):
            return "snp_interaction"

        if (
            cls.MULTI_RSID_PATTERN
            .fullmatch(text)
        ):
            return "multi_rsid"

        if (
            cls.VARIANT_COORDINATE_PATTERN
            .fullmatch(text)
        ):
            return "coordinate_variant"

        if (
            cls.VARIANT_ALLELIC_COLON_PATTERN
            .fullmatch(text)
        ):
            return "allelic_variant"

        if (
            cls.VARIANT_ALLELIC_UNDERSCORE_PATTERN
            .fullmatch(text)
        ):
            return "allelic_variant"

        if (
            cls.COORDINATE_UNDERSCORE_ALLELE_PATTERN
            .fullmatch(text)
        ):
            return "allelic_variant"

        if (
            cls.COORDINATE_INDEL_CODE_PATTERN
            .fullmatch(text)
        ):
            return "coordinate_indel"

        if (
            cls.RSID_ALLELE_PATTERN
            .fullmatch(text)
        ):
            return "rsid_allelic_variant"

        if (
            cls.HLA_PATTERN
            .fullmatch(text)
            or cls.HLA_COMPACT_PATTERN
            .fullmatch(text)
        ):
            return "hla_allele"

        if (
            cls.ARRAY_PROBE_PATTERN
            .fullmatch(text)
        ):
            return "array_probe"

        if (
            cls.STRUCTURAL_VARIANT_PATTERN
            .fullmatch(text)
        ):
            return "structural_variant"

        if (
            cls.GENE_LABEL_PATTERN
            .fullmatch(text)
        ):
            return "gene_or_source_label"

        return "other"

    @classmethod
    def validate_variant_descriptors(
        cls,
        values: list[Any],
        field_name: str = "SNPS",
        *,
        allow_missing: bool = True,
    ) -> list[str]:
        """
        Structurally validate heterogeneous GWAS Catalog descriptors.

        Recognized and unclassified non-empty source descriptors are
        retained. Classification is handled separately from validation.

        Missing descriptors are reported only when allow_missing=False.
        """

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):
            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            (
                                "missing variant "
                                "descriptor"
                            ),
                        )
                    )

                continue

            text = str(
                value
            ).strip()

            if not text:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "variant"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "empty variant "
                            "descriptor"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # Interaction-aware field splitting
    # ==================================================================

    @classmethod
    def _split_interaction_field(
        cls,
        value: Any,
    ) -> list[str]:
        """
        Split single, multi-value, or interaction fields.

        Examples:
            6
            6;7
            6,7
            6 x 6

            12345
            12345;67890
            12345 x 67890
        """

        text = str(
            value
        ).strip()

        return [
            token.strip()
            for token
            in re.split(
                r"\s+[xX]\s+|[;,]+",
                text,
            )
            if token.strip()
        ]

    # ==================================================================
    # Chromosome validation
    # ==================================================================

    @classmethod
    def validate_chromosomes(
        cls,
        values: list[Any],
        field_name: str = "CHR_ID",
        *,
        allow_missing: bool = True,
    ) -> list[str]:
        """
        Validate chromosome identifiers.

        Supported source representations include:
            6
            chr6
            6;7
            6,7
            6 x 6
            chr1 x chr2
        """

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):
            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            "missing chromosome",
                        )
                    )

                continue

            tokens = (
                cls._split_interaction_field(
                    value
                )
            )

            if not tokens:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "chromosome"
                        ],
                        field_name,
                        index,
                        value,
                        "invalid chromosome",
                    )
                )

                continue

            invalid = any(
                cls.CHROMOSOME_PATTERN
                .fullmatch(
                    token
                )
                is None
                for token
                in tokens
            )

            if invalid:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "chromosome"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "invalid chromosome "
                            "identifier"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # Genomic-position validation
    # ==================================================================

    @classmethod
    def validate_positions(
        cls,
        values: list[Any],
        field_name: str = "CHR_POS",
        *,
        allow_missing: bool = True,
    ) -> list[str]:
        """
        Validate positive integer genomic positions.

        Supported examples:
            100
            12345678
            100;200
            100,200
            31268092 x 31463174
        """

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):
            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            (
                                "missing genomic "
                                "position"
                            ),
                        )
                    )

                continue

            tokens = (
                cls._split_interaction_field(
                    value
                )
            )

            if not tokens:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "position"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "invalid genomic "
                            "position"
                        ),
                    )
                )

                continue

            valid = True

            for token in tokens:
                try:
                    numeric = float(
                        token
                    )
                except (
                    TypeError,
                    ValueError,
                ):
                    valid = False
                    break

                if (
                    not math.isfinite(
                        numeric
                    )
                    or numeric <= 0
                    or not numeric.is_integer()
                ):
                    valid = False
                    break

            if not valid:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "position"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "invalid positive integer "
                            "genomic position"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # P-value validation
    # ==================================================================

    @classmethod
    def validate_p_values(
        cls,
        values: list[Any],
        field_name: str = "P-VALUE",
        *,
        allow_zero: bool = True,
        allow_missing: bool = False,
    ) -> list[str]:
        """
        Validate association P-values.

        Parameters
        ----------
        values:
            P-value values to validate.

        field_name:
            Source field name used in validation messages.

        allow_zero:
            If True, P=0 is accepted at field-validation level.

            This is appropriate for GWAS Catalog source QC, where zero
            may represent values below numeric or export precision.

            If False, valid values must satisfy:

                0 < P <= 1

        allow_missing:
            If True, missing P-values are retained without generating a
            field-level error.

            Dataset-level runners should count and report these values
            separately as source-data warnings.

            If False, missing P-values are validation errors.

        Default valid range:
            0 <= P <= 1

        Notes
        -----
        This function never replaces, imputes, or transforms source
        P-values.
        """

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):
            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            "missing P-value",
                        )
                    )

                continue

            try:
                numeric = float(
                    value
                )
            except (
                TypeError,
                ValueError,
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "pvalue"
                        ],
                        field_name,
                        index,
                        value,
                        "non-numeric P-value",
                    )
                )

                continue

            invalid = (
                not math.isfinite(
                    numeric
                )
                or numeric < 0
                or numeric > 1
            )

            if (
                not allow_zero
                and numeric == 0
            ):
                invalid = True

            if invalid:
                valid_range = (
                    "0 <= P <= 1"
                    if allow_zero
                    else "0 < P <= 1"
                )

                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "pvalue"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "P-value outside valid "
                            f"range {valid_range}"
                        ),
                    )
                )

        return errors

    # ==================================================================
    # Generic numeric validation
    # ==================================================================

    @classmethod
    def validate_numeric(
        cls,
        values: list[Any],
        field_name: str,
        *,
        allow_missing: bool = True,
    ) -> list[str]:
        """Validate generic finite numeric values."""

        errors: list[str] = []

        for index, value in enumerate(
            values
        ):
            if cls._is_missing(
                value
            ):
                if not allow_missing:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS[
                                "missing"
                            ],
                            field_name,
                            index,
                            value,
                            (
                                "missing numeric "
                                "value"
                            ),
                        )
                    )

                continue

            try:
                numeric = float(
                    value
                )
            except (
                TypeError,
                ValueError,
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "numeric"
                        ],
                        field_name,
                        index,
                        value,
                        "non-numeric value",
                    )
                )

                continue

            if not math.isfinite(
                numeric
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS[
                            "numeric"
                        ],
                        field_name,
                        index,
                        value,
                        (
                            "non-finite numeric "
                            "value"
                        ),
                    )
                )

        return errors