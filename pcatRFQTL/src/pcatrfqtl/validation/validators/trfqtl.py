"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/trfqtl.py

Description:
    Provides validation utilities for tRNA-derived fragment
    quantitative trait loci (tRFQTL) datasets.

    This module validates identifiers, genomic coordinates,
    alleles, statistical values, false discovery rates, and
    linkage disequilibrium measurements used in the
    PCa-tRFQTL research pipeline.

    Validation anomalies are reported using stable anomaly IDs
    so that validation reports remain machine-readable,
    reproducible, and traceable across pipeline executions.

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
from collections.abc import Iterable
from typing import Any


class TRFQTLValidator:
    """
    Validator for tRFQTL dataset fields.

    The validator returns human-readable anomaly messages while
    assigning deterministic anomaly IDs through the internal
    anomaly registry.

    Anomaly ID format
    -----------------
    TRFQTL-<CATEGORY>-<NUMBER>

    Examples
    --------
    TRFQTL-SNP-001
    TRFQTL-COORD-001
    TRFQTL-ALLELE-001
    TRFQTL-TRF-001
    TRFQTL-PVALUE-001
    TRFQTL-FDR-001
    TRFQTL-LDR2-001
    """

    # ------------------------------------------------------------------
    # Stable anomaly identifiers
    # ------------------------------------------------------------------

    ANOMALY_IDS = {
        "snp": "TRFQTL-SNP-001",
        "coordinate": "TRFQTL-COORD-001",
        "allele": "TRFQTL-ALLELE-001",
        "trf_id": "TRFQTL-TRF-001",
        "p_value": "TRFQTL-PVALUE-001",
        "fdr": "TRFQTL-FDR-001",
        "ld_r2": "TRFQTL-LDR2-001",
        "effect_allele": "TRFQTL-EFFECT-ALLELE-001",
        "effect": "TRFQTL-EFFECT-001",
    }

    # ------------------------------------------------------------------
    # Regular expressions
    # ------------------------------------------------------------------

    SNP_PATTERN = re.compile(
        r"^rs[0-9]+$",
        re.IGNORECASE,
    )

    TRF_PATTERN = re.compile(
        r"^tRF-[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*$",
        re.IGNORECASE,
    )

    COORDINATE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)"
        r":"
        r"[0-9]+$",
        re.IGNORECASE,
    )

    ALLELE_PATTERN = re.compile(
        r"^(?:[ACGT]|[ACGT]+/[ACGT]+)$",
        re.IGNORECASE,
    )

    EFFECT_ALLELE_PATTERN = re.compile(
        r"^[ACGT]+$",
        re.IGNORECASE,
    )

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _is_missing(value: Any) -> bool:
        """
        Determine whether a value should be treated as missing.
        """

        if value is None:
            return True

        try:
            return bool(value != value)
        except Exception:
            return False

    @classmethod
    def _format_error(
        cls,
        anomaly_id: str,
        field: str,
        index: int,
        value: Any,
        message: str,
    ) -> str:
        """
        Construct a standardized anomaly message.
        """

        return (
            f"[{anomaly_id}] "
            f"{field}[{index}] "
            f"{message}: {value!r}."
        )

    # ------------------------------------------------------------------
    # SNP identifiers
    # ------------------------------------------------------------------

    @classmethod
    def validate_snp_ids(
        cls,
        values: Iterable[Any],
        field: str = "SNP ID",
    ) -> list[str]:
        """
        Validate dbSNP identifiers.

        Valid examples
        ---------------
        rs123
        rs456789
        RS123

        Parameters
        ----------
        values:
            Iterable containing SNP identifiers.
        field:
            Field name used in error messages.
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["snp"],
                        field,
                        index,
                        value,
                        "missing SNP identifier",
                    )
                )
                continue

            value_str = str(value).strip()

            if not cls.SNP_PATTERN.fullmatch(value_str):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["snp"],
                        field,
                        index,
                        value,
                        "invalid SNP identifier",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # Genomic coordinates
    # ------------------------------------------------------------------

    @classmethod
    def validate_coordinates(
        cls,
        values: Iterable[Any],
        field: str = "SNP position",
    ) -> list[str]:
        """
        Validate genomic coordinates.

        Supported formats
        -----------------
        chr1:123456
        1:123456
        chrX:123456
        X:123456
        chrMT:123456
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["coordinate"],
                        field,
                        index,
                        value,
                        "missing genomic coordinate",
                    )
                )
                continue

            if not isinstance(value, str):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["coordinate"],
                        field,
                        index,
                        value,
                        "coordinate must be a string",
                    )
                )
                continue

            value_str = value.strip()

            if not cls.COORDINATE_PATTERN.fullmatch(value_str):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["coordinate"],
                        field,
                        index,
                        value,
                        "invalid genomic coordinate",
                    )
                )
                continue

            try:
                position = int(value_str.split(":")[1])

                if position <= 0:
                    errors.append(
                        cls._format_error(
                            cls.ANOMALY_IDS["coordinate"],
                            field,
                            index,
                            value,
                            "genomic position must be greater than zero",
                        )
                    )
            except (IndexError, ValueError):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["coordinate"],
                        field,
                        index,
                        value,
                        "invalid genomic position",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # Alleles
    # ------------------------------------------------------------------

    @classmethod
    def validate_alleles(
        cls,
        values: Iterable[Any],
        field: str = "Alleles",
    ) -> list[str]:
        """
        Validate allele representations.

        Accepted representations
        -------------------------
        A
        C
        G
        T
        A/T
        G/C
        AT/GC

        Single alleles are accepted because some source datasets
        represent effect/reference alleles independently.
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["allele"],
                        field,
                        index,
                        value,
                        "missing allele representation",
                    )
                )
                continue

            value_str = str(value).strip()

            if not cls.ALLELE_PATTERN.fullmatch(value_str):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["allele"],
                        field,
                        index,
                        value,
                        "invalid allele representation",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # tRF identifiers
    # ------------------------------------------------------------------

    @classmethod
    def validate_trf_ids(
        cls,
        values: Iterable[Any],
        field: str = "tRF",
    ) -> list[str]:
        """
        Validate tRF identifiers.

        Examples
        --------
        tRF-30-34HWH3RXSINH
        tRF-23-V47PU9XW0N
        tRF-20-9N15WV7W
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["trf_id"],
                        field,
                        index,
                        value,
                        "missing tRF identifier",
                    )
                )
                continue

            value_str = str(value).strip()

            if not cls.TRF_PATTERN.fullmatch(value_str):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["trf_id"],
                        field,
                        index,
                        value,
                        "invalid tRF identifier",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # P-values
    # ------------------------------------------------------------------

    @classmethod
    def validate_p_values(
        cls,
        values: Iterable[Any],
        field: str = "P-value",
    ) -> list[str]:
        """
        Validate statistical P-values.

        Valid range
        -----------
        0 < p <= 1
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["p_value"],
                        field,
                        index,
                        value,
                        "missing P-value",
                    )
                )
                continue

            try:
                numeric_value = float(value)
            except (TypeError, ValueError):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["p_value"],
                        field,
                        index,
                        value,
                        "P-value must be numeric",
                    )
                )
                continue

            if not 0 < numeric_value <= 1:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["p_value"],
                        field,
                        index,
                        value,
                        "P-value must satisfy 0 < p <= 1",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # FDR
    # ------------------------------------------------------------------

    @classmethod
    def validate_fdr(
        cls,
        values: Iterable[Any],
        field: str = "FDR",
    ) -> list[str]:
        """
        Validate false discovery rate values.

        Valid range
        -----------
        0 <= FDR <= 1
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["fdr"],
                        field,
                        index,
                        value,
                        "missing FDR",
                    )
                )
                continue

            try:
                numeric_value = float(value)
            except (TypeError, ValueError):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["fdr"],
                        field,
                        index,
                        value,
                        "FDR must be numeric",
                    )
                )
                continue

            if not 0 <= numeric_value <= 1:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["fdr"],
                        field,
                        index,
                        value,
                        "FDR must satisfy 0 <= FDR <= 1",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # LD r2
    # ------------------------------------------------------------------

    @classmethod
    def validate_ld_r2(
        cls,
        values: Iterable[Any],
        field: str = "LD (r2)",
    ) -> list[str]:
        """
        Validate linkage disequilibrium r² values.

        Valid range
        -----------
        0 <= r² <= 1
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["ld_r2"],
                        field,
                        index,
                        value,
                        "missing LD r2 value",
                    )
                )
                continue

            try:
                numeric_value = float(value)
            except (TypeError, ValueError):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["ld_r2"],
                        field,
                        index,
                        value,
                        "LD r2 must be numeric",
                    )
                )
                continue

            if not 0 <= numeric_value <= 1:
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["ld_r2"],
                        field,
                        index,
                        value,
                        "LD r2 must satisfy 0 <= r2 <= 1",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # Effect allele
    # ------------------------------------------------------------------

    @classmethod
    def validate_effect_alleles(
        cls,
        values: Iterable[Any],
        field: str = "A1 (effect allele)",
    ) -> list[str]:
        """
        Validate effect allele values.

        This method is intentionally stricter than
        ``validate_alleles`` because the Cancer-tRFQTL S10
        dataset stores a single effect allele.
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["effect_allele"],
                        field,
                        index,
                        value,
                        "missing effect allele",
                    )
                )
                continue

            value_str = str(value).strip()

            if not cls.EFFECT_ALLELE_PATTERN.fullmatch(value_str):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["effect_allele"],
                        field,
                        index,
                        value,
                        "invalid effect allele",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # Effect sizes
    # ------------------------------------------------------------------

    @classmethod
    def validate_effects(
        cls,
        values: Iterable[Any],
        field: str = "Effect",
    ) -> list[str]:
        """
        Validate GWAS effect-size values.

        Effect sizes may be positive or negative.
        """

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["effect"],
                        field,
                        index,
                        value,
                        "missing effect size",
                    )
                )
                continue

            try:
                float(value)
            except (TypeError, ValueError):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["effect"],
                        field,
                        index,
                        value,
                        "effect size must be numeric",
                    )
                )

        return errors

    # ------------------------------------------------------------------
    # Anomaly metadata
    # ------------------------------------------------------------------

    @classmethod
    def anomaly_catalog(cls) -> dict[str, dict[str, str]]:
        """
        Return the stable anomaly catalog.

        This is useful when constructing machine-readable
        validation reports.
        """

        return {
            "TRFQTL-SNP-001": {
                "category": "identifier",
                "field": "SNP ID",
                "description": "Invalid or missing SNP identifier.",
            },
            "TRFQTL-COORD-001": {
                "category": "genomic_coordinate",
                "field": "SNP position",
                "description": "Invalid or missing genomic coordinate.",
            },
            "TRFQTL-ALLELE-001": {
                "category": "allele",
                "field": "Alleles",
                "description": "Invalid or missing allele representation.",
            },
            "TRFQTL-TRF-001": {
                "category": "identifier",
                "field": "tRF",
                "description": "Invalid or missing tRF identifier.",
            },
            "TRFQTL-PVALUE-001": {
                "category": "statistical",
                "field": "P-value",
                "description": "Invalid or missing P-value.",
            },
            "TRFQTL-FDR-001": {
                "category": "statistical",
                "field": "FDR",
                "description": "Invalid or missing false discovery rate.",
            },
            "TRFQTL-LDR2-001": {
                "category": "linkage_disequilibrium",
                "field": "LD (r2)",
                "description": "Invalid or missing LD r² value.",
            },
            "TRFQTL-EFFECT-ALLELE-001": {
                "category": "effect_allele",
                "field": "A1 (effect allele)",
                "description": "Invalid or missing effect allele.",
            },
            "TRFQTL-EFFECT-001": {
                "category": "effect_size",
                "field": "Effect",
                "description": "Invalid or missing GWAS effect size.",
            },
        }