"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/validators/moradi.py

Description:
    Field-level validation utilities for Moradi supplementary
    QTL datasets.

    Supports dbSNP identifiers, genomic coordinates, reference and
    alternate alleles including indels and multi-allelic variants,
    splicing-event identifiers, Ensembl transcript identifiers,
    statistical probabilities, LD values, and generic numerical fields.

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
from typing import Any, Iterable


class MoradiQTLValidator:
    """Field-level validator for Moradi QTL supplementary data."""

    SNP_PATTERN = re.compile(
        r"^rs\d+$",
        re.IGNORECASE,
    )

    COORDINATE_PATTERN = re.compile(
        r"^(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT):[1-9]\d*$",
        re.IGNORECASE,
    )

    ALLELE_PATTERN = re.compile(
        r"^[ACGT]+(?:,[ACGT]+)*$",
        re.IGNORECASE,
    )

    INTRON_EVENT_PATTERN = re.compile(
        r"^INT\d+$",
        re.IGNORECASE,
    )

    EXON_EVENT_PATTERN = re.compile(
        r"^EX\d+$",
        re.IGNORECASE,
    )

    TRANSCRIPT_PATTERN = re.compile(
        r"^ENST\d+(?:\.\d+)?$",
        re.IGNORECASE,
    )

    COMPOSITE_SNP_PATTERN = re.compile(
        r"^rs\d+(?::rs\d+)*$",
        re.IGNORECASE,
    )

    ANOMALY_IDS = {
        "snp": "MORADI-SNP-001",
        "coordinate": "MORADI-COORD-001",
        "allele": "MORADI-ALLELE-001",
        "event": "MORADI-EVENT-001",
        "transcript": "MORADI-TRANSCRIPT-001",
        "p_value": "MORADI-PVALUE-001",
        "fdr": "MORADI-FDR-001",
        "ld": "MORADI-LD-001",
        "maf": "MORADI-MAF-001",
        "call_rate": "MORADI-CALLRATE-001",
        "rsq": "MORADI-RSQ-001",
        "numeric": "MORADI-NUMERIC-001",
        "missing": "MORADI-MISSING-001",
    }

    @staticmethod
    def _is_missing(value: Any) -> bool:
        """Return True when a value should be treated as missing."""

        if value is None:
            return True

        if isinstance(value, str):
            return not value.strip()

        try:
            return bool(math.isnan(value))
        except (TypeError, ValueError):
            return False

    @staticmethod
    def _format_error(
        anomaly_id: str,
        field: str,
        index: int,
        value: Any,
        message: str,
    ) -> str:
        """Return a consistent validation message."""

        return (
            f"[{anomaly_id}] "
            f"{field}[{index}] "
            f"{message}: {value!r}."
        )

    @classmethod
    def validate_snp_ids(
        cls,
        values: Iterable[Any],
        field: str = "SNP",
    ) -> list[str]:
        """Validate dbSNP rs identifiers."""

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

            if not cls.SNP_PATTERN.fullmatch(
                str(value).strip()
            ):
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

    @classmethod
    def validate_composite_snp_ids(
        cls,
        values: list[object],
        field_name: str = "SNP",
    ) -> list[str]:
        """
        Validate single or colon-separated dbSNP identifiers.

        Accepted examples:
            rs12345
            rs12345:rs67890
            rs12345:rs67890:rs11111

        This representation is present in Moradi GWAS-linked QTL tables
        and reflects multiple dbSNP identifiers associated with one
        reported genomic tag position.
        """

        errors: list[str] = []

        for index, value in enumerate(values):

            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        "MORADI-SNP-001",
                        field_name,
                        index,
                        value,
                        "missing SNP identifier",
                    )
                )

                continue

            text = str(value).strip()

            if not cls.COMPOSITE_SNP_PATTERN.fullmatch(text):
                errors.append(
                    cls._format_error(
                        "MORADI-SNP-001",
                        field_name,
                        index,
                        value,
                        (
                            "invalid single or colon-separated "
                            "dbSNP identifier"
                        ),
                    )
                )

        return errors

    @classmethod
    def validate_coordinates(
        cls,
        values: Iterable[Any],
        field: str = "SNP_pos",
    ) -> list[str]:
        """Validate chromosome:position coordinates."""

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

            if not cls.COORDINATE_PATTERN.fullmatch(
                str(value).strip()
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["coordinate"],
                        field,
                        index,
                        value,
                        "invalid genomic coordinate",
                    )
                )

        return errors

    @classmethod
    def validate_alleles(
        cls,
        values: Iterable[Any],
        field: str,
    ) -> list[str]:
        """
        Validate Moradi allele representations.

        Examples accepted:
            A
            C
            CT
            GGT
            A,T
            C,T
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
                        "missing allele",
                    )
                )
                continue

            if not cls.ALLELE_PATTERN.fullmatch(
                str(value).strip()
            ):
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

    @classmethod
    def validate_splicing_events(
        cls,
        values: Iterable[Any],
        field: str = "splicing_event",
    ) -> list[str]:
        """Validate INT/EX splicing-event identifiers."""

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["event"],
                        field,
                        index,
                        value,
                        "missing splicing-event identifier",
                    )
                )
                continue

            value_str = str(value).strip()

            if not (
                cls.INTRON_EVENT_PATTERN.fullmatch(value_str)
                or cls.EXON_EVENT_PATTERN.fullmatch(value_str)
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["event"],
                        field,
                        index,
                        value,
                        "invalid splicing-event identifier",
                    )
                )

        return errors

    @classmethod
    def validate_transcripts(
        cls,
        values: Iterable[Any],
        field: str = "mRNA_isoform",
    ) -> list[str]:
        """Validate Ensembl transcript identifiers."""

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["transcript"],
                        field,
                        index,
                        value,
                        "missing transcript identifier",
                    )
                )
                continue

            if not cls.TRANSCRIPT_PATTERN.fullmatch(
                str(value).strip()
            ):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["transcript"],
                        field,
                        index,
                        value,
                        "invalid Ensembl transcript identifier",
                    )
                )

        return errors

    @classmethod
    def _validate_probability(
        cls,
        values: Iterable[Any],
        field: str,
        anomaly_id: str,
    ) -> list[str]:
        """Validate numerical values in the inclusive range [0, 1]."""

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                errors.append(
                    cls._format_error(
                        anomaly_id,
                        field,
                        index,
                        value,
                        "missing value",
                    )
                )
                continue

            try:
                numeric = float(value)
            except (TypeError, ValueError):
                errors.append(
                    cls._format_error(
                        anomaly_id,
                        field,
                        index,
                        value,
                        "non-numeric value",
                    )
                )
                continue

            if not math.isfinite(numeric):
                errors.append(
                    cls._format_error(
                        anomaly_id,
                        field,
                        index,
                        value,
                        "non-finite value",
                    )
                )
                continue

            if not 0 <= numeric <= 1:
                errors.append(
                    cls._format_error(
                        anomaly_id,
                        field,
                        index,
                        value,
                        "value outside range [0, 1]",
                    )
                )

        return errors

    @classmethod
    def validate_p_values(
        cls,
        values: Iterable[Any],
        field: str = "p_value",
    ) -> list[str]:
        """Validate P-values."""

        return cls._validate_probability(
            values,
            field,
            cls.ANOMALY_IDS["p_value"],
        )

    @classmethod
    def validate_fdr(
        cls,
        values: Iterable[Any],
        field: str = "FDR",
    ) -> list[str]:
        """Validate FDR values."""

        return cls._validate_probability(
            values,
            field,
            cls.ANOMALY_IDS["fdr"],
        )

    @classmethod
    def validate_ld(
        cls,
        values: Iterable[Any],
        field: str = "LD",
    ) -> list[str]:
        """Validate LD r² values."""

        return cls._validate_probability(
            values,
            field,
            cls.ANOMALY_IDS["ld"],
        )

    @classmethod
    def validate_maf(
        cls,
        values: Iterable[Any],
        field: str = "MAF",
    ) -> list[str]:
        """Validate minor-allele frequencies."""

        return cls._validate_probability(
            values,
            field,
            cls.ANOMALY_IDS["maf"],
        )

    @classmethod
    def validate_call_rate(
        cls,
        values: Iterable[Any],
        field: str = "AvgCall",
    ) -> list[str]:
        """Validate average genotype call rates."""

        return cls._validate_probability(
            values,
            field,
            cls.ANOMALY_IDS["call_rate"],
        )

    @classmethod
    def validate_rsq(
        cls,
        values: Iterable[Any],
        field: str = "Rsq",
    ) -> list[str]:
        """Validate imputation-quality Rsq values."""

        return cls._validate_probability(
            values,
            field,
            cls.ANOMALY_IDS["rsq"],
        )

    @classmethod
    def validate_numeric(
        cls,
        values: Iterable[Any],
        field: str,
        *,
        allow_missing: bool = False,
    ) -> list[str]:
        """Validate finite numerical values."""

        errors: list[str] = []

        for index, value in enumerate(values):
            if cls._is_missing(value):
                if allow_missing:
                    continue

                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["missing"],
                        field,
                        index,
                        value,
                        "missing numeric value",
                    )
                )
                continue

            try:
                numeric = float(value)
            except (TypeError, ValueError):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["numeric"],
                        field,
                        index,
                        value,
                        "non-numeric value",
                    )
                )
                continue

            if not math.isfinite(numeric):
                errors.append(
                    cls._format_error(
                        cls.ANOMALY_IDS["numeric"],
                        field,
                        index,
                        value,
                        "non-finite numeric value",
                    )
                )

        return errors