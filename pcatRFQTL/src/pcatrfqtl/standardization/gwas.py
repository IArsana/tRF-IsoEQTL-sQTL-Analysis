"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/standardization/gwas.py

Description:
    Source-specific standardization utilities for GWAS Catalog
    association records.

    The GWAS Catalog "SNPS" field is heterogeneous and may contain:

        - single rsIDs
        - multiple rsIDs
        - SNP-by-SNP interactions
        - coordinate variants
        - allele-resolved variants
        - HLA alleles
        - array/probe identifiers
        - structural variants
        - source-specific descriptors

    This module derives standardized representations while preserving
    raw source values and provenance.

    Responsibilities:

        - variant descriptor classification
        - canonical rsID derivation where unambiguous
        - chromosome standardization
        - genomic-position standardization
        - interaction-aware handling
        - P-value standardization
        - reconstruction from PVALUE_MLOG when P-VALUE is zero
        - risk-allele extraction
        - downstream usability flags
        - row-level standardization status

    This module does not:

        - read source files
        - write output files
        - perform genome-build liftover
        - perform disease filtering
        - deduplicate source associations

    Raw source data are never modified.

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

from pcatrfqtl.standardization.common import (
    StandardizationUtils,
)
from pcatrfqtl.validation.validators.gwas import (
    GWASValidator,
)


class GWASStandardizer:
    """Standardize validated GWAS Catalog association records."""

    INTERACTION_SEPARATOR = re.compile(
        r"\s+[xX]\s+"
    )

    MULTI_RSID_SEPARATOR = re.compile(
        r"[;,\s]+"
    )

    STRONGEST_RISK_ALLELE_PATTERN = re.compile(
        r"^(?P<variant>.+?)-(?P<allele>[A-Za-z0-9?]+)$"
    )

    # ==================================================================
    # Variant classification
    # ==================================================================

    @classmethod
    def classify_variant(
        cls,
        value: Any,
    ) -> str:
        """Return the validated GWAS source variant class."""

        return (
            GWASValidator
            .classify_variant_descriptor(
                value
            )
        )

    # ==================================================================
    # Canonical rsID
    # ==================================================================

    @classmethod
    def standardize_canonical_rsid(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Derive an unambiguous canonical rsID when possible.

        Only a single-rsID source representation produces a canonical
        rsID. Multi-rsID and interaction descriptors are intentionally
        not collapsed into one identifier.
        """

        variant_class = (
            cls.classify_variant(
                value
            )
        )

        if (
            variant_class
            != "single_rsid"
        ):
            return {
                "canonical_rsid":
                    None,
                "canonical_rsid_usable":
                    False,
                "canonical_rsid_source":
                    (
                        "not_single_rsid"
                    ),
            }

        result = (
            StandardizationUtils
            .standardize_rsid(
                value
            )
        )

        return {
            "canonical_rsid":
                result.standardized_value,
            "canonical_rsid_usable":
                result.usable,
            "canonical_rsid_source":
                result.source,
        }

    # ==================================================================
    # Interaction utilities
    # ==================================================================

    @classmethod
    def split_interaction(
        cls,
        value: Any,
    ) -> list[str]:
        """Split source interaction fields such as '6 x 6'."""

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return []

        text = str(
            value
        ).strip()

        return [
            token.strip()
            for token
            in cls.INTERACTION_SEPARATOR
            .split(
                text
            )
            if token.strip()
        ]

    @classmethod
    def split_multi_rsid(
        cls,
        value: Any,
    ) -> list[str]:
        """Extract rsIDs from a validated multi-rsID representation."""

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return []

        text = str(
            value
        ).strip()

        tokens = [
            token
            for token
            in cls.MULTI_RSID_SEPARATOR
            .split(
                text
            )
            if token
        ]

        standardized: list[str] = []

        for token in tokens:

            result = (
                StandardizationUtils
                .standardize_rsid(
                    token
                )
            )

            if result.usable:
                standardized.append(
                    result.standardized_value
                )

        return standardized

    # ==================================================================
    # Chromosome
    # ==================================================================

    @classmethod
    def standardize_chromosome(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Standardize a GWAS chromosome field.

        Single chromosome:
            chr6 -> 6

        Interaction:
            6 x 6 -> ["6", "6"]

        Interaction fields are preserved separately and are not reduced
        to a single coordinate.
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "chromosome":
                    None,
                "interaction_chromosomes":
                    None,
                "chromosome_usable":
                    False,
                "chromosome_source":
                    "missing",
            }

        interaction_parts = (
            cls.split_interaction(
                value
            )
        )

        if (
            len(
                interaction_parts
            )
            > 1
        ):
            standardized_parts: list[
                str
            ] = []

            for part in interaction_parts:

                result = (
                    StandardizationUtils
                    .standardize_chromosome(
                        part
                    )
                )

                if not result.usable:
                    return {
                        "chromosome":
                            None,
                        "interaction_chromosomes":
                            None,
                        "chromosome_usable":
                            False,
                        "chromosome_source":
                            (
                                "unsupported_interaction"
                            ),
                    }

                standardized_parts.append(
                    result.standardized_value
                )

            return {
                "chromosome":
                    None,
                "interaction_chromosomes":
                    standardized_parts,
                "chromosome_usable":
                    False,
                "chromosome_source":
                    "interaction",
            }

        result = (
            StandardizationUtils
            .standardize_chromosome(
                value
            )
        )

        return {
            "chromosome":
                result.standardized_value,
            "interaction_chromosomes":
                None,
            "chromosome_usable":
                result.usable,
            "chromosome_source":
                result.source,
        }

    # ==================================================================
    # Position
    # ==================================================================

    @classmethod
    def standardize_position(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """Standardize single or interaction genomic positions."""

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "position":
                    None,
                "interaction_positions":
                    None,
                "position_usable":
                    False,
                "position_source":
                    "missing",
            }

        interaction_parts = (
            cls.split_interaction(
                value
            )
        )

        if (
            len(
                interaction_parts
            )
            > 1
        ):
            standardized_parts: list[
                int
            ] = []

            for part in interaction_parts:

                result = (
                    StandardizationUtils
                    .standardize_position(
                        part
                    )
                )

                if not result.usable:
                    return {
                        "position":
                            None,
                        "interaction_positions":
                            None,
                        "position_usable":
                            False,
                        "position_source":
                            (
                                "unsupported_interaction"
                            ),
                    }

                standardized_parts.append(
                    result.standardized_value
                )

            return {
                "position":
                    None,
                "interaction_positions":
                    standardized_parts,
                "position_usable":
                    False,
                "position_source":
                    "interaction",
            }

        result = (
            StandardizationUtils
            .standardize_position(
                value
            )
        )

        return {
            "position":
                result.standardized_value,
            "interaction_positions":
                None,
            "position_usable":
                result.usable,
            "position_source":
                result.source,
        }

    # ==================================================================
    # P-value
    # ==================================================================

    @classmethod
    def standardize_p_value(
        cls,
        p_value: Any,
        pvalue_mlog: Any,
    ) -> dict[str, Any]:
        """
        Standardize GWAS association P-values.

        Strategy
        --------
        1. Missing P-VALUE:
             standardized value remains missing.

        2. P-VALUE > 0:
             reported value is used directly.

        3. P-VALUE == 0 and valid PVALUE_MLOG:
             attempt reconstruction using 10 ** (-PVALUE_MLOG).

        4. P-VALUE == 0 but the reconstructed float underflows to zero:
             p_value_standardized remains zero, while the mlog value and
             provenance preserve the available significance information.

        Raw source values are never replaced.
        """

        if (
            StandardizationUtils
            .is_missing(
                p_value
            )
        ):
            return {
                "p_value_standardized":
                    None,
                "p_value_source":
                    "missing",
                "p_value_usable":
                    False,
                "p_value_reconstructed":
                    False,
            }

        probability = (
            StandardizationUtils
            .standardize_probability(
                p_value,
                allow_zero=True,
                allow_missing=False,
            )
        )

        if not probability.usable:
            return {
                "p_value_standardized":
                    None,
                "p_value_source":
                    "unsupported",
                "p_value_usable":
                    False,
                "p_value_reconstructed":
                    False,
            }

        numeric_p = float(
            probability
            .standardized_value
        )

        if numeric_p > 0:
            return {
                "p_value_standardized":
                    numeric_p,
                "p_value_source":
                    "reported_p",
                "p_value_usable":
                    True,
                "p_value_reconstructed":
                    False,
            }

        mlog_result = (
            StandardizationUtils
            .standardize_numeric(
                pvalue_mlog,
                allow_missing=True,
            )
        )

        if (
            not mlog_result.usable
            or mlog_result.standardized_value
            is None
        ):
            return {
                "p_value_standardized":
                    0.0,
                "p_value_source":
                    "reported_zero",
                "p_value_usable":
                    False,
                "p_value_reconstructed":
                    False,
            }

        mlog = float(
            mlog_result
            .standardized_value
        )

        if mlog < 0:
            return {
                "p_value_standardized":
                    0.0,
                "p_value_source":
                    "reported_zero",
                "p_value_usable":
                    False,
                "p_value_reconstructed":
                    False,
            }

        try:
            reconstructed = (
                10.0 ** (
                    -mlog
                )
            )
        except (
            OverflowError,
            ValueError,
        ):
            reconstructed = 0.0

        if (
            reconstructed > 0
            and math.isfinite(
                reconstructed
            )
        ):
            return {
                "p_value_standardized":
                    reconstructed,
                "p_value_source":
                    "reconstructed_from_mlog",
                "p_value_usable":
                    True,
                "p_value_reconstructed":
                    True,
            }

        return {
            "p_value_standardized":
                0.0,
            "p_value_source":
                "mlog_below_float_precision",
            "p_value_usable":
                False,
            "p_value_reconstructed":
                False,
        }

    # ==================================================================
    # Risk allele
    # ==================================================================

    @classmethod
    def standardize_risk_allele(
        cls,
        value: Any,
    ) -> dict[str, Any]:
        """
        Extract risk allele from STRONGEST SNP-RISK ALLELE.

        Example:
            rs123-A
                -> risk_allele = A

        Unknown allele:
            rs123-?
                -> risk_allele = None
        """

        if (
            StandardizationUtils
            .is_missing(
                value
            )
        ):
            return {
                "risk_allele":
                    None,
                "risk_allele_usable":
                    False,
                "risk_allele_source":
                    "missing",
            }

        text = str(
            value
        ).strip()

        match = (
            cls.STRONGEST_RISK_ALLELE_PATTERN
            .fullmatch(
                text
            )
        )

        if match is None:
            return {
                "risk_allele":
                    None,
                "risk_allele_usable":
                    False,
                "risk_allele_source":
                    "unparsed",
            }

        allele = (
            match.group(
                "allele"
            )
            .strip()
            .upper()
        )

        if allele == "?":
            return {
                "risk_allele":
                    None,
                "risk_allele_usable":
                    False,
                "risk_allele_source":
                    "unknown",
            }

        return {
            "risk_allele":
                allele,
            "risk_allele_usable":
                True,
            "risk_allele_source":
                "strongest_snp_risk_allele",
        }

    # ==================================================================
    # Single-row standardization
    # ==================================================================

    @classmethod
    def standardize_record(
        cls,
        record: dict[str, Any],
    ) -> dict[str, Any]:
        """
        Standardize one validated GWAS Catalog association record.

        Raw source fields used for standardization are retained in
        explicit *_raw fields.
        """

        source_variant = (
            record.get(
                "SNPS"
            )
        )

        variant_class = (
            cls.classify_variant(
                source_variant
            )
        )

        rsid = (
            cls
            .standardize_canonical_rsid(
                source_variant
            )
        )

        chromosome = (
            cls
            .standardize_chromosome(
                record.get(
                    "CHR_ID"
                )
            )
        )

        position = (
            cls
            .standardize_position(
                record.get(
                    "CHR_POS"
                )
            )
        )

        pvalue = (
            cls
            .standardize_p_value(
                record.get(
                    "P-VALUE"
                ),
                record.get(
                    "PVALUE_MLOG"
                ),
            )
        )

        risk_allele = (
            cls
            .standardize_risk_allele(
                record.get(
                    "STRONGEST SNP-RISK ALLELE"
                )
            )
        )

        single_variant = (
            variant_class
            not in {
                "multi_rsid",
                "snp_interaction",
            }
        )

        coordinate_usable = (
            chromosome[
                "chromosome_usable"
            ]
            and position[
                "position_usable"
            ]
            and single_variant
        )

        variant_usable = (
            variant_class
            in {
                "single_rsid",
                "coordinate_variant",
                "allelic_variant",
                "coordinate_indel",
                "rsid_allelic_variant",
            }
        )

        if (
            variant_usable
            and coordinate_usable
            and pvalue[
                "p_value_usable"
            ]
        ):
            standardization_status = (
                "STANDARDIZED"
            )

        elif (
            variant_usable
            or coordinate_usable
        ):
            standardization_status = (
                "PARTIAL"
            )

        else:
            standardization_status = (
                "PRESERVED"
            )

        result = {
            "study_accession":
                record.get(
                    "STUDY ACCESSION"
                ),

            "pubmed_id":
                record.get(
                    "PUBMEDID"
                ),

            "disease_trait":
                record.get(
                    "DISEASE/TRAIT"
                ),

            "mapped_trait":
                record.get(
                    "MAPPED_TRAIT"
                ),

            "mapped_trait_uri":
                record.get(
                    "MAPPED_TRAIT_URI"
                ),

            "source_variant":
                source_variant,

            "variant_class":
                variant_class,

            **rsid,

            "chromosome_raw":
                record.get(
                    "CHR_ID"
                ),

            **chromosome,

            "position_raw":
                record.get(
                    "CHR_POS"
                ),

            **position,

            "risk_allele_raw":
                record.get(
                    "STRONGEST SNP-RISK ALLELE"
                ),

            **risk_allele,

            "p_value_raw":
                record.get(
                    "P-VALUE"
                ),

            "pvalue_mlog_raw":
                record.get(
                    "PVALUE_MLOG"
                ),

            **pvalue,

            "mapped_gene":
                record.get(
                    "MAPPED_GENE"
                ),

            "reported_gene":
                record.get(
                    "REPORTED GENE(S)"
                ),

            "coordinate_usable":
                bool(
                    coordinate_usable
                ),

            "variant_usable":
                bool(
                    variant_usable
                ),

            "standardization_status":
                standardization_status,
        }

        return result

    # ==================================================================
    # DataFrame standardization
    # ==================================================================

    @classmethod
    def standardize_dataframe(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Standardize a validated GWAS Catalog DataFrame.

        This method performs no file I/O.
        """

        records = (
            dataframe
            .to_dict(
                orient="records"
            )
        )

        standardized = [
            cls.standardize_record(
                record
            )
            for record
            in records
        ]

        return pd.DataFrame(
            standardized
        )