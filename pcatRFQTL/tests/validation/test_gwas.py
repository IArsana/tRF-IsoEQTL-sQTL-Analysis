"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/validation/test_gwas.py

Description:
    Unit tests for GWAS Catalog field-level validation and variant
    descriptor classification utilities.

    Tests cover:
        - missing-value handling
        - strict dbSNP validation
        - heterogeneous GWAS variant descriptors
        - coordinate-based variants
        - SNP-by-SNP interactions
        - HLA allele descriptors
        - array/probe identifiers
        - structural-variant identifiers
        - chromosome validation
        - interaction-aware chromosome fields
        - genomic-position validation
        - interaction-aware position fields
        - P-value validation
        - zero P-value policy
        - generic numeric validation
        - stable validation-rule identifiers

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

import numpy as np

from pcatrfqtl.validation.validators.gwas import (
    GWASValidator,
)


def create_validator() -> GWASValidator:
    """Create a GWAS Catalog validator."""

    return GWASValidator()


# ===========================================================================
# Missing values
# ===========================================================================


def test_is_missing() -> None:
    """Common missing-value representations should be detected."""

    validator = create_validator()

    assert validator._is_missing(None)
    assert validator._is_missing(float("nan"))
    assert validator._is_missing(np.nan)
    assert validator._is_missing("")
    assert validator._is_missing(" ")
    assert validator._is_missing("NA")
    assert validator._is_missing("nan")
    assert validator._is_missing("None")
    assert validator._is_missing("null")

    assert not validator._is_missing("rs12345")
    assert not validator._is_missing(0)
    assert not validator._is_missing(1)


# ===========================================================================
# Strict dbSNP validation
# ===========================================================================


def test_valid_single_snp_ids() -> None:
    """Canonical dbSNP identifiers should pass strict validation."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "rs12345",
            "rs7412",
            "RS429358",
            "rs100000000",
        ],
        "SNPS",
    )

    assert result == []


def test_valid_multi_snp_ids() -> None:
    """Conventional multi-rsID representations should pass."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "rs12345; rs67890",
            "rs12345,rs67890",
            "rs12345 rs67890",
            "rs12345;rs67890;rs11111",
            "rs12345, rs67890, rs11111",
        ],
        "SNPS",
    )

    assert result == []


def test_invalid_strict_snp_ids() -> None:
    """Non-rsID variant representations should fail strict validation."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "12345",
            "rsABC",
            "chr1:100",
            "rs12345 x rs67890",
            "HLA-DRB1*14:04",
            "exm130158",
        ],
        "SNPS",
    )

    assert len(result) == 6


def test_missing_snp_ids_allowed() -> None:
    """Missing strict SNP identifiers should be allowed when configured."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            None,
            "",
            np.nan,
            "rs12345",
        ],
        "SNPS",
        allow_missing=True,
    )

    assert result == []


def test_missing_snp_ids_not_allowed() -> None:
    """Missing strict SNP identifiers should fail when required."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            None,
            "",
            np.nan,
            "rs12345",
        ],
        "SNPS",
        allow_missing=False,
    )

    assert len(result) == 3


def test_strict_snp_error_contains_rule_id() -> None:
    """Strict SNP validation errors should contain a stable rule ID."""

    validator = create_validator()

    result = validator.validate_snp_ids(
        [
            "invalid",
        ],
        "SNPS",
    )

    assert len(result) == 1
    assert "GWAS-SNP-001" in result[0]
    assert "SNPS[0]" in result[0]


# ===========================================================================
# Variant descriptor classification
# ===========================================================================


def test_classify_single_rsid() -> None:
    """A canonical rsID should classify as single_rsid."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "rs12345"
        )
    )

    assert result == "single_rsid"


def test_classify_multi_rsid() -> None:
    """Multiple rsIDs should classify as multi_rsid."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "rs12345; rs67890"
        )
    )

    assert result == "multi_rsid"


def test_classify_snp_interaction() -> None:
    """SNP-by-SNP interaction descriptors should be recognized."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "rs2853952 x rs45533135"
        )
    )

    assert result == "snp_interaction"


def test_classify_coordinate_variant() -> None:
    """Coordinate-only variant descriptors should be recognized."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "chr12:9098995"
        )
    )

    assert result == "coordinate_variant"


def test_classify_colon_allelic_variant() -> None:
    """Colon-delimited allele-specific variants should be recognized."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "9:139737088:G:A"
        )
    )

    assert result == "allelic_variant"


def test_classify_underscore_allelic_variant() -> None:
    """Underscore-delimited allele-specific variants should be recognized."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "15:79372468:A_ATG"
        )
    )

    assert result == "allelic_variant"


def test_classify_coordinate_indel() -> None:
    """Coordinate-based insertion/deletion descriptors should pass."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "chr2:134471945:D"
        )
    )

    assert result == "coordinate_indel"


def test_classify_rsid_allelic_variant() -> None:
    """rsID plus allele representation should be recognized."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "rs11639303_G_A"
        )
    )

    assert result == "rsid_allelic_variant"


def test_classify_hla_allele() -> None:
    """HLA allele descriptors should be recognized."""

    validator = create_validator()

    values = [
        "HLA-DRB1*14:04",
        "HLA-B*07:02",
        "B*08:01",
        "DQA1*03:01",
        "HLA-A*3101",
    ]

    for value in values:
        assert (
            validator
            .classify_variant_descriptor(
                value
            )
            == "hla_allele"
        )


def test_classify_array_probe() -> None:
    """Array/probe identifiers should be recognized."""

    validator = create_validator()

    values = [
        "exm130158",
        "kgp4136779",
        "i4000416",
    ]

    for value in values:
        assert (
            validator
            .classify_variant_descriptor(
                value
            )
            == "array_probe"
        )


def test_classify_structural_variant() -> None:
    """Structural-variant identifiers should be recognized."""

    validator = create_validator()

    values = [
        "nsv831124",
        "esv3596105",
    ]

    for value in values:
        assert (
            validator
            .classify_variant_descriptor(
                value
            )
            == "structural_variant"
        )


def test_classify_gene_or_source_label() -> None:
    """Simple gene/source labels should remain valid source descriptors."""

    validator = create_validator()

    values = [
        "SLC22A8",
        "LYPD6B",
        "PPP5C",
        "downstreamRASGEF1B",
    ]

    for value in values:
        assert (
            validator
            .classify_variant_descriptor(
                value
            )
            == "gene_or_source_label"
        )


def test_classify_other_variant_descriptor() -> None:
    """Unrecognized non-empty representations should classify as other."""

    result = (
        GWASValidator
        .classify_variant_descriptor(
            "Val86 HLA-DRB1"
        )
    )

    assert result == "other"


def test_classify_missing_variant_descriptor() -> None:
    """Missing values should classify as missing."""

    assert (
        GWASValidator
        .classify_variant_descriptor(
            None
        )
        == "missing"
    )


# ===========================================================================
# Variant descriptor validation
# ===========================================================================


def test_valid_variant_descriptors() -> None:
    """
    Heterogeneous GWAS Catalog source descriptors should pass structural
    validation.
    """

    validator = create_validator()

    result = (
        validator
        .validate_variant_descriptors(
            [
                "rs12345",
                "rs12345; rs67890",
                "rs2853952 x rs45533135",
                "chr12:9098995",
                "9:139737088:G:A",
                "15:79372468:A_ATG",
                "chr2:134471945:D",
                "HLA-DRB1*14:04",
                "exm130158",
                "nsv831124",
                "SLC22A8",
                "Val86 HLA-DRB1",
            ],
            "SNPS",
        )
    )

    assert result == []


def test_missing_variant_descriptors_allowed() -> None:
    """Missing variant descriptors may be retained when allowed."""

    validator = create_validator()

    result = (
        validator
        .validate_variant_descriptors(
            [
                None,
                "",
                np.nan,
                "rs12345",
            ],
            "SNPS",
            allow_missing=True,
        )
    )

    assert result == []


def test_missing_variant_descriptors_not_allowed() -> None:
    """Missing variant descriptors should fail when explicitly required."""

    validator = create_validator()

    result = (
        validator
        .validate_variant_descriptors(
            [
                None,
                "",
                np.nan,
                "rs12345",
            ],
            "SNPS",
            allow_missing=False,
        )
    )

    assert len(result) == 3


# ===========================================================================
# Chromosome validation
# ===========================================================================


def test_valid_chromosomes() -> None:
    """Canonical chromosome identifiers should pass."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            "1",
            "2",
            "10",
            "22",
            "X",
            "Y",
            "MT",
            "M",
            "chr1",
            "chr22",
            "chrX",
            "chrY",
            "chrMT",
        ],
        "CHR_ID",
    )

    assert result == []


def test_valid_multi_chromosome_fields() -> None:
    """Comma- and semicolon-separated chromosome fields should pass."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            "1;2",
            "1,2",
            "X;Y",
            "chr1;chr2",
        ],
        "CHR_ID",
    )

    assert result == []


def test_valid_interaction_chromosomes() -> None:
    """SNP interaction chromosome fields should pass."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            "6 x 6",
            "1 x 2",
            "X x Y",
            "chr1 x chr2",
        ],
        "CHR_ID",
    )

    assert result == []


def test_invalid_chromosomes() -> None:
    """Noncanonical chromosome identifiers should fail."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            "0",
            "23",
            "24",
            "chr23",
            "Z",
            "1A",
            "chr",
        ],
        "CHR_ID",
    )

    assert len(result) == 7


def test_missing_chromosomes_allowed() -> None:
    """Missing chromosomes should be allowed for unmapped associations."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            None,
            "",
            np.nan,
            "1",
        ],
        "CHR_ID",
        allow_missing=True,
    )

    assert result == []


def test_missing_chromosomes_not_allowed() -> None:
    """Missing chromosomes should fail when mapping is required."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            None,
            "",
            np.nan,
            "1",
        ],
        "CHR_ID",
        allow_missing=False,
    )

    assert len(result) == 3


def test_invalid_chromosome_error_contains_rule_id() -> None:
    """Chromosome validation errors should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_chromosomes(
        [
            "23",
        ],
        "CHR_ID",
    )

    assert len(result) == 1
    assert "GWAS-CHR-001" in result[0]
    assert "CHR_ID[0]" in result[0]


# ===========================================================================
# Genomic positions
# ===========================================================================


def test_valid_positions() -> None:
    """Positive integer genomic coordinates should pass."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            1,
            100,
            12345678,
            "1",
            "100",
            "12345678",
            100.0,
        ],
        "CHR_POS",
    )

    assert result == []


def test_valid_multi_position_fields() -> None:
    """Comma- or semicolon-separated positions should pass."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            "100;200",
            "100,200",
            "100;200;300",
            "100, 200, 300",
        ],
        "CHR_POS",
    )

    assert result == []


def test_valid_interaction_positions() -> None:
    """SNP interaction genomic-position fields should pass."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            "31268092 x 31463174",
            "100 x 200",
            "1 x 999999",
        ],
        "CHR_POS",
    )

    assert result == []


def test_invalid_positions() -> None:
    """Invalid genomic coordinates should fail."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            0,
            -1,
            "0",
            "-100",
            "ABC",
            100.5,
            "100.5",
            math.inf,
            -math.inf,
        ],
        "CHR_POS",
    )

    assert len(result) == 9


def test_missing_positions_allowed() -> None:
    """Missing genomic positions should be accepted when permitted."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            None,
            "",
            np.nan,
            100,
        ],
        "CHR_POS",
        allow_missing=True,
    )

    assert result == []


def test_missing_positions_not_allowed() -> None:
    """Missing positions should fail when mapping is required."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            None,
            "",
            np.nan,
            100,
        ],
        "CHR_POS",
        allow_missing=False,
    )

    assert len(result) == 3


def test_invalid_position_error_contains_rule_id() -> None:
    """Position validation errors should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_positions(
        [
            -100,
        ],
        "CHR_POS",
    )

    assert len(result) == 1
    assert "GWAS-POS-001" in result[0]
    assert "CHR_POS[0]" in result[0]


# ===========================================================================
# P-values
# ===========================================================================


def test_valid_p_values() -> None:
    """GWAS Catalog P-values in 0 <= P <= 1 should pass."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0,
            0.0,
            1e-300,
            5e-8,
            0.05,
            0.5,
            1,
            "0.01",
            "5e-8",
        ],
        "P-VALUE",
    )

    assert result == []


def test_zero_p_value_allowed_by_default() -> None:
    """Zero P-values should be accepted in default GWAS Catalog mode."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0,
            0.0,
        ],
        "P-VALUE",
    )

    assert result == []


def test_zero_p_value_invalid_in_strict_mode() -> None:
    """Zero P-values should fail when strict positive P-values are required."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0,
        ],
        "P-VALUE",
        allow_zero=False,
    )

    assert len(result) == 1
    assert "GWAS-PVALUE-001" in result[0]


def test_invalid_p_values() -> None:
    """Out-of-range and non-numeric P-values should fail."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0,
            -0.01,
            1.01,
            "invalid",
            math.inf,
            -math.inf,
        ],
        "P-VALUE",
    )

    # P=0 is valid in default GWAS Catalog mode.
    assert len(result) == 5


def test_missing_p_values_invalid() -> None:
    """Missing P-values should always be invalid."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            None,
            "",
            np.nan,
        ],
        "P-VALUE",
    )

    assert len(result) == 3


def test_pvalue_error_contains_rule_id() -> None:
    """Strict P-value validation errors should contain stable rule IDs."""

    validator = create_validator()

    result = validator.validate_p_values(
        [
            0,
        ],
        "P-VALUE",
        allow_zero=False,
    )

    assert len(result) == 1
    assert "GWAS-PVALUE-001" in result[0]
    assert "P-VALUE[0]" in result[0]


# ===========================================================================
# Generic numeric fields
# ===========================================================================


def test_valid_numeric_values() -> None:
    """Finite positive and negative numeric values should pass."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            -10,
            -0.5,
            0,
            1,
            10.5,
            "3.14",
            "1e-5",
        ],
        "OR-BETA",
    )

    assert result == []


def test_invalid_numeric_values() -> None:
    """Non-numeric and non-finite values should fail."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            "invalid",
            math.inf,
            -math.inf,
        ],
        "OR-BETA",
    )

    assert len(result) == 3


def test_missing_numeric_values_allowed() -> None:
    """Missing generic numeric fields should be allowed by default."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            None,
            "",
            np.nan,
            1.0,
        ],
        "OR-BETA",
        allow_missing=True,
    )

    assert result == []


def test_missing_numeric_values_not_allowed() -> None:
    """Missing numeric values should fail when required."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            None,
            "",
            np.nan,
            1.0,
        ],
        "OR-BETA",
        allow_missing=False,
    )

    assert len(result) == 3


def test_numeric_error_contains_rule_id() -> None:
    """Generic numeric errors should contain stable rule identifiers."""

    validator = create_validator()

    result = validator.validate_numeric(
        [
            "invalid",
        ],
        "OR-BETA",
    )

    assert len(result) == 1
    assert "GWAS-NUMERIC-001" in result[0]
    assert "OR-BETA[0]" in result[0]


# ===========================================================================
# Internal parsing helpers
# ===========================================================================


def test_split_single_snp_field() -> None:
    """Single rsID fields should return one token."""

    result = (
        GWASValidator
        ._split_snp_field(
            "rs12345"
        )
    )

    assert result == [
        "rs12345"
    ]


def test_split_multi_snp_field() -> None:
    """Conventional multi-rsID separators should be parsed."""

    result = (
        GWASValidator
        ._split_snp_field(
            "rs12345; rs67890, rs11111"
        )
    )

    assert result == [
        "rs12345",
        "rs67890",
        "rs11111",
    ]


def test_split_interaction_field() -> None:
    """Interaction-aware fields should split on x separators."""

    result = (
        GWASValidator
        ._split_interaction_field(
            "31268092 x 31463174"
        )
    )

    assert result == [
        "31268092",
        "31463174",
    ]