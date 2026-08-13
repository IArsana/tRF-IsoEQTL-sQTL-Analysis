"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/validation/genomic.py

Description:
    Shared utilities for parsing and validating genomic coordinates
    used throughout the PCa-tRFQTL research pipeline.

    This module provides coordinate parsing functionality for
    chromosome-position representations used by GWAS, tRFQTL,
    iso-eQTL, and sQTL datasets.

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
from dataclasses import dataclass


_COORDINATE_PATTERN = re.compile(
    r"^(?:chr)?(?P<chromosome>[0-9]+|X|Y|MT):(?P<position>[1-9][0-9]*)$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class GenomicCoordinate:
    """
    Represent a genomic coordinate.

    Attributes:
        chromosome:
            Chromosome identifier without the ``chr`` prefix.
        position:
            One-based genomic position.
    """

    chromosome: str
    position: int


def parse_genomic_coordinate(
    value: str,
) -> GenomicCoordinate:
    """
    Parse a genomic coordinate.

    Supported formats include:

        chr1:147308207
        1:147308207
        chrX:123456
        X:123456

    Args:
        value:
            Genomic coordinate string.

    Returns:
        Parsed GenomicCoordinate.

    Raises:
        ValueError:
            If the coordinate format is invalid.
    """

    if not isinstance(value, str):
        raise ValueError(
            "Genomic coordinate must be a string."
        )

    value = value.strip()

    match = _COORDINATE_PATTERN.fullmatch(value)

    if match is None:
        raise ValueError(
            f"Invalid genomic coordinate: {value!r}"
        )

    chromosome = match.group("chromosome").upper()
    position = int(match.group("position"))

    return GenomicCoordinate(
        chromosome=chromosome,
        position=position,
    )


def is_valid_genomic_coordinate(
    value: str,
) -> bool:
    """
    Determine whether a genomic coordinate is valid.

    Args:
        value:
            Genomic coordinate string.

    Returns:
        True if valid, otherwise False.
    """

    try:
        parse_genomic_coordinate(value)
    except ValueError:
        return False

    return True