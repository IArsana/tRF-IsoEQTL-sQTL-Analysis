"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_sumstats_ftp.py

Description:
    GWAS Catalog FTP path utilities for M5.3C.

    Generates official GWAS Catalog summary-statistics accession directory
    URLs using the documented 1000-accession bucket layout.

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


GWAS_SUMSTATS_BASE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
)

GCST_PATTERN = re.compile(
    r"^GCST(\d+)$",
    flags=re.IGNORECASE,
)


def normalize_gcst(
    accession: str,
) -> str:
    """Validate and normalize one GWAS Catalog study accession."""

    text = (
        str(
            accession
        )
        .strip()
        .upper()
    )

    if GCST_PATTERN.fullmatch(
        text
    ) is None:

        raise ValueError(
            f"Invalid GWAS Catalog study accession: {accession}"
        )

    return text


def gcst_bucket(
    accession: str,
) -> str:
    """
    Return GWAS Catalog's 1000-accession FTP bucket.

    Example:
        GCST123456
        ->
        GCST123001-GCST124000
    """

    accession = normalize_gcst(
        accession
    )

    match = GCST_PATTERN.fullmatch(
        accession
    )

    if match is None:
        raise ValueError(
            f"Invalid GCST accession: {accession}"
        )

    number_text = match.group(
        1
    )

    number = int(
        number_text
    )

    width = len(
        number_text
    )

    lower = (
        (
            (number - 1)
            // 1000
        )
        * 1000
        + 1
    )

    upper = (
        lower
        + 999
    )

    lower_gcst = (
        "GCST"
        + str(
            lower
        ).zfill(
            width
        )
    )

    upper_gcst = (
        "GCST"
        + str(
            upper
        ).zfill(
            width
        )
    )

    return (
        f"{lower_gcst}-{upper_gcst}"
    )


def gcst_directory_url(
    accession: str,
) -> str:
    """Return summary-statistics FTP directory URL."""

    accession = normalize_gcst(
        accession
    )

    bucket = gcst_bucket(
        accession
    )

    return (
        f"{GWAS_SUMSTATS_BASE_URL}/"
        f"{bucket}/"
        f"{accession}/"
    )