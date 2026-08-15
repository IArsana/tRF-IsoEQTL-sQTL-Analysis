"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_locus_remote.py

Description:
    Remote GWAS Catalog summary-statistics discovery for M5.3C.3.

    This module resolves GWAS Catalog study directories and discovers
    harmonised summary-statistics files, Tabix indexes, and metadata
    without downloading genome-wide datasets.

    Important safeguards:
        - No genome-wide summary-statistics file is downloaded here.
        - Harmonised files are preferred.
        - File availability is discovered from the official GWAS Catalog.
        - Genome assembly is not assumed from study accession.
        - Coordinate querying is not performed until build compatibility
          has been established.
        - No locus association interpretation is performed.
        - No colocalization is performed.
        - No causal inference is performed.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import html
import re
import urllib.request
from dataclasses import dataclass
from urllib.parse import urljoin


GWAS_FTP_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/"
    "summary_statistics/"
)


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class RemoteGWASFiles:
    """Discovered remote files for one GWAS Catalog study."""

    study_accession: str

    study_directory_url: str

    harmonised_directory_url: str | None

    harmonised_sumstats_url: str | None

    harmonised_index_url: str | None

    harmonised_metadata_url: str | None

    raw_sumstats_url: str | None

    raw_metadata_url: str | None

    harmonised_available: bool

    tabix_available: bool


# ============================================================================
# GCST directory helpers
# ============================================================================


def _normalize_accession(
    accession: str,
) -> str:
    """Normalize and validate a GCST accession."""

    value = (
        str(accession)
        .strip()
        .upper()
    )

    if not re.fullmatch(
        r"GCST\d+",
        value,
    ):

        raise ValueError(
            f"Invalid GWAS Catalog study accession: {accession}"
        )

    return value


def accession_bucket(
    accession: str,
) -> str:
    """
    Convert a GCST accession into its GWAS Catalog thousand-study bucket.

    Example:
        GCST90274713
        ->
        GCST90274001-GCST90275000
    """

    accession = _normalize_accession(
        accession
    )

    numeric = int(
        accession.removeprefix(
            "GCST"
        )
    )

    lower = (
        (
            numeric
            - 1
        )
        // 1000
        * 1000
        + 1
    )

    upper = (
        lower
        + 999
    )

    width = len(
        accession.removeprefix(
            "GCST"
        )
    )

    return (
        f"GCST{lower:0{width}d}"
        "-"
        f"GCST{upper:0{width}d}"
    )


def study_directory_url(
    accession: str,
) -> str:
    """Return official GWAS Catalog study summary-statistics directory."""

    accession = _normalize_accession(
        accession
    )

    bucket = accession_bucket(
        accession
    )

    return (
        f"{GWAS_FTP_BASE}"
        f"{bucket}/"
        f"{accession}/"
    )


# ============================================================================
# HTTP helpers
# ============================================================================


def _fetch_text(
    url: str,
    *,
    timeout: int = 60,
) -> str:
    """Fetch UTF-8 text from one official remote resource."""

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "pcatRFQTL/1.0 "
                "M5.3C.3-locus-discovery",
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=timeout,
    ) as response:

        payload = response.read()

    return payload.decode(
        "utf-8",
        errors="replace",
    )


def _directory_links(
    url: str,
) -> list[str]:
    """Parse file links from an Apache-style directory index."""

    document = _fetch_text(
        url
    )

    hrefs = re.findall(
        r'href=["\']([^"\']+)["\']',
        document,
        flags=re.IGNORECASE,
    )

    result: list[str] = []

    for href in hrefs:

        href = html.unescape(
            href
        )

        if href.startswith(
            "?"
        ):
            continue

        if href in {
            "../",
            "/",
        }:
            continue

        result.append(
            href
        )

    return result


# ============================================================================
# Discovery
# ============================================================================


def _pick_first(
    values: list[str],
    suffixes: tuple[str, ...],
) -> str | None:
    """Return the first filename matching preferred suffix order."""

    lowered = {
        value.lower():
            value
        for value
        in values
    }

    for suffix in suffixes:

        for lower_name, original in lowered.items():

            if lower_name.endswith(
                suffix.lower()
            ):

                return original

    return None


def discover_study_files(
    accession: str,
) -> RemoteGWASFiles:
    """Discover usable summary-statistics resources for one study."""

    accession = _normalize_accession(
        accession
    )

    study_url = study_directory_url(
        accession
    )

    root_links = _directory_links(
        study_url
    )

    harmonised_link = next(
        (
            link
            for link
            in root_links
            if link.rstrip(
                "/"
            ).lower()
            == "harmonised"
        ),
        None,
    )

    raw_sumstats_name = _pick_first(
        root_links,
        (
            ".tsv.gz",
            ".tsv",
            ".txt.gz",
            ".txt",
        ),
    )

    raw_metadata_name = _pick_first(
        root_links,
        (
            "-meta.yaml",
            ".yaml",
            ".yml",
        ),
    )

    harmonised_url: str | None = None

    harmonised_sumstats_name: str | None = None

    harmonised_index_name: str | None = None

    harmonised_metadata_name: str | None = None

    if harmonised_link is not None:

        harmonised_url = urljoin(
            study_url,
            harmonised_link,
        )

        if not harmonised_url.endswith(
            "/"
        ):

            harmonised_url += "/"

        harmonised_links = _directory_links(
            harmonised_url
        )

        harmonised_sumstats_name = _pick_first(
            harmonised_links,
            (
                ".h.tsv.gz",
                ".h.tsv",
            ),
        )

        harmonised_index_name = _pick_first(
            harmonised_links,
            (
                ".h.tsv.gz.tbi",
                ".h.tsv.tbi",
                ".tbi",
            ),
        )

        harmonised_metadata_name = _pick_first(
            harmonised_links,
            (
                "-meta.yaml",
                ".yaml",
                ".yml",
            ),
        )

    harmonised_sumstats_url = (
        urljoin(
            harmonised_url,
            harmonised_sumstats_name,
        )
        if (
            harmonised_url
            and harmonised_sumstats_name
        )
        else None
    )

    harmonised_index_url = (
        urljoin(
            harmonised_url,
            harmonised_index_name,
        )
        if (
            harmonised_url
            and harmonised_index_name
        )
        else None
    )

    harmonised_metadata_url = (
        urljoin(
            harmonised_url,
            harmonised_metadata_name,
        )
        if (
            harmonised_url
            and harmonised_metadata_name
        )
        else None
    )

    return RemoteGWASFiles(
        study_accession=accession,

        study_directory_url=(
            study_url
        ),

        harmonised_directory_url=(
            harmonised_url
        ),

        harmonised_sumstats_url=(
            harmonised_sumstats_url
        ),

        harmonised_index_url=(
            harmonised_index_url
        ),

        harmonised_metadata_url=(
            harmonised_metadata_url
        ),

        raw_sumstats_url=(
            urljoin(
                study_url,
                raw_sumstats_name,
            )
            if raw_sumstats_name
            else None
        ),

        raw_metadata_url=(
            urljoin(
                study_url,
                raw_metadata_name,
            )
            if raw_metadata_name
            else None
        ),

        harmonised_available=(
            harmonised_sumstats_url
            is not None
        ),

        tabix_available=(
            harmonised_index_url
            is not None
        ),
    )