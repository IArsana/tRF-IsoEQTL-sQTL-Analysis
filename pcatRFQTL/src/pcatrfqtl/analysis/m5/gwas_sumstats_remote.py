"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_sumstats_remote.py

Description:
    Remote GWAS Catalog summary-statistics discovery helpers for M5.3C.

    This module checks official GWAS Catalog FTP directories and discovers
    available raw or harmonised summary-statistics files for selected studies.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import html.parser
import urllib.error
import urllib.request
from dataclasses import dataclass

from pcatrfqtl.analysis.m5.gwas_sumstats_ftp import (
    gcst_directory_url,
)


@dataclass(frozen=True)
class GWASSumstatsRemoteStudy:
    """Remote summary-statistics availability record."""

    study_accession: str

    study_directory_url: str

    available: bool

    harmonised_files: tuple[str, ...]

    raw_sumstats_files: tuple[str, ...]

    metadata_files: tuple[str, ...]

    preferred_file_url: str | None

    preferred_file_type: str | None


class _LinkParser(
    html.parser.HTMLParser
):
    """Minimal HTML directory listing parser."""

    def __init__(
        self,
    ) -> None:

        super().__init__()

        self.links: list[str] = []

    def handle_starttag(
        self,
        tag: str,
        attrs: list[
            tuple[
                str,
                str | None,
            ]
        ],
    ) -> None:

        if tag.lower() != "a":
            return

        for key, value in attrs:

            if (
                key.lower()
                == "href"
                and value
            ):

                self.links.append(
                    value
                )


def _list_directory(
    url: str,
    *,
    timeout: int = 60,
) -> list[str]:
    """List hyperlinks from an FTP-over-HTTPS directory index."""

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "pcatRFQTL/1.0 "
                "M5.3C-sumstats-discovery",
        },
    )

    try:

        with urllib.request.urlopen(
            request,
            timeout=timeout,
        ) as response:

            content = (
                response.read()
                .decode(
                    "utf-8",
                    errors="replace",
                )
            )

    except urllib.error.HTTPError as exc:

        if exc.code == 404:
            return []

        raise RuntimeError(
            f"GWAS Catalog HTTP error: "
            f"{exc.code} {url}"
        ) from exc

    except urllib.error.URLError as exc:

        raise RuntimeError(
            f"GWAS Catalog network error for "
            f"{url}: {exc.reason}"
        ) from exc

    parser = _LinkParser()

    parser.feed(
        content
    )

    return parser.links


def inspect_remote_study(
    study_accession: str,
) -> GWASSumstatsRemoteStudy:
    """Inspect one GWAS Catalog sumstats accession directory."""

    directory_url = gcst_directory_url(
        study_accession
    )

    root_links = _list_directory(
        directory_url
    )

    if not root_links:

        return GWASSumstatsRemoteStudy(
            study_accession=study_accession,
            study_directory_url=directory_url,
            available=False,
            harmonised_files=(),
            raw_sumstats_files=(),
            metadata_files=(),
            preferred_file_url=None,
            preferred_file_type=None,
        )

    metadata_files = tuple(
        sorted(
            link
            for link
            in root_links
            if (
                link.endswith(".yaml")
                or link.endswith(".yml")
            )
        )
    )

    raw_sumstats_files = tuple(
        sorted(
            link
            for link
            in root_links
            if (
                link.endswith(".tsv")
                or link.endswith(".tsv.gz")
            )
            and ".h.tsv" not in link
            and ".h.tsv.gz" not in link
        )
    )

    harmonised_directory = (
        directory_url
        + "harmonised/"
    )

    harmonised_links = (
        _list_directory(
            harmonised_directory
        )
        if any(
            link.rstrip("/")
            == "harmonised"
            for link
            in root_links
        )
        else []
    )

    harmonised_files = tuple(
        sorted(
            link
            for link
            in harmonised_links
            if (
                link.endswith(".h.tsv")
                or link.endswith(".h.tsv.gz")
            )
        )
    )

    preferred_file_url: str | None = None
    preferred_file_type: str | None = None

    if harmonised_files:

        preferred_file_url = (
            harmonised_directory
            + harmonised_files[0]
        )

        preferred_file_type = (
            "HARMONISED"
        )

    elif raw_sumstats_files:

        preferred_file_url = (
            directory_url
            + raw_sumstats_files[0]
        )

        preferred_file_type = (
            "RAW_SUMSTATS"
        )

    return GWASSumstatsRemoteStudy(
        study_accession=study_accession,
        study_directory_url=directory_url,
        available=(
            preferred_file_url
            is not None
        ),
        harmonised_files=harmonised_files,
        raw_sumstats_files=raw_sumstats_files,
        metadata_files=metadata_files,
        preferred_file_url=preferred_file_url,
        preferred_file_type=preferred_file_type,
    )