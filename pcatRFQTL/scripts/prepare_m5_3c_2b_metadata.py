"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/prepare_m5_3c_2b_metadata.py

Description:
    Download and cache official GWAS Catalog All Studies and All Ancestry
    metadata for M5.3C.2B verification.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def parse_arguments() -> argparse.Namespace:
    """Parse metadata source URLs."""

    parser = argparse.ArgumentParser(
        description=(
            "Cache official GWAS Catalog study and ancestry metadata."
        )
    )

    parser.add_argument(
        "--studies-url",
        required=True,
        help=(
            "Official GWAS Catalog All Studies v1.0.3.1 TSV URL."
        ),
    )

    parser.add_argument(
        "--ancestry-url",
        required=True,
        help=(
            "Official GWAS Catalog All Ancestry v1.0.3.1 TSV URL."
        ),
    )

    parser.add_argument(
        "--force",
        action="store_true",
    )

    return parser.parse_args()


def download(
    url: str,
    destination: Path,
    *,
    force: bool,
) -> None:
    """Download one metadata artifact."""

    if (
        destination.exists()
        and not force
    ):

        print(
            f"cache exists: {destination}"
        )

        return

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent":
                "pcatRFQTL/1.0 "
                "GWAS-metadata-verification",
        },
    )

    print(
        f"downloading: {url}"
    )

    with (
        urllib.request.urlopen(
            request,
            timeout=120,
        )
        as response,
        destination.open(
            "wb"
        )
        as output,
    ):

        while True:

            chunk = response.read(
                1024
                * 1024
            )

            if not chunk:
                break

            output.write(
                chunk
            )

    print(
        f"saved: {destination}"
    )


def main() -> None:
    """Cache official metadata."""

    args = parse_arguments()

    raw_directory = (
        PROJECT_ROOT
        / "data"
        / "raw"
        / "gwas"
        / "metadata"
    )

    studies_path = (
        raw_directory
        / "gwas_catalog_all_studies.tsv"
    )

    ancestry_path = (
        raw_directory
        / "gwas_catalog_all_ancestry.tsv"
    )

    download(
        args.studies_url,
        studies_path,
        force=args.force,
    )

    download(
        args.ancestry_url,
        ancestry_path,
        force=args.force,
    )

    print()
    print(
        "M5.3C.2B metadata preparation complete."
    )

    print(
        f"Studies : {studies_path}"
    )

    print(
        f"Ancestry: {ancestry_path}"
    )


if __name__ == "__main__":
    main()