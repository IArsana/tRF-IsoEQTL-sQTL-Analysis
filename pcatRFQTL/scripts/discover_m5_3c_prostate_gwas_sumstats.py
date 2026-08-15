"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/discover_m5_3c_prostate_gwas_sumstats.py

Description:
    M5.3C.1 prostate cancer GWAS summary-statistics discovery runner.

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
import json
import time
from pathlib import Path

import pandas as pd

from pcatrfqtl.analysis.m5.gwas_sumstats_discovery import (
    M53CProstateGWASStudyDiscovery,
)
from pcatrfqtl.analysis.m5.gwas_sumstats_remote import (
    inspect_remote_study,
)
from pcatrfqtl.io.parquet import (
    write_parquet,
)


PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Discover prostate cancer GWAS Catalog studies "
            "with downloadable summary statistics."
        )
    )

    parser.add_argument(
        "--gwas-directory",
        default=(
            "data/interim/standardized/"
            "gwas_catalog"
        ),
    )

    parser.add_argument(
        "--delay",
        type=float,
        default=0.5,
    )

    return parser.parse_args()


def main() -> None:
    """Run study discovery and remote availability audit."""

    args = parse_arguments()

    gwas_directory = (
        PROJECT_ROOT
        / args.gwas_directory
    )

    gwas_paths = sorted(
        gwas_directory.glob(
            "*.parquet"
        )
    )

    if not gwas_paths:

        raise RuntimeError(
            "No standardized GWAS Parquet parts found under "
            f"{gwas_directory}"
        )

    studies = (
        M53CProstateGWASStudyDiscovery
        .discover(
            gwas_paths
        )
    )

    print(
        "Candidate prostate cancer studies:",
        len(
            studies
        ),
    )

    remote_records: list[
        dict[str, object]
    ] = []

    for index, row in studies.iterrows():

        accession = str(
            row[
                "study_accession"
            ]
        )

        print(
            f"[{index + 1}/{len(studies)}] "
            f"{accession}"
        )

        remote = inspect_remote_study(
            accession
        )

        remote_records.append(
            {
                "study_accession":
                    accession,

                "sumstats_available":
                    remote.available,

                "preferred_file_type":
                    remote.preferred_file_type,

                "preferred_file_url":
                    remote.preferred_file_url,

                "study_directory_url":
                    remote.study_directory_url,

                "harmonised_files":
                    list(
                        remote.harmonised_files
                    ),

                "raw_sumstats_files":
                    list(
                        remote.raw_sumstats_files
                    ),

                "metadata_files":
                    list(
                        remote.metadata_files
                    ),
            }
        )

        time.sleep(
            args.delay
        )

    remote_df = pd.DataFrame(
        remote_records
    )

    result = studies.merge(
        remote_df,
        on="study_accession",
        how="left",
        validate="one_to_one",
    )

    output_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "sumstats"
    )

    qc_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "qc"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    qc_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = (
        output_directory
        / "prostate_gwas_sumstats_studies.parquet"
    )

    write_parquet(
        result,
        output_path,
        index=False,
    )

    available = (
        result[
            "sumstats_available"
        ]
        .fillna(False)
    )

    harmonised = (
        result[
            "preferred_file_type"
        ]
        .eq(
            "HARMONISED"
        )
        .fillna(False)
    )

    report = {
        "milestone":
            "M5.3C.1",

        "stage":
            "prostate_gwas_sumstats_discovery",

        "gwas_catalog_association_files":
            len(
                gwas_paths
            ),

        "candidate_prostate_studies":
            len(
                result
            ),

        "studies_with_sumstats":
            int(
                available.sum()
            ),

        "studies_with_harmonised_sumstats":
            int(
                harmonised.sum()
            ),

        "study_accessions_with_sumstats":
            sorted(
                result.loc[
                    available,
                    "study_accession",
                ]
                .astype(str)
                .tolist()
            ),

        "output":
            str(
                output_path
            ),
    }

    qc_path = (
        qc_directory
        / "m5_3c_prostate_gwas_sumstats_discovery.json"
    )

    with qc_path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            report,
            handle,
            indent=2,
        )

    print(
        json.dumps(
            report,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()