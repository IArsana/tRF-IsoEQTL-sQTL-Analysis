"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/discover_locus_sumstats.py

Description:
    M5.3C.3A runner for remote GWAS summary-statistics resource discovery.

    The runner reads the locked M5.3C.3 study selection configuration,
    discovers official GWAS Catalog harmonised resources, and writes a
    retrieval manifest.

    No genome-wide summary-statistics file is downloaded.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pcatrfqtl.analysis.m5.gwas_locus_remote import (
    discover_study_files,
)
from pcatrfqtl.io.parquet import (
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


class M53C3RemoteDiscoveryRunner:
    """Build M5.3C.3 locus-retrieval remote manifest."""

    def __init__(
        self,
        *,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.config_path = Path(
            config_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    @property
    def output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_locus_remote_manifest.parquet"
        )

    @property
    def qc_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_3a_remote_sumstats_discovery.json"
        )

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load locked locus-study configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M5.3C.3 config not found: {self.config_path}"
            )

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            config = yaml.safe_load(
                handle
            )

        if not isinstance(
            config,
            dict,
        ):

            raise ValueError(
                "M5.3C.3 configuration must be a mapping."
            )

        return config

    def run(
        self,
    ) -> dict[str, Any]:
        """Discover remote resources for all selected studies."""

        config = self._load_config()

        studies = config.get(
            "studies",
            [],
        )

        enabled_studies = [
            study
            for study
            in studies
            if study.get(
                "enabled",
                True,
            )
        ]

        if not enabled_studies:

            raise RuntimeError(
                "M5.3C.3 configuration contains no enabled studies."
            )

        records: list[
            dict[str, Any]
        ] = []

        for study in enabled_studies:

            accession = str(
                study[
                    "study_accession"
                ]
            ).upper()

            logger.info(
                "Discovering summary-statistics resources for %s.",
                accession,
            )

            try:

                remote = discover_study_files(
                    accession
                )

                status = (
                    "READY_FOR_INDEXED_LOCUS_QUERY"
                    if (
                        remote.harmonised_available
                        and remote.tabix_available
                    )
                    else (
                        "HARMONISED_WITHOUT_INDEX"
                        if remote.harmonised_available
                        else "NO_HARMONISED_FILE"
                    )
                )

                error = None

            except Exception as exc:

                remote = None

                status = (
                    "REMOTE_DISCOVERY_ERROR"
                )

                error = (
                    f"{type(exc).__name__}: {exc}"
                )

                logger.warning(
                    "Remote discovery failed for %s: %s",
                    accession,
                    error,
                )

            records.append(
                {
                    "study_accession":
                        accession,

                    "role":
                        study.get(
                            "role"
                        ),

                    "ancestry_role":
                        study.get(
                            "ancestry_role"
                        ),

                    "phenotype":
                        study.get(
                            "phenotype"
                        ),

                    "study_directory_url":
                        (
                            remote.study_directory_url
                            if remote
                            else None
                        ),

                    "harmonised_directory_url":
                        (
                            remote.harmonised_directory_url
                            if remote
                            else None
                        ),

                    "harmonised_sumstats_url":
                        (
                            remote.harmonised_sumstats_url
                            if remote
                            else None
                        ),

                    "harmonised_index_url":
                        (
                            remote.harmonised_index_url
                            if remote
                            else None
                        ),

                    "harmonised_metadata_url":
                        (
                            remote.harmonised_metadata_url
                            if remote
                            else None
                        ),

                    "raw_sumstats_url":
                        (
                            remote.raw_sumstats_url
                            if remote
                            else None
                        ),

                    "raw_metadata_url":
                        (
                            remote.raw_metadata_url
                            if remote
                            else None
                        ),

                    "harmonised_available":
                        (
                            remote.harmonised_available
                            if remote
                            else False
                        ),

                    "tabix_available":
                        (
                            remote.tabix_available
                            if remote
                            else False
                        ),

                    "retrieval_status":
                        status,

                    "remote_discovery_error":
                        error,
                }
            )

        dataframe = pd.DataFrame(
            records
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        write_parquet(
            dataframe,
            self.output_path,
            index=False,
        )

        status_counts = (
            dataframe[
                "retrieval_status"
            ]
            .value_counts()
            .to_dict()
        )

        report = {
            "milestone":
                "M5.3C.3A",

            "stage":
                "remote_sumstats_discovery",

            "policy": {
                "selected_studies_only":
                    True,

                "harmonised_preferred":
                    True,

                "tabix_index_preferred":
                    True,

                "genome_wide_sumstats_downloaded":
                    False,

                "genome_build_assumed":
                    False,

                "coordinate_query_performed":
                    False,

                "colocalization_performed":
                    False,
            },

            "summary": {
                "configured_studies":
                    len(
                        studies
                    ),

                "enabled_studies":
                    len(
                        enabled_studies
                    ),

                "harmonised_available":
                    int(
                        dataframe[
                            "harmonised_available"
                        ]
                        .fillna(
                            False
                        )
                        .sum()
                    ),

                "tabix_available":
                    int(
                        dataframe[
                            "tabix_available"
                        ]
                        .fillna(
                            False
                        )
                        .sum()
                    ),

                "ready_for_indexed_query":
                    int(
                        (
                            dataframe[
                                "retrieval_status"
                            ]
                            == "READY_FOR_INDEXED_LOCUS_QUERY"
                        )
                        .sum()
                    ),

                "remote_discovery_errors":
                    int(
                        (
                            dataframe[
                                "retrieval_status"
                            ]
                            == "REMOTE_DISCOVERY_ERROR"
                        )
                        .sum()
                    ),
            },

            "status_counts":
                {
                    str(key):
                        int(
                            value
                        )
                    for key, value
                    in status_counts.items()
                },

            "candidate_leads":
                config.get(
                    "candidate_leads",
                    []
                ),

            "locus_window_bp":
                config.get(
                    "locus_window_bp"
                ),

            "output":
                str(
                    self.output_path
                ),
        }

        with self.qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M5.3C.3A complete: %d/%d studies ready for indexed query.",
            report[
                "summary"
            ][
                "ready_for_indexed_query"
            ],
            len(
                enabled_studies
            ),
        )

        return report