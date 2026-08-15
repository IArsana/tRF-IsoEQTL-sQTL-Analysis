"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/audit_whole_catalog_proxies.py

Description:
    Runner for M5.3B whole-GWAS-Catalog proxy audit.

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

from pcatrfqtl.analysis.m5.whole_catalog_proxy_audit import (
    M53BWholeCatalogProxyAudit,
)
from pcatrfqtl.io.parquet import (
    read_parquet,
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


class M53BWholeCatalogProxyAuditRunner:
    """Execute M5.3B for one LD population."""

    def __init__(
        self,
        *,
        normalized_ld_path: str | Path,
        gwas_directory: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
        population: str,
        secondary_r2_threshold: float = 0.5,
    ) -> None:

        self.normalized_ld_path = Path(
            normalized_ld_path
        )

        self.gwas_directory = Path(
            gwas_directory
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

        self.population = (
            str(
                population
            )
            .strip()
            .upper()
        )

        self.secondary_r2_threshold = (
            secondary_r2_threshold
        )

    @property
    def output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / (
                "whole_catalog_proxy_audit_"
                f"{self.population}.parquet"
            )
        )

    @property
    def qc_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / (
                "m5_3b_whole_catalog_proxy_audit_"
                f"{self.population}_summary.json"
            )
        )

    def _discover_gwas_parts(
        self,
    ) -> list[Path]:
        """Discover standardized GWAS Parquet parts."""

        if not self.gwas_directory.exists():

            raise FileNotFoundError(
                "Standardized GWAS directory not found: "
                f"{self.gwas_directory}"
            )

        paths = sorted(
            self.gwas_directory.rglob(
                "*.parquet"
            )
        )

        if not paths:

            raise RuntimeError(
                "No GWAS Parquet files found under: "
                f"{self.gwas_directory}"
            )

        return paths

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute the M5.3B whole-Catalog audit."""

        if not self.normalized_ld_path.exists():

            raise FileNotFoundError(
                "Normalized LD evidence not found: "
                f"{self.normalized_ld_path}"
            )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        gwas_paths = (
            self._discover_gwas_parts()
        )

        logger.info(
            "Starting M5.3B whole-Catalog proxy audit [%s].",
            self.population,
        )

        logger.info(
            "GWAS Parquet parts discovered: %d",
            len(
                gwas_paths
            ),
        )

        ld_evidence = read_parquet(
            self.normalized_ld_path
        )

        (
            result,
            report,
        ) = (
            M53BWholeCatalogProxyAudit
            .build(
                ld_evidence=ld_evidence,
                gwas_paths=gwas_paths,
                population=self.population,
                secondary_r2_threshold=(
                    self.secondary_r2_threshold
                ),
            )
        )

        if not result.empty:

            write_parquet(
                result,
                self.output_path,
                index=False,
            )

        report[
            "output"
        ] = (
            str(
                self.output_path
            )
            if not result.empty
            else None
        )

        report[
            "normalized_ld_input"
        ] = str(
            self.normalized_ld_path
        )

        report[
            "gwas_directory"
        ] = str(
            self.gwas_directory
        )

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
            "M5.3B [%s] complete.",
            self.population,
        )

        logger.info(
            "Whole-Catalog matches: %d",
            report[
                "summary"
            ][
                "whole_catalog_match_rows"
            ],
        )

        logger.info(
            "Unique proxies found: %d",
            report[
                "summary"
            ][
                "unique_proxy_rsids_found_in_catalog"
            ],
        )

        logger.info(
            "Prostate-text matches: %d",
            report[
                "summary"
            ][
                "prostate_text_match_rows"
            ],
        )

        logger.info(
            "Prostate-cancer-text matches: %d",
            report[
                "summary"
            ][
                "prostate_cancer_text_match_rows"
            ],
        )

        return report