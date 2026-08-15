"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/figures.py

Description:
    Figure-generation utilities for M3.6.

    Figures summarize the M3 prostate-cancer GWAS × tRF-QTL analysis
    using descriptive counts only.

    The figures do not imply causal relationships or LD-based overlap.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


class M3FigureBuilder:
    """Build descriptive figures for M3."""

    @staticmethod
    def dataset_flow(
        statistics: dict[str, Any],
        output_directory: Path,
    ) -> list[Path]:
        """Create M3 dataset-reduction flow figure."""

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        gwas_rows = (
            statistics[
                "gwas"
            ][
                "association_rows"
            ]
        )

        gwas_unique = (
            statistics[
                "gwas"
            ][
                "unique_eligible_rsids"
            ]
        )

        trfqtl_rows = (
            statistics[
                "trfqtl"
            ][
                "association_rows"
            ]
        )

        trfqtl_unique = (
            statistics[
                "trfqtl"
            ][
                "unique_eligible_rsids"
            ]
        )

        shared = (
            statistics[
                "direct_overlap"
            ][
                "unique_shared_rsids"
            ]
        )

        figure, axis = plt.subplots(
            figsize=(
                10,
                6,
            )
        )

        axis.axis(
            "off"
        )

        axis.text(
            0.25,
            0.85,
            (
                "Prostate cancer GWAS\n"
                f"{gwas_rows:,} association rows"
            ),
            ha="center",
            va="center",
            fontsize=12,
            bbox={
                "boxstyle": "round,pad=0.5",
            },
        )

        axis.text(
            0.25,
            0.55,
            (
                "Canonical GWAS variants\n"
                f"{gwas_unique:,} unique rsIDs"
            ),
            ha="center",
            va="center",
            fontsize=12,
            bbox={
                "boxstyle": "round,pad=0.5",
            },
        )

        axis.annotate(
            "",
            xy=(
                0.25,
                0.63,
            ),
            xytext=(
                0.25,
                0.77,
            ),
            arrowprops={
                "arrowstyle": "->",
            },
        )

        axis.text(
            0.75,
            0.85,
            (
                "PRAD tRF-QTL\n"
                f"{trfqtl_rows:,} association rows"
            ),
            ha="center",
            va="center",
            fontsize=12,
            bbox={
                "boxstyle": "round,pad=0.5",
            },
        )

        axis.text(
            0.75,
            0.55,
            (
                "tRF-QTL variants\n"
                f"{trfqtl_unique:,} unique rsIDs"
            ),
            ha="center",
            va="center",
            fontsize=12,
            bbox={
                "boxstyle": "round,pad=0.5",
            },
        )

        axis.annotate(
            "",
            xy=(
                0.75,
                0.63,
            ),
            xytext=(
                0.75,
                0.77,
            ),
            arrowprops={
                "arrowstyle": "->",
            },
        )

        axis.text(
            0.50,
            0.20,
            (
                "Exact canonical-rsID overlap\n"
                f"{shared:,} shared variants"
            ),
            ha="center",
            va="center",
            fontsize=13,
            bbox={
                "boxstyle": "round,pad=0.6",
            },
        )

        axis.annotate(
            "",
            xy=(
                0.44,
                0.29,
            ),
            xytext=(
                0.28,
                0.48,
            ),
            arrowprops={
                "arrowstyle": "->",
            },
        )

        axis.annotate(
            "",
            xy=(
                0.56,
                0.29,
            ),
            xytext=(
                0.72,
                0.48,
            ),
            arrowprops={
                "arrowstyle": "->",
            },
        )

        axis.set_title(
            "M3 Prostate Cancer GWAS × tRF-QTL Integration",
            fontsize=14,
        )

        figure.tight_layout()

        png_path = (
            output_directory
            / "m3_dataset_flow.png"
        )

        pdf_path = (
            output_directory
            / "m3_dataset_flow.pdf"
        )

        figure.savefig(
            png_path,
            dpi=300,
            bbox_inches="tight",
        )

        figure.savefig(
            pdf_path,
            bbox_inches="tight",
        )

        plt.close(
            figure
        )

        return [
            png_path,
            pdf_path,
        ]

    @staticmethod
    def variant_intersection(
        statistics: dict[str, Any],
        output_directory: Path,
    ) -> list[Path]:
        """Create variant-set count and intersection figure."""

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        gwas_unique = (
            statistics[
                "gwas"
            ][
                "unique_eligible_rsids"
            ]
        )

        trfqtl_unique = (
            statistics[
                "trfqtl"
            ][
                "unique_eligible_rsids"
            ]
        )

        shared = (
            statistics[
                "direct_overlap"
            ][
                "unique_shared_rsids"
            ]
        )

        labels = [
            "Prostate GWAS\nvariants",
            "PRAD tRF-QTL\nvariants",
            "Exact shared\nvariants",
        ]

        values = [
            gwas_unique,
            trfqtl_unique,
            shared,
        ]

        figure, axis = plt.subplots(
            figsize=(
                8,
                5,
            )
        )

        bars = axis.bar(
            labels,
            values,
        )

        axis.set_ylabel(
            "Number of unique canonical rsIDs"
        )

        axis.set_title(
            "Direct Variant Intersection"
        )

        for (
            bar,
            value,
        ) in zip(
            bars,
            values,
            strict=True,
        ):
            axis.text(
                bar.get_x()
                + (
                    bar.get_width()
                    / 2
                ),
                bar.get_height(),
                f"{value:,}",
                ha="center",
                va="bottom",
            )

        figure.tight_layout()

        png_path = (
            output_directory
            / "m3_variant_intersection.png"
        )

        pdf_path = (
            output_directory
            / "m3_variant_intersection.pdf"
        )

        figure.savefig(
            png_path,
            dpi=300,
            bbox_inches="tight",
        )

        figure.savefig(
            pdf_path,
            bbox_inches="tight",
        )

        plt.close(
            figure
        )

        return [
            png_path,
            pdf_path,
        ]