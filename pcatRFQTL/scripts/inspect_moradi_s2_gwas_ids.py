"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/inspect_moradi_s2_gwas_ids.py

Description:
    Inspect GWAS tag-SNP identifiers in Moradi Supplementary Table S2
    that fail the standard single-rsID validator.

    The script does not modify source data. It separates invalid values
    into likely composite rsID representations and genuinely malformed
    values to support evidence-based validation policy.

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
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

WORKBOOK = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "moradi"
    / "Supplementary Table S2.xlsx"
)

SHEETS = [
    "Cis-intron-sQTL-GWAS-LD_0.5",
    "Cis-exon-sQTL-GWAS-LD_0.5",
    "Cis-iso-eQTl-GWAS-LD_0.5",
]

SINGLE_RSID_PATTERN = re.compile(
    r"^rs\d+$",
    re.IGNORECASE,
)

COMPOSITE_RSID_PATTERN = re.compile(
    r"^rs\d+(?::rs\d+)+$",
    re.IGNORECASE,
)


def classify_rsid(
    value: object,
) -> str:
    """Classify a GWAS tag identifier."""

    if pd.isna(value):
        return "missing"

    text = str(value).strip()

    if SINGLE_RSID_PATTERN.fullmatch(text):
        return "single"

    if COMPOSITE_RSID_PATTERN.fullmatch(text):
        return "composite"

    return "malformed"


def main() -> int:
    """Inspect S2 tag_SNP_pos identifiers."""

    excel = pd.ExcelFile(
        WORKBOOK
    )

    for sheet_name in SHEETS:
        print()
        print("=" * 80)
        print(sheet_name)
        print("=" * 80)

        df = pd.read_excel(
            excel,
            sheet_name=sheet_name,
            header=0,
        )

        classifications = (
            df["tag_SNP_pos"]
            .map(classify_rsid)
        )

        counts = (
            classifications
            .value_counts(
                dropna=False
            )
            .to_dict()
        )

        print(
            "Classification counts:",
            counts,
        )

        non_single = df.loc[
            classifications != "single",
            [
                "sQTL_SNP",
                "sQTL_SNP-pos",
                "tag_SNP",
                "tag_SNP_pos",
                "LD",
                "gwas_cancer",
            ],
        ].copy()

        non_single[
            "classification"
        ] = classifications.loc[
            non_single.index
        ]

        print()
        print(
            "Non-single values:"
        )

        print(
            non_single.head(
                30
            ).to_string()
        )

        malformed = non_single[
            non_single[
                "classification"
            ]
            == "malformed"
        ]

        print()
        print(
            "Malformed count:",
            len(malformed),
        )

        if not malformed.empty:
            print(
                malformed[
                    "tag_SNP_pos"
                ]
                .value_counts()
                .head(50)
                .to_string()
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )