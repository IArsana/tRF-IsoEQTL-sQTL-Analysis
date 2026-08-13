"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/inspect_gwas_validation_errors.py

Description:
    Inspect GWAS Catalog records rejected by the current field-level
    validation rules without modifying the raw source dataset.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import re

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

GWAS_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "gwas"
    / "gwas-catalog-download-associations-alt-full.tsv"
)

RSID = re.compile(
    r"^rs\d+$",
    re.IGNORECASE,
)

MULTI_RSID = re.compile(
    r"^rs\d+(?:[;,\s]+rs\d+)+$",
    re.IGNORECASE,
)

INTERACTION = re.compile(
    r"^rs\d+\s+[xX]\s+rs\d+$",
    re.IGNORECASE,
)

COORDINATE = re.compile(
    r"^(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT):\d+$",
    re.IGNORECASE,
)


def classify_snps(
    value: object,
) -> str:
    """Classify a GWAS Catalog SNPS field."""

    if pd.isna(value):
        return "missing"

    text = str(value).strip()

    if RSID.fullmatch(text):
        return "single_rsid"

    if INTERACTION.fullmatch(text):
        return "snp_interaction"

    if COORDINATE.fullmatch(text):
        return "coordinate_variant"

    if MULTI_RSID.fullmatch(text):
        return "multi_rsid"

    return "other"


def main() -> int:
    """Inspect rejected GWAS representations."""

    snp_counts: Counter[str] = Counter()

    chromosome_examples: list[tuple[object, object, object]] = []
    position_examples: list[tuple[object, object, object]] = []
    pvalue_examples: list[object] = []

    other_snp_examples: list[str] = []

    for chunk in pd.read_csv(
        GWAS_FILE,
        sep="\t",
        chunksize=100_000,
        low_memory=False,
    ):
        # --------------------------------------------------------------
        # SNPS
        # --------------------------------------------------------------

        classifications = (
            chunk["SNPS"]
            .map(classify_snps)
        )

        snp_counts.update(
            classifications
            .value_counts()
            .to_dict()
        )

        other = (
            chunk.loc[
                classifications == "other",
                "SNPS",
            ]
            .dropna()
            .astype(str)
            .unique()
        )

        for value in other:
            if (
                len(other_snp_examples)
                < 100
            ):
                other_snp_examples.append(
                    value
                )

        # --------------------------------------------------------------
        # CHR_ID
        # --------------------------------------------------------------

        chr_text = (
            chunk["CHR_ID"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        valid_chr = chr_text.str.fullmatch(
            r"(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT|M)"
            r"(?:[;,](?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT|M))*",
            case=False,
            na=False,
        )

        invalid_chr_indices = (
            valid_chr[
                ~valid_chr
            ]
            .index
        )

        for index in invalid_chr_indices:
            if (
                len(chromosome_examples)
                < 50
            ):
                chromosome_examples.append(
                    (
                        chunk.at[
                            index,
                            "SNPS",
                        ],
                        chunk.at[
                            index,
                            "CHR_ID",
                        ],
                        chunk.at[
                            index,
                            "CHR_POS",
                        ],
                    )
                )

        # --------------------------------------------------------------
        # P-VALUE
        # --------------------------------------------------------------

        p_values = pd.to_numeric(
            chunk["P-VALUE"],
            errors="coerce",
        )

        invalid_p = (
            p_values.isna()
            | (p_values < 0)
            | (p_values > 1)
            | (p_values == 0)
        )

        values = (
            chunk.loc[
                invalid_p,
                "P-VALUE",
            ]
            .tolist()
        )

        for value in values:
            if len(pvalue_examples) < 100:
                pvalue_examples.append(
                    value
                )

    print("=" * 80)
    print("SNPS classification")
    print("=" * 80)

    for key, value in (
        snp_counts
        .most_common()
    ):
        print(
            f"{key:20s}: {value}"
        )

    print()
    print("=" * 80)
    print("Examples of unclassified SNPS values")
    print("=" * 80)

    for value in other_snp_examples:
        print(
            repr(value)
        )

    print()
    print("=" * 80)
    print("Invalid chromosome/position examples")
    print("=" * 80)

    for snps, chromosome, position in chromosome_examples:
        print(
            f"SNPS={snps!r} | "
            f"CHR_ID={chromosome!r} | "
            f"CHR_POS={position!r}"
        )

    print()
    print("=" * 80)
    print("Invalid P-VALUE examples")
    print("=" * 80)

    for value in pvalue_examples:
        print(
            repr(value)
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )