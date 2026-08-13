"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/inspect_moradi_schema.py

Description:
    Schema profiling utility for the Moradi supplementary datasets.

    The script inspects Excel worksheets and CSV supplementary files
    before formal validation rules are defined. It reports sheet names,
    dimensions, column names, inferred data types, missing-value counts,
    and representative rows without modifying the source data.

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

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

MORADI_DIR = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "moradi"
)


EXCEL_FILES = [
    "Supplementary Table S1.xlsx",
    "Supplementary Table S2.xlsx",
    "Supplementary Table S3.xlsx",
    "Supplementary Table S4.xlsx",
]

CSV_FILES = [
    "Supplementary Table S5_DE_EP_cis.csv",
    "Supplementary Table S6_DE_EP_trans.csv",
    "Supplementary Table S7_DE_INT_cis.csv",
    "Supplementary Table S8_DE_INT_trans.csv",
    "Supplementary Table S9_DE_Isoform_cis.csv",
    "Supplementary Table S10_DE_Isoform_trans_cancer_vs_control.csv",
]


def print_separator(
    character: str = "=",
    width: int = 100,
) -> None:
    """Print a visual separator."""

    print(character * width)


def profile_dataframe(
    dataframe: pd.DataFrame,
    *,
    name: str,
) -> None:
    """
    Display basic schema information for a dataframe.
    """

    print()
    print_separator("-")
    print(name)
    print_separator("-")

    print(f"Rows    : {len(dataframe)}")
    print(f"Columns : {len(dataframe.columns)}")

    print("\nColumn names:")
    for index, column in enumerate(
        dataframe.columns,
        start=1,
    ):
        print(
            f"  {index:02d}. "
            f"{column!r} "
            f"[dtype={dataframe[column].dtype}]"
        )

    print("\nMissing values:")
    for column in dataframe.columns:
        missing = int(
            dataframe[column]
            .isna()
            .sum()
        )

        print(
            f"  {column!r}: {missing}"
        )

    print("\nFirst 5 rows:")
    print(
        dataframe
        .head(5)
        .to_string(index=False)
    )


def inspect_excel_file(
    path: Path,
) -> None:
    """
    Inspect every sheet in a Moradi Excel workbook.
    """

    print()
    print_separator()
    print(f"FILE: {path.name}")
    print_separator()

    if not path.exists():
        print(f"[MISSING] {path}")
        return

    workbook = pd.ExcelFile(path)

    print(
        f"Sheets ({len(workbook.sheet_names)}):"
    )

    for sheet_name in workbook.sheet_names:
        print(f"  - {sheet_name}")

    for sheet_name in workbook.sheet_names:
        try:
            dataframe = pd.read_excel(
                path,
                sheet_name=sheet_name,
                header=0,
            )

            profile_dataframe(
                dataframe,
                name=f"{path.name} :: {sheet_name}",
            )

        except Exception as exc:
            print(
                f"[ERROR] Could not read "
                f"{path.name} :: {sheet_name}: {exc}"
            )


def inspect_csv_file(
    path: Path,
) -> None:
    """
    Inspect one Moradi CSV supplementary table.
    """

    print()
    print_separator()
    print(f"FILE: {path.name}")
    print_separator()

    if not path.exists():
        print(f"[MISSING] {path}")
        return

    try:
        dataframe = pd.read_csv(path)

        profile_dataframe(
            dataframe,
            name=path.name,
        )

    except Exception as exc:
        print(
            f"[ERROR] Could not read "
            f"{path.name}: {exc}"
        )


def main() -> int:
    """Run Moradi supplementary-data schema profiling."""

    print_separator()
    print("pcatRFQTL — Moradi Dataset Schema Profiling")
    print_separator()

    if not MORADI_DIR.exists():
        print(
            f"[FAIL] Moradi directory not found: "
            f"{MORADI_DIR}"
        )
        return 1

    for filename in EXCEL_FILES:
        inspect_excel_file(
            MORADI_DIR / filename
        )

    for filename in CSV_FILES:
        inspect_csv_file(
            MORADI_DIR / filename
        )

    print()
    print_separator()
    print("Moradi schema profiling completed.")
    print_separator()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())