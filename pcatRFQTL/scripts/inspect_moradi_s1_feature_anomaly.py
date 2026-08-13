from pathlib import Path
import re

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

WORKBOOK = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "moradi"
    / "Supplementary Table S1.xlsx"
)

SHEET = "Cis-intron-retention-sQTL-0.05"

INT_PATTERN = re.compile(
    r"^INT\d+$",
    re.IGNORECASE,
)

df = pd.read_excel(
    WORKBOOK,
    sheet_name=SHEET,
    header=0,
)

invalid = df[
    ~df["splicing_event"]
    .astype(str)
    .str.strip()
    .str.match(INT_PATTERN)
]

print(
    invalid[
        [
            "SNP",
            "SNP_pos",
            "splicing_event",
            "Stat",
            "p_value",
            "FDR",
            "beta",
        ]
    ].to_string(index=True)
)