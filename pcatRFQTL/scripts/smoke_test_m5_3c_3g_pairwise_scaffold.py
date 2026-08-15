"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_3c_3g_pairwise_scaffold.py

Description:
    Scaffold-only smoke test for M5.3C.3G.

    This script builds the disease-proxy × regulatory-proxy candidate
    matrix for rs10216902 without performing any remote LD query.

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

import pyarrow.parquet as pq

from pcatrfqtl.analysis.m5.disease_regulatory_pairwise_ld import (
    build_pairwise_candidate_scaffold,
)
from pcatrfqtl.io.parquet import write_parquet


ROOT = Path(__file__).resolve().parents[1]

M5_ROOT = ROOT / "data/processed/m5"

RECURRENCE_PATH = (
    M5_ROOT
    / "sumstats"
    / "gwas_cross_study_variant_recurrence.parquet"
)

OUTPUT_PATH = (
    M5_ROOT
    / "pairwise_ld"
    / "disease_regulatory_pairwise_candidates.parquet"
)

POPULATIONS = (
    "EAS",
    "EUR",
    "SAS",
)


def read_parquet_compat(path: Path):
    """Read Parquet while ignoring pandas reconstruction metadata."""

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    return table.to_pandas(
        ignore_metadata=True,
    )


def main() -> None:
    """Build and inspect the pairwise candidate scaffold."""

    recurrence = read_parquet_compat(
        RECURRENCE_PATH
    )

    regulatory_bridges = {}

    for population in POPULATIONS:

        path = (
            M5_ROOT
            / f"moradi_ld_bridge_{population}.parquet"
        )

        regulatory_bridges[population] = (
            read_parquet_compat(
                path
            )
        )

    result = build_pairwise_candidate_scaffold(
        recurrence=recurrence,
        regulatory_bridges=regulatory_bridges,
        lead_rsid="rs10216902",
    )

    scaffold = result.dataframe

    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_parquet(
        scaffold,
        OUTPUT_PATH,
        index=False,
    )

    print()
    print("=" * 72)
    print("M5.3C.3G PAIRWISE CANDIDATE SCAFFOLD")
    print("=" * 72)

    print(
        f"Disease proxies:      "
        f"{result.disease_variants}"
    )

    print(
        f"Regulatory proxies:   "
        f"{result.regulatory_variants}"
    )

    print(
        f"Candidate pairs:       "
        f"{result.candidate_pairs}"
    )

    print(
        f"Reference populations: "
        f"{len(POPULATIONS)}"
    )

    print(
        f"Expected LD queries:   "
        f"{result.candidate_pairs * len(POPULATIONS)}"
    )

    print()
    print(
        f"Output: {OUTPUT_PATH}"
    )

    print("=" * 72)


if __name__ == "__main__":
    main()