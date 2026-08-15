"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/analysis/m3/test_statistics.py

Description:
    Unit tests for M3.5 descriptive statistical analysis.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import pandas as pd

from pcatrfqtl.analysis.m3.statistics import (
    M3StatisticsAnalyzer,
)


def make_gwas() -> pd.DataFrame:

    return pd.DataFrame(
        {
            "harm_variant_rsid": [
                "rs1",
                "rs1",
                "rs2",
            ],
            "m3_direct_overlap_eligible": [
                True,
                True,
                True,
            ],
        }
    )


def make_trfqtl() -> pd.DataFrame:

    return pd.DataFrame(
        {
            "harm_variant_rsid": [
                "rs3",
                "rs4",
            ],
            "harm_feature_id": [
                "tRF-A",
                "tRF-B",
            ],
            "harm_feature_identity_key": [
                "trf:A",
                "trf:B",
            ],
            "m3_direct_overlap_eligible": [
                True,
                True,
            ],
        }
    )


def make_candidates() -> pd.DataFrame:

    return pd.DataFrame(
        {
            "m3_candidate_id": [
                "M3C000001",
                "M3C000002",
            ],
            "harm_variant_rsid": [
                "rs3",
                "rs4",
            ],
            "harm_variant_identity_key": [
                "rsid:rs3",
                "rsid:rs4",
            ],
            "harm_feature_id": [
                "tRF-A",
                "tRF-B",
            ],
            "harm_feature_identity_key": [
                "trf:A",
                "trf:B",
            ],
            "m3_direct_overlap_eligible": [
                True,
                True,
            ],
            "m3_direct_gwas_match": [
                False,
                False,
            ],
            "m3_gwas_association_count": [
                0,
                0,
            ],
            "m3_candidate_direct_evidence": [
                "NO_DIRECT_GWAS_MATCH",
                "NO_DIRECT_GWAS_MATCH",
            ],
        }
    )


def test_zero_direct_overlap_statistics() -> None:

    result = (
        M3StatisticsAnalyzer
        .summarize(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
            candidates=make_candidates(),
        )
    )

    assert (
        result[
            "direct_overlap"
        ][
            "unique_shared_rsids"
        ]
        == 0
    )

    assert (
        result[
            "direct_overlap"
        ][
            "trfqtl_variant_overlap_rate_percent"
        ]
        == 0.0
    )


def test_unique_gwas_variants() -> None:

    result = (
        M3StatisticsAnalyzer
        .summarize(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
            candidates=make_candidates(),
        )
    )

    assert (
        result[
            "gwas"
        ][
            "unique_eligible_rsids"
        ]
        == 2
    )


def test_candidate_table_is_not_ranked() -> None:

    result = (
        M3StatisticsAnalyzer
        .build_candidate_statistics_table(
            make_candidates()
        )
    )

    assert (
        result[
            "m3_rank"
        ]
        .isna()
        .all()
    )

    assert (
        result[
            "m3_ranking_performed"
        ]
        .eq(False)
        .all()
    )