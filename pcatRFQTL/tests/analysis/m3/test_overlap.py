"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/analysis/m3/test_overlap.py

Description:
    Unit tests for M3.3 direct GWAS × tRF-QTL overlap.

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

from pcatrfqtl.analysis.m3.overlap import (
    DirectOverlapAnalyzer,
)


def make_gwas() -> pd.DataFrame:
    """Create minimal M3.1 GWAS fixture."""

    return pd.DataFrame(
        {
            "harm_variant_rsid": [
                "rs100",
                "rs100",
                "rs200",
                None,
            ],
            "m3_direct_overlap_eligible": [
                True,
                True,
                True,
                False,
            ],
            "study": [
                "A",
                "B",
                "C",
                "D",
            ],
        }
    )


def make_trfqtl() -> pd.DataFrame:
    """Create minimal M3.2 tRF-QTL fixture."""

    return pd.DataFrame(
        {
            "harm_variant_rsid": [
                "rs100",
                "rs300",
            ],
            "harm_feature_id": [
                "tRF-A",
                "tRF-B",
            ],
            "harm_feature_identity_key": [
                "trf:source:tRF-A",
                "trf:source:tRF-B",
            ],
            "m3_direct_overlap_eligible": [
                True,
                True,
            ],
        }
    )


def test_shared_rsids() -> None:
    """Exact canonical rsID intersection should be identified."""

    result = (
        DirectOverlapAnalyzer
        .shared_rsids(
            make_gwas(),
            make_trfqtl(),
        )
    )

    assert result == {
        "rs100"
    }


def test_overlap_preserves_gwas_association_multiplicity() -> None:
    """Multiple GWAS rows for one rsID should remain represented."""

    result = (
        DirectOverlapAnalyzer
        .build_overlap_table(
            make_gwas(),
            make_trfqtl(),
        )
    )

    assert (
        len(
            result
        )
        == 2
    )

    assert (
        result[
            "harm_variant_rsid"
        ]
        .eq(
            "rs100"
        )
        .all()
    )


def test_overlap_method() -> None:
    """Overlap provenance should state exact canonical-rsID matching."""

    result = (
        DirectOverlapAnalyzer
        .build_overlap_table(
            make_gwas(),
            make_trfqtl(),
        )
    )

    assert (
        result[
            "m3_overlap_method"
        ]
        .eq(
            "canonical_rsid_exact"
        )
        .all()
    )


def test_zero_overlap_returns_empty_dataframe() -> None:
    """No shared variant should produce an empty overlap table."""

    trfqtl = (
        make_trfqtl()
        .copy()
    )

    trfqtl[
        "harm_variant_rsid"
    ] = [
        "rs999",
        "rs998",
    ]

    result = (
        DirectOverlapAnalyzer
        .build_overlap_table(
            make_gwas(),
            trfqtl,
        )
    )

    assert result.empty