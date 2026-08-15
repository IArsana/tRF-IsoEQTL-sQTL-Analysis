"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/analysis/m3/test_trfqtl.py

Description:
    Unit tests for M3.2 prostate-cancer tRF-QTL selection.

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
import pytest

from pcatrfqtl.analysis.m3.trfqtl import (
    ProstateTRFQTLSelector,
)


def make_dataframe() -> pd.DataFrame:
    """Create minimal harmonized Cancer-tRFQTL fixture."""

    return pd.DataFrame(
        {
            "harm_is_primary_disease": [
                True,
                True,
                False,
            ],
            "harm_disease_id": [
                "prostate_cancer",
                "prostate_cancer",
                None,
            ],
            "harm_variant_rsid": [
                "rs100",
                None,
                "rs300",
            ],
            "harm_variant_rsid_usable": [
                True,
                False,
                True,
            ],
            "harm_variant_identity_key": [
                "rsid:rs100",
                "coord:hg19:1:200",
                "rsid:rs300",
            ],
            "harm_feature_id": [
                "tRF-A",
                "tRF-B",
                "tRF-C",
            ],
            "harm_feature_identity_key": [
                "trf:source:tRF-A",
                "trf:source:tRF-B",
                "trf:source:tRF-C",
            ],
            "harm_feature_usable": [
                True,
                True,
                True,
            ],
            "harm_feature_type": [
                "TRF",
                "TRF",
                "TRF",
            ],
        }
    )


def test_primary_disease_selection() -> None:
    """Only prostate-cancer tRF-QTL rows should remain."""

    result = (
        ProstateTRFQTLSelector
        .select_primary_disease(
            make_dataframe()
        )
    )

    assert (
        len(
            result
        )
        == 2
    )


def test_overlap_eligibility() -> None:
    """Canonical rsID and usable tRF identity control eligibility."""

    result = (
        ProstateTRFQTLSelector
        .prepare(
            make_dataframe()
        )
    )

    assert (
        result[
            "m3_direct_overlap_eligible"
        ]
        .tolist()
        == [
            True,
            False,
        ]
    )


def test_unresolved_variant_is_preserved() -> None:
    """Coordinate-only prostate tRF-QTL should remain in M3.2."""

    result = (
        ProstateTRFQTLSelector
        .prepare(
            make_dataframe()
        )
    )

    assert (
        len(
            result
        )
        == 2
    )

    assert (
        result.iloc[
            1
        ][
            "harm_feature_id"
        ]
        == "tRF-B"
    )


def test_schema_validation() -> None:
    """Missing harmonized fields should fail loudly."""

    dataframe = pd.DataFrame(
        {
            "harm_disease_id": [
                "prostate_cancer"
            ]
        }
    )

    with pytest.raises(
        ValueError
    ):
        (
            ProstateTRFQTLSelector
            .prepare(
                dataframe
            )
        )