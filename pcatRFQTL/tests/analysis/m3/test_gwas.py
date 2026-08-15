"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/analysis/m3/test_gwas.py

Description:
    Unit tests for M3.1 prostate-cancer GWAS selection.

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

from pcatrfqtl.analysis.m3.gwas import (
    ProstateGWASSelector,
)


def make_dataframe() -> pd.DataFrame:
    """Create minimal harmonized GWAS fixture."""

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
                None,
                "rsid:rs300",
            ],
            "harm_variant_status": [
                "RESOLVED",
                "UNRESOLVED",
                "RESOLVED",
            ],
        }
    )


def test_primary_disease_selection() -> None:
    """Only prostate-cancer association rows should remain."""

    dataframe = (
        make_dataframe()
    )

    result = (
        ProstateGWASSelector
        .select_primary_disease(
            dataframe
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
            "harm_is_primary_disease"
        ]
        .all()
    )


def test_direct_overlap_eligibility() -> None:
    """Single canonical rsID controls initial direct-overlap eligibility."""

    result = (
        ProstateGWASSelector
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
    """
    Prostate-cancer records without a canonical rsID should remain in
    M3.1 and be marked overlap-ineligible.
    """

    result = (
        ProstateGWASSelector
        .prepare(
            make_dataframe()
        )
    )

    unresolved = (
        result[
            result[
                "harm_variant_rsid"
            ]
            .isna()
        ]
    )

    assert (
        len(
            unresolved
        )
        == 1
    )

    assert (
        unresolved.iloc[
            0
        ][
            "m3_direct_overlap_eligible"
        ]
        == False
    )


def test_selection_reason() -> None:
    """Selection provenance should be explicit."""

    result = (
        ProstateGWASSelector
        .prepare(
            make_dataframe()
        )
    )

    assert (
        result[
            "m3_selection_reason"
        ]
        .eq(
            "harmonized_primary_disease"
        )
        .all()
    )


def test_schema_validation() -> None:
    """Missing harmonization fields should fail loudly."""

    dataframe = pd.DataFrame(
        {
            "harm_is_primary_disease": [
                True
            ]
        }
    )

    with pytest.raises(
        ValueError
    ):
        (
            ProstateGWASSelector
            .prepare(
                dataframe
            )
        )