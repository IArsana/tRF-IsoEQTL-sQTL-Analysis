"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/analysis/m3/test_candidates.py

Description:
    Unit tests for M3.4 integrated SNP-tRF candidate generation.

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

from pcatrfqtl.analysis.m3.candidates import (
    CandidateDirectEvidence,
    SNPTRFCandidateBuilder,
)


def make_gwas() -> pd.DataFrame:
    """Create minimal prostate-cancer GWAS fixture."""

    return pd.DataFrame(
        {
            "harm_variant_rsid": [
                "rs100",
                "rs100",
                "rs200",
            ],
            "m3_direct_overlap_eligible": [
                True,
                True,
                True,
            ],
        }
    )


def make_trfqtl() -> pd.DataFrame:
    """Create minimal prostate-cancer tRF-QTL fixture."""

    return pd.DataFrame(
        {
            "harm_variant_rsid": [
                "rs100",
                "rs300",
                None,
            ],
            "harm_variant_identity_key": [
                "rsid:rs100",
                "rsid:rs300",
                "coord:hg19:1:300",
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
            "m3_direct_overlap_eligible": [
                True,
                True,
                False,
            ],
        }
    )


def test_gwas_evidence_collapses_association_rows() -> None:
    """GWAS evidence should summarize association multiplicity per rsID."""

    result = (
        SNPTRFCandidateBuilder
        .build_gwas_evidence_table(
            make_gwas()
        )
    )

    rs100 = (
        result.loc[
            result[
                "harm_variant_rsid"
            ]
            .eq(
                "rs100"
            )
        ]
        .iloc[
            0
        ]
    )

    assert (
        rs100[
            "m3_gwas_association_count"
        ]
        == 2
    )


def test_candidate_cardinality_is_preserved() -> None:
    """Every tRF-QTL row should become one M3.4 candidate row."""

    trfqtl = (
        make_trfqtl()
    )

    result = (
        SNPTRFCandidateBuilder
        .build(
            gwas=make_gwas(),
            trfqtl=trfqtl,
        )
    )

    assert (
        len(
            result
        )
        == len(
            trfqtl
        )
    )


def test_direct_gwas_candidate() -> None:
    """Shared canonical rsID should receive direct GWAS evidence."""

    result = (
        SNPTRFCandidateBuilder
        .build(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
        )
    )

    candidate = (
        result.loc[
            result[
                "harm_variant_rsid"
            ]
            .eq(
                "rs100"
            )
        ]
        .iloc[
            0
        ]
    )

    assert (
        candidate[
            "m3_direct_gwas_match"
        ]
        == True
    )

    assert (
        candidate[
            "m3_gwas_association_count"
        ]
        == 2
    )

    assert (
        candidate[
            "m3_candidate_direct_evidence"
        ]
        == (
            CandidateDirectEvidence
            .DIRECT_GWAS_MATCH
            .value
        )
    )


def test_no_direct_gwas_candidate_is_preserved() -> None:
    """Candidate without direct overlap must not be removed."""

    result = (
        SNPTRFCandidateBuilder
        .build(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
        )
    )

    candidate = (
        result.loc[
            result[
                "harm_variant_rsid"
            ]
            .eq(
                "rs300"
            )
        ]
        .iloc[
            0
        ]
    )

    assert (
        candidate[
            "m3_direct_gwas_match"
        ]
        == False
    )

    assert (
        candidate[
            "m3_gwas_association_count"
        ]
        == 0
    )

    assert (
        candidate[
            "m3_candidate_direct_evidence"
        ]
        == (
            CandidateDirectEvidence
            .NO_DIRECT_GWAS_MATCH
            .value
        )
    )


def test_ineligible_candidate_is_preserved() -> None:
    """Coordinate-only tRF-QTL candidate remains represented."""

    result = (
        SNPTRFCandidateBuilder
        .build(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
        )
    )

    candidate = (
        result.loc[
            result[
                "harm_variant_rsid"
            ]
            .isna()
        ]
        .iloc[
            0
        ]
    )

    assert (
        candidate[
            "m3_candidate_direct_evidence"
        ]
        == (
            CandidateDirectEvidence
            .DIRECT_MATCH_INELIGIBLE
            .value
        )
    )


def test_candidate_ids_are_unique() -> None:
    """Candidate identifiers should be deterministic and unique."""

    result = (
        SNPTRFCandidateBuilder
        .build(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
        )
    )

    assert (
        result[
            "m3_candidate_id"
        ]
        .tolist()
        == [
            "M3C000001",
            "M3C000002",
            "M3C000003",
        ]
    )

    assert (
        result[
            "m3_candidate_id"
        ]
        .is_unique
    )


def test_future_evidence_flags_are_false() -> None:
    """M3.4 must not imply LD or colocalization was performed."""

    result = (
        SNPTRFCandidateBuilder
        .build(
            gwas=make_gwas(),
            trfqtl=make_trfqtl(),
        )
    )

    assert (
        result[
            "m3_ld_evidence_evaluated"
        ]
        .eq(False)
        .all()
    )

    assert (
        result[
            "m3_colocalization_evaluated"
        ]
        .eq(False)
        .all()
    )

    assert (
        result[
            "m3_candidate_ranked"
        ]
        .eq(False)
        .all()
    )