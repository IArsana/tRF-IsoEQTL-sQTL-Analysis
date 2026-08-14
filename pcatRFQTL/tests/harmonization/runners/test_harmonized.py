"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/harmonization/runners/test_harmonize.py

Description:
    Unit tests for M3.6 harmonized dataset execution.

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
import pytest

from pcatrfqtl.harmonization.features import (
    FeatureType,
)
from pcatrfqtl.harmonization.runners.harmonize import (
    HarmonizationInputs,
    HarmonizationRunner,
)


def make_runner(
    tmp_path: Path,
) -> HarmonizationRunner:
    """Create temporary M3 runner."""

    standardized = (
        tmp_path
        / "standardized"
    )

    standardized.mkdir()

    return HarmonizationRunner(
        inputs=HarmonizationInputs(
            standardized_directory=(
                standardized
            )
        ),
        output_directory=(
            tmp_path
            / "harmonized"
        ),
    )


def test_inputs() -> None:
    """Input dataclass should preserve path."""

    inputs = HarmonizationInputs(
        standardized_directory=Path(
            "data/interim/standardized"
        )
    )

    assert (
        inputs.standardized_directory
        == Path(
            "data/interim/standardized"
        )
    )


def test_verify_cardinality() -> None:
    """Equal cardinality should pass."""

    HarmonizationRunner._verify_cardinality(
        dataset_name="test",
        input_rows=10,
        output_rows=10,
    )


def test_verify_cardinality_failure() -> None:
    """Different cardinality should fail."""

    with pytest.raises(
        RuntimeError
    ):
        HarmonizationRunner._verify_cardinality(
            dataset_name="test",
            input_rows=10,
            output_rows=9,
        )


def test_append_records() -> None:
    """Metadata should append without altering rows."""

    dataframe = pd.DataFrame(
        {
            "value": [
                1,
                2,
            ]
        }
    )

    records = [
        {
            "harm_test": "A",
        },
        {
            "harm_test": "B",
        },
    ]

    result = (
        HarmonizationRunner
        ._append_records(
            dataframe,
            records,
        )
    )

    assert len(
        result
    ) == 2

    assert (
        "harm_test"
        in result.columns
    )


def test_append_records_rejects_column_collision() -> None:
    """Harmonization must not overwrite standardized fields."""

    dataframe = pd.DataFrame(
        {
            "harm_test": [
                "old"
            ]
        }
    )

    records = [
        {
            "harm_test": "new",
        }
    ]

    with pytest.raises(
        RuntimeError
    ):
        (
            HarmonizationRunner
            ._append_records(
                dataframe,
                records,
            )
        )


def test_boolean_count() -> None:
    """True rows should be counted."""

    dataframe = pd.DataFrame(
        {
            "flag": [
                True,
                False,
                True,
            ]
        }
    )

    assert (
        HarmonizationRunner
        ._boolean_count(
            dataframe,
            "flag",
        )
        == 2
    )


def test_harmonize_trfqtl_source_anomaly(
    tmp_path: Path,
) -> None:
    """
    Cancer-tRFQTL non-rsID source value should retain coordinate
    identity.
    """

    runner = make_runner(
        tmp_path
    )

    dataframe = pd.DataFrame(
        {
            "snp_id_raw": [
                "6"
            ],
            "canonical_rsid": [
                None
            ],
            "coordinate": [
                "6:28958399"
            ],
            "Cancer type": [
                "PRAD"
            ],
        }
    )

    result = (
        runner
        ._harmonize_trfqtl_dataframe(
            dataframe,
            table="S10",
        )
    )

    assert (
        result.loc[
            0,
            "harm_variant_rsid",
        ]
        is None
    )

    assert (
        result.loc[
            0,
            "harm_build_aware_coordinate_key",
        ]
        == "hg19:6:28958399"
    )

    assert (
        result.loc[
            0,
            "harm_coordinate_join_allowed",
        ]
        == True
    )

    assert (
        result.loc[
            0,
            "harm_is_primary_disease",
        ]
        == True
    )


def test_harmonize_moradi_isoform(
    tmp_path: Path,
) -> None:
    """Moradi ENST feature should harmonize correctly."""

    runner = make_runner(
        tmp_path
    )

    dataframe = pd.DataFrame(
        {
            "SNP": [
                "rs71762629"
            ],
            "SNP_pos": [
                "10:12658508"
            ],
            "splicing_event": [
                "ENST00000264639.9"
            ],
        }
    )

    result = (
        runner
        ._harmonize_moradi_qtl_dataframe(
            dataframe,
            feature_type=(
                FeatureType
                .TRANSCRIPT_ISOFORM
            ),
            implicit_disease="PrCa",
        )
    )

    assert (
        result.loc[
            0,
            "harm_variant_rsid",
        ]
        == "rs71762629"
    )

    assert (
        result.loc[
            0,
            "harm_feature_id",
        ]
        == "ENST00000264639.9"
    )

    assert (
        result.loc[
            0,
            "harm_feature_base_id",
        ]
        == "ENST00000264639"
    )

    assert (
        result.loc[
            0,
            "harm_is_primary_disease",
        ]
        == True
    )

    assert (
        result.loc[
            0,
            "harm_disease_context_source",
        ]
        == "DATASET"
    )


def test_moradi_known_intron_anomaly(
    tmp_path: Path,
) -> None:
    """INT1e+05 must remain a source anomaly."""

    runner = make_runner(
        tmp_path
    )

    dataframe = pd.DataFrame(
        {
            "SNP": [
                "rs12786544"
            ],
            "SNP_pos": [
                "11:69561656"
            ],
            "splicing_event": [
                "INT1e+05"
            ],
        }
    )

    result = (
        runner
        ._harmonize_moradi_qtl_dataframe(
            dataframe,
            feature_type=(
                FeatureType
                .INTRON_EVENT
            ),
            implicit_disease="PrCa",
        )
    )

    assert (
        result.loc[
            0,
            "harm_feature_id",
        ]
        is None
    )

    assert (
        result.loc[
            0,
            "harm_feature_usable",
        ]
        == False
    )

    assert (
        result.loc[
            0,
            "harm_feature_source_anomaly",
        ]
        == True
    )


def test_gwas_coordinate_join_remains_blocked(
    tmp_path: Path,
) -> None:
    """GWAS numeric coordinates must remain unjoinable in M3."""

    runner = make_runner(
        tmp_path
    )

    dataframe = pd.DataFrame(
        {
            "SNPS": [
                "rs123"
            ],
            "CHR_ID": [
                "6"
            ],
            "CHR_POS": [
                100
            ],
            "DISEASE/TRAIT": [
                "Prostate cancer"
            ],
        }
    )

    result = (
        runner
        ._harmonize_gwas_dataframe(
            dataframe
        )
    )

    assert (
        result.loc[
            0,
            "harm_variant_rsid",
        ]
        == "rs123"
    )

    assert (
        result.loc[
            0,
            "harm_coordinate_join_allowed",
        ]
        == False
    )

    assert (
        result.loc[
            0,
            "harm_is_primary_disease",
        ]
        == True
    )

def test_value_skips_missing_candidate() -> None:
    """Missing preferred values should fall back to later source fields."""

    row = pd.Series(
        {
            "canonical_feature_id": None,
            "feature_id": None,
            "feature_id_raw": "INT1e+05",
        }
    )

    result = (
        HarmonizationRunner
        ._value(
            row,
            (
                "canonical_feature_id",
                "feature_id",
                "feature_id_raw",
            ),
        )
    )

    assert (
        result
        == "INT1e+05"
    )

def test_value_prefers_first_non_missing_candidate() -> None:
    """Preferred standardized values should win when usable."""

    row = pd.Series(
        {
            "canonical_rsid": "rs123",
            "snp_id_raw": "RS123",
        }
    )

    result = (
        HarmonizationRunner
        ._value(
            row,
            (
                "canonical_rsid",
                "snp_id_raw",
            ),
        )
    )

    assert (
        result
        == "rs123"
    )

def test_value_skips_missing_preferred_candidate() -> None:
    """Missing canonical values should fall back to source values."""

    row = pd.Series(
        {
            "feature_id": None,
            "feature_raw": "INT1e+05",
        }
    )

    result = (
        HarmonizationRunner
        ._value(
            row,
            (
                "feature_id",
                "feature_raw",
            ),
        )
    )

    assert (
        result
        == "INT1e+05"
    )

def test_harmonize_moradi_ld_variants(
    tmp_path: Path,
) -> None:
    """Moradi LD tables should preserve both QTL and tag variants."""

    runner = make_runner(
        tmp_path
    )

    dataframe = pd.DataFrame(
        {
            "qtl_rsid_raw": [
                "rs200367988"
            ],
            "qtl_rsid": [
                "rs200367988"
            ],
            "qtl_coordinate_raw": [
                "10:122674849"
            ],
            "qtl_chromosome": [
                "10"
            ],
            "qtl_position": [
                122674849
            ],

            "feature_raw": [
                "INT78964"
            ],
            "feature_id": [
                "INT78964"
            ],

            "tag_rsid_raw": [
                "rs1439466"
            ],
            "tag_primary_rsid": [
                "rs1439466"
            ],
            "tag_coordinate_raw": [
                "10:122654969"
            ],
            "tag_chromosome": [
                "10"
            ],
            "tag_position": [
                122654969
            ],

            "gwas_cancer": [
                "prostate_cancer"
            ],
            "ld_r2": [
                0.696408
            ],
        }
    )

    result = (
        runner
        ._harmonize_moradi_qtl_dataframe(
            dataframe,
            feature_type=(
                FeatureType.INTRON_EVENT
            ),
            implicit_disease=None,
        )
    )

    assert (
        result.loc[
            0,
            "harm_qtl_variant_rsid",
        ]
        == "rs200367988"
    )

    assert (
        result.loc[
            0,
            "harm_tag_variant_rsid",
        ]
        == "rs1439466"
    )

    assert (
        result.loc[
            0,
            "harm_qtl_variant_identity_key",
        ]
        == "rsid:rs200367988"
    )

    assert (
        result.loc[
            0,
            "harm_tag_variant_identity_key",
        ]
        == "rsid:rs1439466"
    )

    assert (
        result.loc[
            0,
            "harm_is_primary_disease",
        ]
        == True
    )

def test_moradi_de_source_anomaly_falls_back_to_raw_feature(
    tmp_path: Path,
) -> None:
    """Known S7 anomaly should survive M3 harmonization."""

    runner = make_runner(
        tmp_path
    )

    dataframe = pd.DataFrame(
        {
            "feature_id": [
                None
            ],
            "feature_raw": [
                "INT1e+05"
            ],
            "source_feature_anomaly": [
                True
            ],
        }
    )

    result = (
        runner
        ._harmonize_moradi_de_dataframe(
            dataframe,
            feature_type=(
                FeatureType.INTRON_EVENT
            ),
        )
    )

    assert (
        result.loc[
            0,
            "harm_feature_source_id",
        ]
        == "INT1e+05"
    )

    assert (
        result.loc[
            0,
            "harm_feature_usable",
        ]
        == False
    )

    assert (
        result.loc[
            0,
            "harm_feature_source_anomaly",
        ]
        == True
    )

