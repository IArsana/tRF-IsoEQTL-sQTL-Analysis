"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/standardization/runners/test_standardize.py

Description:
    Unit tests for the M2 dataset-level standardization runner.

    Scientific transformations are tested separately under
    tests/standardization/. These tests focus on orchestration:

        - input configuration
        - directory handling
        - cardinality checks
        - Parquet writing
        - status summaries
        - boolean summaries
        - chunked GWAS execution

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

from pcatrfqtl.standardization.runners.standardize import (
    StandardizationInputs,
    StandardizationRunner,
)


def make_runner(
    tmp_path: Path,
    *,
    chunk_size: int = 2,
) -> StandardizationRunner:
    """Create a temporary standardization runner."""

    inputs = StandardizationInputs(
        gwas_catalog=(
            tmp_path
            / "gwas.tsv"
        ),
        cancer_trfqtl_workbook=(
            tmp_path
            / "trfqtl.xlsx"
        ),
        moradi_directory=(
            tmp_path
            / "moradi"
        ),
    )

    return StandardizationRunner(
        inputs=inputs,
        output_directory=(
            tmp_path
            / "standardized"
        ),
        gwas_chunk_size=chunk_size,
    )


def test_standardization_inputs() -> None:
    """Input paths should remain unchanged."""

    inputs = StandardizationInputs(
        gwas_catalog=Path(
            "gwas.tsv"
        ),
        cancer_trfqtl_workbook=Path(
            "trfqtl.xlsx"
        ),
        moradi_directory=Path(
            "moradi"
        ),
    )

    assert (
        inputs.gwas_catalog
        == Path(
            "gwas.tsv"
        )
    )

    assert (
        inputs.cancer_trfqtl_workbook
        == Path(
            "trfqtl.xlsx"
        )
    )

    assert (
        inputs.moradi_directory
        == Path(
            "moradi"
        )
    )


def test_invalid_gwas_chunk_size(
    tmp_path: Path,
) -> None:
    """Chunk size must be positive."""

    with pytest.raises(
        ValueError
    ):
        make_runner(
            tmp_path,
            chunk_size=0,
        )


def test_prepare_directory(
    tmp_path: Path,
) -> None:
    """Runner should create output directories."""

    output = (
        tmp_path
        / "nested"
        / "output"
    )

    StandardizationRunner._prepare_directory(
        output
    )

    assert output.exists()
    assert output.is_dir()


def test_verify_cardinality_passes() -> None:
    """Matching cardinality should not raise."""

    StandardizationRunner._verify_cardinality(
        dataset_name="test",
        input_rows=10,
        output_rows=10,
    )


def test_verify_cardinality_fails() -> None:
    """Changed row cardinality should raise."""

    with pytest.raises(
        RuntimeError
    ):
        StandardizationRunner._verify_cardinality(
            dataset_name="test",
            input_rows=10,
            output_rows=9,
        )


def test_write_parquet(
    tmp_path: Path,
) -> None:
    """Runner should produce readable Parquet output."""

    dataframe = pd.DataFrame(
        {
            "value": [
                1,
                2,
            ]
        }
    )

    output = (
        tmp_path
        / "test.parquet"
    )

    StandardizationRunner._write_parquet(
        dataframe,
        output,
    )

    assert output.exists()

    recovered = pd.read_parquet(
        output
    )

    assert len(
        recovered
    ) == 2


def test_status_counts() -> None:
    """Status counts should be summarized correctly."""

    dataframe = pd.DataFrame(
        {
            "standardization_status": [
                "STANDARDIZED",
                "STANDARDIZED",
                "PARTIAL",
                "PRESERVED",
            ]
        }
    )

    result = (
        StandardizationRunner
        ._status_counts(
            dataframe
        )
    )

    assert result == {
        "STANDARDIZED": 2,
        "PARTIAL": 1,
        "PRESERVED": 1,
    }


def test_status_counts_without_column() -> None:
    """Missing status column should return an empty summary."""

    dataframe = pd.DataFrame(
        {
            "value": [
                1,
                2,
            ]
        }
    )

    assert (
        StandardizationRunner
        ._status_counts(
            dataframe
        )
        == {}
    )


def test_boolean_count() -> None:
    """Only true values should be counted."""

    dataframe = pd.DataFrame(
        {
            "usable": [
                True,
                False,
                True,
                None,
            ]
        }
    )

    assert (
        StandardizationRunner
        ._boolean_count(
            dataframe,
            "usable",
        )
        == 2
    )


def test_boolean_count_missing_column() -> None:
    """Absent usability field should return zero."""

    dataframe = pd.DataFrame(
        {
            "value": [
                1
            ]
        }
    )

    assert (
        StandardizationRunner
        ._boolean_count(
            dataframe,
            "usable",
        )
        == 0
    )


def test_gwas_chunked_standardization(
    tmp_path: Path,
) -> None:
    """GWAS execution should preserve cardinality across chunks."""

    runner = make_runner(
        tmp_path,
        chunk_size=2,
    )

    dataframe = pd.DataFrame(
        [
            {
                "STUDY ACCESSION":
                    "GCST1",
                "PUBMEDID":
                    "1",
                "DISEASE/TRAIT":
                    "Prostate cancer",
                "MAPPED_TRAIT":
                    "prostate carcinoma",
                "MAPPED_TRAIT_URI":
                    "test",
                "SNPS":
                    "rs1",
                "CHR_ID":
                    "1",
                "CHR_POS":
                    "100",
                "STRONGEST SNP-RISK ALLELE":
                    "rs1-A",
                "P-VALUE":
                    1e-8,
                "PVALUE_MLOG":
                    8,
                "MAPPED_GENE":
                    "GENE1",
                "REPORTED GENE(S)":
                    "GENE1",
            },
            {
                "STUDY ACCESSION":
                    "GCST2",
                "PUBMEDID":
                    "2",
                "DISEASE/TRAIT":
                    "Other trait",
                "MAPPED_TRAIT":
                    "other trait",
                "MAPPED_TRAIT_URI":
                    "test",
                "SNPS":
                    "rs2",
                "CHR_ID":
                    "2",
                "CHR_POS":
                    "200",
                "STRONGEST SNP-RISK ALLELE":
                    "rs2-G",
                "P-VALUE":
                    1e-6,
                "PVALUE_MLOG":
                    6,
                "MAPPED_GENE":
                    "GENE2",
                "REPORTED GENE(S)":
                    "GENE2",
            },
            {
                "STUDY ACCESSION":
                    "GCST3",
                "PUBMEDID":
                    "3",
                "DISEASE/TRAIT":
                    "Other trait",
                "MAPPED_TRAIT":
                    "other trait",
                "MAPPED_TRAIT_URI":
                    "test",
                "SNPS":
                    "rs3",
                "CHR_ID":
                    "3",
                "CHR_POS":
                    "300",
                "STRONGEST SNP-RISK ALLELE":
                    "rs3-C",
                "P-VALUE":
                    1e-5,
                "PVALUE_MLOG":
                    5,
                "MAPPED_GENE":
                    "GENE3",
                "REPORTED GENE(S)":
                    "GENE3",
            },
        ]
    )

    dataframe.to_csv(
        runner.inputs.gwas_catalog,
        sep="\t",
        index=False,
    )

    result = (
        runner.standardize_gwas()
    )

    assert result[
        "input_rows"
    ] == 3

    assert result[
        "output_rows"
    ] == 3

    assert (
        result[
            "cardinality_preserved"
        ]
        is True
    )

    assert result[
        "parts"
    ] == 2

    output_directory = Path(
        result[
            "output"
        ]
    )

    parquet_parts = sorted(
        output_directory.glob(
            "part-*.parquet"
        )
    )

    assert len(
        parquet_parts
    ) == 2

    combined = pd.concat(
        [
            pd.read_parquet(
                part
            )
            for part
            in parquet_parts
        ],
        ignore_index=True,
    )

    assert len(
        combined
    ) == 3


def test_missing_gwas_source(
    tmp_path: Path,
) -> None:
    """Absent GWAS source should fail explicitly."""

    runner = make_runner(
        tmp_path
    )

    with pytest.raises(
        FileNotFoundError
    ):
        runner.standardize_gwas()

def test_s3_trans_iso_mapping_uses_splicing_event() -> None:
    """S3 trans-iso-eQTL stores ENST identifiers in splicing_event."""

    specification = (
        StandardizationRunner
        .MORADI_S3[
            "Trans-iso-eQTL"
        ]
    )

    assert (
        specification[
            "feature_type"
        ]
        == "transcript"
    )

    assert (
        specification[
            "feature_column"
        ]
        == "splicing_event"
    )

    assert (
        specification[
            "variant_column"
        ]
        == "SNP"
    )