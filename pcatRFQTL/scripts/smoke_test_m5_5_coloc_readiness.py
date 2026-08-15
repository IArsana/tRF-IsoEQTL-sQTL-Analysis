"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_5_coloc_readiness.py

Description:
    Smoke test and execution entry point for M5.5 colocalization readiness
    assessment.

    This script evaluates whether the currently available prostate cancer
    GWAS summary statistics and Moradi regulatory-QTL resources are suitable
    for formal colocalization analysis.

    Inputs:
        Disease-side GWAS:
            data/processed/m5/sumstats/standardized/
            gwas_locus_standardized.parquet

        Regulatory-side Moradi QTL index:
            data/processed/m4/moradi_qtl_index/

        Configuration:
            configs/m5_coloc_readiness.yaml

    The Moradi QTL source is intentionally passed as a directory because
    M4.1 stores the standardized QTL index as multiple partitioned Parquet
    artifacts.

    Outputs:
        data/processed/m5/coloc_readiness/
            gwas_coloc_readiness.parquet
            qtl_coloc_readiness.parquet
            coloc_pair_readiness.parquet

        data/processed/m5/qc/
            m5_5_coloc_readiness.json

    Scientific safeguards:
        - This stage assesses readiness only.
        - No formal colocalization is performed.
        - No fine-mapping is performed.
        - No causal inference is performed.
        - Significant-only QTL evidence is not treated as dense regional
          summary statistics.
        - Cross-build coordinate joins are not performed.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import json
from pathlib import Path

from pcatrfqtl.analysis.m5.runners.assess_coloc_readiness import (
    M55ColocReadinessRunner,
)


# ============================================================================
# Project root
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[1]


# ============================================================================
# Input paths
# ============================================================================


GWAS_PATH = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "sumstats"
    / "standardized"
    / "gwas_locus_standardized.parquet"
)


QTL_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m4"
    / "moradi_qtl_index"
)


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m5_coloc_readiness.yaml"
)


# ============================================================================
# Output paths
# ============================================================================


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "coloc_readiness"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "qc"
)


# ============================================================================
# Input validation
# ============================================================================


def validate_inputs() -> None:
    """Validate required M5.5 source artifacts before execution."""

    if not GWAS_PATH.exists():

        raise FileNotFoundError(
            "M5.5 GWAS input was not found:\n"
            f"  {GWAS_PATH}"
        )

    if not GWAS_PATH.is_file():

        raise RuntimeError(
            "M5.5 GWAS input exists but is not a file:\n"
            f"  {GWAS_PATH}"
        )

    if not QTL_DIRECTORY.exists():

        raise FileNotFoundError(
            "M5.5 Moradi QTL directory was not found:\n"
            f"  {QTL_DIRECTORY}"
        )

    if not QTL_DIRECTORY.is_dir():

        raise RuntimeError(
            "M5.5 Moradi QTL input exists but is not a directory:\n"
            f"  {QTL_DIRECTORY}"
        )

    qtl_partitions = sorted(
        QTL_DIRECTORY.glob(
            "*.parquet"
        )
    )

    if not qtl_partitions:

        raise FileNotFoundError(
            "No Moradi QTL Parquet partitions were found under:\n"
            f"  {QTL_DIRECTORY}"
        )

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            "M5.5 configuration file was not found:\n"
            f"  {CONFIG_PATH}"
        )

    if not CONFIG_PATH.is_file():

        raise RuntimeError(
            "M5.5 configuration path exists but is not a file:\n"
            f"  {CONFIG_PATH}"
        )


# ============================================================================
# Input summary
# ============================================================================


def print_input_summary() -> None:
    """Print deterministic input summary before execution."""

    qtl_partitions = sorted(
        QTL_DIRECTORY.glob(
            "*.parquet"
        )
    )

    print()
    print("=" * 72)
    print("M5.5 COLOCALIZATION READINESS ASSESSMENT")
    print("=" * 72)

    print(
        f"Project root:       {ROOT}"
    )

    print(
        f"GWAS input:         {GWAS_PATH}"
    )

    print(
        f"Moradi QTL index:   {QTL_DIRECTORY}"
    )

    print(
        f"QTL partitions:     {len(qtl_partitions)}"
    )

    print(
        f"Config:             {CONFIG_PATH}"
    )

    print(
        f"Output directory:   {OUTPUT_DIRECTORY}"
    )

    print(
        f"QC directory:       {QC_DIRECTORY}"
    )

    print()
    print(
        "Moradi QTL partitions:"
    )

    for path in qtl_partitions:

        print(
            f"  - {path.name}"
        )

    print("=" * 72)


# ============================================================================
# Result summary
# ============================================================================


def print_result_summary(
    report: dict,
) -> None:
    """Print compact M5.5 result summary."""

    summary = report.get(
        "summary",
        {},
    )

    policy = report.get(
        "policy",
        {},
    )

    print()
    print("=" * 72)
    print("M5.5 RESULT SUMMARY")
    print("=" * 72)

    print(
        f"GWAS rows:                         "
        f"{summary.get('gwas_rows')}"
    )

    print(
        f"Moradi QTL rows:                   "
        f"{summary.get('qtl_rows')}"
    )

    print(
        f"GWAS loci assessed:                "
        f"{summary.get('gwas_loci_assessed')}"
    )

    print(
        f"QTL units assessed:                "
        f"{summary.get('qtl_units_assessed')}"
    )

    print(
        f"Study-lead pairs assessed:         "
        f"{summary.get('study_lead_pairs_assessed')}"
    )

    print()

    print(
        f"GWAS basic vector-ready units:     "
        f"{summary.get('gwas_basic_coloc_vector_ready_units')}"
    )

    print(
        f"QTL basic vector-ready units:      "
        f"{summary.get('qtl_basic_coloc_vector_ready_units')}"
    )

    print(
        f"QTL full-locus units:              "
        f"{summary.get('qtl_full_locus_summary_units')}"
    )

    print(
        f"QTL dense-locus units:             "
        f"{summary.get('qtl_dense_locus_units')}"
    )

    print()

    print(
        f"Formal coloc-ready pairs:          "
        f"{summary.get('formal_coloc_ready_pairs')}"
    )

    print(
        f"Formal coloc-blocked pairs:        "
        f"{summary.get('formal_coloc_blocked_pairs')}"
    )

    print(
        f"coloc.abf-ready pairs:             "
        f"{summary.get('coloc_abf_ready_pairs')}"
    )

    print(
        f"coloc.susie-ready pairs:           "
        f"{summary.get('coloc_susie_ready_pairs')}"
    )

    print()

    print(
        f"Total shared-variant observations: "
        f"{summary.get('shared_variant_observations_total')}"
    )

    print(
        f"Maximum shared variants per pair:  "
        f"{summary.get('maximum_shared_variants_per_pair')}"
    )

    print()

    print(
        "Formal colocalization performed:   "
        f"{policy.get('formal_colocalization_performed')}"
    )

    print(
        "Fine-mapping performed:            "
        f"{policy.get('fine_mapping_performed')}"
    )

    print(
        "Causal inference performed:        "
        f"{policy.get('causal_inference_performed')}"
    )

    print("=" * 72)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Run M5.5 colocalization readiness assessment."""

    validate_inputs()

    print_input_summary()

    runner = M55ColocReadinessRunner(
        gwas_path=GWAS_PATH,
        qtl_directory=QTL_DIRECTORY,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print_result_summary(
        report
    )

    print()
    print("Full QC report:")
    print()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        )
    )

    print()
    print("=" * 72)
    print("M5.5 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()