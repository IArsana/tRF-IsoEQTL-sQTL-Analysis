"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_6c_1_cancersplicingqtl.py

Description:
    Smoke-test entry point for M5.6C.1 CancerSplicingQTL PRAD download
    inspection.

    This script validates and inspects the actual downloaded PRAD
    CancerSplicingQTL dataset.

    The current expected raw artifact is:

        data/raw/m5/cancersplicingqtl/prad/PRAD_sQTLs.xlsx

    M5.6C.1 evaluates:
        - source-file discovery;
        - actual dataset row count;
        - canonical rsID coverage;
        - feature coverage;
        - source significance-filtering provenance;
        - availability of a complete tested SNP × feature matrix;
        - availability of dense unfiltered summary statistics;
        - beta and T-stat availability;
        - diagnostic standard-error derivation feasibility;
        - sample-size and allele-frequency metadata;
        - candidate lead presence;
        - formal colocalization readiness.

    Important safeguards:
        - Significant-only QTL records are not treated as dense summary data.
        - Missing candidate variants are not interpreted as null effects.
        - Diagnostic derived SE is not used for formal colocalization.
        - Cross-build coordinate joins are not performed.
        - No formal colocalization is performed.
        - No fine-mapping is performed.
        - No causal inference is performed.

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
from typing import Any

from pcatrfqtl.analysis.m5.runners.inspect_cancersplicingqtl import (
    M56C1CancerSplicingQTLRunner,
)


# ============================================================================
# Project paths
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m5_cancersplicingqtl_inspection.yaml"
)


RAW_DIRECTORY = (
    ROOT
    / "data"
    / "raw"
    / "m5"
    / "cancersplicingqtl"
    / "prad"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "cancersplicingqtl_inspection"
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
    """Validate M5.6C.1 configuration and expected directories."""

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            "M5.6C.1 configuration not found:\n"
            f"  {CONFIG_PATH}"
        )

    if not CONFIG_PATH.is_file():

        raise RuntimeError(
            "M5.6C.1 configuration path is not a file:\n"
            f"  {CONFIG_PATH}"
        )


# ============================================================================
# Formatting helpers
# ============================================================================


def _format_bool(
    value: Any,
) -> str:
    """Format boolean-like output consistently."""

    if value is None:

        return "NA"

    return str(
        bool(
            value
        )
    )


def _format_optional(
    value: Any,
) -> str:
    """Format optional scalar output."""

    if value is None:

        return "NA"

    return str(
        value
    )


def _format_fraction(
    value: Any,
    *,
    digits: int = 4,
) -> str:
    """Format fraction-like output."""

    if value is None:

        return "NA"

    try:

        return (
            f"{float(value):.{digits}f}"
        )

    except (
        TypeError,
        ValueError,
    ):

        return str(
            value
        )


# ============================================================================
# Header
# ============================================================================


def print_header() -> None:
    """Print M5.6C.1 execution header."""

    print()
    print("=" * 72)
    print("M5.6C.1 CANCERSPLICINGQTL PRAD DOWNLOAD INSPECTION")
    print("=" * 72)

    print(
        f"Project root:       {ROOT}"
    )

    print(
        f"Config:             {CONFIG_PATH}"
    )

    print(
        f"Expected raw dir:   {RAW_DIRECTORY}"
    )

    print(
        f"Output directory:   {OUTPUT_DIRECTORY}"
    )

    print(
        f"QC directory:       {QC_DIRECTORY}"
    )

    print("=" * 72)


# ============================================================================
# Candidate summary
# ============================================================================


def print_candidate_summary(
    report: dict[str, Any],
) -> None:
    """Print per-candidate CancerSplicingQTL status."""

    candidates = report.get(
        "candidate_status",
        [],
    )

    if not candidates:

        return

    print()
    print("Candidate lead status:")
    print()

    for candidate in candidates:

        lead = candidate.get(
            "lead_rsid",
            "NA",
        )

        present = candidate.get(
            "direct_rsid_present"
        )

        evidence = candidate.get(
            "candidate_evidence_status",
            "NA",
        )

        formal_ready = candidate.get(
            "candidate_formal_coloc_ready"
        )

        action = candidate.get(
            "recommended_action",
            "NA",
        )

        print(
            f"  {lead}"
        )

        print(
            "    Direct rsID present:       "
            f"{_format_bool(present)}"
        )

        print(
            "    Evidence status:           "
            f"{evidence}"
        )

        print(
            "    Formal coloc ready:        "
            f"{_format_bool(formal_ready)}"
        )

        print(
            "    Recommended action:        "
            f"{action}"
        )


# ============================================================================
# Result summary
# ============================================================================


def print_result_summary(
    report: dict[str, Any],
) -> None:
    """Print current M5.6C.1 QC summary."""

    summary = report.get(
        "summary",
        {},
    )

    print()
    print("=" * 72)
    print("M5.6C.1 RESULT SUMMARY")
    print("=" * 72)

    # ----------------------------------------------------------------------
    # Dataset not present
    # ----------------------------------------------------------------------

    if (
        summary.get(
            "classification"
        )
        ==
        "DATASET_NOT_PRESENT"
    ):

        print(
            "Source files discovered:             "
            f"{summary.get('source_files_discovered', 0)}"
        )

        print(
            "Classification:                       "
            "DATASET_NOT_PRESENT"
        )

        print(
            "Formal coloc ready:                   "
            "False"
        )

        print()
        print(
            "Expected CancerSplicingQTL PRAD input:"
        )

        print(
            f"  {RAW_DIRECTORY}"
        )

        print("=" * 72)

        return

    # ----------------------------------------------------------------------
    # Dataset identity
    # ----------------------------------------------------------------------

    print(
        "Source files discovered:             "
        f"{summary.get('source_files_discovered', 'NA')}"
    )

    print(
        "Source tables loaded:                "
        f"{summary.get('source_tables_loaded', 'NA')}"
    )

    print(
        "Association rows:                    "
        f"{summary.get('association_rows', 'NA')}"
    )

    print(
        "Genome build:                        "
        f"{summary.get('genome_build', 'NA')}"
    )

    print()

    # ----------------------------------------------------------------------
    # Variant / feature coverage
    # ----------------------------------------------------------------------

    print(
        "Canonical rsIDs:                     "
        f"{summary.get('canonical_rsids', 'NA')}"
    )

    print(
        "Canonical rsID fraction:             "
        f"{_format_fraction(summary.get('canonical_rsid_fraction'))}"
    )

    print(
        "Features:                            "
        f"{summary.get('features', 'NA')}"
    )

    print()

    # ----------------------------------------------------------------------
    # Source completeness
    # ----------------------------------------------------------------------

    print(
        "Source significance-filtered:        "
        f"{_format_bool(summary.get('source_significance_filtered'))}"
    )

    print(
        "Full tested variant-feature matrix:  "
        f"{_format_bool(summary.get(
            'complete_tested_variant_feature_matrix_available'
        ))}"
    )

    print(
        "Dense unfiltered statistics:         "
        f"{_format_bool(summary.get(
            'dense_unfiltered_summary_statistics_available'
        ))}"
    )

    print()

    # ----------------------------------------------------------------------
    # P-value diagnostics
    # ----------------------------------------------------------------------

    print(
        "Minimum observed p-value:            "
        f"{_format_optional(summary.get('minimum_p_value'))}"
    )

    print(
        "Maximum observed p-value:            "
        f"{_format_optional(summary.get('maximum_p_value'))}"
    )

    print(
        "P-values above 0.05:                 "
        f"{summary.get('p_values_above_0_05', 'NA')}"
    )

    print()

    # ----------------------------------------------------------------------
    # Effect-statistic availability
    # ----------------------------------------------------------------------

    print(
        "Beta directly reported:              "
        f"{_format_bool(summary.get('beta_directly_reported'))}"
    )

    print(
        "T-stat directly reported:            "
        f"{_format_bool(summary.get('t_statistic_directly_reported'))}"
    )

    print(
        "Standard error directly reported:    "
        f"{_format_bool(summary.get('standard_error_directly_reported'))}"
    )

    print(
        "Derived-SE diagnostic valid fraction:"
        f" {_format_fraction(
            summary.get('derived_se_diagnostic_valid_fraction')
        )}"
    )

    print()

    # ----------------------------------------------------------------------
    # Additional coloc metadata
    # ----------------------------------------------------------------------

    print(
        "Sample-size column available:        "
        f"{_format_bool(summary.get('sample_size_column_available'))}"
    )

    print(
        "Allele frequency available:          "
        f"{_format_bool(summary.get('allele_frequency_available'))}"
    )

    print(
        "Effect allele orientation verified:  "
        f"{_format_bool(summary.get('effect_allele_orientation_verified'))}"
    )

    print()

    # ----------------------------------------------------------------------
    # Final classification
    # ----------------------------------------------------------------------

    print(
        "Classification:                      "
        f"{summary.get('classification', 'NA')}"
    )

    print(
        "Formal coloc ready:                  "
        f"{_format_bool(summary.get('formal_coloc_ready'))}"
    )

    print(
        "coloc.abf ready:                     "
        f"{_format_bool(summary.get('coloc_abf_ready'))}"
    )

    print(
        "coloc.susie ready:                   "
        f"{_format_bool(summary.get('coloc_susie_ready'))}"
    )

    print_candidate_summary(
        report
    )

    print("=" * 72)


# ============================================================================
# Policy summary
# ============================================================================


def print_policy_summary(
    report: dict[str, Any],
) -> None:
    """Print critical scientific safeguards."""

    policy = report.get(
        "policy",
        {},
    )

    if not policy:

        return

    print()
    print("=" * 72)
    print("M5.6C.1 SCIENTIFIC POLICY")
    print("=" * 72)

    keys = [
        (
            "dataset_directly_inspected",
            "Dataset directly inspected",
        ),
        (
            "raw_dataset_modified",
            "Raw dataset modified",
        ),
        (
            "significant_only_qtl_allowed_for_coloc",
            "Significant-only QTL allowed for coloc",
        ),
        (
            "missing_variants_interpreted_as_null",
            "Missing variants interpreted as null",
        ),
        (
            "allele_orientation_inferred",
            "Allele orientation inferred",
        ),
        (
            "diagnostic_se_derivation_allowed",
            "Diagnostic SE derivation allowed",
        ),
        (
            "derived_se_used_for_formal_coloc",
            "Derived SE used for formal coloc",
        ),
        (
            "cross_build_coordinate_join_performed",
            "Cross-build coordinate join performed",
        ),
        (
            "formal_colocalization_performed",
            "Formal colocalization performed",
        ),
        (
            "fine_mapping_performed",
            "Fine-mapping performed",
        ),
        (
            "causal_inference_performed",
            "Causal inference performed",
        ),
    ]

    for key, label in keys:

        print(
            f"{label:<43} "
            f"{_format_bool(policy.get(key))}"
        )

    print("=" * 72)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M5.6C.1 smoke test."""

    validate_inputs()

    print_header()

    runner = M56C1CancerSplicingQTLRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print_result_summary(
        report
    )

    print_policy_summary(
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
    print("M5.6C.1 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()