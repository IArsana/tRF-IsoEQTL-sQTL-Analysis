"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_6d_final_coloc_resolution.py

Description:
    Smoke-test entry point for M5.6D final colocalization resolution.

    M5.6D integrates the locked outputs from:
        M5.5
        M5.6A
        M5.6B
        M5.6C.1

    and determines whether any valid formal-colocalization route remains.

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

from pcatrfqtl.analysis.m5.runners.resolve_final_coloc_status import (
    M56DFinalColocResolutionRunner,
)


# ============================================================================
# Paths
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m5_final_coloc_resolution.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "final_coloc_resolution"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "qc"
)


# ============================================================================
# Validation
# ============================================================================


def validate_inputs() -> None:
    """Validate M5.6D configuration."""

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            "M5.6D configuration not found:\n"
            f"  {CONFIG_PATH}"
        )

    if not CONFIG_PATH.is_file():

        raise RuntimeError(
            "M5.6D config path is not a file:\n"
            f"  {CONFIG_PATH}"
        )


# ============================================================================
# Header
# ============================================================================


def print_header() -> None:
    """Print execution header."""

    print()
    print("=" * 72)
    print("M5.6D FINAL COLOCALIZATION RESOLUTION")
    print("=" * 72)

    print(
        f"Project root:       {ROOT}"
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

    print("=" * 72)


# ============================================================================
# Evidence chain
# ============================================================================


def print_evidence_chain(
    report: dict[str, Any],
) -> None:
    """Print locked colocalization decision chain."""

    evidence = report.get(
        "evidence_chain",
        [],
    )

    if not evidence:

        return

    print()
    print("Colocalization evidence chain:")
    print()

    for row in evidence:

        print(
            f"  {row.get('stage', 'NA')}"
        )

        print(
            "    Source:              "
            f"{row.get('evidence_source', 'NA')}"
        )

        print(
            "    State:               "
            f"{row.get('observed_state', 'NA')}"
        )

        print(
            "    Coloc route:         "
            f"{row.get('formal_coloc_route_available', False)}"
        )

        print(
            "    Biological negative: "
            f"{row.get('biological_negative_evidence', False)}"
        )

        print()


# ============================================================================
# Candidate summary
# ============================================================================


def print_candidate_status(
    report: dict[str, Any],
) -> None:
    """Print final per-candidate status."""

    candidates = report.get(
        "candidate_status",
        [],
    )

    if not candidates:

        return

    print()
    print("Candidate final status:")
    print()

    for candidate in candidates:

        print(
            f"  {candidate.get('lead_rsid', 'NA')}"
        )

        print(
            "    Final status:       "
            f"{candidate.get('final_coloc_status', 'NA')}"
        )

        print(
            "    Formal coloc ready: "
            f"{candidate.get('formal_coloc_ready', False)}"
        )

        print(
            "    Framework:          "
            f"{candidate.get('recommended_analysis_framework', 'NA')}"
        )

        print()


# ============================================================================
# Final summary
# ============================================================================


def print_result_summary(
    report: dict[str, Any],
) -> None:
    """Print final M5.6D result."""

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M5.6D RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidate leads assessed:          "
        f"{summary['candidate_leads_assessed']}"
    )

    print(
        "Formal coloc routes assessed:      "
        f"{summary['formal_coloc_routes_assessed']}"
    )

    print(
        "Formal coloc routes available:     "
        f"{summary['formal_coloc_routes_available']}"
    )

    print(
        "Valid dense QTL route available:   "
        f"{summary['valid_dense_qtl_route_available']}"
    )

    print()

    print(
        "Formal coloc-ready candidates:     "
        f"{summary['formal_coloc_ready_candidates']}"
    )

    print(
        "Formal coloc-blocked candidates:   "
        f"{summary['formal_coloc_blocked_candidates']}"
    )

    print()

    print(
        "FINAL COLOCALIZATION STATUS:"
    )

    print(
        f"  {summary['final_coloc_status']}"
    )

    print_evidence_chain(
        report
    )

    print_candidate_status(
        report
    )

    print("=" * 72)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M5.6D."""

    validate_inputs()

    print_header()

    runner = M56DFinalColocResolutionRunner(
        project_root=ROOT,
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
    print("M5.6D EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()