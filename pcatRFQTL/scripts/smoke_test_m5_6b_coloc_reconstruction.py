"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m5_6b_coloc_reconstruction.py

Description:
    Smoke test and execution entry point for M5.6B colocalization
    reconstruction feasibility audit.

    This stage evaluates whether the full Moradi regulatory-QTL analysis can
    be reconstructed with sufficient fidelity for potential future formal
    colocalization.

    It does not execute the reconstruction itself.

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

from pcatrfqtl.analysis.m5.runners.assess_coloc_reconstruction import (
    M56BColocReconstructionRunner,
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
    / "m5_coloc_reconstruction.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m5"
    / "coloc_reconstruction"
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
    """Validate M5.6B inputs."""

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            "M5.6B config not found:\n"
            f"  {CONFIG_PATH}"
        )

    if not CONFIG_PATH.is_file():

        raise RuntimeError(
            "M5.6B config path is not a file:\n"
            f"  {CONFIG_PATH}"
        )


# ============================================================================
# Printing
# ============================================================================


def print_header() -> None:
    """Print execution header."""

    print()
    print("=" * 72)
    print("M5.6B COLOCALIZATION RECONSTRUCTION FEASIBILITY AUDIT")
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


def print_summary(
    report: dict,
) -> None:
    """Print compact reconstruction readiness summary."""

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M5.6B RESULT SUMMARY")
    print("=" * 72)

    print(
        f"Components assessed:              "
        f"{summary['components_assessed']}"
    )

    print(
        f"Required components:              "
        f"{summary['required_components']}"
    )

    print(
        f"Fully available components:       "
        f"{summary['fully_available_components']}"
    )

    print(
        f"Nonblocking components:           "
        f"{summary['nonblocking_components']}"
    )

    print(
        f"Not-identified components:        "
        f"{summary['not_identified_components']}"
    )

    print()

    print(
        f"Exact availability fraction:      "
        f"{summary['exact_availability_fraction']:.3f}"
    )

    print(
        f"Approx. nonblocking fraction:     "
        f"{summary['approximate_nonblocking_fraction']:.3f}"
    )

    print()

    print(
        f"Reconstruction readiness:         "
        f"{summary['reconstruction_readiness']}"
    )

    print(
        f"Reconstruction can proceed:       "
        f"{summary['reconstruction_can_proceed']}"
    )

    print(
        f"Exact replication claim allowed:  "
        f"{summary['exact_replication_claim_allowed']}"
    )

    print()

    print(
        f"Candidate leads assessed:         "
        f"{summary['candidate_leads_assessed']}"
    )

    print(
        f"Formal coloc-ready candidates:    "
        f"{summary['formal_coloc_ready_candidates']}"
    )

    print(
        f"Formal coloc-blocked candidates:  "
        f"{summary['formal_coloc_blocked_candidates']}"
    )

    print()

    print(
        "Blocking components:"
    )

    blocking = report.get(
        "blocking_components",
        [],
    )

    if blocking:

        for component in blocking:

            print(
                f"  - {component}"
            )

    else:

        print(
            "  - none"
        )

    print("=" * 72)


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Run M5.6B."""

    validate_inputs()

    print_header()

    runner = M56BColocReconstructionRunner(
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print_summary(
        report
    )

    print()
    print(
        "Full QC report:"
    )
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
    print("M5.6B EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()