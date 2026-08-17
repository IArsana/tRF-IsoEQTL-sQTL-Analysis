"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_2c_1_event_reconstruction.py

Description:
    Smoke test for M6.2C.1 FASE/Event Reconstruction Feasibility Audit.

    This smoke test reports whether the unresolved regulatory feature can
    be reconstructed using an exact or source-compatible FASE event
    namespace.

    A critical distinction is displayed explicitly:

        EVENT IDENTIFIER OBSERVED
            !=
        EVENT NAMESPACE REPRODUCIBLE

    Presence of INT98200 in Moradi source tables confirms only that the
    identifier exists in the published source output. It does not establish
    that the original FASE event-generation namespace or event-to-gene
    annotation can be reproduced.

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

from pcatrfqtl.analysis.m6.runners.assess_event_reconstruction import (
    M62C1EventReconstructionRunner,
)


# ============================================================================
# Paths
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[
    1
]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m6_event_reconstruction.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "event_reconstruction"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "qc"
)


# ============================================================================
# Validation
# ============================================================================


def validate_inputs() -> None:
    """Validate required M6.2C.1 inputs."""

    if not CONFIG_PATH.exists():

        raise FileNotFoundError(
            "M6.2C.1 configuration not found: "
            f"{CONFIG_PATH}"
        )

    if not CONFIG_PATH.is_file():

        raise RuntimeError(
            "M6.2C.1 configuration path is not a file: "
            f"{CONFIG_PATH}"
        )


# ============================================================================
# Header
# ============================================================================


def print_header() -> None:
    """Print execution header."""

    print()
    print("=" * 72)
    print("M6.2C.1 FASE/EVENT RECONSTRUCTION FEASIBILITY AUDIT")
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
# Reconstruction components
# ============================================================================


def print_components(
    report: dict[str, Any],
) -> None:
    """Print reconstruction-component audit."""

    components = report.get(
        "reconstruction_components",
        [],
    )

    print()
    print("Reconstruction components:")
    print()

    if not components:

        print(
            "  No reconstruction components were reported."
        )

        return

    for row in components:

        print(
            f"  {row['component_id']}"
        )

        print(
            "    Label:                     "
            f"{row['component_label']}"
        )

        print(
            "    Required:                  "
            f"{row['required']}"
        )

        print(
            "    State:                     "
            f"{row['component_state']}"
        )

        print(
            "    Local artifact available:  "
            f"{row['local_artifact_available']}"
        )

        print(
            "    Local artifact count:      "
            f"{row['local_artifact_count']}"
        )

        # ------------------------------------------------------------------
        # Event-namespace-specific reporting
        # ------------------------------------------------------------------

        if row.get(
            "component_id"
        ) == "event_namespace":

            print()

            print(
                "    Target event:              "
                f"{row.get('target_event_identifier')}"
            )

            print(
                "    Target event observed:     "
                f"{row.get('target_event_observed')}"
            )

            print(
                "    Namespace reproducible:    "
                f"{row.get('namespace_reproducible')}"
            )

            observed_locations = row.get(
                "target_event_observed_locations",
                [],
            )

            if observed_locations:

                print(
                    "    Observed locations:"
                )

                for location in observed_locations:

                    print(
                        f"      - {location}"
                    )

        # ------------------------------------------------------------------
        # Generic artifact paths
        # ------------------------------------------------------------------

        artifact_paths = row.get(
            "local_artifact_paths",
            [],
        )

        if artifact_paths:

            print(
                "    Local artifact paths:"
            )

            for path in artifact_paths:

                print(
                    f"      - {path}"
                )

        print()


# ============================================================================
# Feature readiness
# ============================================================================


def print_feature_readiness(
    report: dict[str, Any],
) -> None:
    """Print feature-level reconstruction readiness."""

    rows = report.get(
        "feature_readiness",
        [],
    )

    print()
    print("Feature reconstruction readiness:")
    print()

    if not rows:

        print(
            "  No features were forwarded to M6.2C.1."
        )

        return

    for row in rows:

        print(
            f"  {row['regulatory_feature_id']}"
        )

        print(
            "    Lead candidate:                  "
            f"{row['lead_rsid']}"
        )

        print(
            "    Priority rank:                   "
            f"{row['priority_rank']}"
        )

        print(
            "    Priority class:                  "
            f"{row['priority_class']}"
        )

        print(
            "    Feature class:                   "
            f"{row['feature_class']}"
        )

        print(
            "    Upstream annotation status:      "
            f"{row['upstream_gene_annotation_status']}"
        )

        print()

        print(
            "    Reconstruction classification:   "
            f"{row['reconstruction_classification']}"
        )

        print(
            "    Exact reconstruction allowed:    "
            f"{row['exact_reconstruction_allowed']}"
        )

        print(
            "    Source-compatible allowed:       "
            f"{row['source_compatible_reconstruction_allowed']}"
        )

        print(
            "    Parent-gene assignment allowed:  "
            f"{row['parent_gene_assignment_from_reconstruction_allowed']}"
        )

        print()

        print(
            "    Parent gene:                     "
            f"{row['parent_gene']}"
        )

        print(
            "    Parent gene assigned:            "
            f"{row['parent_gene_assigned']}"
        )

        print(
            "    Nearest-gene mapping:            "
            f"{row['nearest_gene_assignment_performed']}"
        )

        print(
            "    Approximate mapping used:        "
            f"{row['approximate_mapping_used']}"
        )

        print()

        print(
            "    Next action:                     "
            f"{row['next_action']}"
        )

        print(
            "    Causal claim allowed:            "
            f"{row['causal_claim_allowed']}"
        )

        print()


# ============================================================================
# Summary
# ============================================================================


def print_summary(
    report: dict[str, Any],
) -> None:
    """Print M6.2C.1 result summary."""

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M6.2C.1 RESULT SUMMARY")
    print("=" * 72)

    print(
        "Features assessed:                    "
        f"{summary['features_assessed']}"
    )

    print(
        "Reconstruction components assessed:   "
        f"{summary['reconstruction_components_assessed']}"
    )

    print(
        "Required components:                  "
        f"{summary['required_components']}"
    )

    print(
        "Required components available:        "
        f"{summary['required_components_available']}"
    )

    print(
        "Required components unavailable:      "
        f"{summary['required_components_unavailable']}"
    )

    print()

    print(
        "Target event identifier observed:      "
        f"{summary['target_event_identifier_observed']}"
    )

    print(
        "Event namespace reproducible:          "
        f"{summary['event_namespace_reproducible']}"
    )

    print()

    print(
        "Reconstruction classification:"
    )

    print(
        f"  {summary['reconstruction_classification']}"
    )

    print()

    print(
        "Exact reconstruction allowed:         "
        f"{summary['exact_reconstruction_allowed']}"
    )

    print(
        "Source-compatible reconstruction:     "
        f"{summary['source_compatible_reconstruction_allowed']}"
    )

    print(
        "Features ready for parent gene:        "
        f"{summary['features_ready_for_parent_gene_assignment']}"
    )

    print()

    print(
        "Parent gene assignment performed:      "
        f"{summary['parent_gene_assignment_performed']}"
    )

    print(
        "Nearest-gene mapping performed:        "
        f"{summary['nearest_gene_mapping_performed']}"
    )

    print(
        "Approximate mapping used:              "
        f"{summary['approximate_mapping_used']}"
    )

    print(
        "Causal inference performed:            "
        f"{summary['causal_inference_performed']}"
    )

    print_components(
        report
    )

    print_feature_readiness(
        report
    )

    print()
    print("Next stage:")

    print(
        f"  {summary['next_stage']}"
    )

    print()
    print("=" * 72)


# ============================================================================
# Full QC
# ============================================================================


def print_full_qc(
    report: dict[str, Any],
) -> None:
    """Print standards-compliant full QC report."""

    print()
    print("Full QC report:")
    print()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
    )


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Execute M6.2C.1 smoke test."""

    validate_inputs()

    print_header()

    runner = M62C1EventReconstructionRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    print_summary(
        report
    )

    print_full_qc(
        report
    )

    print()
    print("=" * 72)
    print("M6.2C.1 EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()